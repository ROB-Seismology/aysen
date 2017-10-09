# -*- coding: iso-Latin-1 -*-

import os
import numpy as np
import pylab
import openquake.hazardlib as oqhazlib
import hazard.rshalib as rshalib
from rupture_probabilities import *


#project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")
fig_folder = os.path.join(project_folder, "Figures", "Sensitivity", "v4")


## Thresholds (depends on type of evidence)
## Note: thresholds chosen such that they are not all positive or negative for chosen scenario

## NE: 7.5, 8.5
## PE: 5.5, 6.0, 6.5 7.5
threshold_mmi = 7.5


## Note: it is not possible to apply delta_threshold, as a particular site
## may be both positive and negative if delta_threshold > 0 !!
delta_threshold = 0.


## Parameters
truncation_level = 1.0
soil_params = rshalib.site.REF_SOIL_PARAMS
polygon_discretization = 2.5
imt = oqhazlib.imt.MMI()
strict_intersection = True


## Map parameters
#map_region = (-74, -72, -46, -44.5)
map_region = (-74, -72, -46.25, -44.75)
output_format = "png"


## IPE Logic tree
ipe_models = ["Barrientos1980WithSigma", "BakunWentworth1997WithSigma", "AllenEtAl2012", "AtkinsonWald2007"]
weights = [0.32, 0.26, 0.25, 0.17]
lt_gmpe_system_def = {TRT: rshalib.pmf.GMPEPMF(ipe_models, weights)}
lt_integration_distance_dict = {"AtkinsonWald2007": (None, 30)}
#TODO: apply distance filtering for BakunWentworth1997WithSigma when M larger than
#threshold causing too high values in near field (to be determined)
#lt_integration_distance_dict = {"AtkinsonWald2007": (None, 30), "BakunWentworth1997WithSigma": (20, None)}


## Read site info
all_site_models = []
for geom_type in ["Polygons", "Points"]:
	shapefile = os.path.join(gis_folder, "%s.shp" % geom_type)
	site_models = read_evidence_sites_from_gis(shapefile, polygon_discretization)
	for site_model in site_models:
		## Cuervo opp is part of Centre
		if site_model.name != "Cuervo opp":
			all_site_models.append(site_model)
print len(all_site_models)


## Discretzie faults as network
fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel3.TAB")
dM = 0.2
fault_mags, fault_networks = [], []
for M, source_model in read_fault_source_model_as_network(fault_filespec, dM=dM):
	fault_mags.append(M)
	fault_networks.append(source_model)


## Loop over IPEs
for ipe_name in ipe_models[:1] + ipe_models[2:]:
	print ipe_name
	if ipe_name != "LogicTree":
		trt_gsim_dict = {TRT: ipe_name}
		gmpe_system_def = {TRT: rshalib.pmf.GMPEPMF([ipe_name], [1])}
		integration_distance_dict = {}
	else:
		gmpe_system_def = lt_gmpe_system_def
		integration_distance_dict = lt_integration_distance_dict

	## Loop over scenario magnitudes and corresponding fault networks
	for scen_mag, scen_flt_network in zip(fault_mags, fault_networks):
		#if not (6.2 <= scen_mag <= 6.3):
		#if scen_mag < 6.5:
		#	continue
		print("M=%.2f" % scen_mag)
		respow_dict = {}

		#src_id = "2#01+2#02+2#03+2#04+2#05"
		#sources = [src for src in scen_flt_network if src.source_id == src_id]
		#scen_flt_network = rshalib.source.SourceModel(src_id, sources)

		## Loop over rupture scenarios
		for s, scenflt in enumerate(scen_flt_network):
			scentroid = scenflt.get_surface()._get_top_edge_centroid()

			## Set characteristic magnitude from MSR and fault area
			msr = getattr(oqhazlib.scalerel, MSR)()
			char_mag = msr.get_median_mag(scenflt.get_rupture().surface.get_area(), scenflt.rake)
			scenflt.mfd.min_mag = np.round(char_mag, 1)
			scenflt.mfd.char_mag = scenflt.mfd.min_mag + scenflt.mfd.bin_width/2.
			scenario_src_model = rshalib.source.SourceModel("", [scenflt])

			## Run DSHA model and determine pe_site_models and ne_site_models
			## based on IPE logic tree
			pe_thresholds, pe_site_models, ne_thresholds, ne_site_models = [], [], [], []
			pe_threshold = threshold_mmi - delta_threshold
			ne_threshold = threshold_mmi + delta_threshold

			for site_model in all_site_models:
				dsha_model = rshalib.shamodel.DSHAModel("", scenario_src_model, lt_gmpe_system_def,
							grid_outline=[], grid_spacing=None,
							soil_site_model=site_model, imt_periods={"MMI": [0]},
							truncation_level=0, integration_distance=500)

				uhs_field = dsha_model.calc_gmf_fixed_epsilon()
				hm = uhs_field.getHazardMap()
				if (hm.intensities > pe_threshold).all():
					pe_site_models.append(site_model)
					pe_thresholds.append(pe_threshold)
				elif (hm.intensities <= ne_threshold).all():
					ne_site_models.append(site_model)
					ne_thresholds.append(ne_threshold)

			#print("Positive: n=%d; Negative: n=%d" %(len(pe_site_models), len(ne_site_models)))

			## Compute probability for this scenario and threshold
			## using same IPE (not the logic tree!)
			prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
								scenario_src_model, gmpe_system_def, imt, pe_site_models,
								pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
								integration_distance_dict=integration_distance_dict,
								strict_intersection=strict_intersection)
			[scenario_prob] = prob_dict.values()[0]
			#print("Scenario prob.: %s" % scenario_prob)


			## Perform test
			## Second loop over all possible ruptures
			all_probs = []
			distances, mag_diffs = [], []
			for M, source_model in zip(fault_mags, fault_networks):
				for flt in source_model:
					mag_diffs.append(np.abs(scen_mag - M))
					centroid = flt.get_surface()._get_top_edge_centroid()
					distances.append(scentroid.distance(centroid))

				## Compute rupture probabilities
				prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
									source_model, gmpe_system_def, imt, pe_site_models,
									pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
									integration_distance_dict=integration_distance_dict,
									strict_intersection=strict_intersection)

				probs = np.array(prob_dict.values())
				probs = probs[:, 0]
				all_probs.extend(probs)


			## Compute resolution power
			distances = np.array(distances)
			max_distance = 175. ## use same distance for each case!
			distances /= max_distance
			mag_diffs = np.array(mag_diffs)
			mag_diffs /= (fault_mags[-1] - fault_mags[0])
			dx = np.sqrt(0.5 * ((1-mag_diffs)**2 + (1-distances)**2))
			## Reduce impact of probabilities slightly larger than scenario prob
			## for nearby ruptures with similar magnitude
			prob_diffs = scenario_prob - np.array(all_probs)
			prob_margin = 0.05
			idxs = np.where((-prob_margin <= prob_diffs) & (prob_diffs < 0) & (mag_diffs <= 0.1) & (distances <= 0.06))
			prob_diffs[idxs] = 0
			N = len(prob_diffs)
			## N - 1 because exact scenario (which should have zero prob) is in the results as well
			#res_pow = np.sum(prob_diffs * dx) / (scenario_prob * (N - 1))
			res_pow = np.sum(prob_diffs * dx) / (N - 1)
			if s % 10 == 0:
				print s, res_pow
			respow_dict[scenflt.source_id] = [res_pow]

			#print scenario_prob, res_pow
			#for n in range(N):
			#	print mag_diffs[n], distances[n], prob_diffs[n], dx[n], (prob_diffs * dx)[n]


		## Plot map
		if "WithSigma" in ipe_name:
			ipe_label = ipe_name[:ipe_name.find("WithSigma")]
		else:
			ipe_label = ipe_name
		text_box = "Scenario M:%.2f\nThreshold MMI: %s\nIPE: %s"
		roman_intensity = get_roman_intensity(threshold_mmi)
		text_box %= (scen_mag, roman_intensity, ipe_label)

		title = ""
		fig_filename = "ResPow_%s_M=%.2f_MMI=%s.%s" % (ipe_name, scen_mag, threshold_mmi, output_format)
		#fig_filename = "ResPow_test.png"
		fig_filespec = os.path.join(fig_folder, fig_filename)
		#fig_filespec = None

		## Colormaps: RdBu_r, YlOrRd, BuPu, RdYlBu_r, Greys
		site_model_gis_file = os.path.join(gis_folder, "Polygons.shp")
		plot_rupture_probabilities(scen_flt_network, respow_dict, [], [],
									map_region, plot_point_ruptures=True, colormap="RdYlBu_r",
									title=title, text_box=text_box, site_model_gis_file=site_model_gis_file,
									legend_label="Sensitivity", prob_max=0.75,
									highlight_max_prob_section=True,
									max_prob_mag_precision=0,
									neutral_site_models=all_site_models,
									fig_filespec=fig_filespec)