# -*- coding: iso-Latin-1 -*-

import os, sys
import numpy as np
import pylab
import openquake.hazardlib as oqhazlib
import hazard.rshalib as rshalib
from hazard.rshalib.source_estimation import calc_rupture_probability_from_ground_motion_thresholds
from eqcatalog.macro import get_roman_intensity

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SCRIPT_DIR)
from aysenlib import (project_folder, gis_folder, read_evidence_sites_from_gis,
						read_fault_source_model, read_fault_source_model_as_floating_ruptures,
						read_fault_source_model_as_network,
						plot_rupture_probabilities, TRT, MSR)
from create_animated_gif import create_animated_gif


#fig_folder = os.path.join(project_folder, "Figures", "Sensitivity", "v5", "Scenarios")
fig_folder = r"C:\Temp"


## Scenarios
#scenario, threshold_mmi = "Quitralco", 8
#scenario, threshold_mmi = "Azul Tigre South", 8
#scenario, threshold_mmi = "2007", 8
#scenario, threshold_mmi = "Due East", 7
#scenario, threshold_mmi = "Due West", 7


## Thresholds (depends on type of evidence)
## Note: thresholds chosen such that they are not all positive or negative for chosen scenario

## NE: 7.5, 8.5
## PE: 5.5, 6.0, 6.5 7.5
#threshold_mmi = 7


## Parameters
truncation_level = 1.0
soil_params = rshalib.site.REF_SOIL_PARAMS
polygon_discretization = 2.5
imt = oqhazlib.imt.MMI()
strict_intersection = True


## Map parameters
#map_region = (-74, -72, -46, -44.5)
map_region = (-74, -72, -46.25, -44.8)
output_format = "pdf"


## IPE Logic tree
#models = ["Barrientos1980WithSigma", "BakunWentworth1997WithSigma", "AllenEtAl2012", "AtkinsonWald2007"]
#weights = [0.32, 0.26, 0.25, 0.17]
ipe_models = ["Barrientos1980WithSigma", "BakunWentworth1997WithSigma", "AllenEtAl2012"]
weights = [0.34, 0.33, 0.33]

lt_gmpe_system_def = {TRT: rshalib.pmf.GMPEPMF(ipe_models, weights)}
lt_integration_distance_dict = {"AtkinsonWald2007": (None, 30)}
#TODO: apply distance filtering for BakunWentworth1997WithSigma when M larger than
#threshold causing too high values in near field (to be determined)
#lt_integration_distance_dict = {"AtkinsonWald2007": (None, 30), "BakunWentworth1997WithSigma": (20, None)}


## Read site info
all_site_models = []
for geom_type in ["Polygons_v3", "Points"]:
	shapefile = os.path.join(gis_folder, "%s.shp" % geom_type)
	site_models = read_evidence_sites_from_gis(shapefile, polygon_discretization)
	for site_model in site_models:
		## Cuervo opp is part of Centre
		if site_model.name != "Cuervo opp":
			all_site_models.append(site_model)
print(len(all_site_models))


#scenarios = ["Quitralco", "Azul Tigre South", "2007", "Due East", "Due West"]
#scenarios = ["Azul Tigre South", "Quitralco East", "Due East"]
#scenarios = ["Azul Tigre South", "Quitralco East", "Rio Manihuales", "Quitralco West"]
scenarios = ["Quitralco West"]
for scenario in scenarios:
	## Read rupture scenario
	scenario_filespec = os.path.join(gis_folder, 'LOFZ_rupture_scenarios.TAB')
	flt_src_model = read_fault_source_model(scenario_filespec)
	[scenflt] = [flt for flt in flt_src_model if flt.name == scenario]
	aspect_ratio = scenflt.get_length() / scenflt.width
	scentroid = scenflt.get_surface()._get_top_edge_centroid()


	## Set characteristic magnitude from MSR and fault area
	msr = getattr(oqhazlib.scalerel, MSR)()
	char_mag = msr.get_median_mag(scenflt.get_rupture().surface.get_area(), scenflt.rake)
	print(scenflt.get_length(), scenflt.width, char_mag)
	scenflt.mfd.min_mag = np.round(char_mag, 1)
	scenflt.mfd.char_mag = scenflt.mfd.min_mag + scenflt.mfd.bin_width/2.
	scenario_src_model = rshalib.source.SourceModel(scenario, [scenflt])

	## Note: it is not possible to apply delta_threshold, as a particular site
	## may be both positive and negative if delta_threshold > 0 !!
	#delta_threshold = 0.5
	delta_threshold = 0.

	#thresholds = [7.5, 6.5, 5.5]
	#threshold_sets = [[t] for t in thresholds] + [thresholds]
	#threshold_sets = [[9]]
	threshold_sets = [[7.5]]
	for threshold_set in threshold_sets:
		## Construct ground-motion model
		#ipe_name = "BakunWentworth1997WithSigma"
		#ipe_name = "AtkinsonWald2007"
		#ipe_name = "AllenEtAl2012"
		#ipe_name = "LogicTree"

		#ipe_names = ["LogicTree", "AllenEtAl2012", "AtkinsonWald2007", "BakunWentworth1997WithSigma", "Barrientos1980WithSigma"]
		ipe_names = ["LogicTree"] + ipe_models
		#ipe_names = ipe_models[:1]
		#ipe_names = ipe_names[3:4]
		#ipe_names = ["BakunWentworth1997WithSigma"]
		max_probs, section_probs, scenario_probs = {}, {}, {}
		for ipe_name in ipe_names:
			if ipe_name != "LogicTree":
				trt_gsim_dict = {TRT: ipe_name}
				gmpe_system_def = {TRT: rshalib.pmf.GMPEPMF([ipe_name], [1])}
				integration_distance_dict = {}
			else:
				gmpe_system_def = lt_gmpe_system_def
				integration_distance_dict = lt_integration_distance_dict

			## Run DSHA model and determine pe_site_models and ne_site_models
			## based on IPE logic tree
			pe_thresholds, pe_site_models, ne_thresholds, ne_site_models = [], [], [], []
			for threshold_mmi in threshold_set:
				pe_threshold = threshold_mmi - delta_threshold
				ne_threshold = threshold_mmi + delta_threshold

				for site_model in all_site_models:
					dsha_model = rshalib.shamodel.DSHAModel(scenario, scenario_src_model, lt_gmpe_system_def,
								site_model=site_model, imt_periods={"MMI": [0]},
								truncation_level=0, integration_distance=500)

					uhs_field = dsha_model.calc_gmf_fixed_epsilon()
					hm = uhs_field.get_hazard_map()
					if (hm.intensities > pe_threshold).all():
						pe_site_models.append(site_model)
						pe_thresholds.append(pe_threshold)
					elif (hm.intensities <= ne_threshold).all():
						ne_site_models.append(site_model)
						ne_thresholds.append(ne_threshold)

			print("Positive: n=%d; Negative: n=%d" %(len(pe_site_models), len(ne_site_models)))

			#if ipe_name == "BakunWentworth1997WithSigma" or not 0 in (len(pe_site_models), len(ne_site_models)):
			if True:
				## Compute probability for this scenario and threshold
				prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
									scenario_src_model, lt_gmpe_system_def, imt, pe_site_models,
									pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
									integration_distance_dict=lt_integration_distance_dict,
									strict_intersection=strict_intersection)
				[scenario_prob] = list(prob_dict.values())[0]
				scenario_probs[ipe_name] = scenario_prob
				print("Scenario prob.: %s" % scenario_prob)


				max_probs[ipe_name], section_probs[ipe_name] = [], []

				## Perform test

				## Discretize faults into floating ruptures
				"""
				## Magnitude range to test
				dM = 0.25
				Mrange = np.arange(scenflt.mfd.min_mag - 1, scenflt.mfd.min_mag + 1 + dM/2, dM)
				dM = 0.5
				Mrange = np.linspace(scenflt.mfd.min_mag - dM, scenflt.mfd.min_mag + dM, 3)
				#Mrange = [scenflt.mfd.min_mag]
				Mrange = [6.9]
				print Mrange

				for M in Mrange:
					## Read fault source model
					fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel2.TAB")
					source_model = read_fault_source_model_as_floating_ruptures(fault_filespec,
										M - dM/2, M, dM, depth=0.1, aspect_ratio=aspect_ratio)
					print M, len(source_model)
					if len(source_model) == 0:
						print("Warning: no sources for M=%s!" % M)
						max_probs[ipe_name].append(np.nan)
						section_probs[ipe_name].append(np.nan)
						continue
				"""

				## Discretzie faults as network
				dM = 0.2
				Mrange = []
				fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel4.TAB")
				for M, source_model in read_fault_source_model_as_network(fault_filespec, dM=dM):
					## 6.2, 6.8 and 7.2 crash with logictree!
					#if M <= 7.1:
					#	continue
					Mrange.append(M)

					## Find fault section closest to scenario rupture
					distances = []
					for flt in source_model:
						centroid = flt.get_surface()._get_top_edge_centroid()
						distances.append(scentroid.distance(centroid))
					section_idx = np.argmin(distances)
					section_name = source_model.sources[section_idx].source_id

					## Compute rupture probabilities
					prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
										source_model, gmpe_system_def, imt, pe_site_models,
										pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
										integration_distance_dict=integration_distance_dict,
										strict_intersection=strict_intersection)

					probs = np.array(list(prob_dict.values()))
					probs = probs[:, 0]
					max_prob_idx = probs.argmax()
					max_prob = probs[max_prob_idx]
					rup_name = list(prob_dict.keys())[max_prob_idx]
					print(M, rup_name, max_prob)
					max_probs[ipe_name].append(max_prob)

					for rup_name, prob in prob_dict.items():
						if rup_name == section_name:
							section_probs[ipe_name].append(prob)
						if prob >= scenario_prob:
							print("  ", M, rup_name, prob)


					## Plot map
					#if ipe_name in ("LogicTree", "BakunWentworth1997WithSigma", "AllenEtAl2012"):
					#if ipe_name in ("AllenEtAl2012", "AtkinsonWald2007"):
					if True:
						if "WithSigma" in ipe_name:
							ipe_label = ipe_name[:ipe_name.find("WithSigma")]
						elif ipe_name == "LogicTree":
							ipe_label = "Mix"
						else:
							ipe_label = ipe_name
						text_box = "Scenario: %s\nThreshold MMI: %s\nIPE: %s\nM: %.2f, Pmax: %.2f"
						roman_intensity = ', '.join([get_roman_intensity(mmi) for mmi in threshold_set])
						text_box %= (scenario, roman_intensity, ipe_label, M, max_prob)

						title = ""
						if len(threshold_set) == 1:
							intensity = threshold_set[0]
						else:
							intensity = "mix"
						fig_filename = "%s_MMI=%s_%s_M=%.2f.%s" % (scenario, intensity, ipe_name, M, output_format)
						fig_filespec = os.path.join(fig_folder, fig_filename)
						#fig_filespec = None

						## Colormaps: RdBu_r, YlOrRd, BuPu, RdYlBu_r, Greys
						site_model_gis_file = os.path.join(gis_folder, "Polygons_v3.shp")
						#plot_rupture_probabilities(source_model, prob_dict, pe_site_models, ne_site_models,
						#							map_region, plot_point_ruptures=True, colormap="RdYlBu_r",
						#							title=title, text_box=text_box, site_model_gis_file=site_model_gis_file,
						#							fig_filespec=fig_filespec)

				## Generate animated GIF
				"""
				img_basename = "%s_MMI=%s_%s" % (scenario, intensity, ipe_name)
				try:
					create_animated_gif(fig_folder, img_basename)
				except:
					pass
				"""


		## Plot max_prob vs magnitude
		if section_probs:
			colors = ['r', 'b', 'g', 'm', 'k']
			for ipe_name, color in zip(ipe_names, colors):
				#pylab.plot(Mrange, section_probs, 'kx-', label='Scenario section')
				#pylab.plot(Mrange, max_probs, 'kx--', label='Max. probability')
				if "WithSigma" in ipe_name:
					label = ipe_name[:ipe_name.find("WithSigma")]
				elif ipe_name == "LogicTree":
					label = "IPE Mix"
				else:
					label = ipe_name
				if ipe_name in max_probs:
					pylab.plot(Mrange, max_probs[ipe_name], 'x-', color=color, label=label)
					pylab.plot(Mrange, section_probs[ipe_name], 'x--', color=color, label='_nolegend_')
					#pylab.hlines(scenario_probs[ipe_name], Mrange[0], Mrange[-1], colors=color, linestyles="dotted")
			pylab.vlines(scenflt.mfd.min_mag, 0, 1, linestyles="dotted")
			pylab.xlim(Mrange[0], Mrange[-1])
			pylab.ylim(0, 1)
			pylab.xlabel("Magnitude")
			pylab.ylabel("Max. normalized probability")
			if len(threshold_set) == 1:
				intensity = threshold_set[0]
			else:
				intensity = "mix"
			intensity_label = roman_intensity = ', '.join([get_roman_intensity(mmi) for mmi in threshold_set])
			pylab.title("%s scenario, MMI threshold = %s" % (scenario, intensity_label))
			pylab.legend(loc=3)
			fig_filename = "%s_MMI=%s_M_vs_prob.%s" % (scenario, intensity, output_format)
			fig_filespec = os.path.join(fig_folder, fig_filename)
			#fig_filespec = None
			if fig_filespec:
				pylab.savefig(fig_filespec, dpi=200)
			else:
				pylab.show()

			pylab.clf()
