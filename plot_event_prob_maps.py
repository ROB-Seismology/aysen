import os
import numpy as np
import hazard.rshalib as rshalib
from rupture_probabilities import *


# TODO: check if negative evidence can be more strict (higher or lower probabilities?): event SL-C


project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
#project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")


## Event names
events = ['2007', 'SL-A', 'SL-B', 'SL-C', 'SL-CD', 'SL-D', 'SL-DE', 'SL-EF', 'SL-F', 'SL-G']
#events = ["SL-C"]
events = events[:1]


## IPE names
#ipe_names = ["LogicTree", "AllenEtAl2012", "AtkinsonWald2007", "BakunWentworth1997WithSigma", "Barrientos1980WithSigma"]
#ipe_names = ipe_names[3:4]
ipe_names = ["BakunWentworth1997WithSigma"]


## Parameters
#truncation_level = 2.5
truncation_level = 1
#if "AllenEtAl2012" in ipe_name:
#	truncation_level /= 2
soil_params = rshalib.site.REF_SOIL_PARAMS
polygon_discretization = 2.5
imt = oqhazlib.imt.MMI()
strict_intersection = True


## Map parameters
#map_region = (-74, -72, -46, -44.5)
map_region = (-74, -72, -46.25, -44.75)
output_format = "png"


## Read fault source model

## Faults discretized as floating ruptures
"""
dM = 0.25
min_mag, max_mag = 6.0 - dM/2, 7.0
Mrange = np.arange(min_mag, max_mag, dM) + dM/2
print Mrange
for M in Mrange:
	#fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel.shp")
	fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel2.TAB")
	source_model = read_fault_source_model_as_floating_ruptures(fault_filespec, M, M + dM/2, dM, depth=0.1)
"""

## Fault network
dM = 0.2
Mrange = []
max_prob_dict = {}
section_prob_dict = {}
fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel3.TAB")
for M, source_model in read_fault_source_model_as_network(fault_filespec, dM=dM):
	Mrange.append(M)

	## Loop over events
	for event in events:
		if M == Mrange[0]:
			max_prob_dict[event] = {}
			section_prob_dict[event] = {}
		fig_folder = os.path.join(project_folder, "Figures", event)
		if not os.path.exists(fig_folder):
			os.mkdir(fig_folder)

		## Read MTD evidence
		pe_thresholds, pe_site_models, ne_thresholds, ne_site_models = [], [], [], []
		for geom_type in ["Polygons", "Points"]:
			shapefile = os.path.join(gis_folder, "%s.shp" % geom_type)
			(_pe_thresholds, _pe_site_models,
			_ne_thresholds, _ne_site_models) = read_evidence_site_info_from_gis(shapefile, event, polygon_discretization)
			pe_thresholds.extend(_pe_thresholds)
			pe_site_models.extend(_pe_site_models)
			ne_thresholds.extend(_ne_thresholds)
			ne_site_models.extend(_ne_site_models)

		intensity_correction = 0
		#if ipe_name == "AllenEtAl2012Rrup":
		#	intensity_correction = -0.5
		#elif ipe_name == "AtkinsonWald2007":
		#	intensity_correction = -1
		pe_thresholds = np.array(pe_thresholds) + intensity_correction
		ne_thresholds = np.array(ne_thresholds) + intensity_correction

		print pe_thresholds
		#ne_thresholds = np.minimum(7.5, ne_thresholds)
		print ne_thresholds

		## Construct ground-motion model
		for ipe_name in ipe_names:
			if M == Mrange[0]:
				max_prob_dict[event][ipe_name] = []
				section_prob_dict[event][ipe_name] = {}
			if ipe_name != "LogicTree":
				trt_gsim_dict = {TRT: ipe_name}
				gmpe_system_def = {TRT: rshalib.pmf.GMPEPMF([ipe_name], [1])}
				integration_distance_dict = {}
			else:
				models = ["Barrientos1980WithSigma", "BakunWentworth1997WithSigma", "AllenEtAl2012", "AtkinsonWald2007"]
				weights = [0.32, 0.26, 0.25, 0.17]
				gmpe_system_def = {TRT: rshalib.pmf.GMPEPMF(models, weights)}
				integration_distance_dict = {"AtkinsonWald2007": (None, 30)}


			## Compute rupture probabilities
			prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
								source_model, gmpe_system_def, imt, pe_site_models,
								pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
								integration_distance_dict=integration_distance_dict,
								strict_intersection=strict_intersection)
			#print prob_dict
			probs = np.array(prob_dict.values())
			probs = probs[:, 0]
			max_prob = probs.max()
			max_prob_dict[event][ipe_name].append(max_prob)
			print M, max_prob
			for rup_name, prob in zip(prob_dict.keys(), probs):
				for section in rup_name.split('+'):
					if not section in section_prob_dict[event][ipe_name]:
						section_prob_dict[event][ipe_name][section] = [prob]
					else:
						section_prob_dict[event][ipe_name][section].append(prob)

			## Plot
			if "WithSigma" in ipe_name:
				ipe_label = ipe_name[:ipe_name.find("WithSigma")]
			else:
				ipe_label = ipe_name

			text_box = "Event: %s\nIPE: %s\nM: %.2f, Pmax: %.2f"
			text_box %= (event, ipe_label, M, max_prob)

			#title = "Event: %s, IPE: %s, M=%.2f" % (event, ipe_name, M)
			title = ""

			fig_filename = "%s_%s_M=%.2f.%s" % (event, ipe_label, M, output_format)
			fig_filespec = os.path.join(fig_folder, fig_filename)
			#fig_filespec = None

			## Colormaps: RdBu_r, YlOrRd, BuPu, RdYlBu_r, Greys
			site_model_gis_file = os.path.join(gis_folder, "Polygons.shp")
			plot_rupture_probabilities(source_model, prob_dict, pe_site_models, ne_site_models,
										map_region, plot_point_ruptures=True, colormap="RdYlBu_r",
										title=title, text_box=text_box, site_model_gis_file=site_model_gis_file,
										fig_filespec=fig_filespec)

## Determine which sections have highest probability
for event in events:
	print event
	for ipe_name in ipe_names:
		print ipe_name
		sections = section_prob_dict[event][ipe_name].keys()
		probs = [np.array(l) for l in section_prob_dict[event][ipe_name].values()]
		#print probs[0]
		mean_probs = np.array([p.mean() for p in probs])
		max_probs = np.array([p.max() for p in probs])
		idxs = np.argsort(mean_probs)[::-1]
		for idx in idxs[:10]:
			print("  %s: %.2f, %.2f" % (sections[idx], mean_probs[idx], max_probs[idx]))


## Plot max_prob vs magnitude for different IPEs
colors = ['r', 'b', 'g', 'm', 'k']
for event in events:
	for ipe_name, color in zip(ipe_names, colors):
		if "WithSigma" in ipe_name:
			label = ipe_name[:ipe_name.find("WithSigma")]
		else:
			label = ipe_name
		pylab.plot(Mrange, max_prob_dict[event][ipe_name], 'x-', color=color, label=label)
	pylab.xlim(Mrange[0], Mrange[-1])
	pylab.ylim(0, 1)
	pylab.xlabel("Magnitude")
	pylab.ylabel("Probability")
	pylab.title("Event: %s" % event)
	pylab.legend(loc=3)

	fig_folder = os.path.join(project_folder, "Figures", event)
	fig_filename = "%s_M_vs_prob.%s" % (event, output_format)
	fig_filespec = os.path.join(fig_folder, fig_filename)
	#fig_filespec = None
	if fig_filespec:
		pylab.savefig(fig_filespec, dpi=200)
	else:
		pylab.show()
