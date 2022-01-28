import os, sys
import numpy as np

import hazard.rshalib as rshalib

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SCRIPT_DIR)
from aysenlib import *
from create_animated_gif import create_animated_gif



from aysenlib import (project_folder, gis_folder, data_points, data_polygons, LOFZ_model)
base_fig_folder = os.path.join(project_folder, "Figures", "Events", "v2", "NE=+0.5")

## Event names
events = ['2007', 'SL-G', 'SL-F', 'SL-EF', 'SL-DE', 'SL-D', 'SL-CD', 'SL-C', 'SL-B', 'SL-A']
#events = ['SL-A']
#events = ['SL-EF', 'SL-DE', 'SL-A']
#events = ["SL-A"]
#events = events[:6]
#events = events[6:]

## Selected magnitudes for final figure in paper
events_mags = {'2007': 6.28,
				'SL-G': 5.72,
				'SL-F': 6.48,
				'SL-EF': 6.48,
				'SL-DE': 6.28,
				'SL-D': 5.13,
				'SL-CD': 6.28,
				'SL-C': 6.48,
				'SL-B': 5.13,
				'SL-A': 6.28}


## IPE names
#ipe_names = ["LogicTree", "AllenEtAl2012", "AtkinsonWald2007", "BakunWentworth1997WithSigma", "Barrientos1980WithSigma"]
#ipe_names = ["AtkinsonWald2007"]
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
map_region = (-74, -72, -46.25, -44.8)
output_format = "pdf"


## Read fault source model

## Faults discretized as floating ruptures
"""
dM = 0.25
min_mag, max_mag = 6.0 - dM/2, 7.0
fault_mags = np.arange(min_mag, max_mag, dM) + dM/2
print(fault_mags)
for M in fault_mags:
	#fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel.shp")
	fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel2.TAB")
	source_model = read_fault_source_model_as_floating_ruptures(fault_filespec, M, M + dM/2, dM, depth=0.1)
"""

## Discretize faults as network
fault_filespec = os.path.join(gis_folder, LOFZ_model)
dM = 0.2
fault_mags, fault_networks = [], []
for M, source_model in read_fault_source_model_as_network(fault_filespec, dM=dM):
	fault_mags.append(M)
	fault_networks.append(source_model)


## Loop over events
max_prob_dict = {}
section_prob_dict = {}
for event in events:
	print(event)
	max_prob_dict[event] = {}
	section_prob_dict[event] = {}

	fig_folder = os.path.join(base_fig_folder, event)
	if not os.path.exists(fig_folder):
		os.mkdir(fig_folder)

	## Read MTD evidence
	pe_thresholds, pe_site_models, ne_thresholds, ne_site_models = [], [], [], []
	for geom_type in [data_polygons, data_points]:
		shapefile = os.path.join(gis_folder, "%s.shp" % geom_type)
		(_pe_thresholds, _pe_site_models,
		_ne_thresholds, _ne_site_models) = read_evidence_site_info_from_gis(shapefile, event, polygon_discretization)
		## Remove "West" polygon (only used for 2007 event)
		for _pe_threshold, _pe_site_model in zip(_pe_thresholds, _pe_site_models):
			pe_thresholds.append(_pe_threshold)
			pe_site_models.append(_pe_site_model)
		for _ne_threshold, _ne_site_model in zip(_ne_thresholds, _ne_site_models):
			ne_thresholds.append(_ne_threshold)
			ne_site_models.append(_ne_site_model)
		#pe_thresholds.extend(_pe_thresholds)
		#pe_site_models.extend(_pe_site_models)
		#ne_thresholds.extend(_ne_thresholds)
		#ne_site_models.extend(_ne_site_models)

	intensity_correction = 0
	#if ipe_name == "AllenEtAl2012Rrup":
	#	intensity_correction = -0.5
	#elif ipe_name == "AtkinsonWald2007":
	#	intensity_correction = -1
	pe_thresholds = np.array(pe_thresholds) + intensity_correction
	ne_thresholds = np.array(ne_thresholds) + intensity_correction

	ne_thresholds -= 0.5

	for pe_site_model, pe_threshold in zip(pe_site_models, pe_thresholds):
		print("+%s (n=%d): %s" % (pe_site_model.name.encode(errors='replace'), len(pe_site_model), pe_threshold))
	for ne_site_model, ne_threshold in zip(ne_site_models, ne_thresholds):
		print("-%s (n=%d): %s" % (ne_site_model.name.encode(errors='replace'), len(ne_site_model), ne_threshold))


	## Construct ground-motion model
	for ipe_name in ipe_names:
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

	for M, source_model in zip(fault_mags, fault_networks):
		## Compute rupture probabilities
		prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
							source_model, gmpe_system_def, imt, pe_site_models,
							pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
							integration_distance_dict=integration_distance_dict,
							strict_intersection=strict_intersection)
		#print(prob_dict)
		probs = np.array(list(prob_dict.values()))
		probs = probs[:, 0]
		max_prob = probs.max()
		max_prob_dict[event][ipe_name].append(max_prob)
		print(M, max_prob)
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

		#text_box = "Event: %s\nIPE: %s\nM: %.2f, Pmax: %.2f"
		text_box = "Event: %s\nM: %.2f\nPmax: %.2f"
		text_box %= (event, M, max_prob)

		#title = "Event: %s, IPE: %s, M=%.2f" % (event, ipe_name, M)
		title = ""

		fig_filename = "%s_%s_M=%.2f.%s" % (event, ipe_label, M, output_format)
		fig_filespec = os.path.join(fig_folder, fig_filename)
		fig_filespec = None

		## Colormaps: RdBu_r, YlOrRd, BuPu, RdYlBu_r, Greys
		site_model_gis_file = os.path.join(gis_folder, "Polygons_v3.shp")
		if np.isclose(M, events_mags[event], atol=0.01):
			plot_rupture_probabilities(source_model, prob_dict, pe_site_models, ne_site_models,
										map_region, plot_point_ruptures=True, colormap="RdYlBu_r",
										title=title, text_box=text_box, site_model_gis_file=site_model_gis_file,
										fig_filespec=fig_filespec)


	## Generate animated GIF
	for ipe_name in ipe_names:
		img_basename = "%s_%s" % (event, ipe_label)
		#create_animated_gif(fig_folder, img_basename)
#exit()


## Determine which sections have highest probability
for event in events:
	print(event)
	for ipe_name in ipe_names:
		print(ipe_name)
		sections = list(section_prob_dict[event][ipe_name].keys())
		probs = [np.array(list(l)) for l in section_prob_dict[event][ipe_name].values()]
		#print(probs[0])
		mean_probs = np.array([p.mean() for p in probs])
		max_probs = np.array([p.max() for p in probs])
		idxs = np.argsort(mean_probs)[::-1]
		for idx in idxs[:10]:
			print("  %s: %.2f, %.2f" % (sections[idx], mean_probs[idx], max_probs[idx]))



## Plot max_prob vs magnitude for different IPEs per event
colors = ['r', 'b', 'g', 'm', 'k']
for event in events:
	pylab.cla()
	for ipe_name, color in zip(ipe_names, colors):
		if "WithSigma" in ipe_name:
			label = ipe_name[:ipe_name.find("WithSigma")]
		else:
			label = ipe_name
		pylab.plot(fault_mags, max_prob_dict[event][ipe_name], 'x-', color=color, label=label)
	pylab.xlim(fault_mags[0], fault_mags[-1])
	pylab.ylim(0, 1)
	pylab.xlabel("Magnitude")
	pylab.ylabel("Max. normalized probability")
	pylab.title("Event: %s" % event)
	pylab.legend(loc=3)

	fig_folder = os.path.join(base_fig_folder, event)
	fig_filename = "%s_M_vs_prob.%s" % (event, output_format)
	#fig_filespec = os.path.join(fig_folder, fig_filename)
	fig_filespec = None
	if fig_filespec:
		pylab.savefig(fig_filespec, dpi=200)
	else:
		pylab.show()


## Plot max_prob vs magnitude for different events in one plot

colors = ['r', 'b', 'g', 'm', 'c', 'k']
pylab.cla()
for event, color in zip(events, colors):
	label = event
	ipe_name = ipe_names[0]
	if "WithSigma" in ipe_name:
		ipe_label = ipe_name[:ipe_name.find("WithSigma")]
	else:
		ipe_label = ipe_name
	pylab.plot(fault_mags, max_prob_dict[event][ipe_name], 'x-', color=color, label=label)
	pylab.xlim(fault_mags[0], fault_mags[-1])
	pylab.ylim(0, 1)
	pylab.xlabel("Magnitude")
	pylab.ylabel("Max. normalized probability")
	pylab.legend(loc=4)

	fig_folder = os.path.join(base_fig_folder)
	#fig_filename = "events1-6_M_vs_prob_%s.%s" % (ipe_label, output_format)
	fig_filename = "events7-10_M_vs_prob_%s.%s" % (ipe_label, output_format)
	#fig_filespec = os.path.join(fig_folder, fig_filename)
	fig_filespec = None
	if fig_filespec:
		pylab.savefig(fig_filespec, dpi=200)
	else:
		pylab.show()
