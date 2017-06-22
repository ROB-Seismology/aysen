import os
import numpy as np
import hazard.rshalib as rshalib
from rupture_probabilities import *


project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
#project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")
fig_folder = os.path.join(project_folder, "Figures", "Sensitivity")


## Scenarios
#scenario = "Quitralco"
scenario = "Azul Tigre South"
#scenario = "Due East"


## Thresholds (depends on type of evidence)
## NE: 7.5, 8.5
## PE: 5.5, 6.0, 6.5 7.5
threshold = 8
## Note: increasing delta_threshold increases computed probabilities,
## but decreases discriminating power!
## Note2: it is not possible to test this, as a particular site may
## be both positive and negative if delta_threshold > 0 !!
#delta_threshold = 0.5
delta_threshold = 0.
pe_threshold = threshold - delta_threshold
ne_threshold = threshold + delta_threshold


## Construct ground-motion model
#ipe_name = "BakunWentworth1997WithSigma"
#ipe_name = "AtkinsonWald2007"
#ipe_name = "AllenEtAl2012"
ipe_name = "logictree"

if ipe_name != "logictree":
	trt_gsim_dict = {TRT: ipe_name}
	ground_motion_model = rshalib.gsim.GroundMotionModel(ipe_name, trt_gsim_dict)
	gmpe_system_def = {TRT: rshalib.pmf.GMPEPMF([ipe_name], [1])}
	integration_distance_dict = {}
else:
	gmpe_system_def = {TRT: rshalib.pmf.GMPEPMF(["BakunWentworth1997WithSigma", "AllenEtAl2012", "AtkinsonWald2007"], [0.6, 0.3, 0.1])}
	integration_distance_dict = {"AtkinsonWald2007": (None, 30)}
	#TODO: apply distance filtering for BakunWentworth1997WithSigma when M larger than
	#threshold causing too high values in near field (to be determined)
	#integration_distance_dict = {"AtkinsonWald2007": (None, 30), "BakunWentworth1997WithSigma": (20, None)}
truncation_level = 2.5


## Parameters
soil_params = rshalib.site.REF_SOIL_PARAMS
polygon_discretization = 2.5
imt = oqhazlib.imt.MMI()
strict_intersection = True


## Read site info
all_site_models = []
for geom_type in ["Polygons", "Points"]:
	shapefile = os.path.join(gis_folder, "%s.shp" % geom_type)
	site_models = read_evidence_sites_from_gis(shapefile, polygon_discretization)
	all_site_models.extend(site_models)
print len(all_site_models)


## Read rupture scenario
scenario_filespec = os.path.join(gis_folder, 'LOFZ_rupture_scenarios.TAB')
flt_src_model = read_fault_source_model(scenario_filespec)
[flt] = [flt for flt in flt_src_model if flt.name == scenario]
aspect_ratio = flt.get_length() / flt.width

## Set characteristic magnitude from MSR and fault area
import openquake.hazardlib as oqhazlib
msr = getattr(oqhazlib.scalerel, MSR)()
char_mag = msr.get_median_mag(flt.get_rupture().surface.get_area(), flt.rake)
print flt.get_length(), flt.width, char_mag
flt.mfd.min_mag = np.round(char_mag, 1)
flt.mfd.char_mag = flt.mfd.min_mag + flt.mfd.bin_width/2.

flt_src_model.sources = [flt]


## Run DSHA model and determine pe_site_models and ne_site_models
pe_thresholds, pe_site_models, ne_thresholds, ne_site_models = [], [], [], []
for site_model in all_site_models:
	dsha_model = rshalib.shamodel.DSHAModel(scenario, flt_src_model, gmpe_system_def,
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

print len(pe_site_models), len(ne_site_models)

## Average probability
avg_prob = 0.5 ** (len(pe_site_models) + len(ne_site_models))
print avg_prob


## Compute probability for this scenario
scenario_src_model = rshalib.source.SourceModel(scenario, [flt])
prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
					scenario_src_model, gmpe_system_def, imt, pe_site_models,
					pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
					integration_distance_dict=integration_distance_dict,
					strict_intersection=strict_intersection)
[scenario_prob] = prob_dict.values()[0]
print scenario_prob


## Magnitude range to test
dM = 0.5
min_mag, max_mag = 6.0 - dM/2, 7.0
#Mrange = np.arange(min_mag, max_mag, dM) + dM/2
Mrange = np.linspace(flt.mfd.min_mag - dM, flt.mfd.min_mag + dM, 3)
#Mrange = [flt.mfd.min_mag]
print Mrange


## Map parameters
#map_region = (-74, -72, -46, -44.5)
map_region = (-74, -72, -46.25, -44.75)
output_format = "png"


## Test
M_prob_dict = {}
for M in Mrange:
	## Read fault source model
	fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel2.TAB")
	source_model = read_fault_source_model_as_floating_ruptures(fault_filespec,
						M - dM/2, M, dM, depth=0.1, aspect_ratio=aspect_ratio)

	#for flt in source_model:
	#	if flt.source_id[0] == '4':
	#		print flt.get_length(), flt.get_width()

	## Compute rupture probabilities
	prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
						source_model, gmpe_system_def, imt, pe_site_models,
						pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
						integration_distance_dict=integration_distance_dict,
						strict_intersection=strict_intersection)

	#for rup_name in prob_dict:
	#	if rup_name[0] == '3':
	#		prob_dict[rup_name] *= 0
	probs = np.array(prob_dict.values())
	max_prob_idx = probs.argmax()
	max_prob = probs[max_prob_idx]
	rup_name = prob_dict.keys()[max_prob_idx]
	print M, rup_name, max_prob

	for rup_name, prob in prob_dict.items():
		if prob >= scenario_prob:
			print M, rup_name, prob


	## Plot
	text_box = "Scenario: %s\nM: %.2f\nPmax: %.2f" % (scenario, M, max_prob)

	title = ""
	fig_filename = "%s_%s_M=%.2f.%s" % (scenario, ipe_name, M, output_format)
	fig_filespec = os.path.join(fig_folder, fig_filename)
	#fig_filespec = None

	## Colormaps: RdBu_r, YlOrRd, BuPu, RdYlBu_r, Greys
	plot_rupture_probabilities(source_model, prob_dict, pe_site_models, ne_site_models,
								map_region, plot_point_ruptures=True, colormap="RdYlBu_r",
								title=title, text_box=text_box, fig_filespec=fig_filespec)
