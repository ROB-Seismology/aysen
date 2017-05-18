import os
import numpy as np
import hazard.rshalib as rshalib
from rupture_probabilities import *



project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
#project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")


## Scenarios
scenario = "Due East"


## Thresholds (depends on type of evidence)
threshold = 8
delta_threshold = 0.5
#delta_threshold = 0.
pe_threshold = threshold - delta_threshold
ne_threshold = threshold + delta_threshold


## Construct ground-motion model
ipe_name = "BakunWentworth1997WithSigma"
#ipe_name = "AtkinsonWald2007"
#ipe_name = "AllenEtAl2012Rrup"
trt_gsim_dict = {TRT: ipe_name}
ground_motion_model = rshalib.gsim.GroundMotionModel(ipe_name, trt_gsim_dict)
gmpe_system_def = {TRT: rshalib.pmf.GMPEPMF([ipe_name], [1])}
integration_distance = 300
truncation_level = 2.5
#truncation_level = 0.
if "AllenEtAl2012" in ipe_name:
	truncation_level /= 2


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

## Set characteristic magnitude from MSR and fault area
import openquake.hazardlib as oqhazlib
msr = getattr(oqhazlib.scalerel, MSR)()
char_mag = msr.get_median_mag(flt.get_rupture().surface.get_area(), flt.rake)
flt.mfd.min_mag = np.round(char_mag, 1)
flt.mfd.char_mag = flt.mfd.min_mag + flt.mfd.bin_width/2.

flt_src_model.sources = [flt]


## Run DSHA model and determine pe_site_models and ne_site_models
pe_thresholds, pe_site_models, ne_thresholds, ne_site_models = [], [], [], []
for site_model in all_site_models:
	dsha_model = rshalib.shamodel.DSHAModel(scenario, flt_src_model, gmpe_system_def,
				grid_outline=[], grid_spacing=None,
				soil_site_model=site_model, imt_periods={"MMI": [0]},
				truncation_level=0, integration_distance=integration_distance)

	uhs_field = dsha_model.calc_gmf_fixed_epsilon()
	hm = uhs_field.getHazardMap()
	if (hm.intensities > pe_threshold).all():
		pe_site_models.append(site_model)
		pe_thresholds.append(pe_threshold)
	elif (hm.intensities <= ne_threshold).all():
		ne_site_models.append(site_model)
		ne_thresholds.append(ne_threshold)

print len(pe_site_models), len(ne_site_models)


## Magnitude range to test
dM = 0.5
min_mag, max_mag = 6.0 - dM/2, 7.0
#Mrange = np.arange(min_mag, max_mag, dM) + dM/2
Mrange = np.linspace(flt.mfd.min_mag - dM, flt.mfd.min_mag + dM, 3)
print Mrange


## Map parameters
max_prob_color = 1.0
map_region = (-74, -72, -46, -44.5)
output_format = "PDF"


## Test
M_prob_dict = {}
for M in Mrange:
	## Read fault source model
	fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel2.TAB")
	source_model = read_fault_source_model_as_floating_ruptures(fault_filespec, M - dM/2, M, dM, depth=0.1)

	## Compute rupture probabilities
	prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
						source_model, ground_motion_model, imt, pe_site_models,
						pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
						strict_intersection=strict_intersection)
	M_prob_dict[M] = prob_dict


#for M in sorted(M_prob_dict.keys()):
	prob_dict = M_prob_dict[M]
	probs = np.array(prob_dict.values())
	max_prob_idx = probs.argmax()
	max_prob = probs[max_prob_idx]
	rup_name = prob_dict.keys()[max_prob_idx]
	print M, rup_name, max_prob


	## Plot
	text_box = "Scenario: %s\nM: %.2f\nPmax: %.2f" % (scenario, M, max_prob)

	title = ""
	fig_filename = "%s_%s_M=%.2f.%s" % (scenario, ipe_name, M, output_format)
	#fig_filespec = os.path.join(fig_folder, fig_filename)
	fig_filespec = None

	plot_rupture_probabilities(source_model, prob_dict, pe_site_models, ne_site_models,
								map_region, max_prob_color, plot_point_ruptures=True,
								title=title, text_box=text_box, fig_filespec=fig_filespec)
