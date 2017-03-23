import os
import numpy as np
import hazard.rshalib as rshalib
from rupture_probabilities import *


# TODO: check if negative evidence can be more strict (higher or lower probabilities?): event SL-C
# TODO: forward modelling to check sensitivity


project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
#project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")


## Event names
events = ['2007', 'SL-A', 'SL-B', 'SL-C', 'SL-CD', 'SL-D', 'SL-DE', 'SL-EF', 'SL-F', 'SL-G']


## Magnitude range to test
dM = 0.25
min_mag, max_mag = 6.0 - dM/2, 7.0
Mrange = np.arange(min_mag, max_mag, dM) + dM/2
print Mrange


## Construct ground-motion model
ipe_name = "BakunWentworth1997WithSigma"
#ipe_name = "AtkinsonWald2007"
#ipe_name = "AllenEtAl2012Rrup"
trt_gsim_dict = {TRT: ipe_name}
ground_motion_model = rshalib.gsim.GroundMotionModel(ipe_name, trt_gsim_dict)
truncation_level = 2.5
if "AllenEtAl2012" in ipe_name:
	truncation_level /= 2


## Parameters
soil_params = rshalib.site.REF_SOIL_PARAMS
polygon_discretization = 2.5
imt = oqhazlib.imt.MMI()
strict_intersection = True
max_prob_color = 1.0
map_region = (-74, -72, -46, -44.5)


for M in Mrange:
	## Read fault source model
	#fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel.shp")
	fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel2.TAB")
	#fault_filespec = os.path.join(gis_folder, "EnergiAustral_faults.TAB")
	source_model = read_fault_source_model_as_floating_ruptures(fault_filespec, M, M + dM/2, dM, depth=0.1)

	for event in events:
		fig_folder = os.path.join(project_folder, "Figures", event)
		if not os.path.exists(fig_folder):
			os.mkdir(fig_folder)

		## Read lake evidence
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
		if ipe_name == "AllenEtAl2012Rrup":
			intensity_correction = -0.5
		elif ipe_name == "AtkinsonWald2007":
			intensity_correction = -1
		pe_thresholds = np.array(pe_thresholds) + intensity_correction
		ne_thresholds = np.array(ne_thresholds) + intensity_correction

		print pe_thresholds
		print ne_thresholds


		## Compute rupture probabilities
		prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
							source_model, ground_motion_model, imt, pe_site_models,
							pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
							strict_intersection=strict_intersection)
		#print prob_dict
		probs = np.array(prob_dict.values())
		max_prob = probs.max()


		## Plot
		text_box = "Event: %s\nM: %.2f\nPmax: %.2f" % (event, M, max_prob)

		#itle = "Event: %s, IPE: %s, M=%.2f" % (event, ipe_name, M)
		title = ""
		fig_filename = "%s_%s_M=%.2f.PNG" % (event, ipe_name, M)

		fig_filespec = os.path.join(fig_folder, fig_filename)
		#fig_filespec = None

		plot_rupture_probabilities(source_model, prob_dict, pe_site_models, ne_site_models,
									map_region, max_prob_color, plot_point_ruptures=True,
									title=title, text_box=text_box, fig_filespec=fig_filespec)
