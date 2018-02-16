
import os
import numpy as np
import openquake.hazardlib as oqhazlib
import mapping.layeredbasemap as lbm
import hazard.rshalib as rshalib
from hazard.rshalib.source import SimpleUniformGridSourceModel
from hazard.rshalib.source_estimation import estimate_epicenter_location_and_magnitude_from_intensities

from aysenlib import (create_uniform_grid_source_model, plot_gridsearch_map,
						read_evidence_site_info_from_gis,
						TRT, USD, LSD, RAR, RMS, project_folder, gis_folder)


fig_folder = os.path.join(project_folder, "Figures", "Events", "bw1997", "NE=+1.0")


## Event names
events = ['2007', 'SL-A', 'SL-B', 'SL-C', 'SL-CD', 'SL-D', 'SL-DE', 'SL-EF', 'SL-F', 'SL-G']
#events = ["SL-C"]
#events = events[:1]

## IPE and IMT
ipe_name = "BakunWentworth1997"
#ipe_name = "AtkinsonWald2007"
imt = oqhazlib.imt.MMI()

polygon_discretization = 2.5


## Construct grid source model
grid_outline = (-74, -72, -46.25, -44.75)
grid_spacing = 0.1

min_mag, max_mag, mag_bin_width = 4.5, 8.5, 0.2
depth = 10
strike, dip, rake = 20, 80, 180
point_msr = oqhazlib.scalerel.PointMSR()
wc1994_msr = oqhazlib.scalerel.WC1994()

grd_src_model = SimpleUniformGridSourceModel(grid_outline, grid_spacing,
			min_mag, min_mag + mag_bin_width, mag_bin_width, depth,
			strike, dip, rake, wc1994_msr, USD, LSD, RMS, RAR, TRT)
lon_grid, lat_grid = grd_src_model.lon_grid, grd_src_model.lat_grid


## Loop over events
for event in events:
	## Read MTD evidence
	pe_site_models, ne_site_models = [], []
	pe_thresholds, pe_sites, ne_thresholds, ne_sites = [], [], [], []
	for geom_type in ["Polygons", "Points"]:
		shapefile = os.path.join(gis_folder, "%s.shp" % geom_type)
		(_pe_thresholds, _pe_site_models,
		_ne_thresholds, _ne_site_models) = read_evidence_site_info_from_gis(shapefile,
													event, polygon_discretization)
		pe_site_models.extend(_pe_site_models)
		ne_site_models.extend(_ne_site_models)
		for pesm, pe_threshold in zip(_pe_site_models, _pe_thresholds):
			for pe_site in pesm.get_sites():
				pe_sites.append(pe_site)
				pe_thresholds.append(pe_threshold)
				print("+%s: %s" % (pe_site.name, pe_threshold))
		for nesm, ne_threshold in zip(_ne_site_models, _ne_thresholds):
			for ne_site in nesm.get_sites():
				ne_sites.append(ne_site)
				ne_thresholds.append(ne_threshold)
				print("-%s: %s" % (ne_site.name, ne_threshold))
		print sum([len(pesm) for pesm in _pe_site_models]), sum([len(nesm) for nesm in _ne_site_models])
	print len(pe_sites), len(ne_sites)
	pe_thresholds = np.array(pe_thresholds)
	ne_thresholds = np.array(ne_thresholds)

	#ne_thresholds -= 0.5

	## Compute magnitudes and RMS errors at grid points
	(mag_grid, rms_grid) = (
		estimate_epicenter_location_and_magnitude_from_intensities(
		ipe_name, imt, grd_src_model, pe_sites, pe_thresholds,
		ne_sites, ne_thresholds, method="reverse", mag_bounds=(min_mag, max_mag)))
	idx = np.unravel_index(rms_grid.argmin(), rms_grid.shape)
	print mag_grid[idx], lon_grid[idx], lat_grid[idx]

	#rms_grid[np.isnan(mag_grid)] = np.inf


	## Plot map
	# TODO: blend alpha in function of rms
	# See: https://matplotlib.org/devdocs/gallery/images_contours_and_fields/image_transparency_blend.html
	site_model_gis_file = os.path.join(gis_folder, "Polygons.shp")
	map = plot_gridsearch_map(grd_src_model, mag_grid, rms_grid,
							pe_site_models, ne_site_models,
							site_model_gis_file=site_model_gis_file,
							plot_rms_as_alpha=False, plot_epicenter_as="both")

	fig_filespec = os.path.join(fig_folder, "%s_bw1997.png" % event)
	#fig_filespec = None

	dpi = 200 if fig_filespec else 90
	map.plot(fig_filespec=fig_filespec, dpi=dpi)
