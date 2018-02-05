
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



## Event names
events = ['2007', 'SL-A', 'SL-B', 'SL-C', 'SL-CD', 'SL-D', 'SL-DE', 'SL-EF', 'SL-F', 'SL-G']
events = ["SL-C"]
#events = events[:1]

## IPE and IMT
ipe_name = "BakunWentworth1997"
imt = oqhazlib.imt.MMI()

polygon_discretization = 2.5


## Construct grid source model
grid_outline = (-74, -72, -46.25, -44.75)
grid_spacing = 0.1

min_mag, max_mag, mag_bin_width = 4.5, 7.5, 0.2
depth = 10
strike, dip, rake = 90, 60, 0
point_msr = oqhazlib.scalerel.PointMSR()

grd_src_model = SimpleUniformGridSourceModel(grid_outline, grid_spacing,
			min_mag, min_mag + mag_bin_width, mag_bin_width, depth,
			strike, dip, rake, point_msr, USD, LSD, RMS, RAR, TRT)
lon_grid, lat_grid = grd_src_model.lon_grid, grd_src_model.lat_grid


## Loop over events
for event in events:
	## Read MTD evidence
	pe_thresholds, pe_sites, ne_thresholds, ne_sites = [], [], [], []
	for geom_type in ["Polygons", "Points"]:
		shapefile = os.path.join(gis_folder, "%s.shp" % geom_type)
		(_pe_thresholds, _pe_site_models,
		_ne_thresholds, _ne_site_models) = read_evidence_site_info_from_gis(shapefile, event, polygon_discretization)
		for pesm, pe_threshold in zip(_pe_site_models, _pe_thresholds):
			pe_sites.extend(pesm.get_sites())
			pe_thresholds.extend([pe_threshold] * len(pesm))
		pe_site_model = rshalib.site.SoilSiteModel("Positive evidence", pe_sites)
		for nesm, ne_threshold in zip(_ne_site_models, _ne_thresholds):
			ne_sites.extend(nesm.get_sites())
			ne_thresholds.extend([ne_threshold] * len(nesm))
		ne_site_model = rshalib.site.SoilSiteModel("Negative evidence", ne_sites)
		print len(pe_sites), len(ne_sites)
		print len(pe_thresholds), len(ne_thresholds)

	## Compute magnitudes and RMS errors at grid points
	(mag_grid, rms_grid) = (
		estimate_epicenter_location_and_magnitude_from_intensities(
		ipe_name, imt, grd_src_model, pe_sites, pe_thresholds,
		ne_sites, ne_thresholds, method="reverse", mag_bounds=(min_mag, max_mag)))
	idx = np.unravel_index(rms_grid.argmin(), rms_grid.shape)
	print mag_grid[idx], lon_grid[idx], lat_grid[idx]


	## Plot map
	map = plot_gridsearch_map(grd_src_model, mag_grid, rms_grid,
							[pe_site_model], [ne_site_model])

	point_data = lbm.PointData(lon_grid[idx], lat_grid[idx])
	point_style = lbm.PointStyle(shape='*', fill_color='yellow', size=14)
	layer = lbm.MapLayer(point_data, point_style)
	map.layers.append(layer)

	map.plot()
