
import os, sys
import numpy as np
import openquake.hazardlib as oqhazlib
from plotting.generic_mpl import plot_xy
import hazard.rshalib as rshalib
from hazard.rshalib.source import SimpleUniformGridSourceModel
from hazard.rshalib.source_estimation import estimate_epicenter_location_and_magnitude_from_intensities

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SCRIPT_DIR)
from aysenlib import (create_uniform_grid_source_model, plot_gridsearch_map,
						read_evidence_site_info_from_gis,
						TRT, USD, LSD, RAR, RMS, project_folder, gis_folder)

fig_folder = os.path.join(project_folder, "Figures", "Events", "grid_probabilistic")

output_format = 'png'

## Event names
events = ['2007', 'SL-A', 'SL-B', 'SL-C', 'SL-CD', 'SL-D', 'SL-DE', 'SL-EF', 'SL-F', 'SL-G']
#events = ['SL-A']
#events = ["2007"]
#events = events[3:4]

## Method
method = 'probabilistic_mean'
#method = 'forward'

## IPE and IMT
ipe_name = "BakunWentworth1997WithSigma"
#ipe_name = "AtkinsonWald2007"
imt = oqhazlib.imt.MMI()

polygon_discretization = 2.5


## Construct grid source model
grid_outline = (-74, -72, -46.25, -44.8)
grid_spacing = 0.1

min_mag, max_mag, mag_bin_width = 4.5, 8.5, 0.2
depth = 10
strike, dip, rake = 20, 80, 180
point_msr = oqhazlib.scalerel.PointMSR()
wc1994_msr = oqhazlib.scalerel.WC1994()

grd_src_model = SimpleUniformGridSourceModel(grid_outline, grid_spacing,
			min_mag, max_mag + mag_bin_width, mag_bin_width, depth,
			strike, dip, rake, wc1994_msr, USD, LSD, RMS, RAR, TRT)
lon_grid, lat_grid = grd_src_model.lon_grid, grd_src_model.lat_grid


## Loop over events
for event in events:
	print(event)
	## Read MTD evidence
	pe_site_models, ne_site_models = [], []
	pe_thresholds, pe_sites, ne_thresholds, ne_sites = [], [], [], []
	for geom_type in ["Polygons_v3", "Points"]:
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
		print(sum([len(pesm) for pesm in _pe_site_models]), sum([len(nesm)
						for nesm in _ne_site_models]))
	print(len(pe_sites), len(ne_sites))
	pe_thresholds = np.array(pe_thresholds)
	ne_thresholds = np.array(ne_thresholds)

	## Test NE+0/+0.5/+1
	#ne_thresholds -= 0.5

	## Compute magnitudes and probabilities at grid points
	# TODO: check with Katleen!
	ne_margin = 0 if 'probabilistic' in method else 0.25
	result = estimate_epicenter_location_and_magnitude_from_intensities(
					ipe_name, imt, grd_src_model, pe_sites, pe_thresholds,
					ne_sites, ne_thresholds, method=method, ne_margin=ne_margin,
					mag_pdf_idx='max')
	if 'probabilistic' in method:
		(mag_grid, rms_grid, mag_pdf, pe_curves, ne_curves) = result
	else:
		(mag_grid, rms_grid) = result

	#idx = np.unravel_index(rms_grid.argmin(), rms_grid.shape)
	#print(mag_grid[idx], lon_grid[idx], lat_grid[idx])

	rms_grid[np.isinf(rms_grid)] = 10.0


	## Plot map
	# TODO: blend alpha in function of rms
	# See: https://matplotlib.org/devdocs/gallery/images_contours_and_fields/image_transparency_blend.html
	text_box = "Event: %s\n" % event
	if 'probabilistic' in method:
		text_box += "P: %.2f - %.2f"
	else:
		text_box += "RMSE: %.2f - %.2f"
	try:
		text_box %= (rms_grid.min(), rms_grid[rms_grid < 10].max())
	except:
		pass

	rms_is_prob = ('probabilistic' in method)
	site_model_gis_file = os.path.join(gis_folder, "Polygons_v3.shp")
	map = plot_gridsearch_map(grd_src_model, mag_grid, rms_grid,
							pe_site_models, ne_site_models,
							site_model_gis_file=site_model_gis_file,
							text_box=text_box, rms_is_prob=rms_is_prob,
							plot_rms_as_alpha=False, plot_epicenter_as="area")

	fig_filespec = os.path.join(fig_folder, "%s_grid_%s.%s")
	fig_filespec %= (event, method, output_format)
	#fig_filespec = None

	dpi = 200 if fig_filespec else 90
	map.plot(fig_filespec=fig_filespec, dpi=dpi)

	## Plot probabilities
	if 'probabilistic' in method:
		fig_filespec =  os.path.join(fig_folder, "%s_pdf_%s.%s")
		fig_filespec %= (event, method, output_format)
		#fig_filespec = None
		#mag_pdf.plot(fig_filespec=None, ylabel='Probability')

		datasets = []
		num_pe, num_ne = len(pe_sites), len(ne_sites)
		for p in range(num_pe):
			datasets.append((mag_pdf.values, pe_curves[p]))
		for n in range(num_ne):
			datasets.append((mag_pdf.values, ne_curves[n]))
		prod = np.prod(pe_curves, axis=0) * np.prod(ne_curves, axis=0)
		datasets.append((mag_pdf.values, prod))
		colors = ['m'] * num_pe + ['c'] * num_ne + ['k']
		labels = (['Positive'] + ['_nolegend_'] * (num_pe - 1)
					+ ['Negative'] + ['_nolegend_'] * (num_ne - 1)
					+ ['Product'])
		linestyles = ['-'] * (num_pe + num_ne) + ['--']
		title = 'Event %s' % event

		dpi = 200 if fig_filespec else 90
		plot_xy(datasets, colors=colors, labels=labels, linestyles=linestyles,
					xlabel='Magnitude', ylabel='Probability', title=title,
					fig_filespec=fig_filespec, dpi=dpi)
