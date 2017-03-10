
import os
import numpy as np
import matplotlib
import pylab
import mapping.Basemap as lbm
import openquake.hazardlib as oqhazlib
import hazard.rshalib as rshalib
from hazard.rshalib.calc.calc import calc_rupture_probability_from_ground_motion_thresholds
from hazard.rshalib.source.read_from_gis import import_source_model_from_gis


## Common parameters for area and fault sources
TRT = "ASC"
USD = 0
LSD = 15
RAR = 1
MSR = "WC1994"
RMS = 2.5


def create_uniform_grid_source_model(grid_outline, grid_spacing, min_mag,
									max_mag, dM, depth, strike, dip, rake):
	num_mags = int(round((max_mag - min_mag) / dM))
	mfd = rshalib.mfd.EvenlyDiscretizedMFD(min_mag + dM/2, dM, np.ones(num_mags)/float(num_mags))

	nopl = rshalib.geo.NodalPlane(strike, dip, rake)
	npd = rshalib.pmf.NodalPlaneDistribution([nopl], [1])

	hdd = rshalib.pmf.HypocentralDepthDistribution([depth], [1])

	lons = np.arange(grid_outline[0], grid_outline[1] + grid_spacing, grid_spacing)
	lats = np.arange(grid_outline[2], grid_outline[3] + grid_spacing, grid_spacing)
	sources = []
	i = 0
	for lon in lons:
		for lat in lats:
			point = rshalib.geo.Point(lon, lat)
			name = "%.2f, %.2f" % (lon, lat)
			source = rshalib.source.PointSource(i, name, TRT, mfd, RMS, MSR, RAR,
												USD, LSD, point, npd, hdd)
			sources.append(source)
			i += 1
	return rshalib.source.SourceModel("Grid", sources)


def read_fault_source_model(gis_filespec, characteristic=True):
	column_map = {
		'id': '#',
		'name': 'Name',
		'tectonic_region_type': TRT,
		'rupture_aspect_ratio': RAR,
		'upper_seismogenic_depth': USD,
		'lower_seismogenic_depth': LSD,
		'magnitude_scaling_relationship': MSR,
		'rupture_mesh_spacing': RMS,
		'min_dip': 'Dip',
		'max_dip': 'Dip',
		'min_mag': None,
		'max_mag': None,
		'min_rake': 'Rake',
		'max_rake': 'Rake',
		'min_slip_rate': 1,
		'max_slip_rate': 1,
		'bg_zone': 'Id'
	}

	fault_ids = []

	somo = import_source_model_from_gis(gis_filespec, column_map=column_map,
										source_ids=fault_ids)
	if characteristic:
		for i, flt in enumerate(somo.sources):
			somo.sources[i] = flt.to_characteristic_source(convert_mfd=True)
			somo.sources[i].mfd.set_num_sigma(0)

	return somo


def read_fault_source_model_as_points(gis_filespec, min_mag, max_mag, dM, depth):
	from copy import deepcopy
	from openquake.hazardlib.geo.surface.simple_fault import SimpleFaultSurface

	num_mags = int(round((max_mag - min_mag) / dM))
	mfd = rshalib.mfd.EvenlyDiscretizedMFD(min_mag + dM/2, dM, np.ones(num_mags)/float(num_mags))

	fault_somo = read_fault_source_model(gis_filespec, characteristic=False)
	## Divide fault trace in points
	sources = []
	for f, flt in enumerate(fault_somo.sources):
		fault_surface = SimpleFaultSurface.from_fault_data(flt.fault_trace, USD, LSD, flt.dip, RMS)
		fault_mesh = fault_surface.get_mesh()
		surface_locations = fault_mesh[0:1]
		#print surface_locations.lons
		#print surface_locations.lats
		#print flt.get_length()
		#print (len(surface_locations) - 1) * RMS
		mesh_rows, mesh_cols = fault_mesh.shape
		hypo_depths = []
		for j in range(1, mesh_rows-1):
			mesh = fault_mesh[j-1:j+2,:]
			hypocenter = mesh.get_middle_point()
			hypo_depths.append(hypocenter.depth)
		hypo_idx = np.abs(np.array(hypo_depths) - depth).argmin()
		hypo_idx += 1
		for i in range(1, mesh_cols-1):
			mesh = fault_mesh[hypo_idx-1:hypo_idx+1, i-1:i+1]
			dip, strike = mesh.get_mean_inclination_and_azimuth()
			hypocenter = mesh.get_middle_point()
			distance_to_start = i * RMS - RMS/2.
			distance_to_end = (mesh_cols - i) * RMS - RMS/2.
			nodal_plane = rshalib.geo.NodalPlane(strike, dip, flt.rake)
			npd = rshalib.pmf.NodalPlaneDistribution([nodal_plane], [1])
			hdd = rshalib.pmf.HypocentralDepthDistribution([depth], [1])
			name = "%s #%02d" % (flt.name, i+1)
			ID = flt.source_id + "#%02d" % (i+1)

			point_source = rshalib.source.PointSource(ID, name, TRT, mfd, RMS, MSR, RAR,
												USD, LSD, hypocenter, npd, hdd)
			## Check if rupture stays within fault limits
			for mag in mfd.get_center_magnitudes():
				rup_length, rup_width = point_source._get_rupture_dimensions(mag, nodal_plane)
				if rup_length / 2 <= min(distance_to_start, distance_to_end):
					pt_col = i + 0.5
					rup_col_num = int(round(rup_length / RMS))
					start_col = min(i, int(i - (rup_col_num / 2.)))
					end_col = max(i, int(i + (rup_col_num / 2.)))
					#rup_row_num = int(round(rup_width / RMS))
					subfault_trace = list(surface_locations)[start_col:end_col+1]
					subfault_trace = oqhazlib.geo.Line(subfault_trace)
					subfault_mfd = rshalib.mfd.EvenlyDiscretizedMFD(mag, mfd.bin_width, [1])
					subfault_source_id = ID + "_M=%s" % mag
					subfault_lsd =  rup_width * np.cos(np.radians(90 - flt.dip))
					subfault = rshalib.source.SimpleFaultSource(subfault_source_id, flt.name,
									TRT, subfault_mfd, RMS, MSR, RAR, USD, subfault_lsd,
									subfault_trace, flt.dip, flt.rake)

					#subfault = deepcopy(point_source)
					#subfault.mfd = rshalib.mfd.EvenlyDiscretizedMFD(mag, mfd.bin_width, [1])
					#subfault.source_id = ID + "_M=%s" % mag

					sources.append(subfault)
				#print mag, rup_length/2, distance_to_start, distance_to_end

	somo_name = fault_somo.name + "_pts"
	return rshalib.source.SourceModel(somo_name, sources)


def read_evidence_site_info(filespec, polygon_discretization=5):
	pe_thresholds, ne_thresholds = [], []
	pe_polygons, ne_polygons = [], []

	with open(filespec) as f:
		lons, lats = [], []
		for line in f:
			line = line.strip()
			if line and line[0] in ('<', '>'):
				intensity = float(line[2:])
				points = [oqhazlib.geo.Point(lon, lat) for (lon, lat) in zip(lons, lats)]
				if len(points) > 1:
					polygon = oqhazlib.geo.Polygon(points)
				else:
					[polygon] = points
				if line[0] == '<':
					ne_thresholds.append(intensity)
					ne_polygons.append(polygon)
				else:
					pe_thresholds.append(intensity)
					pe_polygons.append(polygon)
				lons, lats = [], []
			elif line:
				if line[-1] in ('W', 'E', 'N', 'S'):
					deg = int(line[:2])
					min = int(line[3:5])
					sec = int(line[6:8])

					deg = deg + (min / 60.) + (sec / 3600.)
					if line[-1] in ('W', 'S'):
						deg = -deg
					if line[-1] in ('W', 'E'):
						lons.append(deg)
					else:
						lats.append(deg)

	def polygon_to_site_model(polygon, name, polygon_discretization):
		if isinstance(polygon, oqhazlib.geo.Polygon):
			try:
				site_model = rshalib.site.SoilSiteModel.from_polygon(polygon,
											polygon_discretization, name=name)
			except:
				polygon = lbm.PolygonData(polygon.lons, polygon.lats)
				centroid = polygon.get_centroid()
				site = rshalib.site.SoilSite(centroid.lon, centroid.lat)
				site_model = rshalib.site.SoilSiteModel(name, [site])
		else:
			point = polygon
			site = rshalib.site.SoilSite(point.longitude, point.latitude)
			site_model = rshalib.site.SoilSiteModel(name, [site])
		return site_model

	pe_site_models = []
	for p, pe_polygon in enumerate(pe_polygons):
		name = "Positive evidence #%d (I>%.1f)" % (p+1, pe_thresholds[p])
		site_model = polygon_to_site_model(pe_polygon, name, polygon_discretization)
		pe_site_models.append(site_model)

	ne_site_models = []
	for n, ne_polygon in enumerate(ne_polygons):
		name = "Negative evidence #%d (I<%.1f)" % (n+1, ne_thresholds[n])
		site_model = polygon_to_site_model(ne_polygon, name, polygon_discretization)
		ne_site_models.append(site_model)

	pe_thresholds, ne_thresholds = np.array(pe_thresholds), np.array(ne_thresholds)
	return pe_thresholds, pe_site_models, ne_thresholds, ne_site_models


def plot_rupture_probabilities(source_model, prob_dict, pe_site_models, ne_site_models,
								region, max_prob_color=0.5, plot_point_ruptures=True, fig_filespec=None):

	## Extract source locations
	x, y = [], []
	values = {'mag': [], 'prob': []}
	PROB_MIN = 1E-5
	if source_model.get_point_sources():
		## Point sources and discretized fault sources
		for source_id, probs in prob_dict.items():
			source = source_model[source_id]
			center_magnitudes = source.mfd.get_center_magnitudes()
			idx = probs.argmax()
			prob_max = probs[idx]
			## Select non-zero probability rupture locations to be plotted
			if prob_max > PROB_MIN:
				values['prob'].append(prob_max)
				mag = source.mfd.get_center_magnitudes()[idx]
				values['mag'].append(mag)

				if not plot_point_ruptures:
					x.append(source.location.longitude)
					y.append(source.location.latitude)
					print x[-1], y[-1], values['mag'][-1], prob_max
				else:
					## Not sure this is correct if fault is not vertical
					## Point source ruptures
					[nodal_plane] = source.nodal_plane_distribution.nodal_planes
					hypocenter = source.location
					#hypocenter.depth = 0.1
					rup_surface = source._get_rupture_surface(mag, nodal_plane, hypocenter)
					top_left = rup_surface.top_left
					top_right = rup_surface.top_right
					# TODO: extract all top coordinates!
					x.append([top_left.longitude, top_right.longitude])
					y.append([top_left.latitude, top_right.latitude])

		max_prob = np.max(values['prob'])
		print("Max. probability: %.3f" % max_prob)
		"""
		if not plot_point_ruptures:
			source_data = lbm.MultiPointData(x, y, values=values)
		else:
			## Reorder from lowest to highest probability
			idxs = np.argsort(values['prob'])
			values['prob'] = [values['prob'][idx] for idx in idxs]
			values['mag'] = [values['mag'][idx] for idx in idxs]
			x = [x[idx] for idx in idxs]
			y = [y[idx] for idx in idxs]
			source_data = lbm.MultiLineData(x, y, values=values)
		"""

	else:
		## Fault sources
		for fault_id, [prob] in prob_dict.items():
			fault = source_model[fault_id]
			lons = np.array([pt.longitude for pt in fault.fault_trace.points])
			lats = np.array([pt.latitude for pt in fault.fault_trace.points])
			values['mag'].append(fault.mfd.get_center_magnitudes()[0])
			values['prob'].append(prob)
			x.append(lons)
			y.append(lats)
		#source_data = lbm.MultiLineData(x, y, values=values)

	## Reorder from lowest to highest probability
	idxs = np.argsort(values['prob'])
	values['prob'] = [values['prob'][idx] for idx in idxs]
	values['mag'] = [values['mag'][idx] for idx in idxs]
	x = [x[idx] for idx in idxs]
	y = [y[idx] for idx in idxs]
	source_data = lbm.MultiLineData(x, y, values=values)
	if source_model.get_point_sources() and not plot_point_ruptures:
		source_data = lbm.MultiPointData(x, y, values=values)
	else:
		source_data = lbm.MultiLineData(x, y, values=values)

	## Plot histogram of probabilities
	"""
	bin_width = 0.02
	xmax = np.ceil(max_prob / bin_width) * bin_width
	num_bins = xmax / bin_width + 1
	bins = np.linspace(0, xmax, num_bins)
	pylab.hist(values['prob'], bins=bins, log=False)
	pylab.show()
	"""


	layers = []

	## Coastlines
	data = lbm.BuiltinData("coastlines")
	style = lbm.LineStyle()
	layer = lbm.MapLayer(data, style)
	layers.append(layer)

	## Sources
	#colors = matplotlib.cm.jet(np.arange(num_mags))
	#colors = ["purple", "blue", "green", "yellow", "orange", "red"]
	#colors = "random_color"
	colorbar_style = lbm.ColorbarStyle("Probability")
	#thematic_color = lbm.ThematicStyleGradient([1E-3, 1E-2, 1E-1, 1], "RdBu_r", value_key='prob', colorbar_style=colorbar_style)
	#thematic_color = lbm.ThematicStyleGradient([0.01, 0.05, 0.125, 0.25, 0.5, 1], "RdBu_r", value_key='prob', colorbar_style=colorbar_style)
	thematic_color = lbm.ThematicStyleColormap("Reds", vmin=0.001, vmax=max_prob_color, value_key='prob', colorbar_style=colorbar_style, alpha=1)
	#thematic_color.color_map.set_under('w')

	if source_model.get_point_sources() and not plot_point_ruptures:
		edge_magnitudes = np.concatenate([source.mfd.get_magnitude_bin_edges(), [center_magnitudes[-1]+dM/2]])
		mag_sizes = (center_magnitudes - 4) ** 2
		thematic_size = lbm.ThematicStyleRanges(edge_magnitudes, mag_sizes, value_key='mag', labels=center_magnitudes)
		#thematic_size = 8
		thematic_legend_style = lbm.LegendStyle("Magnitude", location="upper left")
		style = lbm.PointStyle(fill_color=thematic_color, size=thematic_size, thematic_legend_style=thematic_legend_style)
	else:
		style = lbm.LineStyle(line_color=thematic_color, line_width=2, thematic_legend_style="main")
		layer = lbm.MapLayer(source_data, style)
	layers.append(layer)

	## Positive evidence
	for pe_site_model in pe_site_models:
		pe_style = lbm.PointStyle('+', size=8, line_width=1, line_color='m')
		pe_data = lbm.MultiPointData(pe_site_model.lons, pe_site_model.lats)
		layer = lbm.MapLayer(pe_data, pe_style)
		layers.append(layer)

	## Negative evidence
	for ne_site_model in ne_site_models:
		ne_style = lbm.PointStyle('_', size=8, line_width=1, line_color='c')
		ne_data = lbm.MultiPointData(ne_site_model.lons, ne_site_model.lats)
		layer = lbm.MapLayer(ne_data, ne_style)
		layers.append(layer)

	map = lbm.LayeredBasemap(layers, title, "merc", region=region,
							graticule_interval=(1, 0.5), resolution='h')
	if fig_filespec:
		dpi = 200
	else:
		dpi = 90
	map.plot(fig_filespec=fig_filespec, dpi=dpi)



if __name__ == "__main__":
	## Construct uniform grid source model
	grid_outline = (-74, -72, -46, -44.5)
	#grid_outline = (-73, -73, -45.4, -45.4)
	grid_spacing = 0.1
	strike, dip, rake = 20, 90, 180
	depth = 5

	dM = 0.25
	min_mag, max_mag = 6.0 - dM/2, 6.75
	print np.arange(min_mag, max_mag, dM) + dM/2
	#min_mag, max_mag = 6.25 - dM/2, 6.25 + dM/2

	#source_model = create_uniform_grid_source_model(grid_outline, grid_spacing,
	#							min_mag, max_mag, dM, depth, strike, dip, rake)


	## Read fault source model
	#gis_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
	gis_folder = r"E:\Home\_kris\Publications\2017 - Aysen"

	#gis_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel.shp")
	gis_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel2.TAB")
	#gis_filespec = os.path.join(gis_folder, "EnergiAustral_faults.TAB")
	source_model = read_fault_source_model(gis_filespec)
	#for flt in source_model:
	#	print flt.source_id, flt.name, flt.mfd.char_mag, flt.mfd.max_mag

	#for mag in np.arange(6., 7.5, 0.1):
	#	rup_cols, rup_rows = flt._get_rupture_dimensions(flt.get_length(), flt.get_width(), mag)
	#	print mag, rup_cols, rup_rows

	#source_model = read_fault_source_model_as_points(gis_filespec, min_mag, max_mag, dM, depth=0.1)

	## Construct ground-motion model
	#ipe_name = "AtkinsonWald2007"
	#ipe_name = "BakunWentworth1997WithSigma"
	ipe_name = "AllenEtAl2012Rrup"
	trt_gsim_dict = {TRT: ipe_name}
	ground_motion_model = rshalib.gsim.GroundMotionModel(ipe_name, trt_gsim_dict)
	truncation_level = 2.5
	if "AllenEtAl2012" in ipe_name:
		truncation_level /= 2

	soil_params = rshalib.site.REF_SOIL_PARAMS


	## Read lake evidence
	event, version = "2007", "3b"
	filespec = "%s_polygons_v%s.txt" % (event, version)
	polygon_discretization = 2.5
	(pe_thresholds, pe_site_models,
	ne_thresholds, ne_site_models) = read_evidence_site_info(filespec, polygon_discretization)


	imt = oqhazlib.imt.MMI()

	## Compute rupture probabilities
	strict_intersection = True
	prob_dict = calc_rupture_probability_from_ground_motion_thresholds(
						source_model, ground_motion_model, imt, pe_site_models,
						pe_thresholds, ne_site_models, ne_thresholds, truncation_level,
						strict_intersection=strict_intersection)
	#print prob_dict


	## Plot
	title = "Event: %s v.%s, IPE: %s, %s sigma" % (event, version, ipe_name, truncation_level)
	fig_filename = "%s_%s_%ssigma" % (event, ipe_name, truncation_level)
	#fig_filename += "_M=6.25"
	fig_filename += "_Mrange"
	if strict_intersection:
		title += ", strict"
		fig_filename += "_strict"
	fig_filename += "_v%s" % version
	fig_filename += "_%s" % source_model.name
	fig_filename += ".PNG"

	#fig_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
	fig_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
	fig_filespec = os.path.join(fig_folder, fig_filename)
	#fig_filespec = None

	max_prob_color = 0.5
	#max_prob_color = 1.0
	plot_rupture_probabilities(source_model, prob_dict, pe_site_models, ne_site_models,
								grid_outline, max_prob_color, plot_point_ruptures=True, fig_filespec=fig_filespec)
