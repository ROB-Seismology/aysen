
import os
import numpy as np
import matplotlib
import mapping.Basemap as lbm
import openquake.hazardlib as oqhazlib
import hazard.rshalib as rshalib
from hazard.rshalib.calc.calc import calc_rupture_probability_from_ground_motion_thresholds



def create_uniform_grid_source_model(grid_outline, grid_spacing, min_mag,
									max_mag, dM, depth, strike, dip, rake,
									trt, usd, lsd, rar, msr, rms):
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
			source = rshalib.source.PointSource(i, name, trt, mfd, rms, msr, rar,
												usd, lsd, point, npd, hdd)
			sources.append(source)
			i += 1
	return rshalib.source.SourceModel("Point source model", sources)


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

	return pe_thresholds, pe_site_models, ne_thresholds, ne_site_models


if __name__ == "__main__":
	## Construct uniform grid source model
	grid_outline = (-74, -72, -46, -44.5)
	grid_spacing = 0.1
	dM = 0.25
	min_mag, max_mag = 6.0 - dM/2, 7.5
	#min_mag, max_mag = 6.25 - dM/2, 6.25 + dM/2
	trt = "ASC"
	strike, dip, rake = 20, 90, 180
	depth = 5
	usd = 0
	lsd = 20
	rar = 1
	msr = "WC1994"
	rms = 2.5

	source_model = create_uniform_grid_source_model(grid_outline, grid_spacing,
								min_mag, max_mag, dM, depth, strike, dip, rake,
								trt, usd, lsd, rar, msr, rms)

	## Construct ground-motion model
	ipe_name = "AtkinsonWald2007"
	#ipe_name = "BakunWentworth1997WithSigma"
	trt_gsim_dict = {trt: ipe_name}
	ground_motion_model = rshalib.gsim.GroundMotionModel(ipe_name, trt_gsim_dict)
	truncation_level = 2.5

	soil_params = rshalib.site.REF_SOIL_PARAMS

	## Define positive evidence sites
	"""
	pe_site = rshalib.site.SoilSite(-73, -45.4, soil_params=soil_params)
	pe_sites = [rshalib.site.SoilSiteModel("Positive evidence", [pe_site])]
	pe_thresholds = [7.5]

	## Define negative evidence sites
	ne_site1 = rshalib.site.SoilSite(-73.2, -45.3, soil_params=soil_params)
	ne_site2 = rshalib.site.SoilSite(-72.7, -45.4, soil_params=soil_params)
	ne_sites = [rshalib.site.SoilSiteModel("Negative evidence", [ne_site1, ne_site2])]
	ne_thresholds = [7.0]
	"""

	event, version = "2007", 2
	filespec = "%s_polygons_v%d.txt" % (event, version)
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

	## Plot

	## Select non-zero probability rupture locations to be plotted
	x, y = [], []
	values = {'mag': [], 'prob': []}
	for source_id, probs in prob_dict.items():
		source = source_model[source_id]
		center_magnitudes = source.mfd.get_center_magnitudes()
		idx = probs.argmax()
		prob_max = probs[idx]
		if prob_max > 1E-5:
			values['prob'].append(prob_max)
			values['mag'].append(center_magnitudes[idx])
			x.append(source.location.longitude)
			y.append(source.location.latitude)

			print x[-1], y[-1], values['mag'][-1], prob_max

	layers = []

	## Coastlines
	data = lbm.BuiltinData("coastlines")
	style = lbm.LineStyle()
	layer = lbm.MapLayer(data, style)
	layers.append(layer)

	## Sources
	data = lbm.MultiPointData(x, y, values=values)
	#colors = matplotlib.cm.jet(np.arange(num_mags))
	#colors = ["purple", "blue", "green", "yellow", "orange", "red"]
	#colors = "random_color"
	colorbar_style = lbm.ColorbarStyle("Probability")
	#thematic_color = lbm.ThematicStyleGradient([1E-3, 1E-2, 1E-1, 1], "RdBu_r", value_key='prob', colorbar_style=colorbar_style)
	#thematic_color = lbm.ThematicStyleGradient([0.01, 0.05, 0.125, 0.25, 0.5, 1], "RdBu_r", value_key='prob', colorbar_style=colorbar_style)
	thematic_color = lbm.ThematicStyleColormap("Reds", vmin=0.01, vmax=0.3, value_key='prob', colorbar_style=colorbar_style)

	edge_magnitudes = np.concatenate([source.mfd.get_magnitude_bin_edges(), [center_magnitudes[-1]+dM/2]])
	mag_sizes = (center_magnitudes - 4) ** 2
	thematic_size = lbm.ThematicStyleRanges(edge_magnitudes, mag_sizes, value_key='mag', labels=center_magnitudes)
	#thematic_size = 8
	thematic_legend_style = lbm.LegendStyle("Magnitude", location="upper left")
	style = lbm.PointStyle(fill_color=thematic_color, size=thematic_size, thematic_legend_style=thematic_legend_style)
	layer = lbm.MapLayer(data, style)
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

	title = "Event: %s, IPE: %s, %s sigma" % (event, ipe_name, truncation_level)
	fig_filename = "%s_%s_%ssigma" % (event, ipe_name, truncation_level)
	#fig_filename += "_M=6.25"
	fig_filename += "_Mrange"
	if strict_intersection:
		title += ", strict"
		fig_filename += "_strict"
	fig_filename += "_v%d.PNG" % version

	#fig_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
	fig_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
	fig_filespec = os.path.join(fig_folder, fig_filename)
	#fig_filespec = None

	map = lbm.LayeredBasemap(layers, title, "merc", region=grid_outline,
							graticule_interval=(1, 0.5), resolution='h')
	map.plot(fig_filespec=fig_filespec, dpi=200)
