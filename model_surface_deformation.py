# -*- coding: iso-Latin-1 -*-

# Todo: have a look at:
# https://github.com/simondaout/Flower2d
# https://github.com/dmelgarm/MudPy
# https://github.com/scottyhq/roipy/blob/master/plot.py

# It is simply a matter of multiplying the unwrapped phase (in the natural radians units)
# by (radar wavelength)/(4 * pi) to get the apparent surface displacement in meters or whatever units you use for the radar wavelength


import os
import numpy as np

import mapping.layeredbasemap as lbm
from eqgeology.faultlib import okada

from aysenlib import (read_fault_source_model_as_network, read_fault_source_model,
					gis_folder)


## Folders
#project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")
fig_folder = os.path.join(project_folder, "Figures", "insar")

fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel3.TAB")


def read_fault_info(fault_id, num_sections):
	for M, flt_network in read_fault_source_model_as_network(fault_filespec, dM=None,
									num_sections=num_sections, characteristic=False):
		print("M=%.2f" % M)
		for flt in flt_network:
			if flt.source_id == fault_id:
				return flt


def plot_okada_map(U, component, wavelength, output="displacement", map_region=[],
					fault_gis_filespec="", elastic_fault=None, fig_filespec=None):
	"""
	Plot map of surface deformation modelled with Okada

	:param U:
		instance of :class:`eqgeology.faultlib.okada.OkadaDisplacementField`,
		modeled displacement field
	:param component:
		str or tuple, component specification
	:param output:
		str, "displacement", "phase" or "vector"
		(default: "displacement")
	:param map_region:
		(lonmin, lonmax, latmin, latmax) tuple
		(default: [], will use bounding box of :param:`U`)
	:param fault_gis_filespec:
		str, full path to GIS file containing background faults
		(default: "")
	:param elastic_fault:
		instance of :class:`eqgeology.faultlib.okada.ElasticFault`,
		fault rupture used to calculate displacement field
		(default: None)
	:param fig_filespec:
		str, full path to output file
		(default: None)
	"""

	from clawpack.visclaw import colormaps
	#import cmocean
	import matplotlib

	dZ = U.get_component(component)
	dzmax = dZ.range()
	print U.U.max(), dzmax

	layers = []

	## DEM
	"""
	url = 'http://seishaz.oma.be:8080/geoserver/wcs'
	layer_name, grid_resolution = "ngdc:etopo1_bedrock", 1./60
	wcs_data = lbm.WCSData(url, layer_name, resolution=grid_resolution, region=map_region)
	colorbar_style = lbm.ColorbarStyle("Elevation (m)", label_size=12)
	hillshade_style = None
	cmap = "gist_earth"
	tsc = lbm.ThematicStyleColormap(color_map=cmap, vmin=0, vmax=2000)
	style = lbm.GridStyle(color_map_theme=tsc, colorbar_style=colorbar_style, line_style=None, pixelated=False, hillshade_style=hillshade_style)
	layer = lbm.MapLayer(wcs_data, style)
	#layers.append(layer)
	"""


	## Displacement / Phase
	if output in ("displacement", "phase"):
		if output == "displacement":
			#contour_interval = dzmax / 10.
			contour_interval = 0.1
			contour_levels = np.arange(contour_interval, dzmax+contour_interval, contour_interval)
			contour_levels = list(-contour_levels[::-1]) + [0.] + list(contour_levels)

			grid_data = lbm.MeshGridData(dZ.X, dZ.Y, dZ.Z)
			color_map = colormaps.blue_white_red
			color_map_theme = lbm.ThematicStyleColormap(color_map=color_map, vmin=-0.5, vmax=0.5)
			colorbar_title = "%s displacement (m)" % dZ.comp_label
			#contour_levels = None
			contour_line_style = lbm.LineStyle(label_style=lbm.TextStyle())

		elif output == "phase":
			from matplotlib.colors import LinearSegmentedColormap
			delta_phi = dZ.get_wrapped_phase_difference(wavelength=wavelength,
											unit="degrees")
			grid_data = lbm.MeshGridData(dZ.X, dZ.Y, delta_phi)
			#color_map = cmocean.cm.phase
			#color_map = matplotlib.cm.hsv
			colors = ['b', 'r', 'yellow', 'c', 'b']
			color_map = LinearSegmentedColormap.from_list("InSAR phase", colors, N=256)
			color_map = lbm.cm.adjust_cmap_saturation(color_map, -0.40, colorspace='HSV')
			color_map_theme = lbm.ThematicStyleColormap(color_map=color_map, vmin=-180, vmax=180)
			colorbar_title = "%s phase (degrees)" % dZ.comp_label
			contour_levels = [-180.,180.]
			contour_line_style = lbm.LineStyle(label_style=None)

		## Mask oceans
		#from mpl_toolkits.basemap import maskoceans
		#grid_data.values = maskoceans(grid_data.lons, grid_data.lats, grid_data.values, resolution='h')

		colorbar_style = lbm.ColorbarStyle(colorbar_title, format="%.2f")
		#colorbar_style = None
		grid_style = lbm.GridStyle(color_map_theme, color_gradient="continuous",
						line_style=contour_line_style, contour_levels=contour_levels,
						colorbar_style=colorbar_style)

		## Add topographic hillshading
		"""
		#elevation_grid = lbm.GdalRasterData(r"D:\GIS-data\DEM\Etopo2.bin", region=map.region)
		url = 'http://seishaz.oma.be:8080/geoserver/wcs'
		layer_name, grid_resolution = "ngdc:etopo1_bedrock", 1./60
		elevation_grid = lbm.WCSData(url, layer_name, resolution=grid_resolution, region=map_region)
		blend_mode = "soft"
		hillshade_style = lbm.HillshadeStyle(0, 45, 1, blend_mode=blend_mode,
												elevation_grid=elevation_grid)
		grid_style.hillshade_style = hillshade_style
		grid_style.pixelated = True
		"""

		layer = lbm.MapLayer(grid_data, grid_style)
		layers.append(layer)

	elif output == "vector":
		grid_data1 = lbm.MeshGridData(X[::4,::4], Y[::4,::4], U.E[::4,::4])
		grid_data2 = lbm.MeshGridData(X[::4,::4], Y[::4,::4], U.N[::4,::4])
		vector_data = lbm.MeshGridVectorData(grid_data1, grid_data2, unit='m')
		colorbar_title = "Displacement (m)"
		contour_levels = None
		vector_style = lbm.VectorStyle(color='b', scale=0.005, width=1, pivot="tail",
									thematic_legend_style=None)
		layer = lbm.MapLayer(vector_data, vector_style)
		layers.append(layer)

		grid_data1 = lbm.MeshGridData(X, Y, np.zeros_like(U.E))
		grid_data2 = lbm.MeshGridData(X, Y, U.Z)
		vector_data = lbm.MeshGridVectorData(grid_data1, grid_data2, unit='m')
		colorbar_title = "Displacement (m)"
		contour_levels = None
		vector_style = lbm.VectorStyle(color='r', scale=0.005, width=1, pivot="tail",
									thematic_legend_style=None)
		layer = lbm.MapLayer(vector_data, vector_style)
		#layers.append(layer)

	text_box = ""
	if output and elastic_fault:
		subflt = elastic_fault.subfaults[0]
		text_box = u"Length: %.1f km\n" % (elastic_fault.get_length() / 1000)
		#text_box += "Depth: %.1f km\n" % flt.lower_seismogenic_depth
		text_box += u"Dip: %d°\n" % subflt.dip
		text_box += u"Rake: %d°\n" % subflt.rake
		#text_box += "Slip: %s m\n" % '/'.join(["%.2f" % subflt.slip for subflt in elastic_fault.subfaults])
		text_box += "Mean slip: %.2f m\n" % elastic_fault.calc_average_slip()
		text_box += "Magnitude: %.1f\n" % elastic_fault.calc_magnitude()
		text_box += "Displacement: %.2f/%.2f m" % (dZ.min(), dZ.max())

	## Coastlines
	coastline_style = countries_style = lbm.LineStyle(line_width=0.75, line_color="k")
	data = lbm.BuiltinData("coastlines")
	layer = lbm.MapLayer(data, coastline_style)
	layers.append(layer)

	#continent_style = lbm.PolygonStyle(fill_color="none", line_width=0)
	#continent_style.bg_color = "white"
	#data = lbm.BuiltinData("continents")
	#layer = lbm.MapLayer(data, continent_style)
	#layers.append(layer)

	## Faults
	data = lbm.GisData(fault_gis_filespec)
	style = lbm.LineStyle(line_color='purple', line_width=1.5)
	layer = lbm.MapLayer(data, style, legend_label="Faults")
	layers.append(layer)

	if elastic_fault:
		top_subfaults = elastic_fault.get_top_subfaults()
		lons, lats = [], []
		for subflt in top_subfaults:
			A, B, C, D = subflt.corners
			lons.append([D[0], A[0]])
			lats.append([D[1], A[1]])
		data = lbm.MultiLineData(lons, lats)
		#data = flt.to_lbm_data()
		style = lbm.LineStyle(line_color='m', line_width=3)
		layer = lbm.MapLayer(data, style, legend_label="Faults")
		layers.append(layer)

	legend_style = None
	title = ""
	graticule_interval = (0.2, 0.1)
	if not map_region:
		map_region = [U.X.min(), U.X.max(), U.Y.min(), U.Y.max()]
	map = lbm.LayeredBasemap(layers, title, projection="merc", region=map_region,
			title_style=lbm.DefaultTitleTextStyle, graticule_style=lbm.GraticuleStyle(),
			graticule_interval=graticule_interval, resolution='f',
			legend_style=legend_style)

	## Add text box
	if text_box:
		pos = (0.965, 0.965)
		text_style = lbm.TextStyle(font_size=10, horizontal_alignment='right',
							vertical_alignment='top', multi_alignment='center',
							background_color='white', border_pad=0.4, border_color='k')
		map.draw_text_box(pos, text_box, text_style)


	if fig_filespec:
		dpi = 200
	else:
		dpi = 90
	map.plot(fig_filespec=fig_filespec, dpi=dpi)



if __name__ == "__main__":
	## Map parameters
	#map_region = (-74, -72, -46.25, -44.75)
	#graticule_interval = (1, 0.5)
	#map_region = (-73.5, -72.5, -45.675, -45.175)
	#graticule_interval = (0.5, 0.25)
	#map_region = (-73.25, -72.75, -45.575, -45.2)
	map_region = (-73.22, -72.77, -45.56, -45.22)

	## Read fault model
	#for flt in read_fault_source_model(fault_filespec):
	#	print flt.source_id, flt.name
	#exit()
	"""
	#flt_id = "2#01+2#02+2#03+2#04+2#05"
	flt_id = "0#05+0#06+0#07+0#08"
	num_sections = len(flt_id.split('+'))
	flt = read_fault_info(flt_id, num_sections)
	"""

	gis_filespec = os.path.join(gis_folder, "InSAR_fault_rupture.TAB")
	src_model = read_fault_source_model(gis_filespec, characteristic=False)

	slip_distribution = {0: [0.1, 0.2, 0.3, 0.4, 0.5, 0.5],
						1: [0.01, 0.1, 0.1, 0.15, 0.2, 0.25],
						2: [0.75, 0.75, 0.5, 0.25, 0.15, 0.1],
						3: [0.2]*5
					}

	#slip_distribution = {0: [0.4]*6, 1:[0.2]*6, 2: [0.6]*6, 3: [0.1]*5}

	elastic_faults = []
	#flt_idxs = [0, 1, 2, 3]
	flt_idxs = [0, 2]
	for flt_idx in flt_idxs:
		flt = src_model.sources[flt_idx]
		#flt.reverse_trace()

		## Override geometry (dip, depth, rake)
		#flt.dip = 88
		flt.dip = 90
		flt.lower_seismogenic_depth = 14
		flt.rake = 176
		#if flt_idx == 3:
		#	flt.rake = -4
		print flt.get_dip_direction()

		#flt.mfd.min_mag = 5.0
		#flt.plot_rupture_bounds_3d(flt.get_ruptures_Poisson())
		#exit()

		section_length = 1.5
		num_sections = int(round(flt.get_length() / section_length))
		subfaults = flt.get_subfaults(num_sections, 1, rigidity=3E+10)
		print subfaults.shape

		## Override kinematics (rake, mu)
		## Note: dip would only work if fault is not subdivided downdip
		for subflt in subfaults.flatten():
			#subflt.rake = 176
			subflt.rake = 180 - (28 - subflt.strike)

		elastic_fault = okada.create_fault(subfaults)

		## Set slip distribution
		for s, slip in enumerate(slip_distribution[flt_idx]):
			elastic_fault.subfaults[s].slip = slip

		elastic_fault.subfaults = list(elastic_fault.subfaults)
		elastic_faults.append(elastic_fault)

	elastic_fault = elastic_faults[0]
	for ef in elastic_faults[1:]:
		elastic_fault.subfaults.extend(ef.subfaults)

	## Set slip distribution
	mag = 6.3
	#elastic_fault.set_slip_from_magnitude(mag)
	#elastic_fault.taper_slip()
	for s, subflt in enumerate(elastic_fault.subfaults):
		print s, subflt.slip, subflt.calc_magnitude(), subflt.strike, subflt.rake
	print elastic_fault.calc_magnitude()

	## Compute moment tensor
	mt = elastic_fault.get_moment_tensor()
	mt.plot_text()
	print(mt.getMW())
	exit()

	elastic_fault.plot_3D()
	#exit()

	## Compute surface deformation
	num_pts = 101
	#num_pts = 25
	x = np.linspace(map_region[0], map_region[1], num_pts)
	y = np.linspace(map_region[-2], map_region[-1], num_pts)
	X, Y = np.meshgrid(x, y)
	U = elastic_fault.okada(X, Y)


	## InSAR parameters
	## wavelength (in m)
	wavelength = 0.2360571
	## Azimuth (satellite looks right)
	insar_az = -11.7 + 90
	## Off-nadir angle
	off_nadir = 34.3

	## Select component
	#component = 'E'
	## Along-strike
	#azimuth = np.round(elastic_fault.subfaults[0].strike)
	#elevation_angle = 0
	#comp_string = "AS"
	## LOS (satellite looks right)
	azimuth = insar_az + 180
	elevation_angle = 90 - off_nadir
	#azimuth, elevation_angle = 0, 75
	component = (azimuth, elevation_angle, "up")


	## Plot
	#output = "displacement"
	output = "phase"
	#output = "vector"
	#output = None

	fault_gis_filename = "LOFZ_breukenmodel.shp"
	fault_gis_filespec = os.path.join(gis_folder, fault_gis_filename)

	fig_filespec = None
	dZ = U.get_component(component)
	fig_filename = "okada_test_%s_%s.png" % (output, dZ.comp_label)
	#fig_filespec = os.path.join(fig_folder, fig_filename)

	plot_okada_map(U, component, wavelength, output=output, map_region=map_region,
					fault_gis_filespec=fault_gis_filespec, elastic_fault=elastic_fault,
					fig_filespec=fig_filespec)

	exit()

	## Test inverting for slip distribution
	original_slip_distribution = [subflt.slip for subflt in elastic_fault.subfaults]
	#elastic_fault.set_slip_from_magnitude(mag-0.5)
	for subflt in elastic_fault.subfaults:
		subflt.slip = 0
	#phase_info = None
	#D_obs = dZ.Z
	phase_info = (wavelength, "degrees")
	D_obs = dZ.get_wrapped_phase_difference(*phase_info)

	## Apply mask
	grid_data = lbm.MeshGridData(X, Y, D_obs)
	pg_gis_file = os.path.join(fig_folder, "insar_reliable_outline.TAB")
	gis_data = lbm.GisData(pg_gis_file)
	_, _, pg_data = gis_data.get_data()
	grid_data.mask_polygons(pg_data, inside=False)
	D_obs = grid_data.values

	inverted_slip_distribution = elastic_fault.invert_slip_distribution(X, Y, D_obs,
									component, max_slip=2, phase_info=phase_info)
	import pylab
	pylab.plot(original_slip_distribution, label="Original")
	pylab.plot(inverted_slip_distribution, label="Inverted")
	pylab.legend()
	pylab.show()


