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

from aysenlib import (read_fault_source_model_as_network, read_fault_source_model)


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


if __name__ == "__main__":
	#for flt in read_fault_source_model(fault_filespec):
	#	print flt.source_id, flt.name
	#exit()
	#flt_id = "2#01+2#02+2#03+2#04+2#05"
	flt_id = "0#05+0#06+0#07"
	num_sections = 3
	flt = read_fault_info(flt_id, num_sections)

	## Override geometry (dip, depth)
	flt.dip = 88
	flt.lower_seismogenic_depth = 14
	flt.rake = 176

	subfaults = flt.get_subfaults(num_sections, 1)

	## Override kinematics (rake, mu)
	## Note: dip would only work if fault is not subdivided downdip
	for subflt in subfaults:
		subflt.mu = 3E+10
		subflt.rake = 176
		#subflt.calculate_geometry()

	## Map parameters
	#map_region = (-74, -72, -46.25, -44.75)
	#graticule_interval = (1, 0.5)
	#map_region = (-73.5, -72.5, -45.675, -45.175)
	#graticule_interval = (0.5, 0.25)
	#map_region = (-73.25, -72.75, -45.575, -45.2)
	map_region = (-73.22, -72.77, -45.56, -45.22)
	graticule_interval = (0.2, 0.1)

	## Compute surface deformation
	elastic_fault = okada.create_fault(subfaults)
	elastic_fault.set_slip_from_magnitude(6.2)
	elastic_fault.taper_slip()
	for s, subflt in enumerate(elastic_fault.subfaults):
		print s, subflt.slip, subflt.calc_magnitude()
	print elastic_fault.calc_magnitude()

	num_pts = 101
	#num_pts = 25
	x = np.linspace(map_region[0], map_region[1], num_pts)
	y = np.linspace(map_region[-2], map_region[-1], num_pts)
	X, Y = np.meshgrid(x, y)
	U = elastic_fault.okada(X, Y)

	## Take component
	#component = 'E'
	#dZ = getattr(U, component)
	## Along-strike
	#azimuth = np.round(elastic_fault.subfaults[0].strike)
	#elevation_angle = 0
	#comp_string = "AS"
	## LOS (satellite looks right)
	azimuth = -11.7 + 90 + 180
	elevation_angle = 90 - 34.3
	component = (azimuth, elevation_angle)
	dZ = U.get_los_component(*component)
	comp_string = component if not isinstance(component, tuple) else "LOS"
	dzmax = np.abs(dZ.max() - dZ.min())
	print U.U.max(), dzmax


	## Plot
	from clawpack.visclaw import colormaps
	#import cmocean
	import matplotlib

	#output = "displacement"
	output = "phase"
	#output = "vector"

	layers = []

	## InSAR image
	img_file = os.path.join(fig_folder, "insar.jpg")
	data = lbm.ImageData(img_file, 0, 0, coord_frame="display")
	style = lbm.ImageStyle(horizontal_alignment='stretch', vertical_alignment='stretch', on_top=False)
	layer = lbm.MapLayer(data, style)
	#layers.append(layer)

	## Displacement / Phase
	if output in ("displacement", "phase"):
		if output == "displacement":
			#contour_interval = dzmax / 10.
			contour_interval = 0.1
			contour_levels = np.arange(contour_interval, dzmax+contour_interval, contour_interval)
			contour_levels = list(-contour_levels[::-1]) + [0.] + list(contour_levels)

			grid_data = lbm.MeshGridData(X, Y, dZ)
			color_map = colormaps.blue_white_red
			color_map_theme = lbm.ThematicStyleColormap(color_map=color_map, vmin=-0.5, vmax=0.5)
			colorbar_title = "%s displacement (m)" % comp_string
			#contour_levels = None
			contour_line_style = lbm.LineStyle(label_style=lbm.TextStyle())

		elif output == "phase":
			wavelength = 0.2362
			#grid_data = lbm.MeshGridData(X, Y, (dZ * 4 * np.pi / wavelength) % (2 * np.pi)
			grid_data = lbm.MeshGridData(X, Y, (dZ * 4 * 180. / wavelength) % (360.))
			#color_map = cmocean.cm.phase
			color_map = matplotlib.cm.hsv
			color_map_theme = lbm.ThematicStyleColormap(color_map=color_map, vmin=0, vmax=360, alpha=0.5)
			colorbar_title = "%s phase (degrees)" % comp_string
			contour_levels = [0.,360.]
			contour_line_style = lbm.LineStyle(label_style=None)

		colorbar_style = lbm.ColorbarStyle(colorbar_title, format="%.2f")
		#colorbar_style = None
		grid_style = lbm.GridStyle(color_map_theme, color_gradient="continuous",
						line_style=contour_line_style, contour_levels=contour_levels,
						colorbar_style=colorbar_style)
		layer = lbm.MapLayer(grid_data, grid_style)
		layers.append(layer)

	elif output == "vector":
		grid_data1 = lbm.MeshGridData(X, Y, U.E)
		grid_data2 = lbm.MeshGridData(X, Y, U.N)
		vector_data = lbm.MeshGridVectorData(grid_data1, grid_data2, unit='m')
		color_map = colormaps.blue_white_red
		color_map_theme = lbm.ThematicStyleColormap(color_map=color_map, vmin=-0.5, vmax=0.5)
		colorbar_title = "Displacement (m)"
		contour_levels = None
		vector_style = lbm.VectorStyle(color='b', scale=0.005, width=1, pivot="tail",
									thematic_legend_style=color_map_theme)
		layer = lbm.MapLayer(vector_data, vector_style)
		layers.append(layer)

	## Coastlines
	coastline_style = countries_style = lbm.LineStyle(line_width=0.75, line_color="dimgrey")
	data = lbm.BuiltinData("coastlines")
	layer = lbm.MapLayer(data, coastline_style)
	layers.append(layer)

	## Faults
	gis_filename = "LOFZ_breukenmodel.shp"
	gis_filespec = os.path.join(gis_folder, gis_filename)
	data = lbm.GisData(gis_filespec)
	style = lbm.LineStyle(line_color='purple', line_width=1.5)
	layer = lbm.MapLayer(data, style, legend_label="Faults")
	layers.append(layer)

	data = flt.to_lbm_data()
	style = lbm.LineStyle(line_color='m', line_width=3)
	layer = lbm.MapLayer(data, style, legend_label="Faults")
	layers.append(layer)

	legend_style = None
	title = ""
	map = lbm.LayeredBasemap(layers, title, projection="merc", region=map_region,
			title_style=lbm.DefaultTitleTextStyle, graticule_style=lbm.GraticuleStyle(),
			graticule_interval=graticule_interval, resolution='f',
			legend_style=legend_style)

	## Add text box
	pos = (0.965, 0.965)
	text_style = lbm.TextStyle(font_size=10, horizontal_alignment='right',
						vertical_alignment='top', multi_alignment='center',
						background_color='white', border_pad=0.4, border_color='k')
	text = u"Length: %.1f km\n" % flt.get_length()
	text += "Depth: %.1f km\n" % flt.lower_seismogenic_depth
	text += u"Dip: %d°\n" % flt.dip
	text += u"Rake: %d°\n" % flt.rake
	text += "Slip: %s m\n" % '/'.join(["%.2f" % subflt.slip for subflt in elastic_fault.subfaults])
	text += "Magnitude: %.1f\n" % elastic_fault.calc_magnitude()
	text += "Displacement: %.2f/%.2f m" % (dZ.min(), dZ.max())
	map.draw_text_box(pos, text, text_style)

	fig_filespec = None
	fig_filename = "okada_test_%s_%s.png" % (output, comp_string)
	#fig_filespec = os.path.join(fig_folder, fig_filename)

	if fig_filespec:
		dpi = 200
	else:
		dpi = 90
	map.plot(fig_filespec=fig_filespec, dpi=dpi)
