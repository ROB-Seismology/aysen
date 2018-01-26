

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

from rupture_probabilities import read_fault_source_model_as_network, read_fault_source_model


## Folders
#project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")
fig_folder = os.path.join(project_folder, "Figures", "insar")

fault_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel3.TAB")


def read_fault_info(fault_id):
	dM = 0.2
	for M, flt_network in read_fault_source_model_as_network(fault_filespec, dM=dM, characteristic=False):
		print("M=%.2f" % M)
		for flt in flt_network:
			if flt.source_id == fault_id:
				return flt


if __name__ == "__main__":
	#for flt in read_fault_source_model(fault_filespec):
	#	print flt.source_id, flt.name
	#exit()
	#flt_id = "2#01+2#02+2#03+2#04+2#05"
	flt_id = "0#04+0#05+0#06+0#07+0#08"
	flt = read_fault_info(flt_id)
	subfaults = flt.get_subfaults(5, 1)

	## Override geometry (dip, rake, depth
	for subflt in subfaults:
		subflt.mu = 3E+10
		#subflt.rake = -165
		subflt.dip = 85
		# Note: only works if fault is not subdivided downdip!
		subflt.calculate_geometry()

	## Map parameters
	#map_region = (-74, -72, -46.25, -44.75)
	#graticule_interval = (1, 0.5)
	#map_region = (-73.5, -72.5, -45.675, -45.175)
	#graticule_interval = (0.5, 0.25)
	map_region = (-73.25, -72.75, -45.575, -45.2)
	graticule_interval = (0.2, 0.1)

	## Compute surface deformation
	elastic_fault = okada.create_fault(subfaults)
	elastic_fault.set_slip_from_magnitude(6.3)
	elastic_fault.taper_slip()
	for s, subflt in enumerate(elastic_fault.subfaults):
		print s, subflt.slip
		#print subflt.calc_magnitude()
	print elastic_fault.calc_magnitude()

	x = np.linspace(map_region[0], map_region[1], 101)
	y = np.linspace(map_region[-2], map_region[-1], 101)
	X, Y = np.meshgrid(x, y)
	U = elastic_fault.okada(X, Y)

	## Plot
	from clawpack.visclaw import colormaps
	import cmocean
	import matplotlib

	layers = []

	component = 'Z'
	dZ = getattr(U, component)
	#azimuth = np.round(elastic_fault.subfaults[0].strike)
	azimuth = 0
	elevation_angle = 0
	#component = (azimuth, elevation_angle)
	#dZ = U.get_los_component(*component)
	dzmax = np.abs(dZ.max() - dZ.min())
	print U.U.max(), dzmax
	#contour_interval = dzmax / 10.
	#contour_levels = np.arange(contour_interval, dzmax+contour_interval, contour_interval)
	#contour_levels = list(-contour_levels[::-1]) + [0.] + list(contour_levels)
	contour_levels = None
	contour_line_style = lbm.LineStyle(label_style=lbm.TextStyle())

	grid_data = lbm.MeshGridData(X, Y, dZ)
	color_map = colormaps.blue_white_red
	color_map_theme = lbm.ThematicStyleColormap(color_map=color_map, vmin=-0.5, vmax=0.5)
	colorbar_title = "%s displacement (m)" % str(component)

	"""
	wavelength = 0.057
	#grid_data = lbm.MeshGridData(X, Y, (dZ * 4 * np.pi / wavelength) % (2 * np.pi)
	grid_data = lbm.MeshGridData(X, Y, (dZ * 4 * 180. / wavelength) % (360.))
	#color_map = cmocean.cm.phase
	color_map = matplotlib.cm.hsv
	color_map_theme = lbm.ThematicStyleColormap(color_map=color_map, vmin=0, vmax=360)
	colorbar_title = "%s phase (degrees)" % str(component)
	"""

	img_file = os.path.join(fig_folder, "insar.jpg")
	data = lbm.ImageData(img_file, 0, 0, coord_frame="display")
	style = lbm.ImageStyle(horizontal_alignment='stretch', vertical_alignment='stretch', on_top=False)
	layer = lbm.MapLayer(data, style)
	#layers.append(layer)

	colorbar_style = lbm.ColorbarStyle(colorbar_title, format="%.2f")
	#colorbar_style = None
	grid_style = lbm.GridStyle(color_map_theme, color_gradient="continuous",
					line_style=contour_line_style, contour_levels=contour_levels,
					colorbar_style=colorbar_style)
	layer = lbm.MapLayer(grid_data, grid_style)
	layers.append(layer)

	coastline_style = countries_style = lbm.LineStyle(line_width=0.75, line_color="dimgrey")
	data = lbm.BuiltinData("coastlines")
	layer = lbm.MapLayer(data, coastline_style)
	layers.append(layer)

	## Add faults
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

	fig_filespec = None

	if fig_filespec:
		dpi = 200
	else:
		dpi = 90
	map.plot(fig_filespec=fig_filespec, dpi=dpi)
