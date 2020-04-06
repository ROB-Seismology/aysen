# -*- coding: iso-Latin-1 -*-

# Todo: have a look at:
# https://github.com/simondaout/Flower2d
# https://github.com/dmelgarm/MudPy
# https://github.com/scottyhq/roipy/blob/master/plot.py

# It is simply a matter of multiplying the unwrapped phase (in the natural radians units)
# by (radar wavelength)/(4 * pi) to get the apparent surface displacement in meters or whatever units you use for the radar wavelength


import os
import numpy as np
import matplotlib
import cmocean

import mapping.layeredbasemap as lbm

from aysenlib import (read_fault_source_model_as_network, read_fault_source_model,
					gis_folder)


## What to plot
output = "phase"
#output = "displacement"


## InSAR parameters
## wavelength (in m)
wavelength = 0.2360571


## Folders
#project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")
fig_folder = os.path.join(project_folder, "Figures", "insar")


## Map parameters
#map_region = (-74, -72, -46.25, -44.75)
#graticule_interval = (1, 0.5)
#map_region = (-73.5, -72.5, -45.675, -45.175)
#graticule_interval = (0.5, 0.25)
#map_region = (-73.25, -72.75, -45.575, -45.2)
map_region = (-73.22, -72.77, -45.56, -45.22)
graticule_interval = (0.2, 0.1)


layers = []

## InSAR image
img_file = os.path.join(fig_folder, "insar.jpg")
data = lbm.ImageData(img_file, 0, 0, coord_frame="display")
style = lbm.ImageStyle(horizontal_alignment='stretch', vertical_alignment='stretch', on_top=False)
layer = lbm.MapLayer(data, style)
#layers.append(layer)


## Wrapped phase
if output == "phase":
	grd_file = os.path.join(fig_folder, "filt_topophase.mph.vrt")
	grid_data = lbm.GdalRasterData(grd_file, band_nr=2, down_sampling=2)
	grid_data.apply_bbox((-73.2, -45.5, -73.0, -45.2))
	grid_data = grid_data.to_mesh_grid()

	## Apply mask
	pg_gis_file = os.path.join(fig_folder, "insar_reliable_outline.TAB")
	gis_data = lbm.GisData(pg_gis_file)
	_, _, pg_data = gis_data.get_data()
	grid_data.mask_polygons(pg_data, inside=False)

	#color_map = matplotlib.cm.jet
	#color_map_theme = lbm.ThematicStyleColormap(color_map=color_map, vmin=0, vmax=2e+6)
	#colorbar_title = "Amplitude"
	color_map = matplotlib.cm.hsv
	color_map_theme = lbm.ThematicStyleColormap(color_map=color_map, vmin=-np.pi, vmax=np.pi)
	colorbar_title = "Wrapped phase"
	contour_levels = None
	#contour_line_style = lbm.LineStyle(label_style=lbm.TextStyle())
	contour_line_style = None

## Unwrapped phase / displacement
elif output == "displacement":
	grd_file = os.path.join(fig_folder, "2007040220070703-absolute.tif")
	grid_data = lbm.GdalRasterData(grd_file, band_nr=1,
						value_conversion=lambda x: x * wavelength / (4 * np.pi))

	grid_data.apply_bbox((map_region[0], map_region[2], map_region[1], map_region[3]))
	grid_data = grid_data.interpolate_grid(X, Y)

	color_map = matplotlib.cm.jet
	color_map_theme = lbm.ThematicStyleColormap(color_map=color_map, vmin=-0.7, vmax=0.3)
	colorbar_title = "LOS Displacement (m)"
	contour_levels = np.arange(-0.6, 0.30, 0.1)
	contour_line_style = lbm.LineStyle(label_style=lbm.TextStyle())

colorbar_style = lbm.ColorbarStyle(colorbar_title, format="%.2f")
grid_style = lbm.GridStyle(color_map_theme, color_gradient="continuous",
				line_style=contour_line_style, contour_levels=contour_levels,
				colorbar_style=colorbar_style)
layer = lbm.MapLayer(grid_data, grid_style)
layers.append(layer)


## Mask oceans
#from mpl_toolkits.basemap import maskoceans
#grid_data.values = maskoceans(grid_data.lons, grid_data.lats, grid_data.values, resolution='h')


## Add topographic hillshading
#elevation_grid = lbm.GdalRasterData(r"D:\GIS-data\DEM\Etopo2.bin", region=map.region)
url = 'http://seishaz.oma.be:8080/geoserver/wcs'
layer_name, grid_resolution = "ngdc:etopo1_bedrock", 1./60
elevation_grid = lbm.WCSData(url, layer_name, resolution=grid_resolution, region=map_region)
blend_mode = "soft"
hillshade_style = lbm.HillshadeStyle(0, 45, 1, blend_mode=blend_mode,
										elevation_grid=elevation_grid)
layer.style.hillshade_style = hillshade_style
layer.style.pixelated = True


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
gis_filename = "LOFZ_breukenmodel.shp"
gis_filespec = os.path.join(gis_folder, gis_filename)
data = lbm.GisData(gis_filespec)
style = lbm.LineStyle(line_color='purple', line_width=1.5)
layer = lbm.MapLayer(data, style, legend_label="Faults")
layers.append(layer)


legend_style = None
title = ""
map = lbm.LayeredBasemap(layers, title, projection="merc", region=map_region,
		title_style=lbm.DefaultTitleTextStyle, graticule_style=lbm.GraticuleStyle(),
		graticule_interval=graticule_interval, resolution='f',
		legend_style=legend_style)


fig_filespec = None
fig_filename = "insar_%s.png" % output
#fig_filespec = os.path.join(fig_folder, fig_filename)

if fig_filespec:
	dpi = 200
else:
	dpi = 90
map.plot(fig_filespec=fig_filespec, dpi=dpi)
