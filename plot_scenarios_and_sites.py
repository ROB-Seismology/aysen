
import os
import mapping.layeredbasemap as lbm
from rupture_probabilities import *


project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
#project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")
fig_folder = os.path.join(project_folder, "Figures")


## Map parameters
#map_region = (-74, -72, -46.25, -44.75)
map_region = (-74, -72, -46.0, -44.75)
output_format = "png"

fig_filename = "Scenario_and_site_map.%s" % output_format
#fig_filespec = os.path.join(fig_folder, fig_filename)
fig_filespec = None


layers = []

## Coastlines
data = lbm.BuiltinData("coastlines")
style = lbm.LineStyle()
layer = lbm.MapLayer(data, style)
layers.append(layer)


## Faults
gis_filespec = os.path.join(gis_folder, "LOFZ_breukenmodel.shp")
data = lbm.GisData(gis_filespec)
style = lbm.LineStyle(line_color='k', line_width=1)
layer = lbm.MapLayer(data, style, legend_label='Main faults')
layers.append(layer)


## Rupture scenarios
gis_filespec = os.path.join(gis_folder, 'LOFZ_rupture_scenarios.TAB')
data = lbm.GisData(gis_filespec, label_colname='Name')
data.style_params['label_anchor'] = [0., 0., 1., None, None]
label_style = lbm.TextStyle(font_size=10, color='m', rotation='auto', vertical_alignment='center', horizontal_alignment='center', offset=(5,-5))
style = lbm.LineStyle(line_color='m', line_width=2, label_style=label_style, label_anchor=0.5)
layer = lbm.MapLayer(data, style, legend_label='Test scenarios')
layers.append(layer)


## Polygon sites
gis_filespec = os.path.join(gis_folder, 'Polygons.shp')
data = lbm.GisData(gis_filespec, label_colname='Name')
multipoint_data, multiline_data, multipolygon_data = data.get_data()
site_names, polygons = [], []
for polygon_data in multipolygon_data:
	if not polygon_data.label in site_names:
		site_names.append(polygon_data.label)
		polygons.append(polygon_data)
		if polygon_data.label == 'Cuervo opp':
			polygon_data.label = " "
data = polygons[0].to_multi_polygon()
for polygon in polygons[1:]:
	data.append(polygon)
print site_names
data.style_params['label_anchor'] = ['east', None, None, 'south', 'south', None]
data.style_params['horizontal_alignment'] = ['left', None, None, 'right', 'left', None]
data.style_params['offset'] = [(3,0), None, None, (0,0), (0,0), None]
label_style = lbm.TextStyle(font_size=10, color='b', horizontal_alignment='right', offset=(-3,0))
style = lbm.PolygonStyle(line_color=None, line_width=0, fill_color='b', label_style=label_style, label_anchor="west", alpha=0.5)
layer = lbm.MapLayer(data, style, legend_label='Polygon sites')
layers.append(layer)


## Point sites
gis_filespec = os.path.join(gis_folder, 'Points.shp')
data = lbm.GisData(gis_filespec, label_colname='Name')
multipoint_data, multiline_data, multipolygon_data = data.get_data()
site_names, points = [], []
for point_data in multipoint_data:
	if not point_data.label in site_names:
		site_names.append(point_data.label)
		points.append(point_data)
data = points[0].to_multi_point()
for point in points[1:]:
	data.append(point)
label_style = lbm.TextStyle(font_size=10, color='b', horizontal_alignment='left', offset=(7,0))
style = lbm.PointStyle(shape='D', size=6, fill_color='b', line_color='k', label_style=label_style)
layer = lbm.MapLayer(data, style, legend_label='Point sites')
layers.append(layer)


title = ""
legend_style = lbm.LegendStyle(location=4)
map = lbm.LayeredBasemap(layers, title, "merc", region=map_region,
						graticule_interval=(1, 0.5), resolution='h',
						legend_style=legend_style)

if fig_filespec:
	dpi = 200
else:
	dpi = 90
map.plot(fig_filespec=fig_filespec, dpi=dpi)
