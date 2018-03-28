"""
Ground-motion field due to fault source
Aysen Fjord
"""

import os


## Folder locations
from aysenlib import (project_folder, gis_folder)
fig_folder = os.path.join(project_folder, "Figures")


## Event to plot
event_ID = "20070421"
#event_ID = "20070402"


if __name__ == "__main__":
	import datetime
	import numpy as np
	import openquake.hazardlib as oqhazlib
	import hazard.rshalib as rshalib
	import mapping.layeredbasemap as lbm
	from mapping.layeredbasemap.cm.norm import PiecewiseLinearNorm, LinearNorm
	import eqcatalog
	from aysenlib import get_roman_intensity


	## Common parameters
	trt = "LOFZ"
	rms = 2.5
	msr = "WC1994"
	rar = 1.0
	usd, lsd = 0, 12.5


	## 2007 earthquake(s)
	## 21/04/2007, epicenter from Legrand et al. (2011)
	if event_ID == "20070421":
		MW = 6.2
		lat, lon, depth = -45.374, -73.045, 4.
		#lat, lon, depth = -45.3293, -73.2073, 4.  ## Russo
		date = datetime.date(2007, 4, 21)
		time = datetime.time(17, 53, int(round(40.80)))
		#strike, dip, rake = 354, 88, 176
		strike, dip, rake = 20, 89, 180	## Strike from LOFZ faults

	## 02/04/2007, epicenter from Russo et al. (2011)
	elif event_ID == "20070402":
		MW = 6.1
		lat, lon, depth = -45.4472, -73.6762, 5.
		date = datetime.date(2007, 4, 2)
		time = datetime.time(2, 49, 31)
		strike, dip, rake = 53, 43, -86


	## Construct source
	eq = eqcatalog.LocalEarthquake(event_ID, date, time, lon, lat, depth, mag={'MW': MW})
	nopl = rshalib.geo.NodalPlane(strike, dip, rake)
	npd = rshalib.pmf.NodalPlaneDistribution([nopl], [1])

	pt_source = rshalib.source.PointSource.from_eq_record(eq, upper_seismogenic_depth=usd,
				lower_seismogenic_depth=lsd, nodal_plane_distribution=npd, synthetic=True,
				tectonic_region_type=trt)

	pt_src_model = rshalib.source.SourceModel("", [pt_source])


	#imt_periods = {'PGA': [0], 'SA': [0.25, 1.]}
	#period_list = sorted(np.sum(imt_periods.values()))
	imt_periods = {'MMI': [0]}
	period_list = [0]

	truncation_level = 0
	integration_distance = 1000


	## Define site model
	#grid_outline = [-74, -71, -46, -44.5]
	grid_outline = [-74, -71, -46.35, -44.85]
	#grid_spacing = (0.1, 0.1)
	grid_spacing = (1./60, 1./60)
	soil_site_model = None


	## Read observed intensities
	csv_filename = "Aysen_intensities.txt"
	csv_filespec = os.path.join(gis_folder, csv_filename)
	observation_sites = []
	observed_intensities = []
	with open(csv_filespec) as f:
		for l, line in enumerate(f):
			if l >= 1 and line[0] != '#':
				name, lat, lon, mmi1, mmi2 = line.split(',')
				lat, lon = float(lat), float(lon)
				mmi = {"20070421": mmi1, "20070402": mmi2}[event_ID]
				try:
					mmi = float(mmi)
				except:
					pass
				else:
					site = rshalib.site.SHASite(lon, lat, name=name)
					observation_sites.append(site)
					observed_intensities.append(mmi)


	# Define GMPEs or IPEs
	gmpe_names = ["AllenEtAl2012", "AtkinsonWald2007", "Barrientos1980", "BakunWentworth1997"]
	#gmpe_names = ["BakunWentworth1997"]
	for gmpe_name in gmpe_names:
		gmpe_system_def = {}
		gmpe_pmf = rshalib.pmf.GMPEPMF([gmpe_name], [1])
		gmpe_system_def[trt] = gmpe_pmf

		## Compute ground_motion field
		print("Computing ground-motion maps...")
		model_name = "Aysen Fjord"
		dsha_model = rshalib.shamodel.DSHAModel(model_name, pt_src_model, gmpe_system_def,
						grid_outline=grid_outline, grid_spacing=grid_spacing,
						soil_site_model=soil_site_model, imt_periods=imt_periods,
						truncation_level=truncation_level, integration_distance=integration_distance)

		#correlation_model = oqhazlib.correlation.JB2009CorrelationModel(vs30_clustering=True)
		#uhs_field = dsha_model.calc_gmf_fixed_epsilon_mp(num_cores=3, stddev_type="total")
		uhs_field = dsha_model.calc_gmf_fixed_epsilon_mp(num_cores=3,
							stddev_type="total", np_aggregation="avg")
		num_sites = uhs_field.num_sites

		## Plot map
		print("Plotting map")
		for T in period_list:
			#norm = None
			contour_interval = 0.5
			breakpoints = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
			norm = PiecewiseLinearNorm(breakpoints)
			#norm = LinearNorm(vmin=0, vmax=12)
			#title = "%s, Mw=%s, %s" % (date, MW, gmpe_name)
			title = ""
			hm = uhs_field.getHazardMap(period_spec=T)
			#map zonder kleuren
			#hm.export_GeoTiff("raster.tiff", num_cells = 100)
			colorbar_style = lbm.ColorbarStyle(format="%.1f", title="MMI", ticks=range(3, 10, 1))
			#colorbar_style = None
			#site_style = lbm.PointStyle(shape=".", line_color="k", size=0.5)
			site_style = None
			coastline_style = countries_style = lbm.LineStyle(line_width=0.75, line_color="dimgrey")
			show_legend = False
			map = hm.get_plot(graticule_interval=(1, 1), cmap="usgs", norm=norm,
							contour_interval=contour_interval, num_grid_cells=None,
							title=title, projection="merc", site_style=site_style,
							coastline_style=coastline_style, countries_style=countries_style,
							source_model=pt_src_model, resolution="h",
							colorbar_style=colorbar_style, show_legend=show_legend)


			## Add topographic hillshading
			idx = map.get_named_layer_index("intensity_grid")
			layer = map.layers[idx]
			#elevation_grid = lbm.GdalRasterData(r"D:\GIS-data\DEM\Etopo2.bin", region=map.region)
			elevation_grid = lbm.WCSData("http://seishaz.oma.be:8080/geoserver/wcs", "ngdc:etopo1_bedrock", region=map.region)
			blend_mode = "soft"
			hillshade_style = lbm.HillshadeStyle(0, 45, 1, blend_mode=blend_mode,
													elevation_grid=elevation_grid)
			layer.style.hillshade_style = hillshade_style
			layer.style.pixelated = True

			## Add faults
			#gis_filespec = r"C:\Users\Katleen\OneDrive\UGent\2de Master\Thesis\Global Mapper\coastline.kml"
			gis_filename = "LOFZ_breukenmodel4.TAB"
			gis_filespec = os.path.join(gis_folder, gis_filename)
			data = lbm.GisData(gis_filespec)
			style = lbm.LineStyle(line_color='purple', line_width=2)
			layer = lbm.MapLayer(data, style, legend_label="Faults")
			map.layers.insert(idx+1, layer)

			## Plot sites with observed intensities
			lons = [site.lon for site in observation_sites]
			lats = [site.lat for site in observation_sites]
			labels = ["%s (%s)" % (site.name, get_roman_intensity(mmi)) for (site, mmi) in zip(observation_sites, observed_intensities)]
			data = lbm.MultiPointData(lons, lats, labels=labels)
			label_style = lbm.TextStyle(font_size=10, horizontal_alignment="left", offset=(8,0))
			style = lbm.PointStyle(shape='s', size=7, fill_color='k', label_style=label_style)
			layer = lbm.MapLayer(data, style, legend_label="Intensity observations")
			map.layers.insert(idx+2, layer)

			## Text box
			pos = (0.965, 0.965)
			text = "%s IPE" % gmpe_name
			text_style = lbm.TextStyle(font_size=12, horizontal_alignment='right',
								vertical_alignment='top', multi_alignment='center',
								background_color='white', border_pad=0.4, border_color='k')
			map.draw_text_box(pos, text, text_style)

			#print map.map.proj4string
			#print map.get_srs().ExportToWkt()

			# colorbar aanpassen, niet weergeven = None
			#layer = map.get_layer_by_name("intensity_grid")
			#	layer.style.color_map_theme.colorbar_style = None


			#if len(gmpe_names) == 1:
			#	gmpe_name = gmpe_names[0]
			#else:
			#	gmpe_name = "AverageGMPE"

			fig_filename = "%s_%s.PNG" % (event_ID, gmpe_name)
			fig_filespec = os.path.join(fig_folder, fig_filename)
			#fig_filespec = r"C:\Temp\Aysen_test.png"
			#fig_filespec = None
			#map.export_geotiff(out_filespec=out_filespec)
			map.plot(fig_filespec=fig_filespec, dpi=200)
