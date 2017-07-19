"""
Ground-motion field due to fault source
Aysen Fjord
"""


## Folder locations
#project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")
fig_folder = os.path.join(project_folder, "Figures")


## Event to plot
event_ID = "20070421"

roman_intensity_dict = {1: 'I', 2: 'II', 3: 'III', 4: 'IV', 5: 'V', 6: 'VI',
						7: 'VII', 8: 'VIII', 9: 'IX', 10: 'X', 11: 'XI', 12: 'XII'}


if __name__ == "__main__":
	import os
	import datetime
	import numpy as np
	import openquake.hazardlib as oqhazlib
	import hazard.rshalib as rshalib
	import mapping.layeredbasemap as lbm
	from mapping.layeredbasemap.cm.norm import PiecewiseLinearNorm, LinearNorm
	import eqcatalog


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
		strike, dip, rake = 20, 90, 180	## Strike from LOFZ faults

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


	# Define GMPEs or IPEs
	#gmpe_names = ["AtkinsonWald2007", "Barrientos2007", "BakunWentworth1997"]
	gmpe_names = ["BakunWentworth1997"]
	for gmpe_name in gmpe_names:
		gmpe_system_def = {}
		gmpe_pmf = rshalib.pmf.GMPEPMF([gmpe_name], [1])
		gmpe_system_def[trt] = gmpe_pmf

	#imt_periods = {'PGA': [0], 'SA': [0.25, 1.]}
	#period_list = sorted(np.sum(imt_periods.values()))
	imt_periods = {'MMI': [0]}
	period_list = [0]

	truncation_level = 0
	integration_distance = 1000


	## Define site model
	#grid_outline = [-74, -71, -46, -44.5]
	grid_outline = [-74, -71, -46.35, -44.85]
	grid_spacing = (0.1, 0.1)
	soil_site_model = None


	## Read observed intensities
	csv_filename = "Aysen_intensities.txt"
	csv_filespec = os.path.join(gis_folder, csv_filename)
	observation_sites = []
	observed_intensities = []
	with open(csv_filespec) as f:
		for l, line in enumerate(f):
			if l >= 1:
				name, lat, lon, mmi1, mmi2 = line.split(',')
				lat, lon = float(lat), float(lon)
				mmi = {"20070421": mmi1, "20070402": mmi2}[event_ID]
				try:
					mmi = int(mmi)
				except:
					pass
				else:
					site = rshalib.site.SHASite(lon, lat, name=name)
					observation_sites.append(site)
					observed_intensities.append(mmi)


	## Compare observed with predicted intensities
	from prettytable import PrettyTable
	header = ["IPE", "RMSE", "MAE", "MBE"]
	tab = PrettyTable(header)
	tab_rows = []
	for gmpe_name in ["AllenEtAl2012", "AtkinsonWald2007", "BakunWentworth1997", "Barrientos2007"]:
		gmpe_system_def = {}
		gmpe_pmf = rshalib.pmf.GMPEPMF([gmpe_name], [1])
		gmpe_system_def[trt] = gmpe_pmf

		model_name = "Aysen Fjord"
		dsha_model = rshalib.shamodel.DSHAModel(model_name, pt_src_model, gmpe_system_def,
						grid_outline=[], grid_spacing=[], soil_site_model=None,
						sites=observation_sites, imt_periods=imt_periods,
						truncation_level=truncation_level, integration_distance=integration_distance)

		uhs_field = dsha_model.calc_gmf_fixed_epsilon()
		predicted = uhs_field.intensities[:,0]
		observed = np.array(observed_intensities)
		misfit = predicted - observed
		#print misfit
		rmse = np.sqrt(np.sum(misfit**2) / len(observed))
		## Mean Absolute Error (see https://medium.com/human-in-a-machine-world/mae-and-rmse-which-metric-is-better-e60ac3bde13d)
		mae = np.sum(np.abs(misfit)) / len(observed)
		## Mean Bias Error
		mbe = np.sum(misfit) / len(observed)
		tab.add_row([gmpe_name, "%.2f" % rmse, "%.2f" % mae, "%.2f" % mbe])


		## Compute ground_motion field
		print("Computing ground-motion maps...")
		dsha_model = rshalib.shamodel.DSHAModel(model_name, pt_src_model, gmpe_system_def,
						grid_outline=grid_outline, grid_spacing=grid_spacing,
						soil_site_model=soil_site_model, imt_periods=imt_periods,
						truncation_level=truncation_level, integration_distance=integration_distance)

		#correlation_model = oqhazlib.correlation.JB2009CorrelationModel(vs30_clustering=True)
		#uhs_field = dsha_model.calc_gmf_fixed_epsilon_mp(num_cores=4, stddev_type="total")
		uhs_field = dsha_model.calc_gmf_fixed_epsilon()
		num_sites = uhs_field.num_sites


		## Plot map
		for T in period_list:
			#norm = None
			contour_interval = 0.5
			breakpoints = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
			norm = PiecewiseLinearNorm(breakpoints)
			#norm = LinearNorm(vmin=0, vmax=12)
			title = "%s, Mw=%s, %s" % (date, MW, gmpe_name)
			hm = uhs_field.getHazardMap(period_spec=T)
			#map zonder kleuren
			#hm.export_GeoTiff("raster.tiff", num_cells = 100)
			#site_style = lbm.PointStyle(shape=".", line_color="k", size=0.5)
			site_style = None
			coastline_style = countries_style = lbm.LineStyle(line_width=2, line_color="w")
			show_legend = False
			map = hm.get_plot(graticule_interval=(1, 1), cmap="usgs", norm=norm,
							contour_interval=contour_interval, num_grid_cells=num_sites,
							title=title, projection="merc", site_style=site_style,
							coastline_style=coastline_style, countries_style=countries_style,
							source_model=pt_src_model, resolution="h", show_legend=show_legend)
			#, contour_format="%.1f", colorbar_interval=1, gridlabel_format="%.1f"

			## Add faults
			#gis_filespec = r"C:\Users\Katleen\OneDrive\UGent\2de Master\Thesis\Global Mapper\coastline.kml"
			gis_filename = "LOFZ_breukenmodel.shp"
			gis_filespec = os.path.join(gis_folder, gis_filename)
			data = lbm.GisData(gis_filespec)
			style = lbm.LineStyle(line_color='purple', line_width=2)
			layer = lbm.MapLayer(data, style, legend_label="Faults")
			map.layers.append(layer)

			## Plot sites with observed intensities
			lons = [site.lon for site in observation_sites]
			lats = [site.lat for site in observation_sites]
			labels = ["%s (%s)" % (site.name, roman_intensity_dict[mmi]) for (site, mmi) in zip(observation_sites, observed_intensities)]
			data = lbm.MultiPointData(lons, lats, labels=labels)
			label_style = lbm.TextStyle(font_size=10, horizontal_alignment="left", offset=(8,0))
			style = lbm.PointStyle(shape='s', size=7, fill_color='k', label_style=label_style)
			layer = lbm.MapLayer(data, style, legend_label="Intensity observations")
			map.layers.append(layer)

			#print map.map.proj4string
			#print map.get_srs().ExportToWkt()

			# colorbar aanpassen, niet weergeven = None
			#for layer in map.layers:
			#    if isinstance(layer.data, lbm.GridData):
			#	layer.style.color_map_theme.colorbar_style = None


			#if len(gmpe_names) == 1:
			#	gmpe_name = gmpe_names[0]
			#else:
			#	gmpe_name = "AverageGMPE"

			fig_filename = "%s_%s.PNG" % (event_ID, gmpe_name)
			fig_filespec = os.path.join(fig_folder, fig_filename)
			#fig_filespec = None
			#map.export_geotiff(out_filespec=out_filespec)
			map.plot(fig_filespec=fig_filespec, dpi=200)
			#exit()

	## Print misfit measures
	print tab
