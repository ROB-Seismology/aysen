"""
Ground-motion field due to fault source
Aysen Fjord
"""

if __name__ == "__main__":
	import os
	import datetime
	import numpy as np
	import openquake.hazardlib as oqhazlib
	import hazard.rshalib as rshalib
	import mapping.layeredbasemap as lbm
	from mapping.layeredbasemap.cm.norm import PiecewiseLinearNorm
	import eqcatalog


	## Common parameters
	trt = "LOFZ"
	rms = 2.5
	msr = "WC1994"
	rar = 1.0

	## CMT example
	ID = "Swarm 2007"
	MW = 6.2
	lat, lon, depth = -45.374, -73.045, 4
	date = datetime.date(2007, 4, 21)
	time = datetime.time(17, 53, int(round(40.80)))
	# Updated (KVN)
	#eq = eqcatalog.LocalEarthquake(ID, date, time, lon, lat, depth, MW=MW, ML=0, MS=0)
	eq = eqcatalog.LocalEarthquake(ID, date, time, lon, lat, depth, mag={'MW': MW})
	print(eq.date)

	strike, dip, rake = 354, 88, 176
	nopl = rshalib.geo.NodalPlane(strike, dip, rake)
	npd = rshalib.pmf.NodalPlaneDistribution([nopl], [1])
	usd, lsd = 0, 14
	trt = "LOFZ"

	pt_source = rshalib.source.PointSource.from_eq_record(eq, upper_seismogenic_depth=usd,
				lower_seismogenic_depth=lsd, nodal_plane_distribution=npd, synthetic=True,
				tectonic_region_type=trt)

	pt_src_model = rshalib.source.SourceModel("", [pt_source])




	# Define GMPEs
	#gmpe_names = ["AtkinsonWald2007", "AllenEtAl2012Rrup", "BakunWentworth1997"]
	gsim = oqhazlib.gsim.get_available_gsims()
	gmpe_names = ["AtkinsonWald2007"]
	for gmpe_name in gmpe_names:
		gmpe_system_def = {}
		gmpe_pmf = rshalib.pmf.GMPEPMF([gmpe_name], [1])
		gmpe_system_def[trt] = gmpe_pmf


	## Define site model
	grid_outline = [-74, -71, -46, -44.5]
	grid_spacing = (0.1, 0.1)
	site_model = rshalib.site.GenericSiteModel.from_grid_spec(grid_outline,
																grid_spacing)

	#imt_periods = {'PGA': [0], 'SA': [0.25, 1.]}
	#period_list = sorted(np.sum(imt_periods.values()))
	imt_periods = {'MMI': [0]}
	period_list = [0]

	truncation_level = 0
	integration_distance = 1000
	model_name = "Aysen Fjord"


	## Compute ground_motion field
	print("Computing ground-motion fields...")

	if len(gmpe_system_def[trt]) == 1:
		gmpe_name = gmpe_system_def[trt].gmpe_names[0]
	else:
		gmpe_name = "AverageGMPE"

	dsha_model = rshalib.shamodel.DSHAModel(model_name, pt_src_model, gmpe_system_def,
										site_model, imt_periods=imt_periods,
										truncation_level=truncation_level,
										integration_distance=integration_distance)

	#uhs_field = dsha_model.calc_gmf_fixed_epsilon_mp(num_cores=4, stddev_type="total")
	#correlation_model = oqhazlib.correlation.JB2009CorrelationModel(vs30_clustering=True)
	correlation_model = None
	uhs_field = dsha_model.calc_gmf_fixed_epsilon()
	num_sites = uhs_field.num_sites


	## Plot map
	for T in period_list:
		#contour_interval = 0.05
		#norm = None
		contour_interval = 0.5
		breakpoints = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
		norm = PiecewiseLinearNorm(breakpoints)
		title = "%s, Mw=%s" % (date, MW)
		hm = uhs_field.get_hazard_map(period_spec=T)
		#map zonder kleuren
		#hm.export_GeoTiff("raster.tiff", num_cells = 100)
		# contour_format parameter is replaced with format property of ColorbarStyle (KVN)
		# colorbar_interval can be replaced by ticks property of ColorbarStyle (KVN)
		colorbar_style = lbm.ColorbarStyle(format="%.1f", title="Intensity", ticks=breakpoints)
		#colorbar_style = None
		site_style = lbm.PointStyle(shape=".", line_color="k", size=0.5)
		map = hm.get_plot(graticule_interval=(1, 1), cmap="usgs", norm=norm,
						contour_interval=contour_interval, num_grid_cells=num_sites,
						title=title, projection="merc", site_style=site_style,
						source_model=pt_src_model, resolution="h",
						colorbar_style=colorbar_style, show_legend=False)
		# These parameters are not supported anymore (KVN)
		#, contour_format="%.1f", colorbar_interval=1, gridlabel_format="%.1f")

		# graticule_style replaces gridlabel_format in get_plot method (KVN)
		map.graticule_style.annot_format = "%.1f"

#		gis_filespec = r"C:\Users\Katleen\OneDrive\UGent\2de Master\Thesis\Global Mapper\coastline.kml"
#		data = lbm.GisData(gis_filespec)
#		line_style = lbm.LineStyle()
#		style = lbm.CompositeStyle(line_style=line_style)
#		layer = lbm.MapLayer(data, style)
#		map.layers.append(layer)

		#print map.map.proj4string
		#print map.get_srs().ExportToWkt()

		# Following line replaced by show_legend parameter in get_plot method above (KVN)
		#map.legend_style = None
		fig_filespec = None

		# Following lines can be replaced by specifying colorbar_style=None in get_plot method above (KVN)
		# colorbar aanpassen, niet weergeven = None
#		for layer in map.layers:
#		    if isinstance(layer.data, lbm.GridData):
#			layer.style.color_map_theme.colorbar_style = None


#		fig_filename = "%s_%s_T=%ss.PNG" % (src_model.name, gmpe_name, T)
#		fig_filespec = os.path.join(out_folder, fig_filename)
#		out_filespec = "%s_%s_%s_%s.tiff" % (lat, lon, depth, MW)
#		map.export_geotiff(out_filespec=out_filespec)
		map.plot(fig_filespec=fig_filespec, dpi=100)

		#exit()
