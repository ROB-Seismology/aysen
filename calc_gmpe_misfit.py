"""
Ground-motion field due to fault source
Aysen Fjord
"""

import os, sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SCRIPT_DIR)


## Folder locations
#project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")
fig_folder = os.path.join(project_folder, "Figures")



import datetime
import numpy as np
import openquake.hazardlib as oqhazlib
import hazard.rshalib as rshalib
import mapping.layeredbasemap as lbm
from mapping.layeredbasemap.cm.norm import PiecewiseLinearNorm, LinearNorm
import eqcatalog
from aysenlib import roman_intensity_dict


## Common parameters
trt = "LOFZ"
rms = 2.5
msr = "WC1994"
rar = 1.0
usd, lsd = 0, 12.5

imt_periods = {'MMI': [0]}
period_list = [0]

truncation_level = 0
integration_distance = 1000


## Define site model
#grid_outline = [-74, -71, -46, -44.5]
grid_outline = [-74, -71, -46.35, -44.85]
grid_spacing = (0.1, 0.1)
soil_site_model = None


from prettytable import PrettyTable
header = ["IPE", "RMSE", "MAE", "MBE"]
tab = PrettyTable(header)
tab_rows = []


for gmpe_name in ["AllenEtAl2012", "AtkinsonWald2007", "BakunWentworth1997", "Barrientos1980"]:
	gmpe_system_def = {}
	gmpe_pmf = rshalib.pmf.GMPEPMF([gmpe_name], [1])
	gmpe_system_def[trt] = gmpe_pmf

	misfits = []
	num_observed = 0

	for event_ID in ["20070421", "20070402"]:
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


		## Compare observed with predicted intensities
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
		misfits.extend(misfit)
		num_observed += len(observed)

	misfits = np.array(misfits)
	print(misfits)
	rmse = np.sqrt(np.sum(misfits**2) / num_observed)
	## Mean Absolute Error (see https://medium.com/human-in-a-machine-world/mae-and-rmse-which-metric-is-better-e60ac3bde13d)
	mae = np.sum(np.abs(misfits)) / num_observed
	## Mean Bias Error
	mbe = np.sum(misfits) / num_observed
	tab.add_row([gmpe_name, "%.2f" % rmse, "%.2f" % mae, "%.2f" % mbe])


## Print misfit measures
print("n=%d" % num_observed)
print(tab)
