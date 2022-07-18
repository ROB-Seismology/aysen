import os, sys
import numpy as np

import mapping.layeredbasemap as lbm
import openquake.hazardlib as oqhazlib
import hazard.rshalib as rshalib
from hazard.rshalib.source.grid_source_model import SimpleUniformGridSourceModel
from hazard.rshalib.source_estimation import estimate_epicenter_location_and_magnitude_from_intensities

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(SCRIPT_DIR)
from aysenlib import (create_uniform_grid_source_model, create_point_source,
						plot_gridsearch_map, project_folder, gis_folder,
						TRT, USD, LSD, RAR, RMS)


## Construct grid source model
lonmin, lonmax, latmin, latmax = (4.0, 7.0, 49.0, 51.0)
grid_outline = (lonmin, lonmax, latmin, latmax)
grid_spacing = 0.1

min_mag, max_mag, mag_bin_width = 4.5, 7.5, 0.2
depth = 10
strike, dip, rake = 90, 60, 0
point_msr = oqhazlib.scalerel.PointMSR()

method = "reverse"

grd_src_model = SimpleUniformGridSourceModel(grid_outline, grid_spacing,
		min_mag, max_mag, mag_bin_width, depth,
		strike, dip, rake, point_msr, USD, LSD, RMS, RAR, TRT)

## IPE and IMT
ipe_name = "BakunWentworth1997"
imt = oqhazlib.imt.MMI()

## Draw random observation locations
num_obs = 20
lon_obs = np.random.choice(grd_src_model.lons, num_obs)
lat_obs = np.random.choice(grd_src_model.lats, num_obs)

## Compute intensities at observation locations for given epicenter and magnitude
mag, lon, lat = 7.2, 5.0, 49.5
epicenter = grd_src_model.create_point_source(lon, lat)
epicenter.mfd.min_mag = mag
epicenter.mfd.modify_set_occurrence_rates([1])
tom = oqhazlib.tom.PoissonTOM(1)
[rupture] = epicenter.iter_ruptures(tom)

all_sites = [rshalib.site.SoilSite(lon_obs[i], lat_obs[i]) for i in range(num_obs)]
site_model = rshalib.site.SoilSiteModel(all_sites, "All sites")
ipe = rshalib.gsim.get_oq_gsim(ipe_name)
print(type(ipe))
sctx, rctx, dctx = ipe.make_contexts(site_model, rupture)
site_intensities, _ = ipe.get_mean_and_stddevs(sctx, rctx, dctx, imt, [oqhazlib.const.StdDev.TOTAL])

for i in range(num_obs):
	print(lon_obs[i], lat_obs[i], site_intensities[i])

## Categorize sites as positive or negative evidence
"""
threshold_mmi = 5.5
pe_sites, ne_sites = [], []
#pe_intensities, ne_intensities = [], []
for s, site in enumerate(all_sites):
	if site_intensities[s] >= threshold_mmi:
		pe_sites.append(site)
	else:
		ne_sites.append(site)
pe_site_model = rshalib.site.SoilSiteModel("Positive evidence", pe_sites)
ne_site_model = rshalib.site.SoilSiteModel("Negative evidence", ne_sites)
pe_intensities = [threshold] * len(pe_sites)
ne_intensities = [threshold] * len(ne_sites)
print len(pe_sites), len(ne_sites)
"""


## Use all intensities as positive evidence and add arbitrary negative evidence
pe_site_model, pe_sites, pe_intensities = site_model, all_sites, site_intensities
#lon_obs_neg, lat_obs_neg, ne_intensities = [], [], []
lon_obs_neg, lat_obs_neg, ne_intensities = [5.], [50.5], [5]
num_obs = len(ne_intensities)
ne_sites = [rshalib.site.SoilSite(lon_obs_neg[i], lat_obs_neg[i]) for i in range(num_obs)]
ne_site_model = rshalib.site.SoilSiteModel(ne_sites, "Negative evidence")


## Grid search
lon_grid, lat_grid = grd_src_model.lon_grid, grd_src_model.lat_grid
print("BakunWentworth %s" % method)
(mag_grid, rms_grid) = (
	estimate_epicenter_location_and_magnitude_from_intensities(
	ipe_name, imt, grd_src_model, pe_sites, pe_intensities,
	ne_sites, ne_intensities, method=method, ne_margin=0.25))
idx = np.unravel_index(rms_grid.argmin(), rms_grid.shape)
print(mag_grid[idx], lon_grid[idx], lat_grid[idx])


map = plot_gridsearch_map(grd_src_model, mag_grid, rms_grid,
						[pe_site_model], [ne_site_model], colormap="RdYlGn_r",
						plot_rms_as_alpha=False, plot_epicenter_as="area")
map.plot()
