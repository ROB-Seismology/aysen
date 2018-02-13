
import numpy as np

import mapping.layeredbasemap as lbm
import openquake.hazardlib as oqhazlib
import hazard.rshalib as rshalib
from hazard.rshalib.source.grid_source_model import SimpleUniformGridSourceModel
from hazard.rshalib.source_estimation import estimate_epicenter_location_and_magnitude_from_intensities
from aysenlib import (create_uniform_grid_source_model, create_point_source,
						plot_gridsearch_map, project_folder, gis_folder,
						TRT, USD, LSD, RAR, RMS)


## Construct grid source model
lonmin, lonmax, latmin, latmax = (4.0, 6.0, 49.0, 51.0)
grid_outline = (lonmin, lonmax, latmin, latmax)
grid_spacing = 0.1

min_mag, max_mag, mag_bin_width = 4.5, 7.5, 0.2
depth = 10
strike, dip, rake = 90, 60, 0
point_msr = oqhazlib.scalerel.PointMSR()

method = "reverse"

if method == "reverse":
	grd_src_model = SimpleUniformGridSourceModel(grid_outline, grid_spacing,
			min_mag, min_mag + mag_bin_width, mag_bin_width, depth,
			strike, dip, rake, point_msr, USD, LSD, RMS, RAR, TRT)

elif method == "forward":
	grd_src_model = SimpleUniformGridSourceModel(grid_outline, grid_spacing,
			min_mag, max_mag, mag_bin_width, depth,
			strike, dip, rake, point_msr, USD, LSD, RMS, RAR, TRT)

## IPE and IMT
ipe_name = "BakunWentworth1997"
imt = oqhazlib.imt.MMI()

## Draw random observation locations
num_obs = 10
lon_obs = np.random.choice(grd_src_model.lons, num_obs)
lat_obs = np.random.choice(grd_src_model.lats, num_obs)

## Compute intensities at observation locations for given epicenter and magnitude
mag, lon, lat = 5.7, 5.0, 49.5
epicenter = grd_src_model.create_point_source(lon, lat)
epicenter.mfd.min_mag = mag
epicenter.mfd.modify_set_occurrence_rates([1])
tom = oqhazlib.tom.PoissonTOM(1)
[rupture] = epicenter.iter_ruptures(tom)
pe_sites = [rshalib.site.SoilSite(lon_obs[i], lat_obs[i]) for i in range(num_obs)]
pe_site_model = rshalib.site.SoilSiteModel("Positive evidence", pe_sites)
ipe = oqhazlib.gsim.get_available_gsims()[ipe_name]()
sctx, rctx, dctx = ipe.make_contexts(pe_site_model, rupture)
pe_intensities, _ = ipe.get_mean_and_stddevs(sctx, rctx, dctx, imt, [oqhazlib.const.StdDev.TOTAL])

for i in range(num_obs):
	print lon_obs[i], lat_obs[i], pe_intensities[i]

## Add negative evidence
#lon_obs_neg, lat_obs_neg, ne_intensities = [], [], []
lon_obs_neg, lat_obs_neg, ne_intensities = [5.], [50.5], [5]
num_obs = len(ne_intensities)
ne_sites = [rshalib.site.SoilSite(lon_obs_neg[i], lat_obs_neg[i]) for i in range(num_obs)]
ne_site_model = rshalib.site.SoilSiteModel("Negative evidence", ne_sites)

## Grid search
lon_grid, lat_grid = grd_src_model.lon_grid, grd_src_model.lat_grid
print("BakunWentworth %s" % method)
(mag_grid, rms_grid) = (
	estimate_epicenter_location_and_magnitude_from_intensities(
	ipe_name, imt, grd_src_model, pe_sites, pe_intensities,
	ne_sites, ne_intensities, method=method, mag_bounds=(min_mag, max_mag)))
idx = np.unravel_index(rms_grid.argmin(), rms_grid.shape)
print mag_grid[idx], lon_grid[idx], lat_grid[idx]


map = plot_gridsearch_map(grd_src_model, mag_grid, rms_grid,
						[pe_site_model], [ne_site_model], colormap="RdYlGn_r",
						plot_rms_as_alpha=False, plot_epicenter_as="both")
map.plot()
