"""
"""

import os
import numpy as np

import mapping.layeredbasemap as lbm
from eqgeology.faultlib import okada

from aysenlib import (read_fault_source_model_as_network, read_fault_source_model,
					gis_folder)
from model_surface_deformation import plot_okada_map


## Folders
#project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")
fig_folder = os.path.join(project_folder, "Figures", "insar")

## Read InSAR ruptures
gis_filespec = os.path.join(gis_folder, "InSAR_fault_rupture.TAB")
src_model = read_fault_source_model(gis_filespec, characteristic=False)

all_subfaults = []
for flt in src_model:
	## Override geometry (dip, depth, rake)
	flt.dip = 88
	flt.lower_seismogenic_depth = 14
	flt.rake = 176
	print flt.get_dip_direction()

	as_num = int(flt.get_length() // 2)
	ad_num = 2
	subfaults = flt.get_subfaults(as_num, ad_num, rigidity=3E+10)
	print subfaults.shape

	## Override kinematics (rake, mu)
	## Note: dip would only work if fault is not subdivided downdip
	for subflt in subfaults.flatten():
		#subflt.rake = 176
		subflt.rake = 145
		subflt.slip = 0.1

	all_subfaults.extend(subfaults.flatten())

elastic_fault = okada.create_fault(all_subfaults)
#elastic_fault.plot_3D()


## InSAR parameters
## wavelength (in m)
wavelength = 0.2360571
## Azimuth (satellite looks right)
insar_az = -11.7 + 90
## Off-nadir angle
off_nadir = 34.3

## LOS (satellite looks right)
azimuth = insar_az + 180
elevation_angle = 90 - off_nadir
#azimuth, elevation_angle = 0, 75
component = (azimuth, elevation_angle, "up")


## Read observed phase difference
grd_file = os.path.join(fig_folder, "filt_topophase.mph.vrt")
grid_data = lbm.GdalRasterData(grd_file, band_nr=2, down_sampling=4, value_conversion=np.degrees)
grid_data.apply_bbox((-73.2, -45.5, -73.0, -45.2))
grid_data = grid_data.to_mesh_grid()

## Apply mask
pg_gis_file = os.path.join(fig_folder, "insar_reliable_outline.TAB")
gis_data = lbm.GisData(pg_gis_file)
_, _, pg_data = gis_data.get_data()
grid_data.mask_polygons(pg_data, inside=False)
D_obs = grid_data.values
X, Y = grid_data.lons, grid_data.lats

## Invert slip distribution
phase_info = (wavelength, "degrees")
slip_distribution = elastic_fault.invert_slip_distribution(X, Y, D_obs,
								component, max_slip=2, phase_info=phase_info)

#slip_distribution = [ 1.21515289,  0.50608087,  0.06913836,  0.67256767,  0.,          0.82989577,
#  0.,          1.1523853,   1.12765895,  1.0780759,   0.94483067,  0.91949458,
#  0.57986581,  1.24608897,  0.83176787,  0.98021408,  0.449108,    0.96729787,
#  0.13963023,  0.80576593,  0.38963646,  0.76551757,  1.07283133,  0.61103275,
#  0.40497023,  0.34862825,  0.37856665,  0.48748454,  0.35868898,  0.40661756]


slip_file = os.path.join(fig_folder, "slip_distribution.npy")
np.save(slip_file, slip_distribution)
#slip_distribution = np.load(slip_file)
print slip_distribution
slip_distribution = [0.25] * len(elastic_fault)


for subflt, slip in zip(elastic_fault.subfaults, slip_distribution):
	subflt.slip = slip
print elastic_fault.calc_magnitude()

elastic_fault.plot_3D()


U = elastic_fault.okada(X, Y)
output = "phase"
map_region = (-73.22, -72.77, -45.56, -45.22)

fig_filespec = None
plot_okada_map(U, component, wavelength, output=output, map_region=map_region,
				fault_gis_filespec=gis_filespec, elastic_fault=elastic_fault,
				fig_filespec=fig_filespec)
