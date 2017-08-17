import os
import hazard.rshalib as rshalib
from rupture_probabilities import *


project_folder = r"C:\Users\kris\Documents\Publications\2017 - Aysen"
#project_folder = r"E:\Home\_kris\Publications\2017 - Aysen"
gis_folder = os.path.join(project_folder, "GIS")

gis_filename = "LOFZ_breukenmodel3.TAB"
gis_filespec = os.path.join(gis_folder, gis_filename)

allow_triple_junctions = False
rms = 5.7
src_model = read_fault_source_model(gis_filespec, characteristic=False)
for flt in src_model:
	flt.rupture_mesh_spacing = rms
	print flt.source_id, flt.name, flt.get_mean_strike(), flt.get_length(), flt.get_length() / rms

flt_network = src_model.get_fault_network(allow_triple_junctions=allow_triple_junctions)
flt_network.check_consistency()
print flt_network.fault_links['0#09']
#exit()
#connections = flt_network.get_connections('0#01', 25, allow_triple_junctions=allow_triple_junctions)
connections = flt_network.get_all_connections(35, allow_triple_junctions=allow_triple_junctions)
print len(connections), max([len(conn) for conn in connections])
for conn in connections:
	if len(conn) == 32:
	#if '8' in [subflt_id.split('#')[0] for subflt_id in conn]:
		print conn

		faults = src_model.get_linked_subfaults([conn])

map = faults.get_plot()
map.plot()

