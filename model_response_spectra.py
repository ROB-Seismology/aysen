"""
Model response spectra for intraslab and megathrust events
"""

import os
import numpy as np
import openquake.hazardlib as oqhazlib
import hazard.rshalib as rshalib

FOLDER = r"E:\Home\_kris\Publications\2018 - Chile_Intraslab-Megathrust\accelerograms Chile 2017"


#base_gmpe_name = 'AtkinsonBoore2003'
#base_gmpe_name = 'ZhaoEtAl2006'
#base_gmpe_name = 'LinLee2008'
#base_gmpe_name = 'YoungsEtAl1997'
base_gmpe_name = 'IdiniEtAl2017'

## These don't work
#base_gmpe_name = 'ZhaoEtAl2016'
#base_gmpe_name = 'SiMidorikawa1999'  ## Only for PGV

gmpe_names = {'intraslab': '%sSSlab' % base_gmpe_name,
				'megathrust': '%sSInter' % base_gmpe_name}

imt = 'SA'
intensity_unit = 'ms2'

vs30 = 761


def get_gmpe(gmpe_name, trt):
	oq_gmpe = oqhazlib.gsim.get_available_gsims()[gmpe_name]
	distance_metric = list(oq_gmpe.REQUIRES_DISTANCES)[0]
	imt_periods = {'PGA': [0.],
			'SA': [0.075, 0.1, 0.13, 0.2, 0.25, 0.33, 0.4, 0.5, 0.67, 1., 2.]}
	gmpe = rshalib.gsim.NhlibGMPE(gmpe_name, gmpe_name, distance_metric,
								4.0, 9.6, 0, 300, "MW", imt_periods=imt_periods)
	try:
		gmpe.gmpe.COEFFS = {'intraslab': getattr(gmpe.gmpe, 'COEFFS_SSLAB'),
						'megathrust': getattr(gmpe.gmpe, 'COEFFS_SINTER')}[trt]
	except:
		pass
	if gmpe.gmpe.COEFFS:
		gmpe.imt_periods = gmpe._get_imt_periods_from_nhlib_coeffs(gmpe_name)
	return gmpe


## Test AtkinsonBoore2003
"""
import scipy
import pylab
gmpe = get_gmpe('AtkinsonBoore2003SInter', 'megathrust')
M = 6.9
dmin, dmax = 20, 300
h = 28.
gmpe.plot_distance([M], dmin, dmax, h=h, imt='PGA', imt_unit='cms2', vs30=vs30,
					plot_style='loglog', amax=100)
G = 10 ** (1.2 - 0.18 * M)
rrup = np.linspace(20, 300, 10)
imt = oqhazlib.imt.PGA()
C = gmpe.gmpe.COEFFS_SINTER[imt]
pga =  10**gmpe.gmpe()._compute_mean(C, G, M, h, rrup, vs30, 0, 0, 0, 1)
#ax = pylab.gca()
pylab.loglog(rrup, pga, 'b', linewidth=1.5)

## Compare with Contreras & Boroschek (2012)
c1 = -1.8559
c2 = 0.2549
c3 = 0.0111
c4 = -0.0013
c5 = 0.3061
c6 = 0.0734
c7 = 0.3552
c8 = 1.5149
c9 = -0.103
g = c8 + c9*M
delta = c6 * 10**(c7*M)
R = np.sqrt(rrup**2 + delta**2)
pga = 10**(c1 + c2*M + c3*h + c4*R - g*np.log10(R))
pga = pga * scipy.constants.g * 100
pylab.loglog(rrup, pga, 'r', linewidth=1.5)

pylab.xlim(dmin, dmax)
pylab.grid(which='both')
pylab.show()
exit()
"""

## Similar to accelerometer stations
mags = {'intraslab': 5.5, 'megathrust': 6.9}
distances = {'FAR1': {'intraslab': 34.4, 'megathrust': 166.4},
			'MT02': {'intraslab': 47, 'megathrust': 88},
			'MT05': {'intraslab': 22, 'megathrust': 128},
			'VA05': {'intraslab': 104, 'megathrust': 80}}
depths = {'intraslab': 88, 'megathrust': 28}

#for station in sorted(distances.keys()):
for station in ('FAR1',):
	rs_list, labels = [], []
	for trt in sorted(gmpe_names.keys()):
		gmpe_name = gmpe_names[trt]
		gmpe = get_gmpe(gmpe_name, trt)
		print(gmpe_name)
		#gmpe.check_imt_unit_scaling()
		#exit()
		mag = mags[trt]
		distance = distances[station][trt]
		depth = depths[trt]
		R = np.sqrt(distance**2 + depth**2)

		periods, intensities = gmpe.get_spectrum(mag, R, depth, imt, include_pgm=True,
									imt_unit=intensity_unit, vs30=vs30, soil_type=None)
		periods, intensities = periods[2:], intensities[2:]
		spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
													intensity_unit)
		rs_list.append(spectrum)
		label = '%s (M=%.1f, d=%.0f km, z=%.0f km)' % (trt.title(), mag, distance, depth)
		labels.append(label)

		csv_filename = "%s_%s_%s.csv" % (station, trt, base_gmpe_name)
		csv_filespec = os.path.join(FOLDER, csv_filename)
		spectrum.export_csv(csv_filespec)

	colors = ['b', 'r']
	rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors)
	title = "RS Comparison, station=%s, GMPE=%s" % (station, base_gmpe_name)
	fig_filename = "%s_RS_%s.pdf" % (station, base_gmpe_name)
	#fig_filespec = os.path.join(FOLDER, fig_filename)
	fig_filespec = None
	rsc.plot(plot_freq=False, intensity_unit='ms2', amax=2.0, legend_location=1,
			legend_label_size='medium', title=title, fig_filespec=fig_filespec)
exit()


## Comparison for different distances (station FAR1), same magnitude
"""
mag = 7.0
distances = {'intraslab': [10, 50, 100, 200], 'megathrust': [170]}
depths = {'intraslab': 100, 'megathrust': 35}

rs_list, labels = [], []
for trt in sorted(gmpe_names.keys()):
	gmpe_name = gmpe_names[trt]
	gmpe = get_gmpe(gmpe_name, trt)
	depth = depths[trt]
	for distance in distances[trt]:
		R = np.sqrt(distance**2 + depth**2)
		periods, intensities = gmpe.get_spectrum(mag, R, depth, imt, include_pgm=True,
								imt_unit=intensity_unit, vs30=vs30, soil_type=None)
		periods, intensities = periods[2:], intensities[2:]
		spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
													intensity_unit)
		rs_list.append(spectrum)
		label = '%s (d=%.0f km, z=%.0f km)' % (trt.title(), distance, depth)
		labels.append(label)

colors = ['b'] * len(distances['intraslab']) + ['r'] * len(distances['megathrust'])
line_styles = ['-', '--', '-.', ':'] * 2
rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors, linestyles=line_styles)
title = "RS Comparison, M=%.1f, GMPE=%s" % (mag, base_gmpe_name)
fig_filename = "FAR1_RS_M=%.1f_vs_distance.PNG" % mag
fig_filespec = os.path.join(FOLDER, fig_filename)
#fig_filespec = None
amax = {7: 8, 8: 15}.get(mag)
rsc.plot(plot_freq=False, intensity_unit='ms2', amax=amax, legend_location=1,
		legend_label_size='medium', title=title, fig_filespec=fig_filespec)
exit()
"""

## Maule vs Tarapaca / Chiapas
"""
#intraslab_event = "Tarapaca"
intraslab_event = "Chiapas"

if intraslab_event == "Tarapaca":
	mags = {'intraslab': 7.8, 'megathrust': 8.8}
elif intraslab_event == "Chiapas":
	mags = {'intraslab': 8.2, 'megathrust': 8.8}
distances = {'intraslab': [10, 50, 100, 200], 'megathrust': [170]}
depths = {'intraslab': 100, 'megathrust': 35}

rs_list, labels = [], []
for trt in sorted(gmpe_names.keys()):
	gmpe_name = gmpe_names[trt]
	gmpe = get_gmpe(gmpe_name, trt)
	mag = mags[trt]
	depth = depths[trt]
	for distance in distances[trt]:
		R = np.sqrt(distance**2 + depth**2)
		periods, intensities = gmpe.get_spectrum(mag, R, depth, imt, include_pgm=True,
								imt_unit=intensity_unit, vs30=vs30, soil_type=None)
		periods, intensities = periods[2:], intensities[2:]
		spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
													intensity_unit)
		rs_list.append(spectrum)
		#label = '%s (M=%.1f, d=%.0f km, z=%.0f km)' % (trt.title(), mag, distance, depth)
		label = '%s (M=%.1f, d=%.0f km)' % (trt.title(), mag, distance)
		labels.append(label)

colors = ['b'] * len(distances['intraslab']) + ['r'] * len(distances['megathrust'])
line_styles = ['-', '--', '-.', ':'] * 2
rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors, linestyles=line_styles)
title = "RS Comparison, Maule vs %s, GMPE=%s" % (intraslab_event, base_gmpe_name)
fig_filename = "FAR1_RS_Maule_vs_%s.PNG" % intraslab_event
fig_filespec = os.path.join(FOLDER, fig_filename)
#fig_filespec = None
amax = 25
rsc.plot(plot_freq=False, intensity_unit='ms2', amax=amax, legend_location=1,
		legend_label_size='medium', title=title, fig_filespec=fig_filespec)
exit()
"""


## Maule vs Chiapas vs historical events
"""
maule = {'trt': 'megathrust', 'M': 8.8, 'd': [230], 'z': 30, 'label': 'Maule 2010', 'color': 'r'}
is1945 = {'trt': 'intraslab', 'M': 7.1, 'd': [35], 'z': 90, 'label': '1945 intraslab', 'color': 'b'}
mt1985 = {'trt': 'megathrust', 'M': 8.0, 'd': [170], 'z': 30, 'label': '1985 megathrust', 'color': 'darkred'}
chiapas = {'trt': 'intraslab', 'M': 8.2, 'd': [10, 50, 100, 200], 'z': 100, 'label': 'Chiapas-like', 'color': 'g'}

all_line_styles = ['-', '--', '-.', ':']
rs_list, labels, colors, line_styles = [], [], [], []
for event in [maule, is1945, mt1985, chiapas]:
	trt = event['trt']
	gmpe_name = gmpe_names[trt]
	gmpe = get_gmpe(gmpe_name, trt)
	mag = event['M']
	depth = event['z']
	distances = event['d']
	color = event['color']
	for d, distance in enumerate(distances):
		R = np.sqrt(distance**2 + depth**2)
		periods, intensities = gmpe.get_spectrum(mag, R, depth, imt, include_pgm=True,
								imt_unit=intensity_unit, vs30=vs30, soil_type=None)
		periods, intensities = periods[2:], intensities[2:]
		spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
													intensity_unit)
		rs_list.append(spectrum)
		#label = '%s (M=%.1f, d=%.0f km, z=%.0f km)' % (trt.title(), mag, distance, depth)
		label = '%s (M=%.1f, d=%.0f km)' % (event['label'], mag, distance)
		labels.append(label)
		colors.append(color)
		line_styles.append(all_line_styles[d])

rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors, linestyles=line_styles)
title = "RS Comparison, GMPE=%s" % base_gmpe_name
fig_filename = "FAR1_RS_Maule_vs_Chiapas_vs_historical.pdf"
fig_filespec = os.path.join(FOLDER, fig_filename)
#fig_filespec = None
amax = 25
rsc.plot(plot_freq=False, intensity_unit='ms2', amax=amax, legend_location=1,
		legend_label_size='medium', title=title, fig_filespec=fig_filespec)
exit()
"""


## Intraslab M=8, d=300 / megathrust M=7, d=150
## Reviewer #3: "My concern here is that a high magnitude intra-slab or even
## crustal earthquake at significant distance from a lake site may cause a
## similar response to a lower magnitude but more proximal megathrust event
## because the long distance from source to site would attenuates the
## high frequency content of the shaking."
"""
mags = {'intraslab': 8.0, 'megathrust': 7.0}
distances = {'intraslab': 300, 'megathrust': 150}
depths = {'intraslab': 100, 'megathrust': 35}

rs_list, labels = [], []
for trt in sorted(gmpe_names.keys()):
	gmpe_name = gmpe_names[trt]
	gmpe = get_gmpe(gmpe_name, trt)
	mag = mags[trt]
	depth = depths[trt]
	distance = distances[trt]
	R = np.sqrt(distance**2 + depth**2)
	periods, intensities = gmpe.get_spectrum(mag, R, depth, imt, include_pgm=True,
							imt_unit=intensity_unit, vs30=vs30, soil_type=None)
	periods, intensities = periods[2:], intensities[2:]
	spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
												intensity_unit)
	rs_list.append(spectrum)
	label = '%s (M=%.1f, d=%.0f km, z=%.0f km)' % (trt.title(), mag, distance, depth)
	labels.append(label)

colors = ['b', 'r']
line_styles = ['-', '-']
rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors, linestyles=line_styles)
title = "RS Comparison, Reviewer3, GMPE=%s" % (base_gmpe_name,)
fig_filename = "FAR1_RS_Reviewer3.PNG"
fig_filespec = os.path.join(FOLDER, fig_filename)
#fig_filespec = None
amax = 2.0
rsc.plot(plot_freq=False, intensity_unit='ms2', amax=amax, legend_location=1,
		legend_label_size='medium', title=title, fig_filespec=fig_filespec)
exit()
"""

## Comparison for different magnitudes, same distance (FAR1)
"""
#mags = {'intraslab': [6.5, 6.0, 5.5, 5.0], 'megathrust': [7.0]}
#distances = {'intraslab': 34, 'megathrust': 170}
mags = {'intraslab': [7.0], 'megathrust': [8.5, 8.0, 7.5, 7.0]}
distances = {'intraslab': 170, 'megathrust': 170}
depths = {'intraslab': 100, 'megathrust': 35}

rs_list, labels = [], []
for trt in sorted(gmpe_names.keys()):
	gmpe_name = gmpe_names[trt]
	gmpe = get_gmpe(gmpe_name, trt)
	depth = depths[trt]
	distance = distances[trt]
	R = np.sqrt(distance**2 + depth**2)
	for mag in mags[trt]:
		periods, intensities = gmpe.get_spectrum(mag, R, depth, imt, include_pgm=True,
								imt_unit=intensity_unit, vs30=vs30, soil_type=None)
		periods, intensities = periods[2:], intensities[2:]
		spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
													intensity_unit)
		rs_list.append(spectrum)
		label = '%s (M=%.1f, d=%.0f km, z=%.0f km)' % (trt.title(), mag, distance, depth)
		labels.append(label)

colors = ['b'] * len(mags['intraslab']) + ['r'] * len(mags['megathrust'])
line_styles = ['-', '--', '-.', ':']
line_styles = line_styles[:len(mags['intraslab'])] + line_styles[:len(mags['megathrust'])]

rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors, linestyles=line_styles)
title = "RS Comparison, GMPE=%s" % (base_gmpe_name,)
#fig_filename = "FAR1_RS_vs_intraslab_magnitude.PNG"
fig_filename = "FAR1_RS_vs_megathrust_magnitude.PNG"
fig_filespec = os.path.join(FOLDER, fig_filename)
#fig_filespec = None
rsc.plot(plot_freq=False, intensity_unit='ms2', amax=4, legend_location=1,
		legend_label_size='medium', title=title, fig_filespec=fig_filespec)
exit()
"""


## Comparison of different GMPEs
## AtkinsonBoore2003 is the only GMPE that peaks at different frequencies
## for intraslab and megathrust
## It also predicts the highest accelerations for intraslab
"""
mag = 7.0
distance = 100
depths = {'intraslab': 100, 'megathrust': 35}
base_gmpe_names = ['AtkinsonBoore2003', 'LinLee2008', 'YoungsEtAl1997', 'ZhaoEtAl2006']

rs_list, labels = [], []
for trt in sorted(gmpe_names.keys()):
	for base_gmpe_name in base_gmpe_names:
		gmpe_name = {'intraslab': '%sSSlab' % base_gmpe_name,
				'megathrust': '%sSInter' % base_gmpe_name}[trt]
		gmpe = get_gmpe(gmpe_name, trt)
		depth = depths[trt]
		R = np.sqrt(distance**2 + depth**2)
		periods, intensities = gmpe.get_spectrum(mag, R, depth, imt, include_pgm=True,
								imt_unit=intensity_unit, vs30=vs30, soil_type=None)
		spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
													intensity_unit)
		rs_list.append(spectrum)
		label = '%s (%s, z=%.0f km)' % (base_gmpe_name, trt.title(), depth)
		labels.append(label)

colors = ['r'] * len(base_gmpe_names) + ['g'] * len(base_gmpe_names)
line_styles = ['-', '--', '-.', ':'][:len(base_gmpe_names)] * 2
rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors, linestyles=line_styles)
title = "RS Comparison, M=%.1f, d=%.0f km" % (mag, distance)
fig_filespec = None
rsc.plot(plot_freq=False, intensity_unit='ms2', title=title, fig_filespec=fig_filespec)
"""
