"""
Model response spectra for intraslab and megathrust events
"""

import os
import openquake.hazardlib as oqhazlib
import hazard.rshalib as rshalib

FOLDER = r"E:\Home\_kris\Publications\2018 - Chile_Intraslab-Megathrust\accelerograms Chile 2017"


base_gmpe_name = 'AtkinsonBoore2003'
#base_gmpe_name = 'ZhaoEtAl2006'
#base_gmpe_name = 'LinLee2008'
#base_gmpe_name = 'YoungsEtAl1997'

## These don't work
#base_gmpe_name = 'ZhaoEtAl2016'
#base_gmpe_name = 'SiMidorikawa1999'  ## Only for PGV

gmpe_names = {'intraslab': '%sSSlab' % base_gmpe_name,
				'megathrust': '%sSInter' % base_gmpe_name}

imt = 'SA'
intensity_unit = 'ms2'


def get_gmpe(gmpe_name, trt):
	oq_gmpe = oqhazlib.gsim.get_available_gsims()[gmpe_name]
	distance_metric = list(oq_gmpe.REQUIRES_DISTANCES)[0]
	imt_periods = {'PGA': [0.],
			'SA': [0.075, 0.1, 0.13, 0.2, 0.25, 0.33, 0.4, 0.5, 0.67, 1., 2.]}
	gmpe = rshalib.gsim.NhlibGMPE(gmpe_name, gmpe_name, distance_metric,
								4.0, 9.6, 0, 300, "MW", imt_periods=imt_periods)
	try:
		gmpe.gmpe.COEFFS = {'intraslab': getattr(gmpe.gmpe, 'COEFFS_SSLAB', None),
						'megathrust': getattr(gmpe.gmpe, 'COEFFS_SINTER', None)}[trt]
	except:
		pass
	else:
		if gmpe.gmpe.COEFFS:
			gmpe.imt_periods = gmpe._get_imt_periods_from_nhlib_coeffs(gmpe_name)
	return gmpe


## Similar to accelerometer stations
"""
mags = {'intraslab': 5.5, 'megathrust': 6.9}
distances = {'FAR1': {'intraslab': 34.4, 'megathrust': 166.4},
			'MT02': {'intraslab': 47, 'megathrust': 88},
			'MT05': {'intraslab': 22, 'megathrust': 128},
			'VA05': {'intraslab': 104, 'megathrust': 80}}
depths = {'intraslab': 88, 'megathrust': 28}

for station in sorted(distances.keys()):
	rs_list, labels = [], []
	for trt in sorted(gmpe_names.keys()):
		gmpe_name = gmpe_names[trt]
		gmpe = get_gmpe(gmpe_name, trt)
		mag = mags[trt]
		distance = distances[station][trt]
		depth = depths[trt]

		periods, intensities = gmpe.get_spectrum(mag, distance, depth, imt, include_pgm=True,
									imt_unit=intensity_unit, vs30=760, soil_type=None)
		spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
													intensity_unit)
		rs_list.append(spectrum)
		label = '%s (M=%.1f, d=%.0f km, z=%.0f km)' % (trt.title(), mag, distance, depth)
		labels.append(label)

	rsc = rshalib.result.UHSCollection(rs_list, labels=labels)
	title = "RS Comparison, station=%s, GMPE=%s" % (station, base_gmpe_name)
	fig_filename = "%s_RS_%s.PNG" % (station, base_gmpe_name)
	fig_filespec = os.path.join(FOLDER, fig_filename)
	#fig_filespec = None
	rsc.plot(plot_freq=True, intensity_unit='ms2', amax=1.0, legend_location=2,
			title=title, fig_filespec=fig_filespec)
exit()
"""

## Comparison for different distances (station FAR1), same magnitude
"""
mag = 7.0
distances = {'intraslab': [10, 50, 100, 200], 'megathrust': [166]}
depths = {'intraslab': 100, 'megathrust': 35}

rs_list, labels = [], []
for trt in sorted(gmpe_names.keys()):
	gmpe_name = gmpe_names[trt]
	gmpe = get_gmpe(gmpe_name, trt)
	depth = depths[trt]
	for distance in distances[trt]:
		periods, intensities = gmpe.get_spectrum(mag, distance, depth, imt, include_pgm=True,
								imt_unit=intensity_unit, vs30=760, soil_type=None)
		spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
													intensity_unit)
		rs_list.append(spectrum)
		label = '%s (d=%.0f km, z=%.0f km)' % (trt.title(), distance, depth)
		labels.append(label)

colors = ['r'] * len(distances['intraslab']) + ['g'] * len(distances['megathrust'])
line_styles = ['-', '--', '-.', ':'] * 2
rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors, linestyles=line_styles)
title = "RS Comparison, M=%.1f, GMPE=%s" % (mag, base_gmpe_name)
fig_filename = "FAR1_RS_M=%.1f_vs_distance.PNG" % mag
fig_filespec = os.path.join(FOLDER, fig_filename)
#fig_filespec = None
amax = {7: 8, 8: 15}.get(mag)
rsc.plot(plot_freq=True, intensity_unit='ms2', amax=amax, legend_location=2,
		title=title, fig_filespec=fig_filespec)
exit()
"""

## Maule vs Chiapas
"""
mags = {'intraslab': 8.2, 'megathrust': 8.8}
distances = {'intraslab': [10, 50, 100, 200], 'megathrust': [166]}
depths = {'intraslab': 100, 'megathrust': 35}

rs_list, labels = [], []
for trt in sorted(gmpe_names.keys()):
	gmpe_name = gmpe_names[trt]
	gmpe = get_gmpe(gmpe_name, trt)
	mag = mags[trt]
	depth = depths[trt]
	for distance in distances[trt]:
		periods, intensities = gmpe.get_spectrum(mag, distance, depth, imt, include_pgm=True,
								imt_unit=intensity_unit, vs30=760, soil_type=None)
		spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
													intensity_unit)
		rs_list.append(spectrum)
		label = '%s (M=%.1f, d=%.0f km, z=%.0f km)' % (trt.title(), mag, distance, depth)
		labels.append(label)

colors = ['r'] * len(distances['intraslab']) + ['g'] * len(distances['megathrust'])
line_styles = ['-', '--', '-.', ':'] * 2
rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors, linestyles=line_styles)
title = "RS Comparison, Maule vs Chiapas, GMPE=%s" % (base_gmpe_name,)
fig_filename = "FAR1_RS_Maule_vs_Chiapas.PNG"
fig_filespec = os.path.join(FOLDER, fig_filename)
#fig_filespec = None
amax = 15
rsc.plot(plot_freq=True, intensity_unit='ms2', amax=amax, legend_location=2,
		title=title, fig_filespec=fig_filespec)
exit()
"""

## Intraslab M=8, d=300 / megathrust M=7, d=150
## Reviewer #3: "My concern here is that a high magnitude intra-slab or even
## crustal earthquake at significant distance from a lake site may cause a
## similar response to a lower magnitude but more proximal megathrust event
## because the long distance from source to site would attenuates the
## high frequency content of the shaking."
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
	periods, intensities = gmpe.get_spectrum(mag, distance, depth, imt, include_pgm=True,
							imt_unit=intensity_unit, vs30=760, soil_type=None)
	spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
												intensity_unit)
	rs_list.append(spectrum)
	label = '%s (M=%.1f, d=%.0f km, z=%.0f km)' % (trt.title(), mag, distance, depth)
	labels.append(label)

colors = ['r', 'g']
line_styles = ['-', '-']
rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors, linestyles=line_styles)
title = "RS Comparison, Reviewer3, GMPE=%s" % (base_gmpe_name,)
fig_filename = "FAR1_RS_Reviewer3.PNG"
fig_filespec = os.path.join(FOLDER, fig_filename)
#fig_filespec = None
amax = 1.5
rsc.plot(plot_freq=True, intensity_unit='ms2', amax=amax, legend_location=2,
		title=title, fig_filespec=fig_filespec)
exit()



## Comparison for different magnitudes, same distance (FAR1)
"""
mags = {'intraslab': [6.5, 6.0, 5.5, 5.0], 'megathrust': [7.0]}
distances = {'intraslab': 34, 'megathrust': 166}
depths = {'intraslab': 100, 'megathrust': 35}

rs_list, labels = [], []
for trt in sorted(gmpe_names.keys()):
	gmpe_name = gmpe_names[trt]
	gmpe = get_gmpe(gmpe_name, trt)
	depth = depths[trt]
	distance = distances[trt]
	for mag in mags[trt]:
		periods, intensities = gmpe.get_spectrum(mag, distance, depth, imt, include_pgm=True,
								imt_unit=intensity_unit, vs30=760, soil_type=None)
		spectrum = rshalib.result.ResponseSpectrum(trt, periods, imt, intensities,
													intensity_unit)
		rs_list.append(spectrum)
		label = '%s (M=%.1f, d=%.0f km, z=%.0f km)' % (trt.title(), mag, distance, depth)
		labels.append(label)

colors = ['r'] * len(mags['intraslab']) + ['g'] * len(mags['megathrust'])
line_styles = ['-', '--', '-.', ':'][:len(mags['intraslab'])] * 2
rsc = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors, linestyles=line_styles)
title = "RS Comparison, GMPE=%s" % (base_gmpe_name,)
fig_filename = "FAR1_RS_vs_magnitude.PNG"
fig_filespec = os.path.join(FOLDER, fig_filename)
#fig_filespec = None
rsc.plot(plot_freq=True, intensity_unit='ms2', amax=3, legend_location=2,
		title=title, fig_filespec=fig_filespec)
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
		periods, intensities = gmpe.get_spectrum(mag, distance, depth, imt, include_pgm=True,
								imt_unit=intensity_unit, vs30=760, soil_type=None)
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
rsc.plot(plot_freq=True, intensity_unit='ms2', title=title, fig_filespec=fig_filespec)
"""
