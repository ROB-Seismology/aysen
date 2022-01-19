import os
import xlrd
import numpy as np
import obspy

import seismology.obspy_tools as ot
import thirdparty.pyrotd as pyrotd


FOLDER = r"E:\Home\_kris\Publications\2018 - Chile_Intraslab-Megathrust"


def read_stream_from_excel(event, station_id):
	stream = obspy.Stream()

	# TODO: read sample_rate and timestamp from excel as well
	sample_rate = 100
	timestamp = "2017-08-02T07:15:06.809900Z"

	basename = "%s-%s" % (event, station_id)
	xl_filename = "%s.xlsx" % basename
	xl_filespec = os.path.join(FOLDER, xl_filename)
	wb = xlrd.open_workbook(xl_filespec)
	for sheet in wb.sheets():
		print(sheet.name)
		if station_id in sheet.name:
			comp = sheet.name[-3:]
			print(comp)
			#time = sheet.col(0, start_rowx=7)
			acc = sheet.col(1, start_rowx=7)
			acc = np.array([cell.value for cell in acc])
			tr = obspy.Trace(data=acc)
			tr.stats.sampling_rate = sample_rate
			tr.stats.starttime = obspy.UTCDateTime(timestamp)
			tr.stats.network = 'C'
			tr.stats.station = station_id
			tr.stats.channel = comp
			stream.append(tr)

	return stream


def get_stream_folder(timestamp, station_id):
	time = obspy.UTCDateTime(timestamp)
	date = time.date
	folder = os.path.join(FOLDER, "accelerograms Chile 2017", date.isoformat(), station_id)
	return folder


def read_stream_from_txt(timestamp, station_id):
	folder = get_stream_folder(timestamp, station_id)
	stream = obspy.Stream()
	for filename in os.listdir(folder):
		if os.path.splitext(filename)[-1] == ".txt":
			filespec = os.path.join(folder, filename)
			tr = read_trace_from_txt(filespec)
			stream.append(tr)
	return stream


def read_trace_from_txt(filespec):
	amplitudes = []
	with open(filespec) as f:
		for line in f:
			if line[0] == '#':
				if "Tiempo de Origen" in line:
					timestamp = line[line.index(':')+1:].strip()
				elif "Tasa de muestreo" in line:
					sample_rate = float(line[line.index(':')+1:].split()[0])
				elif "Estacion" in line:
					cols = line.split(':')
					station_id = cols[1].split()[0].strip()
					component = cols[2].strip()
			elif line:
				amplitudes.append(float(line))
	tr = obspy.Trace(data=np.array(amplitudes))
	tr.stats.sampling_rate = sample_rate
	tr.stats.starttime = obspy.UTCDateTime(timestamp)
	tr.stats.network = get_station_network(station_id)
	tr.stats.station = station_id
	tr.stats.channel = component

	return tr


def get_station_network(station_id):
	if station_id in ("FAR1", "LMEL"):
		network = "C"
	else:
		network = "C1"
	return network



if __name__ == "__main__":
	import pylab
	import hazard.rshalib as rshalib

	"""
	station_id = "FAR1"
	events = ["20170802-071506", "20170424-213820"]
	intensity_unit = 'ms2'

	fas_spectra, rs_spectra = [], []
	for event in events:
		stream = read_stream_from_excel(event, station_id)
		#stream.plot()

		#TODO: rotate to radial / tangential components : stream.rotate
		#TODO: theoretical arrival time of S-wave

		## Isolate S-wave on N component
		stream = stream.select(component='N')
		tr = stream[0]
		tr.detrend()
		if event == "20170802-071506":
			tr.trim(tr.stats.starttime + 31, tr.stats.starttime + 41)
		elif event == "20170424-213820":
			tr.trim(tr.stats.starttime + 57, tr.stats.starttime + 90)
		tr.plot()
		#tr.spectrogram()
		exit()

		## Compute FAS and RS
		smooth = False
		freqs, fas = ot.get_fas(tr, smooth=smooth)
		freqs, rs = ot.get_rs_rvt(tr, smooth=smooth)

		## Plot FAS vs RS
		pylab.semilogx(freqs, fas, label="FAS")
		pylab.semilogx(freqs, rs, label="RS")
		pylab.legend()
		pylab.xlabel("Frequency (Hz)")
		pylab.ylabel("Acceleration")
		#pylab.show()

		fas = rshalib.result.ResponseSpectrum("", 1./freqs, "SA", fas, intensity_unit)
		fas_spectra.append(fas)


		## Export or import to/from CSV
		#channel = tr.stats.channel
		channel = "HNN"
		csv_filename = "%s-%s-%s_RS.csv" % (event, station_id, channel)
		csv_filespec = os.path.join(FOLDER, csv_filename)

		rs = rshalib.result.ResponseSpectrum("", 1./freqs, "SA", rs, intensity_unit)
		#rs = rshalib.result.ResponseSpectrum.from_csv(csv_filespec, intensity_unit=intensity_unit)
		rs_spectra.append(rs)
		#rs.plot(plot_freq=True)
		#rs.export_csv(csv_filespec)

	## Compare RS and FAS for two events
	labels = ["Intraslab", "Megathrust"]
	#labels = events

	coll_fas = rshalib.result.UHSCollection(fas_spectra, labels=labels)
	fig_filespec = os.path.join(FOLDER, "FAS_Comparison.PNG")
	#fig_filespec = None
	#coll_fas.plot(plot_freq=True, title="FAS Comparison", fig_filespec=fig_filespec)

	coll_rs = rshalib.result.UHSCollection(rs_spectra, labels=labels)
	#fig_filespec = os.path.join(FOLDER, "RS_Comparison.PNG")
	fig_filespec = None
	coll_rs.plot(plot_freq=True, title="RS Comparison", fig_filespec=fig_filespec)
	"""

	event_timestamps = ["2017-08-02 07:15:13", "2017-04-24 21:38:28"]
	event_durations = [30, 50]
	#phases = ["Sn", "Sb"]
	#phases = ["S*", "S*"]
	phases = ["s", "S*"]
	## Time difference between IRIS and Chilean catalog
	time_corrections = [-0.98, -2.84]

	#station_ids = ["VA05", "MT02", "MT05", "FAR1"]
	station_ids = ["FAR1"]
	intensity_unit = 'ms2'

	rs_spectra = {}
	distances = {}
	depths = {}
	overwrite = True

	for e, timestamp in enumerate(event_timestamps):
		print(timestamp)
		date = timestamp[:10]
		rs_spectra[timestamp] = {}
		distances[timestamp] = {}
		depths[timestamp] = {}
		event_info = ot.download_event_info(timestamp)
		#print(event_info)
		depths[timestamp] = event_info.preferred_origin().depth / 1000.0

		for station_id in station_ids:
			print(station_id)

			## CSV file
			folder = get_stream_folder(timestamp, station_id)
			#csv_filename = "%s-%s-%s_RS.csv" % (event, station_id, channel)
			rs_filename = "%s-%s-RS.csv" % (date.replace('-', ''), station_id)
			rs_filespec = os.path.join(folder, rs_filename)

			## Read station info
			network = get_station_network(station_id)
			inventory = ot.download_station_metadata(station_id, network)
			station_info = inventory.networks[0].stations[0]
			#print station_info

			## Compute (horizontal) distance and back-azimuth
			dist, az, backaz = ot.calc_geodetics(event_info, station_info)
			distances[timestamp][station_id] = dist

			if overwrite or not os.path.exists(rs_filespec):
				## Read Z, E, N traces
				stream = read_stream_from_txt(timestamp, station_id)
				## Attach response
				## units is available as trace.stats.response.instrument_sensitivity.input_units
				stream.attach_response(inventory)
				fig_filename = "%s_raw.png" % station_id
				#fig_filespec = os.path.join(folder, fig_filename)
				fig_filespec = None
				#stream.plot(outfile=fig_filespec, size=(1680, 1024), dpi=100, number_of_ticks=10, type="normal")

				## Isolate S-wave based on theoretical arrival time
				arrival_time, inclination = ot.calc_arrival_time_and_inclination(event_info, station_info, phase=phases[e], verbose=True)
				arrival_time += time_corrections[e]
				print(arrival_time, inclination)
				fig_filename = "%s_S_phase.png" % station_id
				#fig_filespec = os.path.join(folder, fig_filename)
				fig_filespec = None
				#stream.plot(outfile=fig_filespec, size=(1680, 1024), dpi=100, number_of_ticks=10, type="normal", starttime=arrival_time)
				for tr in stream:
					tr.trim(arrival_time, arrival_time + event_durations[e])

				## Rotate to LQT components to separate P, SV and SH waves
				## Note: RT and LQT rotation give same results for T component
				print(backaz)
				#stream = stream.rotate('NE->RT', back_azimuth=backaz)
				stream = stream.rotate('ZNE->LQT', back_azimuth=backaz, inclination=inclination)
				fig_filename = "%s_S_rotated.png" % station_id
				#fig_filespec = os.path.join(folder, fig_filename)
				fig_filespec = None
				#stream.plot(outfile=fig_filespec, size=(1680, 1024), dpi=100, number_of_ticks=10, type="normal")
				stream = stream.select(component='T')

				## Compute (FAS and) RS
				tr = stream[0]
				tr.detrend()
				print("PGA: %.2f" % ot.get_pgm(tr))
				smooth = False
				#freqs, fas = ot.get_fas(tr, smooth=smooth)
				#freqs, rs = ot.get_rs_rvt(tr, smooth=smooth)
				#tr2 = stream[1]
				#tr2.detrend()
				#freqs, rs = ot.get_geometric_mean_rs(tr, tr2)
				freqs, rs = ot.get_rs(tr)
				#freqs = freqs[::5]
				#rs = pyrotd.calc_spec_accels(1. / tr.stats.sampling_rate, tr.data, freqs).spec_accel

				## Remove zero frequency
				freqs = freqs[1:]
				rs = rs[1:]

				## Export RS
				periods = 1./freqs
				rs = rshalib.result.ResponseSpectrum("", periods, "SA", rs, intensity_unit)
				#rs.export_csv(rs_filespec)

			else:
				## Import RS from CSV file
				rs = rshalib.result.ResponseSpectrum.from_csv(rs_filespec, intensity_unit=intensity_unit)

			#rs.plot(plot_freq=True)
			rs_spectra[timestamp][station_id] = rs

	## Compare RS for two events
	for station_id in station_ids:
		labels = ["Intraslab (M=5.5, d=%.0f km, z=%.0f km)",
				"Megathrust (M=6.9, d=%.0f km, z=%.0f km)"]
		colors = ['b', 'r']
		rs_list = []
		for e, timestamp in enumerate(event_timestamps):
			rs_list.append(rs_spectra[timestamp][station_id])
			labels[e] %= (distances[timestamp][station_id] / 1000, depths[timestamp])

		rs_coll = rshalib.result.UHSCollection(rs_list, labels=labels, colors=colors)
		fig_filename = "%s_RS_comparison.png" % station_id
		fig_filespec = os.path.join(FOLDER, "accelerograms Chile 2017", fig_filename)
		#fig_filespec = None
		title = "RS Comparison, station=%s" % station_id
		rs_coll.plot(plot_freq=False, title=title, fig_filespec=fig_filespec,
					amax=2.5, legend_location=1, intensity_unit=intensity_unit,
					Tmin=0.02, Tmax=10)

	exit()

	## Compare RS for different stations
	for timestamp in event_timestamps:
		labels = []
		rs_list = []
		for station_id in station_ids:
			rs_list.append(rs_spectra[timestamp][station_id])
			d = distances[timestamp][station_id] / 1000
			labels.append("%s (d=%.1f km)" % (station_id, d))

		rs_coll = rshalib.result.UHSCollection(rs_list, labels=labels)
		date = obspy.UTCDateTime(timestamp).date
		fig_filename = "%s_RSmean_comparison.png" % date.isoformat()
		fig_folder = os.path.split(get_stream_folder(timestamp, station_id))[0]
		fig_filespec = os.path.join(fig_folder, fig_filename)
		#fig_filespec = None
		rs_coll.plot(plot_freq=False, title="RS Comparison, event=%s" % timestamp, fig_filespec=fig_filespec, amax=2.5, legend_location=2, intensity_unit=intensity_unit)
