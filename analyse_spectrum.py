import os
import xlrd
import numpy as np
import obspy


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
		print sheet.name
		if station_id in sheet.name:
			comp = sheet.name[-3:]
			print comp
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


def get_fas(tr, smooth=False):
	## Perform FFT
	fas = np.fft.fft(tr.data.astype(np.complex))

	## Construct corresponding frequencies
	n = len(tr.data)
	d = 1. / tr.stats.sampling_rate
	freqs = np.fft.fftfreq(n, d)

	## Limit FFT to positive frequencies
	idxs = freqs > 0
	freqs = freqs[idxs]

	## Take absolute value of complex numbers
	fas = np.abs(fas[idxs])

	## Multiply by the sample interval to get amplitude per frequency unit
	## See https://nl.mathworks.com/matlabcentral/answers/15770-scaling-the-fft-and-the-ifft
	fas *= d

	## Smoothing
	if smooth:
		from obspy.signal.konnoohmachismoothing import konnoOhmachiSmoothing
		fas = konnoOhmachiSmoothing(fas, freqs, normalize=False)

	return (freqs, fas)


def get_rs_rvt(tr, smooth=False, damping=0.05):
	import pyrvt

	freqs, fas = get_fas(tr, smooth=smooth)
	duration = tr.stats.delta * len(tr.data)
	#duration = None
	rvt = pyrvt.motions.RvtMotion(freqs, fas, duration=duration)
	rs = rvt.compute_osc_resp(freqs, damping=damping)
	return (freqs, rs)


def get_max_peak_to_peak_height(tr):
	## Find indexes where signal changes in polarity
	diff = np.diff(tr.data)
	sign = np.sign(diff)
	idxs = np.where(sign[:-1] != sign[1:])[0] + 1

	## Find amplitude differences between successive polarity changes
	amp_diff = tr.data[idxs[1:]] - tr.data[idxs[:-1]]
	abs_amp_diff = np.abs(amp_diff)
	idx = abs_amp_diff.argmax()
	idx0 = idxs[idx]
	max_peak_height = amp_diff[idx]
	if max_peak_height > 0:
		idx1 = idxs[idx-1]
	else:
		idx1 = idxs[idx+1]

	return (np.abs(max_peak_height), sorted([idx0, idx1]))


def get_energy(tr):
	return np.sum(tr.data**2) * tr.stats.delta


def check_parseval_theorem(tr):
	## According to Parseval's theorem, the energy in the time domain
	## must equal the energy in the frequency domain
	e_time = get_energy(tr)
	freqs, fas = get_fas(tr)
	[df] = np.diff(freqs[:2])
	e_freq = np.sum(fas**2) * df * 2  # multiply by 2 to account for negative freqs
	return np.isclose(e_time, e_freq)


def download_station_metadata(station_id, network='C'):
	from obspy.fdsn.client import Client
	client = Client("IRIS")
	inv = client.get_stations(network=network, station=station_id, level="channel")
	return inv


def download_event_info(timestamp, twin=60, min_mag=5):
	from obspy.fdsn.client import Client
	client = Client("IRIS")
	t = obspy.UTCDateTime(timestamp)
	catalog = client.get_events(starttime=t-twin, endtime=t+twin, minmagnitude=min_mag)
	return catalog[0]


def calc_geodetics(eq, station):
	eq_origin = eq.preferred_origin()
	eq_lon, eq_lat = eq_origin.longitude.real, eq_origin.latitude.real
	stat_lon, stat_lat = station.longitude.real, station.latitude.real
	dist, az, backaz = obspy.core.util.geodetics.gps2DistAzimuth(eq_lat, eq_lon,
															stat_lat, stat_lon)
	return (dist, az, backaz)


def calc_distance_degrees(eq, station):
	#from obspy.geodetics import locations2degrees
	from obspy.core.util.geodetics import locations2degrees

	eq_origin = eq.preferred_origin()
	eq_lon, eq_lat = eq_origin.longitude, eq_origin.latitude
	stat_lon, stat_lat = station.longitude, station.latitude
	distance = locations2degrees(eq_lat, eq_lon, stat_lat, stat_lon)
	return distance


def calc_arrival_time_and_inclination(eq, station, phase='S*', model="ak135"):
	import warnings
	from obspy.taup import getTravelTimes
	#from obspy.taup import TauPyModel

	with warnings.catch_warnings():
		warnings.filterwarnings("ignore",category=DeprecationWarning)

		#m = TauPyModel(model="ak135")
		dist_degrees = calc_distance_degrees(eq, station)
		eq_origin = eq.preferred_origin()
		source_depth = eq_origin.depth / 1000.0
		#arrivals = m.get_ray_paths(distance_in_degree=dist_degrees,
		#	source_depth_in_km=source_depth)
		arrivals = getTravelTimes(dist_degrees, source_depth, model=model)
		#for arr in arrivals:
		#	if 'S' in arr["phase_name"] or 's' in arr["phase_name"]:
		#		print arr["phase_name"], arr["time"]

		#[arrival] = [arr for arr in arrivals if arr.name == phase]
		if phase[-1] == '*':
			phase = phase[:-1]
			fastest_arrival = None
			for arr in arrivals:
				if arr["phase_name"][:len(phase)] == phase:
					if not fastest_arrival:
						fastest_arrival = arr
					elif arr["time"] < fastest_arrival["time"]:
						fastest_arrival = arr
			arrival = fastest_arrival
		else:
			[arrival] = [arr for arr in arrivals if arr["phase_name"] == phase]
		print arrival["time"]
		#arrival_time = eq_origin.time + arrival.time
		arrival_time = eq_origin.time + arrival["time"]
		return (arrival_time, arrival["take-off angle"])



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
		freqs, fas = get_fas(tr, smooth=smooth)
		freqs, rs = get_rs_rvt(tr, smooth=smooth)

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
	phases = ["S*", "S*"]
	## Time difference between IRIS and Chilean catalog
	time_corrections = [-0.98, -2.84]

	station_ids = ["VA05", "MT02", "MT05", "FAR1"]
	#station_ids = ["MT02"]
	intensity_unit = 'ms2'

	rs_spectra = {}
	distances = {}
	overwrite = False

	for e, timestamp in enumerate(event_timestamps):
		print timestamp
		date = timestamp[:10]
		rs_spectra[timestamp] = {}
		distances[timestamp] = {}
		event_info = download_event_info(timestamp)
		#print event_info

		for station_id in station_ids:
			print station_id

			## CSV file
			folder = get_stream_folder(timestamp, station_id)
			#csv_filename = "%s-%s-%s_RS.csv" % (event, station_id, channel)
			rs_filename = "%s-%s-RS.csv" % (date, station_id)
			rs_filespec = os.path.join(folder, rs_filename)

			## Read station info
			network = get_station_network(station_id)
			inventory = download_station_metadata(station_id, network)
			station_info = inventory.networks[0].stations[0]
			#print station_info

			## Compute distance and back-azimuth
			dist, az, backaz = calc_geodetics(event_info, station_info)
			distances[timestamp][station_id] = dist

			if overwrite or not os.path.exists(rs_filespec):
				## Read Z, E, N traces
				stream = read_stream_from_txt(timestamp, station_id)
				## Attach response
				## units is available as trace.stats.response.instrument_sensitivity.input_units
				stream.attach_response(inventory)
				fig_filename = "%s_raw.png" % station_id
				fig_filespec = os.path.join(folder, fig_filename)
				stream.plot(outfile=fig_filespec, size=(1680, 1024), dpi=100, number_of_ticks=10, type="normal")

				## Isolate S-wave based on theoretical arrival time
				arrival_time, inclination = calc_arrival_time_and_inclination(event_info, station_info, phase=phases[e])
				arrival_time += time_corrections[e]
				print arrival_time, inclination
				fig_filename = "%s_S_phase.png" % station_id
				fig_filespec = os.path.join(folder, fig_filename)
				stream.plot(outfile=fig_filespec, size=(1680, 1024), dpi=100, number_of_ticks=10, type="normal", starttime=arrival_time)
				for tr in stream:
					tr.trim(arrival_time, arrival_time + event_durations[e])

				## Rotate to LQT components to separate P, SV and SH waves
				## Note: RT and LQT rotation give same results for T component
				print backaz
				#stream = stream.rotate('NE->RT', back_azimuth=backaz)
				stream = stream.rotate('ZNE->LQT', back_azimuth=backaz, inclination=inclination)
				fig_filename = "%s_S_rotated.png" % station_id
				fig_filespec = os.path.join(folder, fig_filename)
				stream.plot(outfile=fig_filespec, size=(1680, 1024), dpi=100, number_of_ticks=10, type="normal")
				stream = stream.select(component='T')

				## Compute (FAS and) RS
				tr = stream[0]
				tr.detrend()
				smooth = False
				#freqs, fas = get_fas(tr, smooth=smooth)
				freqs, rs = get_rs_rvt(tr, smooth=smooth)

				## Export RS
				rs = rshalib.result.ResponseSpectrum("", 1./freqs, "SA", rs, intensity_unit)
				rs.export_csv(rs_filespec)

			else:
				## Import RS from CSV file
				rs = rshalib.result.ResponseSpectrum.from_csv(rs_filespec, intensity_unit=intensity_unit)

			#rs.plot(plot_freq=True)
			rs_spectra[timestamp][station_id] = rs

	## Compare RS for two events
	for station_id in station_ids:
		labels = ["Intraslab (M=5.5, 2017-08-02, d=%.1f km)",
				"Megathrust (M=6.9, 2017-04-24, d=%.1f km)"]
		rs_list = []
		for e, timestamp in enumerate(event_timestamps):
			rs_list.append(rs_spectra[timestamp][station_id])
			labels[e] %= (distances[timestamp][station_id] / 1000)

		rs_coll = rshalib.result.UHSCollection(rs_list, labels=labels)
		fig_filename = "%s_RS_comparison.png" % station_id
		#fig_filespec = os.path.join(FOLDER, "accelerograms Chile 2017", fig_filename)
		fig_filespec = None
		rs_coll.plot(plot_freq=True, title="RS Comparison, station=%s" % station_id, fig_filespec=fig_filespec, amax=1.4, legend_location=2, intensity_unit=intensity_unit)

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
		fig_filename = "%s_RS_comparison.png" % date.isoformat()
		fig_folder = os.path.split(get_stream_folder(timestamp, station_id))[0]
		#fig_filespec = os.path.join(fig_folder, fig_filename)
		fig_filespec = None
		rs_coll.plot(plot_freq=True, title="RS Comparison, event=%s" % timestamp, fig_filespec=fig_filespec, amax=1.4, legend_location=2, intensity_unit=intensity_unit)
