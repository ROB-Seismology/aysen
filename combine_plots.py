"""
Python script to combine different figures into one.
"""

import os
import pylab
import matplotlib
import matplotlib.gridspec as gridspec
from thirdparty.recipes.autocrop import getAutoCroppedImage

from aysenlib import project_folder

base_fig_folder = os.path.join(project_folder, "Figures")


## Fault network illustration
def plot_fault_network(color=False):
	width, height = 18.4 / 2.54, 13.1 / 2.54
	num_cols, num_rows = 3, 1
	FIGSIZE = (width * num_cols, height * num_rows)
	pylab.rcParams['figure.figsize'] = FIGSIZE

	fig = pylab.figure()
	gs = gridspec.GridSpec(num_rows, num_cols)
	gs.update(wspace=0.1)

	sublabels = ['a', 'b', 'c']
	row = 0
	for col, sublabel in enumerate(sublabels):
		fig_filename = "FaultNetwork_%s.PNG" % sublabel
		fig_filespec = os.path.join(base_fig_folder, fig_filename)
		img = pylab.imread(fig_filespec)
		img = getAutoCroppedImage(img)

		ax = pylab.subplot(gs[row, col])
		ax.imshow(img)
		if sublabel == 'c':
			fig_filename = "not_allowed.png"
			fig_filespec = os.path.join(base_fig_folder, fig_filename)
			img = pylab.imread(fig_filespec)
			if not color:
				idxs = img[:,:,0] > 0
				for i in range(3):
					img[:,:,i][idxs] = 0.5
			left, right, bottom, top = pylab.axis()
			h, w = bottom - top, right - left
			x1 = left + (w-h)/2.
			x2 = h + x1
			ax.imshow(img, extent=(x1, x2, bottom, top), aspect=None)
			ax.axis((left, right, bottom, top))
		ax.set_axis_off()
		ax.set_title(sublabel + ')', loc='left', fontsize=26)

	out_filename = "FaultNetwork.png"
	if color:
		out_filename = "FaultNetwork_color.png"
	out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
	pylab.savefig(out_filespec, dpi=120)
	#pylab.show()


## 2007 intensity maps
def plot_2007_intensity_maps():
	#event_ID = "20070421"
	event_ID = "20070402"
	ipe_names = ["AllenEtAl2012", "AtkinsonWald2007", "BakunWentworth1997", "Barrientos1980"]

	height, width = 8.8 / 2.54, 11 / 2.54
	num_cols, num_rows = 2, 2
	FIGSIZE = (width * num_cols, height * num_rows)
	pylab.rcParams['figure.figsize'] = FIGSIZE

	fig = pylab.figure()
	gs = gridspec.GridSpec(num_rows, num_cols)
	gs.update(wspace=0, hspace=0.1)

	sublabels = ['a) ', 'b) ', 'c) ', 'd) ']
	for i, ipe_name in enumerate(ipe_names):
		row, col = divmod(i, 2)
		fig_filename = "%s_%s.PNG" % (event_ID, ipe_name)
		fig_filespec = os.path.join(base_fig_folder, fig_filename)
		img = pylab.imread(fig_filespec)
		img = getAutoCroppedImage(img)

		ax = pylab.subplot(gs[row, col])
		ax.imshow(img)
		ax.set_axis_off()
		ax.set_title(sublabels[i], loc='left', fontsize=13)

	out_filename = "%s_intensity_maps.png" % event_ID
	out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
	pylab.savefig(out_filespec, dpi=300)


## Event results (magnitude vs. max. probability)
def plot_event_mag_vs_maxprob():
	height, width = 9 / 2.54, 11.5 / 2.54
	num_cols, num_rows = 2, 1
	FIGSIZE = (width * num_cols, height * num_rows)
	pylab.rcParams['figure.figsize'] = FIGSIZE

	fig = pylab.figure()
	gs = gridspec.GridSpec(num_rows, num_cols)
	gs.update(wspace=0.05, hspace=0.1)

	sublabels = ['a) ', 'b) ']
	for col, event_range in enumerate(["1-6", "7-10"]):
		row = 0
		fig_filename = "events%s_M_vs_prob_BakunWentworth1997.png" % event_range
		fig_filespec = os.path.join(base_fig_folder, "Events", "v2", "NE=+0.5", fig_filename)
		img = pylab.imread(fig_filespec)
		img = getAutoCroppedImage(img)

		ax = pylab.subplot(gs[row, col])
		ax.imshow(img)
		ax.set_axis_off()
		ax.set_title(sublabels[col], loc='left', fontsize=13)

	out_filename = "all_events_M_vs_prob.png"
	out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
	pylab.savefig(out_filespec, dpi=300)


## Event results

def plot_event_results(event_set="single"):
	event_magnitudes = {"2007": 6.28,
						"SL-G": 5.72,
						"SL-F": 6.48,
						"SL-EF": 6.48,
						"SL-DE": 6.28,
						"SL-D": 5.13,
						"SL-CD": 6.28,
						"SL-C": 6.48,
						"SL-B": 5.13,
						"SL-A": 6.28}

	height, width = 13.4 / 2.54, 12.3 / 2.54

	## 2007
	if event_set == "single":
		num_cols, num_rows = 2, 1
		FIGSIZE = (width * num_cols, height * num_rows)
		pylab.rcParams['figure.figsize'] = FIGSIZE

		event = "2007"
		sublabels = ['a) ', 'b) ']
		fig = pylab.figure()
		gs = gridspec.GridSpec(num_rows, num_cols)
		gs.update(wspace=0, hspace=0)

		row = 0

		## Grid-search map
		col = 0
		fig_folder = os.path.join(base_fig_folder, "Events", "bw1997", "NE=+0.5")
		fig_filename = "%s_bw1997.png" % event
		fig_filespec = os.path.join(fig_folder, fig_filename)
		img = pylab.imread(fig_filespec)
		img = getAutoCroppedImage(img)
		ax = pylab.subplot(gs[row, col])
		ax.imshow(img)
		ax.set_axis_off()
		ax.set_title(sublabels[col], loc='left', fontsize=13)

		## Probabilistic map
		col = 1
		fig_folder = os.path.join(base_fig_folder, "Events", "v2", "NE=+0.5", event)
		mag = event_magnitudes[event]
		fig_filename = "%s_BakunWentworth1997_M=%.2f.png" % (event, mag)
		fig_filespec = os.path.join(fig_folder, fig_filename)
		img = pylab.imread(fig_filespec)
		img = getAutoCroppedImage(img)
		ax = pylab.subplot(gs[row, col])
		ax.imshow(img)
		ax.set_axis_off()
		ax.set_title(sublabels[col], loc='left', fontsize=13)

		out_filename = "%s_results.png" % event
		out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
		pylab.savefig(out_filespec, dpi=200)

	elif event_set == "all":
		## 6 columns, 3 rows
		events = ['SL-G', 'SL-F', 'SL-EF', 'SL-DE', 'SL-D', 'SL-CD', 'SL-C', 'SL-B', 'SL-A']
		num_cols, num_rows = 6, 3
		FIGSIZE = (width * num_cols, height * num_rows)
		pylab.rcParams['figure.figsize'] = FIGSIZE

		sublabels = ['a) ', 'b) ', 'c) ', 'd) ', 'e) ', 'f) ', 'g) ', 'h) ', 'i) ']

		fig = pylab.figure()
		gs = gridspec.GridSpec(num_rows, num_cols)
		gs.update(wspace=0, hspace=0.1)

		for event_nr, event in enumerate(events):

			row = event_nr // num_rows

			## Grid-search map
			col = (event_nr % 3) * 2
			print col, row
			fig_folder = os.path.join(base_fig_folder, "Events", "bw1997", "NE=+0.5")
			fig_filename = "%s_bw1997.png" % event
			fig_filespec = os.path.join(fig_folder, fig_filename)
			img = pylab.imread(fig_filespec)
			img = getAutoCroppedImage(img)
			ax = pylab.subplot(gs[row, col])
			ax.imshow(img)
			ax.set_axis_off()
			ax.set_title(sublabels[event_nr] + event, loc='right', horizontalalignment='center', fontsize=13)

			## Probabilistic map
			col = (event_nr % 3) * 2 + 1
			print col, row
			fig_folder = os.path.join(base_fig_folder, "Events", "v2", "NE=+0.5", event)
			mag = event_magnitudes[event]
			fig_filename = "%s_BakunWentworth1997_M=%.2f.png" % (event, mag)
			fig_filespec = os.path.join(fig_folder, fig_filename)
			img = pylab.imread(fig_filespec)
			img = getAutoCroppedImage(img)
			ax = pylab.subplot(gs[row, col])
			ax.imshow(img)
			ax.set_axis_off()

		out_filename = "Events_6x3.png"
		out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
		pylab.savefig(out_filespec, dpi=300)


		## 4 columns, 2 rows
		"""
		height, width = 9 / 2.54, 8 / 2.54
		#events = [['2007', 'SL-G', 'SL-F', 'SL-EF'], ['SL-D', 'SL-C', 'SL-B', 'SL-A']]
		events = [['2007', 'SL-G', 'SL-F', 'SL-DE'], ['SL-D', 'SL-C', 'SL-B', 'SL-A']]
		num_cols, num_rows = 4, 2
		FIGSIZE = (width * num_cols, height * num_rows)
		pylab.rcParams['figure.figsize'] = FIGSIZE

		sublabels = [['a) ', 'b) ', 'c) ', 'd) '], ['e) ', 'f) ', 'g) ', 'h) ']]
		for page, event_list in enumerate(events):
			fig = pylab.figure()
			gs = gridspec.GridSpec(num_rows, num_cols)
			gs.update(wspace=0, hspace=0)

			for col, event in enumerate(event_list):
				## Grid-search map
				row = 0
				fig_folder = os.path.join(base_fig_folder, "Events", "bw1997", "NE=+0.5")
				fig_filename = "%s_bw1997.png" % event
				fig_filespec = os.path.join(fig_folder, fig_filename)
				img = pylab.imread(fig_filespec)
				img = getAutoCroppedImage(img)
				ax = pylab.subplot(gs[row, col])
				ax.imshow(img)
				ax.set_axis_off()
				ax.set_title(sublabels[page][col] + event, loc='center', fontsize=13)

				## Probabilistic map
				row = 1
				fig_folder = os.path.join(base_fig_folder, "Events", "v2", "NE=+0.5", event)
				mag = event_magnitudes[event]
				fig_filename = "%s_BakunWentworth1997_M=%.2f.png" % (event, mag)
				fig_filespec = os.path.join(fig_folder, fig_filename)
				img = pylab.imread(fig_filespec)
				img = getAutoCroppedImage(img)
				ax = pylab.subplot(gs[row, col])
				ax.imshow(img)
				ax.set_axis_off()

			out_filename = "Events_4x2_%d.png" % page
			out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
			pylab.savefig(out_filespec, dpi=300)
		"""


		## 2 columns, 4 rows
		"""
		height, width = 8.5 / 2.54, 8 / 2.54
		events = [['2007', 'SL-G', 'SL-F', 'SL-EF'], ['SL-D', 'SL-C', 'SL-B', 'SL-A']]
		num_cols, num_rows = 2, 4
		FIGSIZE = (width * num_cols, height * num_rows)
		pylab.rcParams['figure.figsize'] = FIGSIZE

		sublabels = ['a) ', 'b) ', 'c) ', 'd) ']
		for page, event_list in enumerate(events):
			fig = pylab.figure()
			gs = gridspec.GridSpec(num_rows, num_cols)
			gs.update(wspace=0, hspace=0)

			for row, event in enumerate(event_list):
				## Grid-search map
				col = 0
				fig_folder = os.path.join(base_fig_folder, "Events", "bw1997", "NE=+0.5")
				fig_filename = "%s_bw1997.png" % event
				fig_filespec = os.path.join(fig_folder, fig_filename)
				img = pylab.imread(fig_filespec)
				img = getAutoCroppedImage(img)
				if row < 3:
					img = img[:-140]
				ax = pylab.subplot(gs[row, col])
				ax.imshow(img)
				ax.set_axis_off()
				ax.set_title(sublabels[row] + event, loc='center', fontsize=13)

				## Probabilistic map
				col = 1
				fig_folder = os.path.join(base_fig_folder, "Events", "v2", "NE=+0.5", event)
				mag = event_magnitudes[event]
				fig_filename = "%s_BakunWentworth1997_M=%.2f.png" % (event, mag)
				fig_filespec = os.path.join(fig_folder, fig_filename)
				img = pylab.imread(fig_filespec)
				if row < 3:
					img = img[:-140]
				img = getAutoCroppedImage(img)
				ax = pylab.subplot(gs[row, col])
				ax.imshow(img)
				ax.set_axis_off()

			out_filename = "Events_4x2_%d.png" % page
			out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
			pylab.savefig(out_filespec, dpi=200)
		"""

## Sensitivity scenarios
def plot_scenario_sensitivity_maps():
	height, width = 13.4 / 2.54, 12.3 / 2.54
	num_cols, num_rows = 2, 2
	FIGSIZE = (width * num_cols, height * num_rows)
	pylab.rcParams['figure.figsize'] = FIGSIZE

	sublabels = ['a) ', 'b) ', 'c) ', 'd) ']

	fig = pylab.figure()
	gs = gridspec.GridSpec(num_rows, num_cols)
	gs.update(wspace=0, hspace=0.1)

	#scenario = "Quitralco East"
	#mags = [6.06, 6.48, 6.68, 7.09]
	scenario = "Azul Tigre South"
	mags = [5.72, 6.28, 6.48, 7.09]
	ipe_name = "BakunWentworth1997WithSigma"
	threshold_mmi = 7.5

	for m, mag in enumerate(mags):
		row = m // num_cols
		col = (m % num_cols)
		print col, row
		fig_folder = os.path.join(base_fig_folder, "Sensitivity", "v5", "Scenarios")
		fig_filename = "%s_MMI=%.1f_%s_M=%.2f.png"
		fig_filename %= (scenario, threshold_mmi, ipe_name, mag)
		fig_filespec = os.path.join(fig_folder, fig_filename)
		img = pylab.imread(fig_filespec)
		img = getAutoCroppedImage(img)
		ax = pylab.subplot(gs[row, col])
		ax.imshow(img)
		ax.set_title(sublabels[m], loc='left', fontsize=13)
		ax.set_axis_off()

	out_filename = "Sensitivity_scenarios_maps.png"
	out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
	pylab.savefig(out_filespec, dpi=200)


def plot_scenario_mag_vs_maxprob():
	height, width = 13.9 / 2.54, 17.2 / 2.54
	num_cols, num_rows = 2, 2
	FIGSIZE = (width * num_cols, height * num_rows)
	pylab.rcParams['figure.figsize'] = FIGSIZE

	fig = pylab.figure()
	gs = gridspec.GridSpec(num_rows, num_cols)
	gs.update(wspace=0.05, hspace=0.1)

	sublabels = ['a) ', 'b) ', 'c) ', 'd) ']
	scenarios = ["Azul Tigre South", "Quitralco East", "Quitralco West", "Rio Manihuales"]
	threshold_mmi = 7.5

	for s, scenario in enumerate(scenarios):
		row = s // num_cols
		col = (s % num_cols)
		print col, row

		fig_folder = os.path.join(base_fig_folder, "Sensitivity", "v5", "Scenarios")
		fig_filename = "%s_MMI=%.1f_M_vs_prob.png" % (scenario, threshold_mmi)
		fig_filespec = os.path.join(fig_folder, fig_filename)
		img = pylab.imread(fig_filespec)
		img = getAutoCroppedImage(img)

		ax = pylab.subplot(gs[row, col])
		ax.imshow(img)
		ax.set_axis_off()
		ax.set_title(sublabels[s], loc='left', fontsize=13)

	out_filename = "Sensitivity_scenarios_M_vs_prob.png"
	out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
	pylab.savefig(out_filespec, dpi=200)


def plot_resolution_power_maps_by_mag():
	height, width = 13.4 / 2.54, 12.3 / 2.54
	num_cols, num_rows = 2, 3
	#num_cols, num_rows = 3, 2
	FIGSIZE = (width * num_cols, height * num_rows)
	pylab.rcParams['figure.figsize'] = FIGSIZE

	mags = [5.13, 5.72, 6.06, 6.48, 6.68, 7.09]
	sublabels = ['a) ', 'b) ', 'c) ', 'd) ', 'e) ', 'f) ']

	fig = pylab.figure()
	gs = gridspec.GridSpec(num_rows, num_cols)
	gs.update(wspace=0, hspace=0.1)

	ipe_name = "BakunWentworth1997WithSigma"
	threshold_mmi = 7.5

	for m, mag in enumerate(mags):
		row = m // num_cols
		col = (m % num_cols)
		print col, row
		fig_folder = os.path.join(base_fig_folder, "Sensitivity", "v5", "ResolutionPower")
		fig_filename = "ResPow_%s_MMI=%.1f_M=%.2f.png"
		fig_filename %= (ipe_name, threshold_mmi, mag)
		fig_filespec = os.path.join(fig_folder, fig_filename)
		img = pylab.imread(fig_filespec)
		img = getAutoCroppedImage(img)
		ax = pylab.subplot(gs[row, col])
		ax.imshow(img)
		ax.set_title(sublabels[m], loc='left', fontsize=13)

		ax.text(0.46, 0, "Resolving Power", transform=ax.transAxes,
			ha="center", va="bottom", size=11.5,
			bbox=dict(facecolor='w', edgecolor='none', boxstyle='square, pad=0'))

		ax.set_axis_off()

	out_filename = "Sensitivity_respow_maps_by_mag.png"
	out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
	pylab.savefig(out_filespec, dpi=200)


def plot_resolution_power_maps_by_ipe_and_lower_mmi():
	height, width = 13.4 / 2.54, 12.3 / 2.54
	num_cols, num_rows = 2, 2
	FIGSIZE = (width * num_cols, height * num_rows)
	pylab.rcParams['figure.figsize'] = FIGSIZE

	sublabels = ['a) ', 'b) ', 'c) ', 'd) ']

	fig = pylab.figure()
	gs = gridspec.GridSpec(num_rows, num_cols)
	gs.update(wspace=0, hspace=0.1)

	ipe_names = ["AllenEtAl2012", "Barrientos1980WithSigma"]
	mag = 6.48
	threshold_mmi = 7.5

	row = 0
	for i, ipe_name in enumerate(ipe_names):
		col = (i % num_cols)
		print col, row
		fig_folder = os.path.join(base_fig_folder, "Sensitivity", "v5", "ResolutionPower")
		fig_filename = "ResPow_%s_MMI=%.1f_M=%.2f.png"
		fig_filename %= (ipe_name, threshold_mmi, mag)
		fig_filespec = os.path.join(fig_folder, fig_filename)
		img = pylab.imread(fig_filespec)
		img = getAutoCroppedImage(img)
		ax = pylab.subplot(gs[row, col])
		ax.imshow(img)
		ax.set_title(sublabels[i], loc='left', fontsize=13)
		ax.text(0.46, 0, "Resolving Power", transform=ax.transAxes,
			ha="center", va="bottom", size=11.5,
			bbox=dict(facecolor='w', edgecolor='none', boxstyle='square, pad=0'))
		ax.set_axis_off()

	mags = [5.72, 6.48]
	ipe_name = "BakunWentworth1997WithSigma"
	threshold_mmi = 6.5

	row = 1
	for m, mag in enumerate(mags):
		col = (m % num_cols)
		print col, row
		fig_folder = os.path.join(base_fig_folder, "Sensitivity", "v5", "ResolutionPower")
		fig_filename = "ResPow_%s_MMI=%.1f_M=%.2f.png"
		fig_filename %= (ipe_name, threshold_mmi, mag)
		fig_filespec = os.path.join(fig_folder, fig_filename)
		img = pylab.imread(fig_filespec)
		img = getAutoCroppedImage(img)
		ax = pylab.subplot(gs[row, col])
		ax.imshow(img)
		ax.set_title(sublabels[m+2], loc='left', fontsize=13)
		ax.text(0.46, 0, "Resolving Power", transform=ax.transAxes,
			ha="center", va="bottom", size=11.5,
			bbox=dict(facecolor='w', edgecolor='none', boxstyle='square, pad=0'))
		ax.set_axis_off()

	out_filename = "Sensitivity_respow_maps_by_ipe_and_lower_mmi.png"
	out_filespec = os.path.join(base_fig_folder, "Paper", out_filename)
	pylab.savefig(out_filespec, dpi=200)



if __name__ == "__main__":
	plot_resolution_power_maps_by_mag()
