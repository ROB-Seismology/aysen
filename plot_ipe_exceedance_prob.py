import os
import numpy as np
import scipy.stats
import openquake.hazardlib as oqhazlib
import hazard.rshalib as rshalib
from hazard.rshalib.gsim import OqhazlibGMPE
import pylab

from aysenlib import (get_roman_intensity, project_folder)


ipe_names = ["AtkinsonWald2007", "AllenEtAl2012"]
#colors = ['r', 'b']
colors = ['0.25', '0.5']
M, d = 6.5, 10
#M, d = 8.5, 100
Imin = 7
intensities = np.linspace(Imin-4, Imin+4, 1001)
dI = intensities[1] - intensities[0]
fig_filename = "IPE_exceedance_prob.png"
#fig_filename = None

for ipe_name, color in zip(ipe_names, colors):
	oq_ipe = oqhazlib.gsim.get_available_gsims()[ipe_name]()
	[dist_metric] = oq_ipe.REQUIRES_DISTANCES
	ipe = OqhazlibGMPE(ipe_name, "", dist_metric, 4, 9, 0, 300, "MW")
	sctx, rctx, dctx, imt = ipe._get_contexts_and_imt(M, d, 0.1, 800, None, None, None, None, "strike-slip", "MMI", 0, 0.5)
	median, sigma = oq_ipe.get_mean_and_stddevs(sctx, rctx, dctx, imt, [oqhazlib.const.StdDev.TOTAL])
	#print sigma
	dist = scipy.stats.norm(median, sigma)
	[probs] = dist.pdf(intensities)
	probs *= dI

	pylab.plot(intensities, probs, color, lw=2, label=ipe_name)
	idx = (len(intensities) - 1) / 2

	exc_intensities = intensities[intensities > Imin]
	exc_probs = probs[intensities > Imin]
	if len(exc_probs):
		pylab.fill_between(exc_intensities, exc_probs, np.zeros_like(exc_probs),
					np.ones_like(exc_probs), color=color, alpha=0.5)
		x = exc_intensities[exc_probs>1E-4].mean()
		y = np.interp(x, exc_intensities, exc_probs)
		label = "P(I>%s)=%.3f" % (get_roman_intensity(Imin), np.sum(exc_probs)/np.sum(probs))
		pylab.annotate(label, xy=(x, y), xycoords="data", color=color,
				xytext=(15, 15), textcoords="offset points",
				arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color=color))

	xmin, xmax, ymin, ymax = pylab.axis()
	pylab.vlines(intensities[idx], 0, ymax, 'k', lw=2, linestyle='--')

	## Convert X-axis tick labels (intensities) to Roman numerals
	locs, xtick_labels = pylab.xticks()
	num_ticks = len(xtick_labels)
	xtick_values = np.linspace(xmin, xmax, num_ticks)
	for i in range(num_ticks):
		lbl = get_roman_intensity(int(xtick_values[i]))
		xtick_labels[i].set_text(lbl)
	pylab.xticks(locs, xtick_labels)

	pylab.xlabel("Intensity (MMI)")
	pylab.ylabel("Probability Density")
	pylab.legend()

if fig_filename:
	fig_filespec = os.path.join(project_folder, "Figures", "Paper", fig_filename)
	dpi = 300
	pylab.savefig(fig_filespec, dpi=dpi)
else:
	pylab.show()
