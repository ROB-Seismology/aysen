import numpy as np

scenario_prob = 0.8
scenario_mag = 6.2

cases = []

title = "Scenario"
prob = scenario_prob
mag = scenario_mag
dist = 0
cases.append((title, prob, mag, dist))

title = "Maximum"
prob = 0
mag = scenario_mag + 1
dist = 1
cases.append((title, prob, mag, dist))

title = "Different magnitude at scenario location with same probability"
prob = scenario_prob
mag = scenario_mag + 1
dist = 0
cases.append((title, prob, mag, dist))

title = "Different magnitude at scenario location with higher probability"
prob = 1
mag = scenario_mag + 1
dist = 0
cases.append((title, prob, mag, dist))

title = "Different magnitude at scenario location with lower probability"
prob = 0.2
mag = scenario_mag + 1
dist = 0
cases.append((title, prob, mag, dist))

title = "Same magnitude at different location with same probability"
prob = scenario_prob
mag = scenario_mag
dist = 1
cases.append((title, prob, mag, dist))

title = "Same magnitude at different location with higher probability"
prob = 1
mag = scenario_mag
dist = 1
cases.append((title, prob, mag, dist))

title = "Same magnitude at different location with lower probability"
prob = 0.2
mag = scenario_mag
dist = 1
cases.append((title, prob, mag, dist))

title = "Different magnitude, different location, lower probability"
prob = 0.2
mag = scenario_mag + 0.5
dist = 0.5
cases.append((title, prob, mag, dist))

title = "Very different magnitude, very different location, same low probability"
prob = 0.2
mag = scenario_mag + 1
dist = 1
cases.append((title, prob, mag, dist))

title = "Slightly different from scenario"
prob = scenario_prob - 0.1
mag = scenario_mag - 0.2
dist = 0.8
cases.append((title, prob, mag, dist))

title = "Slightly different from scenario, lower probability"
prob = scenario_prob - 0.4
mag = scenario_mag - 0.2
dist = 0.8
cases.append((title, prob, mag, dist))



"""
probs = np.array([scenario_prob, scenario_prob, 0.5, 1.])
mags = np.array([scenario_mag, scenario_mag, scenario_mag, 5.6])
distances = np.array([0, 1, 1, 0])
#probs = np.array([0.1, 0.2, 0.3, 0.2, 0.1])
probs = np.ones(5) * scenario_prob / 2
#probs = np.zeros(5)
#probs = np.ones(5)
mags = np.array([5.8, 6.0, 6.2, 6.4, 6.6])
mags = np.ones(5) * scenario_mag
distances = np.array([10., 20, 40, 75, 100]) / 100.
#distances = np.zeros(5)
#distances = np.ones(1)

#probs = np.array([0.8])
#mags = np.array([scenario_mag])
#distances = np.array([0.])
"""

for title, prob, mag, dist in cases:
	print title
	## prob_diffs: negative when prob > scenario_prob, resulting in stronger penalization
	prob_diff = scenario_prob - prob
	mag_diff = mag - scenario_mag
	print("  prob_diff=%s" % prob_diff)
	print("  mag_diff=%s" % mag_diff)
	## dx: max. = 1
	dx = np.sqrt(0.5 * ((1-mag_diff)**2 + (1-dist)**2))
	print("  dx=%s" % dx)
	## if dx is 1, res_pow varies between 0 and 1:
	## - if all probs are zero, res_pow = 1
	## - if all probs are equal to scenario_prob, res_pow = 0
	#res_pow = np.sum(prob_diffs * dx) / (scenario_prob * (len(prob_diffs)-1))
	#res_pow = np.mean(prob_diffs * dx) / scenario_prob
	res_pow = prob_diff * dx
	print("  res_pow=%s" % res_pow)
