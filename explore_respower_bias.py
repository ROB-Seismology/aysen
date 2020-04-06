"""
Explore bias in resolving power

Consider 4 faults:
- one in the center (1)
- two at the extremities (2 and 3 at d=0.5 from center)
- one intermediate between center and one extremity (4 at d=0.25 from
either of those, and at d=0.75 from the other extremity)

We compute resolving power for three cases:
1/ center fault
2/ fault at extremity close to intermediate fault
3/ fault at other extremity

|2         |4         |1                   |3

We look only at distances, magnitudes will be the same
"""

import numpy as np

def calc_resolving_power(delta_p, distances):
	N = len(distances)
	dx = np.sqrt(0.5 * ((1-distances)**2))
	res_pow = np.sum(delta_p * dx) / N
	return res_pow
	#return np.sum(delta_p * np.sqrt(((1-dists)**2) / 2))


## Normalized distances
dists1 = np.array([0.25, 0.5, 0.5])
dists2 = np.array([0.25, 0.5, 1.])
dists3 = np.array([0.5, 0.75, 1.])
dists4 = np.array([0.25, 0.25, 0.75])

## Constant delta_p
delta_p = 0.25
#print("Constant delta_p (%s)" % delta_p)
rp1 = calc_resolving_power(delta_p, dists1)
rp2 = calc_resolving_power(delta_p, dists2)
rp3 = calc_resolving_power(delta_p, dists3)
rp4 = calc_resolving_power(delta_p, dists4)
#print("  1: %.3f, 2: %.3f, 3: %.3f, 4: %.3f" % (rp1, rp2, rp3, rp4))

## Constant delta_p, renormalize dists to 1
dists1b = dists1 / dists1.max()
dists2b = dists2 / dists2.max()
dists3b = dists3 / dists3.max()
dists4b = dists4 / dists4.max()
#print("Constant delta_p, renormalize distances to 1")
rp1 = calc_resolving_power(delta_p, dists1b)
rp2 = calc_resolving_power(delta_p, dists2b)
rp3 = calc_resolving_power(delta_p, dists3b)
rp4 = calc_resolving_power(delta_p, dists4b)
#print("  1: %.3f, 2: %.3f, 3: %.3f, 4: %.3f" % (rp1, rp2, rp3, rp4))

## Constant delta_p, distances truncated at 0.5
truncation = 0.5
#print("Constant delta_p, distance truncation at %s" % truncation)
rp1 = calc_resolving_power(delta_p, np.minimum(dists1, truncation))
rp2 = calc_resolving_power(delta_p, np.minimum(dists2, truncation))
rp3 = calc_resolving_power(delta_p, np.minimum(dists3, truncation))
rp4 = calc_resolving_power(delta_p, np.minimum(dists4, truncation))
#print("  1: %.3f, 2: %.3f, 3: %.3f, 4: %.3f" % (rp1, rp2, rp3, rp4))


## Increasing delta_p with distance
print("Increasing delta_p with distance (dists/2)")
delta_p1 = dists1 / 2
delta_p2 = dists2 / 2
delta_p3 = dists3 / 2
delta_p4 = dists4 / 2
rp1 = calc_resolving_power(delta_p1, dists1)
rp2 = calc_resolving_power(delta_p2, dists2)
rp3 = calc_resolving_power(delta_p3, dists3)
rp4 = calc_resolving_power(delta_p4, dists4)
print("  1: %.3f, 2: %.3f, 3: %.3f, 4: %.3f" % (rp1, rp2, rp3, rp4))

## Increasing delta_p with distance, renormalize dists to 1
print("Increasing delta_p with distance (dists/2), renormalize distances to 1")
rp1 = calc_resolving_power(delta_p1, dists1b)
rp2 = calc_resolving_power(delta_p2, dists2b)
rp3 = calc_resolving_power(delta_p3, dists3b)
rp4 = calc_resolving_power(delta_p4, dists4b)
print("  1: %.3f, 2: %.3f, 3: %.3f, 4: %.3f" % (rp1, rp2, rp3, rp4))

## Increasing delta_p with distance, distances truncated at 0.5
print("Increasing delta_p with distance (dists/2), distance truncation at %s" % truncation)
rp1 = calc_resolving_power(delta_p1, np.minimum(dists1, truncation))
rp2 = calc_resolving_power(delta_p2, np.minimum(dists2, truncation))
rp3 = calc_resolving_power(delta_p3, np.minimum(dists3, truncation))
rp4 = calc_resolving_power(delta_p4, np.minimum(dists4, truncation))
print("  1: %.3f, 2: %.3f, 3: %.3f, 4: %.3f" % (rp1, rp2, rp3, rp4))

## delta_p divided by 2, truncation = 0.75
truncation = 0.75
print("Increasing delta_p with distance (dists/2), distance truncation at %s" % truncation)
rp1 = calc_resolving_power(delta_p1, np.minimum(dists1, truncation))
rp2 = calc_resolving_power(delta_p2, np.minimum(dists2, truncation))
rp3 = calc_resolving_power(delta_p3, np.minimum(dists3, truncation))
rp4 = calc_resolving_power(delta_p4, np.minimum(dists4, truncation))
print("  1: %.3f, 2: %.3f, 3: %.3f, 4: %.3f" % (rp1, rp2, rp3, rp4))
