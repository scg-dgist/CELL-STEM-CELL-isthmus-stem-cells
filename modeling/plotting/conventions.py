'''
	Conventions for the plotting routines.

	Author: David Jorg
'''

import numpy as np

# ------------------------------------------------------------------------------
# data storage

data_path = '../Data/'
sim_path = '../Simulations/results/'
results_path = './results/'

# ------------------------------------------------------------------------------
# reference time scale to be compared with experiments

timescale = 19.0

# ------------------------------------------------------------------------------
# experimental data specifications

bins = 8

roi = 'I2'

datasets = [ 'control', 'dmp777' ]

barrier_numbers = [ 2, 3 ]

times = {
	'control': [ 0.5, 3, 6, 12 ],
	'dmp777': [ 1, 3, 6 ]
}

time_prefix = { 0.5: '2w', 1: '1m', 3: '3m', 6: '6m', 12: '12m' }

time_name = { 0.5: '2 weeks', 1: '1 month', 3: '3 months', 6: '6 months', 12: '12 months' }

ds_prefix = { 'control': '', 'dmp777': 'DMP' }

ds_name = { 'control': 'Control', 'dmp777': 'DMP777' }

# ------------------------------------------------------------------------------
# appearance

screen_mode = False

colors1 = { 'control': '#00FFC0', 'dmp777': '#FF0080' }
colors2 = { 'control': '#008262', 'dmp777': '#7D003F' }

if screen_mode == False:
	for each in colors1.keys():
		temp = colors1[each]
		colors1[each] = colors2[each]
		colors2[each] = temp

colors = colors1 if screen_mode else colors2
dark_colors = colors2 if screen_mode else colors1

color_barrier = '#3DB4F1'
dark_color_barrier = '#215067'

def lighter(hexcolor, percent):
	value = hexcolor.lstrip('#')
	lv = len(value)
	color = tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

	color = np.array(color)
	white = np.array([255, 255, 255])
	vector = white - color

	rcolor = color + vector * percent
	return '#%02x%02x%02x' % tuple(rcolor)

# ------------------------------------------------------------------------------

def prob(hst):
	h = hst.astype(float)
	norm = np.sum(h)
	return np.array(h/norm) if norm != 0 else np.array(h)

def prob_ci(p, f, n_clones, n_samples=1000, alpha=2.5):
	c = np.cumsum(p)

	hst_resample = []
	mean_resample = []
	for _ in range(n_samples):
		sample = []
		for _ in range(n_clones):
			die = np.random.random()
			for s, each in enumerate(c):
				if die < each:
					sample.append(s+1)
					break
		hst = np.histogram(sample, bins=range(1,8+2))
		hst_resample.append(hst[0])
		mean_resample.append(np.mean(sample))

	p_resample = np.array(map(f, map(prob, hst_resample)))

	lq = np.percentile(p_resample, alpha, axis=0)
	uq = np.percentile(p_resample, 100. - alpha, axis=0)

	lq_mean = np.percentile(mean_resample, alpha)
	uq_mean = np.percentile(mean_resample, 100. - alpha)

	return lq, uq, lq_mean, uq_mean
