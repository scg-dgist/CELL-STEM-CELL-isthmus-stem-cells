'''
	Plots experimental data and simulation results for
	the ratio of labelled glands.

	Author: David Jorg
'''

import numpy as np
import pylab as plt

import csv, sys
import pandas as pd

import conventions as conv
if conv.screen_mode:
	import plotstyle_screen
else:
	import plotstyle_print

flatten = lambda l: [item for sublist in l for item in sublist]

def prob(hst):
	h = hst.astype(float)
	norm = np.sum(h)
	return np.array(h/norm) if norm != 0 else np.array(h)

hstmean_ex = lambda hst: np.sum((hst[:,0] - 0.5) * hst[:,1]) / float(conv.bins)
hstmean_sim = lambda hsttraj: np.sum(map(prob, hsttraj) * (np.arange(1,conv.bins+1) - 0.5), axis=1) / float(conv.bins)

cdf = lambda pdf: [ np.sum(pdf[:i]) for i in range(len(pdf)+1) ]

def get_quantile(cdf, q):
	for i, c in enumerate(cdf):
		if c > q:
			return float(i) + (q - cdf[i-1]) / (cdf[i] - cdf[i-1])

def main(save_average=True):
	sim_id = 'fit'

	# --------------------------------------------------------------------------
	# import experimental data

	data_ex = pd.read_csv(conv.data_path + 'labelled_gland_ratio.csv', header=None).values

	# --------------------------------------------------------------------------
	# import simulation data

	# get parameters
	params = pd.read_csv(conv.sim_path + 'queue_clonal_' + sim_id + '/parameters.csv', header=None, delim_whitespace=True).values.tolist()
	get_param = lambda name: [ e[1] for e in params if e[0] == name ][0]
	R = float(get_param('num_runs'))

	# get simulation data
	data_sim = (1. - np.transpose(pd.read_csv(conv.sim_path + 'queue_clonal_' + sim_id + '/global_unlabelled_glands.csv', header=None, delim_whitespace=True).values)[0] / R)

	# --------------------------------------------------------------------------
	# plot

	dt = float(get_param('dt'))

	fig, ax = plt.subplots(figsize=(2.5,2.5), nrows=1, ncols=1, sharey=True)

	y_sim = data_sim
	x_sim = np.arange(len(y_sim)) * ( dt / conv.timescale )
	ax.plot(x_sim, y_sim, linestyle='-', linewidth=1.5, color=conv.colors['control'])

	y_ex = data_ex[:,1]
	x_ex = data_ex[:,0]
	ax.plot(x_ex, y_ex, marker='o', linestyle='', color=conv.colors['control'])

	ax.set_xlim([-0.5, 12.5])
	ax.set_xlabel('time (months)')

	yticks = np.arange(0,4+1) * 0.25
	ax.set_yticks(yticks)
	ax.set_yticklabels(yticks*100)
	ax.set_ylim([0., 0.5])

	ax.set_ylabel('relative frequency of\nlabelled glands (%)')

	fig.set_tight_layout(True)
	fig.savefig(conv.results_path + 'labelled_gland_ratio.pdf', facecolor=fig.get_facecolor(), edgecolor='none')

if __name__=='__main__':
	main()
