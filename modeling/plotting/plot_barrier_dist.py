'''
	Plots experimental data and simulation results for
	the distribution of barrier/parietal cell numbers and interbarrier angles.

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

def stdprop(hst):
	p = conv.prob(hst)
	return np.sqrt(p * (1. - p) / np.sum(hst))

def main():
	sim_id = 'fit'

	# --------------------------------------------------------------------------
	# import experimental data

	number_dist_ex = pd.read_csv(conv.data_path + 'barrier_number_dist.csv', header=None).values

	angle_dist_ex = { n: np.transpose(pd.read_csv(conv.data_path + 'barrier_angle_dist_' + str(n) + '.csv', header=None).values)[0] for n in conv.barrier_numbers }

	# --------------------------------------------------------------------------
	# import simulation data

	# get simulation data
	number_dist_sim = np.sum(pd.read_csv(conv.sim_path + 'queue_clonal_' + sim_id + '/global_barrier_number_dist.csv', header=None, delim_whitespace=True).values[:,1:], axis=0)

	angle_dist_sim = {
		n: np.sum(pd.read_csv(conv.sim_path + 'queue_clonal_' + sim_id + '/global_barrier_angle_dist_' + str(n).zfill(2) + '.csv', header=None, delim_whitespace=True).values[:,1:], axis=0) for n in conv.barrier_numbers
	}

	# --------------------------------------------------------------------------
	# plot

	fig, axs = plt.subplots(figsize=(6,2.7), nrows=1, ncols=len(conv.barrier_numbers)+1, sharey=True)

	bw = 0.8

	h_ex = number_dist_ex[:,1]
	n_ex = int(np.sum(h_ex))

	y_ex = conv.prob(number_dist_ex[:,1])
	x_ex = number_dist_ex[:,0] + 0.5 - 0.5 * bw
	# get standard error of the proportion
	e_ex = stdprop(number_dist_ex[:,1])
	axs[0].bar(x_ex, y_ex, width=bw, color=conv.dark_color_barrier, ecolor=conv.color_barrier)

	y_sim = conv.prob(number_dist_sim)
	redef = lambda x: x
	lq, uq, _, _ = conv.prob_ci(y_sim, redef, n_ex)
	y_err_sim = [ y_sim - lq, uq - y_sim ]

	x_sim = np.arange(1,conv.bins+1) + 0.5
	axs[0].errorbar(x_sim, y_sim, yerr=y_err_sim, marker='o', markersize=6, markeredgewidth=0.0, linestyle='none', linewidth=1.5, color=conv.color_barrier)

	axs[0].set_xticks(np.arange(0,6) + 0.5)
	axs[0].set_xticklabels(np.arange(0,6))
	axs[0].set_xlim([1, 6])

	axs[0].set_xlabel('number of parietal cells/\nbarrier segments')
	axs[0].set_ylabel('rel. frequency of glands')

	for i, n in enumerate(conv.barrier_numbers):
		x_ex = np.arange(1,conv.bins+1) - 0.5 - 0.5 * bw
		y_ex = conv.prob(angle_dist_ex[n])
		e_ex = stdprop(angle_dist_ex[n])
		axs[i+1].bar(x_ex, y_ex, width=bw, color=conv.dark_color_barrier, ecolor=conv.color_barrier)

		y_sim = conv.prob(angle_dist_sim[n])
		redef = lambda x: x
		lq, uq, _, _ = conv.prob_ci(y_sim, redef, n_ex)
		y_err_sim = [ y_sim - lq, uq - y_sim ]

		x_sim = x_ex + 0.5
		#axs[i+1].bar(x_sim, y_sim, width=0.5, color=conv.color_barrier)
		axs[i+1].errorbar(x_sim, y_sim, yerr=y_err_sim, marker='o', markersize=6, markeredgewidth=0.0, linestyle='none', linewidth=1.5, color=conv.color_barrier)

		axs[i+1].set_xticks(np.array([0, 1./8., 2./8., 3/8., 4./8. ]) * conv.bins)
		axs[i+1].set_xticklabels([ '0', r'$\frac{1}{8}$', r'$\frac{1}{4}$', r'$\frac{3}{8}$', r'$\frac{1}{2}$' ])
		#axs[i+1].set_xlim([-0, conv.bins+1])
		axs[i+1].set_xlim([0, 4])
		axs[i+1].set_xlabel('barrier distance (1/gland)')

		axs[i+1].set_yticks(np.arange(10) * 0.2)
		axs[i+1].set_ylim([0., 1.])

		axs[i+1].set_title('glands with ' + str(n) + ' parietal cells/\nbarrier segments', fontsize=10)

	axs[1].set_ylabel('rel. frequency of distance')

	fig.set_tight_layout(True)
	fig.savefig(conv.results_path + 'angle_dists.pdf', facecolor=fig.get_facecolor(), edgecolor='none')

if __name__=='__main__':
	main()
