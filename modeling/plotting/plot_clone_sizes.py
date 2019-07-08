'''
	Plots experimental data and simulation results for
	the clone size distributions.

	Author: David Jorg
'''

import numpy as np
import scipy.stats as stats
import pylab as plt

import csv, sys
import pandas as pd

import conventions as conv
if conv.screen_mode:
	import plotstyle_screen
else:
	import plotstyle_print

flatten = lambda l: [item for sublist in l for item in sublist]

id = lambda x: x

def pvalue(f_obs, p_exp, modifier=id):
	chi2, p = stats.power_divergence(f_obs=modifier(f_obs), f_exp=modifier(p_exp * np.sum(f_obs)), lambda_='log-likelihood')
	return p

def main():
	sim_id = 'fit'
	mode = 'cdf' if len(sys.argv) < 3 else sys.argv[2]

	# --------------------------------------------------------------------------
	# import experimental data

	data_ex = {
		ds: {
			t: pd.read_csv(conv.data_path + 'clone_sizes_abs_' + conv.ds_prefix[ds] + conv.time_prefix[t] + '_' + conv.roi + '.csv', header=None).values for t in conv.times[ds]
		} for ds in conv.datasets
	}

	# --------------------------------------------------------------------------
	# import simulation data

	# get parameters
	params = pd.read_csv(conv.sim_path + 'queue_clonal_' + sim_id + '/parameters.csv', header=None, delim_whitespace=True).values.tolist()
	get_param = lambda name: [ e[1] for e in params if e[0] == name ][0]

	# get simulation data
	data_sim = {
		ds: pd.read_csv(conv.sim_path + 'queue_clonal' + conv.ds_prefix[ds] + '_' + sim_id + '/global_relative_clone_sizes.csv', header=None, delim_whitespace=True).values[:,1:] for ds in conv.datasets
	}

	# only keep relevant time points
	dt = float(get_param('dt'))
	time_indices = { ds: (np.array(conv.times[ds]) / ( dt / conv.timescale )).astype(int) for ds in conv.datasets }
	data_sim = {
		ds: {
			t: data_sim[ds][time_indices[ds][i]] for i, t in enumerate(conv.times[ds])
		} for ds in conv.datasets
	}

	# --------------------------------------------------------------------------
	# goodness-of-fit test

	def fill_obs(u_obs):
		u = np.array(u_obs)
		uf = []
		for i in range(1,conv.bins+1):
			if i in u[:,0]:
				uf.append([ e[1] for e in u if e[0] == i ][0])
			else:
				uf.append(0)
		return np.array(uf)

	pv = {
		ds: {
			t: pvalue(f_obs=fill_obs(data_ex[ds][t]), p_exp=conv.prob(data_sim[ds][t])) for i, t in enumerate(conv.times[ds])
		} for ds in conv.datasets
	}

	# --------------------------------------------------------------------------
	# plot

	if mode == 'cdf':
		redef = lambda h: np.cumsum(h)
	elif mode == 'cdf_inv':
		redef = lambda h: 1. - np.cumsum(h)
	else:
		redef = lambda h: h

	bw = 0.8
	for ds in conv.datasets:
		ncols = len(conv.times[ds])
		fig, axs = plt.subplots(figsize=(2*ncols,2.5), nrows=1, ncols=ncols, sharey=True)

		for i, t in enumerate(conv.times[ds]):
			h_ex = data_ex[ds][t][:,1]
			n_clones = int(np.sum(h_ex))

			y_ex = redef(conv.prob(h_ex))
			x_ex = data_ex[ds][t][:,0] - 0.5 - 0.5 * bw

			if mode == 'pdf':
				yerr_ex = np.sqrt(y_ex * (1. - y_ex) / np.sum(h_ex)) # standard error of the proportion
			else:
				yerr_ex = np.full(len(x_ex), 0.)

			axs[i].plot(x_ex + 0.5 * bw, y_ex, linestyle='', marker='o', markersize=6, color=conv.dark_colors[ds])

			p_sim = conv.prob(data_sim[ds][t])
			lq, uq, _, _ = conv.prob_ci(p_sim, redef, n_clones)

			y_sim = np.array(redef(p_sim))
			y_err_sim = np.array([y_sim - lq, uq - y_sim])
			x_sim = np.arange(len(y_sim)) + 0.5

			axs[i].plot(x_sim, y_sim, marker='', markersize=6, markeredgewidth=0.0, linestyle='-', linewidth=1.5, color=conv.dark_colors[ds])
			axs[i].fill_between(x_sim, lq, uq, color=conv.lighter(conv.colors[ds], 0.7))

			axs[i].set_xticks(np.arange(0,4+1) * 0.25 * conv.bins)
			axs[i].set_xticklabels([ '0', '1/4', '1/2', '3/4', '1' ])
			axs[i].set_xlim([0., conv.bins])
			axs[i].set_xlabel('clone size (1/gland)')

			axs[i].set_yticks(np.arange(10) * 0.2)
			axs[i].set_ylim([0., 0.85 if mode == 'pdf' else 1.1])

			axs[i].set_title(conv.time_name[t] + ' (' + str(n_clones) + ' clones)', fontsize=10)

		if mode == 'cdf':
			axs[0].set_ylabel('cumulative distribution')
		else:
			axs[0].set_ylabel('rel. frequency')

		fig.set_tight_layout(True)
		fig.savefig(conv.results_path + 'clone_size_dist_' + ds + '.pdf', facecolor=fig.get_facecolor(), edgecolor='none')

if __name__=='__main__':
	main()
