'''
	Plots experimental data and simulation results for
	the average clone size.

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

def main(save_average=False):
	sim_id = 'fit'

	# --------------------------------------------------------------------------
	# import experimental data

	avg_ex = {
		ds: [
			hstmean_ex(pd.read_csv(conv.data_path + 'clone_sizes_' + conv.ds_prefix[ds] + conv.time_prefix[t] + '_' + conv.roi + '.csv', header=None).values) for t in conv.times[ds]
		] for ds in conv.datasets
	}

	if save_average:
		for ds in conv.datasets:
			with open('../Simulations/fit_data/mean_clone_sizes_' + ds + '.csv','wb') as csvfile:
				writer = csv.writer(csvfile)
				for i, t in enumerate(conv.times[ds]):
					writer.writerow([ t, avg_ex[ds][i] ])

	clones = {
		ds: [
			pd.read_csv(conv.data_path + 'clones_' + conv.ds_prefix[ds] + conv.time_prefix[t] + '_' + conv.roi + '.csv', header=None).values / float(conv.bins) for t in conv.times[ds]
		] for ds in conv.datasets
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

	cdf_sim = { ds: np.array(map(cdf, map(prob, data_sim[ds]))) for ds in conv.datasets }

	lb_sim = { ds: np.array(map(lambda x: get_quantile(x, 0.025) / float(conv.bins), cdf_sim[ds])) for ds in conv.datasets }
	ub_sim = { ds: np.array(map(lambda x: get_quantile(x, 0.975) / float(conv.bins), cdf_sim[ds])) for ds in conv.datasets }
	med_sim = { ds: np.array(map(lambda x: get_quantile(x, 0.5) / float(conv.bins), cdf_sim[ds])) for ds in conv.datasets }

	avg_sim = { ds: hstmean_sim(data_sim[ds]) for ds in conv.datasets }

	# --------------------------------------------------------------------------
	# plot

	dt = float(get_param('dt'))

	for i, ds in enumerate(conv.datasets):
		fig, ax = plt.subplots(figsize=(2.5,2.5), nrows=1, ncols=1, sharey=True)

		ax.boxplot(
			clones[ds],
			boxprops=dict(color=conv.colors[ds]),
			medianprops=dict(linewidth=2.5, color=conv.colors[ds]),
			whiskerprops=dict(linestyle='-', color=conv.colors[ds]),
			capprops=dict(linestyle='-', color=conv.colors[ds]),
			flierprops=dict(marker='o', markersize=4, markerfacecolor=conv.colors[ds]),
			positions=conv.times[ds],
			usermedians=avg_ex[ds]
		)

		y_sim = avg_sim[ds]
		x_sim = np.arange(len(y_sim)) * ( dt / conv.timescale )
		ax.plot(x_sim, y_sim, linestyle='-', linewidth=1.5, color=conv.dark_colors[ds])

		ax.set_yticks(np.arange(0,4+1) * 0.25)
		ax.set_yticklabels([ '0', '1/4', '1/2', '3/4', '1' ])

		ax.set_xlim([-0.25, 13.])
		ax.set_xlabel('time (months)')

		ax.set_ylim([0., 1.05])

		ax.set_title(conv.ds_name[ds], fontsize=11)

		ax.set_ylabel('average clone size (1/gland)')

		fig.set_tight_layout(True)
		fig.savefig(conv.results_path + 'avg_clone_size_' + ds + '.pdf', facecolor=fig.get_facecolor(), edgecolor='none')

if __name__=='__main__':
	main()
