
import numpy as np
import pandas as pd
import h5py
from matplotlib import pyplot as plt
from scipy.stats import spearmanr

from pre_process import remove_duplicates_by_umi
from utilities import load_configs, hasOL, load_mc4c, get_nreads_per_bin

# initialization
cfg_fname = 'LVR-BMaj-96x'
# config_file = 'WPL-KOC'
# config_file = 'WPL-KOD2'
# config_file = 'NPC-PCDHa4-96x'
# config_file = 'LVR-HS2-96x'
# cfg_fname = 'BRN-BMaj-96x,BRN-BMaj-96x2'
min_mq = 20
cvg_tresh = 0.02
y_lim = [0, 0.1]


def load_data(config_lst, vp_crd, roi_crd):

    # load mc4c data
    mc4c_pd = load_mc4c(config_lst, unique_only=False, valid_only=True, min_mq=min_mq, reindex_reads=True, verbose=True)
    MAX_ReadID = np.max(mc4c_pd['ReadID'])
    print('There are {:d} reads in the dataset.'.format(len(np.unique(mc4c_pd['ReadID']))))

    # filtering reads according to their MQ
    header_lst = ['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'MQ', 'Flag']
    read_all = mc4c_pd[header_lst].values
    is_mapped = read_all[:, 4] >= min_mq
    is_valid = np.bitwise_and(read_all[:, 5], 1) == 0
    read_all = read_all[is_mapped & is_valid, :4]
    print('Selected non-overlapping fragments with MQ >= {:d}: {:d} reads are left.'.format(min_mq, len(np.unique(read_all[:, 0]))))
    # del is_mapped, is_valid

    # select informative reads (#frg > 1), ignoring VP fragments
    is_vp = hasOL(vp_crd, read_all[:, 1:4], offset=0)
    is_roi = hasOL(roi_crd, read_all[:, 1:4], offset=0)
    frg_roi = read_all[~is_vp & is_roi, :].copy()
    read_n_roi = np.bincount(frg_roi[:, 0], minlength=MAX_ReadID + 1)
    is_inf = np.isin(read_all[:, 0], frg_roi[read_n_roi[frg_roi[:, 0]] > 1, 0])
    read_inf = np.hstack([read_all[is_inf, :].copy(), read_n_roi[read_all[is_inf, 0]].reshape(-1, 1)])
    print('Selected reads #cis fragment > 1: {:d} reads are selected.'.format(len(np.unique(read_inf[:, 0]))))

    return read_all, read_inf


# load data
config_lst = load_configs(cfg_fname)
run_id = ','.join([config['run_id'] for config in config_lst])
configs = config_lst[0]
vp_crd = np.array([configs['vp_cnum'], configs['vp_start'], configs['vp_end']])
roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
read_all, read_inf = load_data(config_lst, vp_crd=vp_crd, roi_crd=roi_crd)

# select reads with #traceable fragment > 1
roi_size = configs['roi_end'] - configs['roi_start']
lcl_crd = np.array([configs['vp_cnum'], configs['roi_start'] - roi_size, configs['roi_end'] + roi_size])
is_lcl = hasOL(lcl_crd, read_inf[:, 1:4], offset=0)
umi_org = read_inf[~is_lcl, :].copy()
print('Selected reads with #trans fragment > 0: {:d} reads are selected.'.format(len(np.unique(umi_org[:, 0]))))
unq_set, lcl_info = remove_duplicates_by_umi(umi_org)
is_unq = np.isin(read_all[:, 0], unq_set[:, 0])
frg_org = read_all[is_unq, :].copy()
n_read = [len(np.unique(frg_org[:, 0]))]
print('\t#reads: {:,d} --> {:,d}'.format(len(np.unique(read_all[:, 0])), n_read[0]))

# use cis chromosome
is_cis = lcl_crd[0] == read_inf[:, 1]
umi_trs = read_inf[~is_cis, :].copy()
print('Selected reads with #trans fragment > 0: {:d} reads are selected.'.format(len(np.unique(umi_trs[:, 0]))))
unq_set, trs_info = remove_duplicates_by_umi(umi_trs)
is_unq = np.isin(read_all[:, 0], unq_set[:, 0])
frg_trs = read_all[is_unq, :].copy()
n_read.append(len(np.unique(frg_trs[:, 0])))
print('\t#reads: {:,d} --> {:,d}'.format(len(np.unique(read_all[:, 0])), n_read[1]))

# make profile
edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
n_bin = bin_bnd.shape[0]
bin_w = bin_bnd[0, 1] - bin_bnd[0, 0]
del edge_lst

# compute profiles
bin_frq = np.zeros([2, n_bin], dtype=np.int)
bin_frq[0, :] = get_nreads_per_bin(frg_org[:, :4], n_bin=200, boundary=roi_crd)
bin_frq[1, :] = get_nreads_per_bin(frg_trs[:, :4], n_bin=200, boundary=roi_crd)

# compare profiles
plt.figure(figsize=[8, 7])
plt.plot(bin_frq[0] / float(n_read[0]), bin_frq[1] / float(n_read[1]), 'o')
plt.xlim([0, 0.1])
plt.ylim([0, 0.1])
plt.title(run_id + '\n' +
          'ROI coverage, far-cis + trans vs. trans only duplicate filter\n' +
          'spr-corr: {:0.5f}'.format(spearmanr(bin_frq.T).correlation))
out_fname = './plots/AdjustingROI_{:s}_'.format(run_id)
plt.savefig(out_fname + 'binCvg_Cmp.pdf', bbox_inches='tight')

# check coverage across the chromosome
from utilities import get_chr_info
chr_size = get_chr_info(genome_str=configs['genome_build'], property='chr_size')
blk_lst = np.linspace(0, chr_size[configs['vp_cnum'] - 1], num=5001, dtype=np.int64).reshape(-1, 1)
blk_bnd = np.hstack([blk_lst[:-1], blk_lst[1:] - 1])
n_blk = blk_bnd.shape[0]
blk_w = blk_bnd[0, 1] - blk_bnd[0, 0]

# looping over bins
blk_cvg = np.zeros([2, n_blk], dtype=np.int)
for bi in range(n_blk):
    is_in = hasOL(blk_bnd[bi, :], frg_org[:, 2:4])
    blk_cvg[0, bi] = len(np.unique(frg_org[is_in, 0]))

    is_in = hasOL(blk_bnd[bi, :], frg_trs[:, 2:4])
    blk_cvg[1, bi] = len(np.unique(frg_trs[is_in, 0]))

# check similarity
plt.figure(figsize=[8, 7])
plt.plot(blk_cvg[0] / float(n_read[0]), blk_cvg[1] / float(n_read[1]), 'o')
plt.xlim([0, cvg_tresh+0.02])
plt.ylim([0, cvg_tresh+0.02])
plt.title('Chromosome coverage, far-cis + trans vs. trans only duplicate filter\n' +
          'spr-corr: {:0.5f}'.format(spearmanr(blk_cvg.T).correlation))
plt.savefig(out_fname + 'ChrCvg_Cmp.pdf', bbox_inches='tight')

# plot chr coverage
plt_h = [None] * 3
plt.figure(figsize=[15, 5])
# plt.plot(bin_bnd[:, 0], chr_cvg[0] / float(n_read[0]), '+', color='green')
p0_pd = pd.Series(blk_cvg[0] / float(n_read[0]))
plt_h[0] = plt.plot(blk_bnd[:, 0], p0_pd.rolling(5).mean().values, ':', color='green')[0]
# plt.plot(bin_bnd[:, 0], chr_cvg[1] / float(n_read[1]), 'o', color='orange', markerfacecolor='None')
p1_pd = pd.Series(blk_cvg[1] / float(n_read[1]))
plt_h[1] = plt.plot(blk_bnd[:, 0], p1_pd.rolling(5).mean().values, ':', color='orange')[0]
plt.plot([roi_crd[1], roi_crd[1]], [0, 1], color='k')
plt.plot([roi_crd[2], roi_crd[2]], [0, 1], color='k')
plt.plot([lcl_crd[1], lcl_crd[1]], [0, 1], color='gray')
plt.plot([lcl_crd[2], lcl_crd[2]], [0, 1], color='gray')
plt.text(roi_crd[1], y_lim[1] * 0.7, '{:,d}> '.format(roi_crd[1]), horizontalalignment='right')
plt.text(roi_crd[2], y_lim[1] * 0.7, ' <{:,d}'.format(roi_crd[2]), horizontalalignment='left')
plt.plot([blk_lst[0], blk_lst[-1]], [cvg_tresh, cvg_tresh], '-.', color='lightgray')
plt.ylim([0, 0.1])
plt.legend(plt_h, [
    'Far-cis (n={:0.0f})'.format(n_read[0]),
    'Trans only (n={:0.0f})'.format(n_read[1])
])
plt.savefig(out_fname + 'ProfileChrom.pdf', bbox_inches='tight')

# use adjusted area
prf_rol = p1_pd.rolling(5).mean().values
high_cvg_blks = np.where(prf_rol > cvg_tresh)[0]
adj_crd = np.array([configs['vp_cnum'], blk_bnd[high_cvg_blks[0], 0], blk_bnd[high_cvg_blks[-1], 1]])
vp_crd = np.array([configs['vp_cnum'], configs['vp_start'], configs['vp_end']])
roi_cen = np.mean([configs['vp_start'], configs['vp_end']])
edge_lst_adj = np.linspace(adj_crd[1], adj_crd[2], num=201, dtype=np.int64).reshape(-1, 1)
bin_w_adj = edge_lst_adj[1] - edge_lst_adj[0]
vp_crd_adj = np.array([configs['vp_cnum'], roi_cen - int(bin_w_adj * 1.5), roi_cen + int(bin_w_adj * 1.5)])
read_all, read_inf = load_data(config_lst, vp_crd=vp_crd_adj, roi_crd=adj_crd)

adj_w = adj_crd[2] - adj_crd[1]
is_adj = hasOL([adj_crd[0], adj_crd[1]-adj_w, adj_crd[2]+adj_w], read_inf[:, 1:4], offset=0)
umi_set = read_inf[~is_adj, :].copy()
print('Selected reads with #trans fragment > 0: {:d} reads are selected.'.format(len(np.unique(umi_set[:, 0]))))
adj_set, adj_info = remove_duplicates_by_umi(umi_set)
is_adj = np.isin(read_all[:, 0], adj_set[:, 0])
frg_adj = read_all[is_adj, :].copy()
n_read.append(len(np.unique(frg_adj[:, 0])))
print('\t#reads: {:,d} --> {:,d}'.format(len(np.unique(read_all[:, 0])), n_read[2]))
plt.xlim([adj_crd[1]-adj_w*3, adj_crd[2]+adj_w*3])
plt.ylim(y_lim)
plt.savefig(out_fname + 'ProfileLocal.pdf', bbox_inches='tight')

# plots for adjust
plt_h[2] = plt.plot([adj_crd[1], adj_crd[1]], [0, 1], color='red')[0]
plt_h[2] = plt.plot([adj_crd[2], adj_crd[2]], [0, 1], color='red')[0]
plt.plot([adj_crd[1]-adj_w, adj_crd[1]-adj_w], [0, 1], color='#fe9090')
plt.plot([adj_crd[2]+adj_w, adj_crd[2]+adj_w], [0, 1], color='#fe9090')
plt.text(adj_crd[1], y_lim[1] * 0.8, '{:,d}> '.format(adj_crd[1]), horizontalalignment='right', color='red')
plt.text(adj_crd[2], y_lim[1] * 0.8, ' <{:,d}'.format(adj_crd[2]), horizontalalignment='left', color='red')
plt.xlim([adj_crd[1]-adj_w*3, adj_crd[2]+adj_w*3])
plt.ylim([0, 0.1])
plt.legend(plt_h, [
    'Far-cis (n={:0.0f})'.format(n_read[0]),
    'Trans only (n={:0.0f})'.format(n_read[1]),
    'Adjusted (n={:0.0f})'.format(n_read[2])
])
plt.title('{:s}\n'.format(run_id) +
          'bin_w={:0.0f}, block_w={:0.0f}k\n'.format(bin_w, blk_w / 1e3))
plt.savefig(out_fname + 'ProfileLocal_withAdj.pdf', bbox_inches='tight')

# assert np.array_equal(read_all, mc4c_pd.loc[is_mapped & is_valid, header_lst[:-2]].values)
print('Done')
# plt.show()