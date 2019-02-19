
import numpy as np
import pandas as pd
import h5py
from matplotlib import pyplot as plt
from scipy.stats import spearmanr

from pre_process import remove_duplicates_by_umi
from utilities import load_configs, hasOL

# initialization
# config_file = 'LVR-BMaj-96x'
# config_file = 'WPL-KOC'
# config_file = 'WPL-KOD2'
config_file = 'NPC-PCDHa4-96x'
# config_file = 'LVR-HS2-96x'
configs = load_configs(config_file, max_n_configs=1)[0]
input_file = './datasets/mc4c_' + configs['run_id'] + '_all.hdf5'
min_mq = 20

# load mc4c data
h5_fid = h5py.File(input_file, 'r')
target_field = 'frg_np'
data_np = h5_fid[target_field].value
header_lst = list(h5_fid[target_field + '_header_lst'].value)
mc4c_pd = pd.DataFrame(data_np, columns=header_lst)
chr_lst = list(h5_fid['chr_lst'].value)
h5_fid.close()
MAX_ReadID = np.max(mc4c_pd['ReadID'])
print 'There are {:d} reads in the dataset.'.format(len(np.unique(mc4c_pd['ReadID'])))

# filtering reads according to their MQ
header_lst = ['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'MQ', 'ErrFlag']
read_all = mc4c_pd[header_lst].values
is_mapped = read_all[:, 4] >= min_mq
is_valid = read_all[:, 5] == 0
read_all = read_all[is_mapped & is_valid, :4]
print 'Selected non-overlapping fragments with MQ >= {:d}: {:d} reads are left.'.format(
    min_mq, len(np.unique(read_all[:, 0])))
# del is_mapped, is_valid

# select informative reads (#frg > 1), ignoring VP fragments
vp_crd = np.array([configs['vp_cnum'], configs['vp_start'], configs['vp_end']])
roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
is_vp = hasOL(vp_crd, read_all[:, 1:4], offset=0)
is_roi = hasOL(roi_crd, read_all[:, 1:4], offset=0)
frg_roi = read_all[~is_vp & is_roi, :].copy()
read_n_roi = np.bincount(frg_roi[:, 0], minlength=MAX_ReadID + 1)
is_inf = np.isin(read_all[:, 0], frg_roi[read_n_roi[frg_roi[:, 0]] > 1, 0])
read_inf = np.hstack([read_all[is_inf, :].copy(), read_n_roi[read_all[is_inf, 0]].reshape(-1, 1)])
print 'Selected reads #cis fragment > 1: {:d} reads are selected.'.format(len(np.unique(read_inf[:, 0])))

# select reads with #traceable fragment > 1
roi_size = configs['roi_end'] - configs['roi_start']
lcl_crd = np.array([configs['vp_cnum'], configs['roi_start'] - roi_size, configs['roi_end'] + roi_size])
is_lcl = hasOL(lcl_crd, read_inf[:, 1:4], offset=0)
umi_org = read_inf[~is_lcl, :].copy()
print 'Selected reads with #trans fragment > 0: {:d} reads are selected.'.format(len(np.unique(umi_org[:, 0])))
unq_set, lcl_info = remove_duplicates_by_umi(umi_org)
is_unq = np.isin(read_all[:, 0], unq_set[:, 0])
frg_org = read_all[is_unq, :].copy()
n_read = [len(np.unique(frg_org[:, 0]))]
print '\t#reads: {:,d} --> {:,d}'.format(len(np.unique(read_all[:, 0])), n_read[0])

# use cis chromosome
is_cis = lcl_crd[0] == read_inf[:, 1]
umi_trs = read_inf[~is_cis, :].copy()
print 'Selected reads with #trans fragment > 0: {:d} reads are selected.'.format(len(np.unique(umi_trs[:, 0])))
unq_set, trs_info = remove_duplicates_by_umi(umi_trs)
is_unq = np.isin(read_all[:, 0], unq_set[:, 0])
frg_trs = read_all[is_unq, :].copy()
n_read.append(len(np.unique(frg_trs[:, 0])))
print '\t#reads: {:,d} --> {:,d}'.format(len(np.unique(read_all[:, 0])), n_read[1])

# make profile
edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
n_bin = bin_bnd.shape[0]
bin_w = bin_bnd[0, 1] - bin_bnd[0, 0]
del edge_lst

# looping over bins
bin_frq = np.zeros([2, n_bin], dtype=np.int)
for bi in range(n_bin):
    is_in = hasOL(bin_bnd[bi, :], frg_org[:, 2:4])
    bin_frq[0, bi] = len(np.unique(frg_org[is_in, 0]))

    is_in = hasOL(bin_bnd[bi, :], frg_trs[:, 2:4])
    bin_frq[1, bi] = len(np.unique(frg_trs[is_in, 0]))

# compare profiles
plt.figure(figsize=[7, 7])
plt.plot(bin_frq[0] / float(n_read[0]), bin_frq[1] / float(n_read[1]), 'o')
plt.xlim([0, 0.1])
plt.ylim([0, 0.1])
plt.title('ROI coverage, far-cis + trans vs. trans only duplicate filter\n' +
          'spr-corr: {:0.5f}'.format(spearmanr(bin_frq.T).correlation))

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
plt.figure(figsize=[7, 7])
plt.plot(blk_cvg[0] / float(n_read[0]), blk_cvg[1] / float(n_read[1]), 'o')
plt.xlim([0, 0.1])
plt.ylim([0, 0.1])
plt.title('Chromosome coverage, far-cis + trans vs. trans only duplicate filter\n' +
          'spr-corr: {:0.5f}'.format(spearmanr(blk_cvg.T).correlation))

# plot chr coverage
plt.figure(figsize=[15, 5])
plt.cla()
# plt.plot(bin_bnd[:, 0], chr_cvg[0] / float(n_read[0]), '+', color='green')
p0_pd = pd.Series(blk_cvg[0] / float(n_read[0]))
plt.plot(blk_bnd[:, 0], p0_pd.rolling(5).mean().values, ':', color='green')
# plt.plot(bin_bnd[:, 0], chr_cvg[1] / float(n_read[1]), 'o', color='orange', markerfacecolor='None')
p1_pd = pd.Series(blk_cvg[1] / float(n_read[1]))
plt.plot(blk_bnd[:, 0], p1_pd.rolling(5).mean().values, ':', color='orange')
plt.plot([roi_crd[1], roi_crd[1]], [0, 1], color='k')
plt.plot([roi_crd[2], roi_crd[2]], [0, 1], color='k')
plt.plot([lcl_crd[1], lcl_crd[1]], [0, 1], color='gray')
plt.plot([lcl_crd[2], lcl_crd[2]], [0, 1], color='gray')
plt.plot([edge_lst[0], edge_lst[-1]], [0.05, 0.05], '-.', color='lightgray')
# plt.xlim(lcl_crd[1:])
plt.ylim([0, 0.1])
plt.legend([
    'Far-cis method (n={:0.0f})'.format(n_read[0]),
    'Trans only (n={:0.0f})'.format(n_read[1])
])
plt.title('{:s}\n'.format(configs['run_id']) +
          'bin_w={:0.0f}, block_w={:0.0f}k\n'.format(bin_w, blk_w / 1e3))


# use adjusted area
prf_rol = p1_pd.rolling(5).mean().values
high_cvg_blks = np.where(prf_rol > 0.05)[0]
adj_crd = np.array([configs['vp_cnum'], blk_bnd[high_cvg_blks[0], 0], blk_bnd[high_cvg_blks[-1], 1]])
adj_w = adj_crd[2] - adj_crd[1]
is_adj = hasOL([adj_crd[0], adj_crd[1]-adj_w, adj_crd[2]+adj_w], read_inf[:, 1:4], offset=0)
umi_set = read_inf[~is_adj, :].copy()
print 'Selected reads with #trans fragment > 0: {:d} reads are selected.'.format(len(np.unique(umi_set[:, 0])))
adj_set, adj_info = remove_duplicates_by_umi(umi_set)
is_adj = np.isin(read_all[:, 0], adj_set[:, 0])
frg_adj = read_all[is_adj, :].copy()
n_read.append(len(np.unique(frg_adj[:, 0])))
print '\t#reads: {:,d} --> {:,d}'.format(len(np.unique(read_all[:, 0])), n_read[2])

# plots for adjust
plt.plot([adj_crd[1], adj_crd[1]], [0, 1], color='red')
plt.plot([adj_crd[2], adj_crd[2]], [0, 1], color='red')
plt.plot([adj_crd[1]-adj_w, adj_crd[1]-adj_w], [0, 1], color='#fe9090')
plt.plot([adj_crd[2]+adj_w, adj_crd[2]+adj_w], [0, 1], color='#fe9090')

assert np.array_equal(read_all, mc4c_pd.loc[is_mapped & is_valid, header_lst[:-2]].values)
print 'Done'
plt.show()