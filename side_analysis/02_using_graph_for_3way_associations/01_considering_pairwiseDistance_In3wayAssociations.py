#! /usr/bin/env python2

import argparse
import sys

import numpy as np
from matplotlib import pyplot as plt, patches

sys.path.insert(0, '../../')
from utilities import load_configs, load_mc4c, hasOL, flatten
from utilities import get_gauss_kernel
from analysis import compute_mc_associations

np.set_printoptions(linewidth=180, threshold=5000, formatter={'float_kind': '{:6.3f}'.format})  # , suppress=True,
# pd.set_option('display.width', 250)
# pd.set_option('display.max_rows', 200)
# pd.set_option('display.max_columns', 25)

# creating argument parser
cfg_fname = ','.join(['../../configs/cfg_Prdm14_Slc_WT.cfg',
                      '../../configs/cfg_Prdm14_Slc_WT2.cfg',
                      '../../configs/cfg_Prdm14_Slc_WT3.cfg'])
parser = argparse.ArgumentParser()
parser.add_argument('--config_file', default=cfg_fname)
parser.add_argument('--min_n_frg', default=2, type=int)
parser.add_argument('--sigma', default=2.0, type=float)
args = parser.parse_args(sys.argv[1:])

# load config file
config_lst = load_configs(args.config_file)
configs = config_lst[0]
vp_bnd = [configs['vp_start'], configs['vp_end']]
roi_w = configs['roi_end'] - configs['roi_start']
roi_bnd = [configs['roi_start'] - roi_w / 2, configs['roi_end'] + roi_w / 2]
soi_bnd = [13122000, 13127000]
# soi_bnd = [12830000, 12835000]
soi_crd = [configs['vp_cnum']] + soi_bnd

# make bins
edge_lst = np.linspace(roi_bnd[0], roi_bnd[1], num=301, dtype=np.int64).reshape(-1, 1)
bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
bin_cen = np.mean(bin_bnd, axis=1)
bin_w = bin_bnd[0, 1] - bin_bnd[0, 0]
n_bin = bin_bnd.shape[0]
y_lim = [0, 10]

# load MC-4C data
frg_dp = load_mc4c(config_lst, unique_only=True, valid_only=True, min_mq=20, reindex_reads=True, verbose=True,
                   data_path='../../datasets/')
frg_np = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'Strand']].values
del frg_dp

# select within roi fragments
vp_crd = [configs['vp_cnum'], configs['vp_start'], configs['vp_end']]
roi_crd = [configs['vp_cnum']] + roi_bnd
is_vp = hasOL(vp_crd, frg_np[:, 1:4])
is_roi = hasOL(roi_crd, frg_np[:, 1:4])
frg_roi = frg_np[~is_vp & is_roi, :]
del frg_np

# filter small circles (>1 roi-frg, ex.)
cir_size = np.bincount(frg_roi[:, 0])[frg_roi[:, 0]]
frg_inf = frg_roi[cir_size >= args.min_n_frg, :]
frg_inf[:, 0] = np.unique(frg_inf[:, 0], return_inverse=True)[1]

# convert fragments to bin-coverage
n_read = len(np.unique(frg_inf[:, 0]))
rd2bins = [list() for i in range(n_read)]
n_frg = frg_inf.shape[0]
for fi in range(n_frg):
    bin_idx = np.where(hasOL(frg_inf[fi, 2:4], bin_bnd))[0]
    if frg_inf[fi, 4] == 1:
        rd2bins[frg_inf[fi, 0]].append(bin_idx[0])
    else:
        rd2bins[frg_inf[fi, 0]].append(bin_idx[-1])

# filter circles for (>1 bin cvg)
valid_lst = []
for rd_idx in range(n_read):
    hit_bins = np.unique(rd2bins[rd_idx])
    if len(hit_bins) >= args.min_n_frg:
        valid_lst.append(rd_idx)
frg_inf = frg_inf[np.isin(frg_inf[:, 0], valid_lst), :]
frg_inf[:, 0] = np.unique(frg_inf[:, 0], return_inverse=True)[1]

# convert fragments to bin-coverage
n_read = len(np.unique(frg_inf[:, 0]))
rd2bins = [list() for i in range(n_read)]
for fi in range(frg_inf.shape[0]):
    bin_idx = np.where(hasOL(frg_inf[fi, 2:4], bin_bnd))[0]
    if frg_inf[fi, 4] == 1:
        rd2bins[frg_inf[fi, 0]].append(bin_idx[0])
    else:
        rd2bins[frg_inf[fi, 0]].append(bin_idx[-1])

# fill in coverage matrix
cvg_mat = np.zeros([n_bin, n_bin])
for bi in range(n_bin):
    if not hasOL(bin_bnd[bi], [configs['roi_start'], configs['roi_end']])[0]:
        continue
    hit_lst = np.unique(frg_inf[hasOL(bin_bnd[bi], frg_inf[:, 2:4]), 0])
    for rd_idx in hit_lst:
        cvg_mat[bi, rd2bins[rd_idx]] += 1
val_idxs = np.where(np.sum(cvg_mat, axis=1) != 0)[0]
cvg_mat = cvg_mat[val_idxs, :]
n_row = cvg_mat.shape[0]

cvg_st1 = cvg_mat / np.sum(cvg_mat, axis=1).reshape(-1, 1)
cvg_med = cvg_st1 - np.median(cvg_st1, axis=0).reshape(1, -1)
cvg_avg = cvg_st1 - np.mean(cvg_st1, axis=0).reshape(1, -1)

cvg_rol = np.zeros([n_row, n_bin])
for bi in range(n_row):
    cvg_rol[bi, :] = np.roll(cvg_avg[bi, :], n_bin // 2 - val_idxs[bi])
rol_med = cvg_rol - np.median(cvg_rol, axis=0).reshape(1, -1)
rol_avg = cvg_rol - np.mean(cvg_rol, axis=0).reshape(1, -1)

# smoothen
rol_smt = np.zeros([n_row, n_bin])
kernel = get_gauss_kernel(size=11, sigma=args.sigma, ndim=1)
print(kernel)
for ri in range(n_row):
    rol_smt[ri, :] = np.convolve(cvg_rol[ri, :], kernel, mode='same')

for var_idx, var_name in enumerate(['cvg_mat', 'cvg_st1', 'cvg_med', 'cvg_avg', 'cvg_rol', 'rol_med', 'rol_avg', 'rol_smt']):
    plt.close()
    fig = plt.figure(figsize=(12, 11))
    img_h = plt.imshow(eval(var_name), vmin=0, vmax=0.005)
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title(var_name)
    # plt.show()
    plt.savefig('./Eval_{:,d}-{:,d}_{:,d}-{:,d}_'.format(*vp_crd + soi_bnd) +
                '{:02d}-varName_Mat_{:s}_'.format(var_idx, var_name) +
                'sig{:0.2f}.pdf'.format(args.sigma), bbox_inches='tight')

rol_avg = np.mean(cvg_rol, axis=0)
rol_med = np.median(cvg_rol, axis=0)
smt_avg = np.mean(rol_smt, axis=0)
smt_med = np.median(rol_smt, axis=0)
plt.close()
fig = plt.figure(figsize=(12, 11))
for var_idx, var_name in enumerate(['rol_avg', 'rol_med', 'smt_avg', 'smt_med']):
    plt.plot(bin_cen, eval(var_name), label=var_name)
plt.vlines([configs['roi_start'], configs['roi_end']], ymin=0, ymax=1, colors='k', linestyles='--', linewidth=1)
plt.hlines([0], xmin=roi_bnd[0], xmax=roi_bnd[1], colors='#bbbbbb', linestyles=':', linewidth=0.5)
plt.xlim(roi_bnd)
plt.ylim([-0.005, 0.006])
plt.legend()
# plt.show()
plt.savefig('./Eval_{:,d}-{:,d}_{:,d}-{:,d}_'.format(*vp_crd + soi_bnd) +
            '{:02d}-varName_Vecs_'.format(10) +
            'sig{:0.2f}.pdf'.format(args.sigma), bbox_inches='tight')
plt.close()
exit(0)

# measure significance
# prf_frq, prf_rnd, frg_pos, frg_neg = compute_mc_associations(frg_inf, soi_crd, bin_bnd, n_perm=1000, sigma=args.sigma)
# n_pos = len(np.unique(frg_pos[:, 0]))
# prf_obs = prf_frq * 100.0 / n_pos
# record bin coverage in adj matrix
# read_cvg = np.zeros([n_read, n_bin])
# for rd_idx, read in enumerate(rd2bins):
#     assert len(read) >= args.min_n_frg
#     print(read)

# compute scores
# nrm_rnd = prf_rnd * 100.0 / n_pos
# prf_exp = np.mean(nrm_rnd, axis=0)
# prf_std = np.std(nrm_rnd, axis=0, ddof=0)
# np.seterr(all='ignore')
# bin_scr = np.divide(prf_obs - prf_exp, prf_std)
# np.seterr(all=None)

# set vp bins to nan
# is_vp = hasOL(vp_bnd, bin_bnd)
# bin_scr[is_vp] = np.nan

fig = plt.figure(figsize=(15, 3))
ax = plt.gca()
ax.plot(bin_cen, prf_obs, color='#17171f', linewidth=1, zorder=4, label='pos profile')
ax.plot(bin_cen, bin_scr, color='#d7171f', linewidth=1, zorder=4, label='z-score')

ax.add_patch(patches.Rectangle([vp_bnd[0], y_lim[0]], vp_bnd[1] - vp_bnd[0], y_lim[1] - y_lim[0], edgecolor='None', facecolor='orange', zorder=100))
ax.add_patch(patches.Rectangle([soi_bnd[0], y_lim[0]], soi_bnd[1] - soi_bnd[0], y_lim[1] - y_lim[0], edgecolor='None', facecolor='green', zorder=100))
ax.vlines(roi_bnd, ymin=y_lim[0], ymax=y_lim[1], colors='#909090', linestyles=':', linewidth=0.3)
ax.vlines([configs['roi_start'], configs['roi_end']], ymin=y_lim[0], ymax=y_lim[1], colors='k', linestyles='--', linewidth=1)
ax.set_title('#read (#roiFrg>{:d}, ex. vp)={:,d}, #pos={:d}\n'.format(args.min_n_frg - 1, n_read, n_pos) +
             'bin-w={:0.0f}; soi-w={:0.0f}; '.format(bin_w, soi_bnd[1] - soi_bnd[0]))

# final adjustments
ax.set_xlim(roi_bnd)
ax.set_ylim(y_lim)
x_ticks = np.linspace(roi_bnd[0], roi_bnd[1], 15, dtype=np.int64)
x_tick_label = ['{:0.2f}m'.format(x / 1e6) for x in x_ticks]
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_tick_label)
ax.legend()

plt.savefig('./plt_{:d}-{:,d}-{:,d}_{:,d}-{:,d}_'.format(*vp_crd + soi_bnd) +
            'sig{:0.2f}.pdf'.format(args.sigma), bbox_inches='tight')
