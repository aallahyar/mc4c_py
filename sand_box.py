
import numpy as np
import scipy.ndimage as ndimage
from copy import copy

from utilities import get_gauss_kernel
from utilities import OnlineStats
from analysis import estimate_decay_effect

np.set_printoptions(linewidth=230, threshold=300, edgeitems=30, formatter={'float_kind': lambda x: "%8.3f" % x})


def contact_test_2d(frg_inf, bin_bnd, n_perm=1000, sigma=1.0):

    # find covered bins
    n_frg = frg_inf.shape[0]
    frg_bdx = np.zeros(n_frg, dtype=np.int64)
    is_fw = frg_inf[:, 4] == 1
    frg_bdx[ is_fw] = np.searchsorted(bin_bnd[:, 0], frg_inf[ is_fw, 2], side='left') - 1
    frg_bdx[~is_fw] = np.searchsorted(bin_bnd[:, 0], frg_inf[~is_fw, 3], side='left') - 1
    frg_bdx[frg_bdx == -1] = 0  # TODO: Better correction for out of bounds fragments

    # reset read indices
    read_idxs = np.unique(frg_inf[:, 0], return_inverse=True)[1]
    n_read = np.max(read_idxs) + 1
    n_bin = bin_bnd.shape[0]
    assert len(np.unique(frg_inf[:, 1])) == 1

    # convert fragments to bin-coverage
    print('Mapping fragments to bins, and bins to reads ...')
    rds2bin = [list() for _ in range(n_read)]
    bin2rds = [list() for _ in range(n_bin)]
    for fi in range(n_frg):
        rds2bin[read_idxs[fi]].append(frg_bdx[fi])
        bin2rds[frg_bdx[fi]].append(read_idxs[fi])
    for ri in range(n_read):
        rds2bin[ri] = np.unique(rds2bin[ri])
    for bi in range(n_bin):
        bin2rds[bi] = np.unique(bin2rds[bi])
    # del frg_inf

    # compute observed coverage
    print('Calculating observed coverage ...')
    obs_cvg = np.zeros([n_bin, n_bin])
    for soi_bdx in range(n_bin):
        for rd_idx in bin2rds[soi_bdx]:
            obs_cvg[soi_bdx, rds2bin[rd_idx]] += 1
    kernel_2d = get_gauss_kernel(size=11, sigma=sigma, ndim=2)
    obs_smt = ndimage.convolve(obs_cvg, kernel_2d, mode='reflect')

    # estimate decay profile
    print('Estimating decay profile ...')
    decay_prob = estimate_decay_effect(rds2bin, n_bin, sigma=sigma)

    # estimate expected distributions
    print('Estimating expected distributions:')
    all_rids = np.arange(n_read)
    exp_obj = []
    for bi in range(n_bin):
        exp_obj.append([OnlineStats() for _ in range(n_bin)])
    for ei in range(n_perm):
        if ei % 50 == 0:
            print('\tEpoch #{:04d}/{:04d}: '.format(ei + 1, n_perm))

        bkg_cvg = np.zeros([n_bin, n_bin])
        for soi_bdx in range(n_bin):

            # select pos/neg reads
            neg_rids = all_rids[~np.isin(all_rids, bin2rds[soi_bdx])]
            n_pos = len(bin2rds[soi_bdx])
            n_neg = len(neg_rids)

            # assign probability to neg fragments
            frg_prob = decay_prob[np.abs(soi_bdx - frg_bdx)]
            frg_prob = frg_prob / np.sum(frg_prob)

            # make background coverage
            neg_fbdx = np.random.choice(frg_bdx, p=frg_prob, size=n_pos)  # TODO: selection from neg, instead of all?
            rnd_rdxs = np.random.randint(n_neg, size=n_pos)
            rnd_fdxs = np.random.randint(30, size=n_pos)
            for ni in range(n_pos):
                rnd_read = copy(rds2bin[neg_rids[rnd_rdxs[ni]]])
                rnd_read[rnd_fdxs[ni] % len(rnd_read)] = neg_fbdx[ni]
                for bin_i in rnd_read:
                    bkg_cvg[soi_bdx, bin_i] += 1
        bkg_smt = ndimage.convolve(bkg_cvg, kernel_2d, mode='reflect')

        # store the current epoch
        for bi in range(n_bin):
            for bj in range(n_bin):
                exp_obj[bi][bj].include(bkg_smt[bi, bj])

    # compute expected values
    blk_scr = np.zeros([n_bin, n_bin])
    for bi in range(n_bin):
        for bj in range(n_bin):
            if np.abs(bi - bj) <= 2:
                continue
            if exp_obj[bi][bj].std != 0:
                blk_scr[bi, bj] = (obs_smt[bi, bj] - exp_obj[bi][bj].mean) / exp_obj[bi][bj].std

    # plot
    exp_avg = np.zeros([n_bin, n_bin])
    exp_std = np.zeros([n_bin, n_bin])
    for bi in range(n_bin):
        for bj in range(n_bin):
            exp_avg[bi, bj] = exp_obj[bi][bj].mean
            exp_std[bi, bj] = exp_obj[bi][bj].std
    # from matplotlib import pyplot as plt, cm
    # from matplotlib.colors import LinearSegmentedColormap
    # plt.close('all')
    # plt.figure(figsize=[15, 13])
    # cmap_lst = [cm.get_cmap('hot_r', 20), cm.get_cmap('hot_r', 20), cm.get_cmap('hot_r', 20),
    #             LinearSegmentedColormap.from_list('rwwb', ['#ff0000', '#ffffff', '#ffffff', '#0000ff']),
    #             ]
    # for pi, mat in enumerate([obs_smt, exp_avg, exp_std, zscr_mat]):
    #     ax = plt.subplot(2, 2, pi + 1)
    #     img_h = ax.imshow(mat, cmap=cmap_lst[pi])
    #     if pi == 3:
    #         img_h.set_clim([-6, 6])
    #     elif pi < 2:
    #         img_h.set_clim([0, np.percentile(obs_smt, 96)])
    #     plt.colorbar(ax=ax, mappable=img_h)
    # out_fname = './plt_sig{:0.2f}_'.format(sigma) + \
    #             'prm{:d}_rnd{:03d}.pdf'.format(n_perm, np.random.randint(1000))
    # plt.savefig(out_fname, bbox_inches='tight')

    return obs_smt, exp_avg, exp_std, blk_scr


