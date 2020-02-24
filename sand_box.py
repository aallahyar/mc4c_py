
import numpy as np
import scipy.ndimage as ndimage
from copy import copy

from utilities import hasOL, flatten, get_gauss_kernel
from utilities import OnlineStats

np.set_printoptions(linewidth=250, threshold=300, edgeitems=30, formatter={'float_kind': lambda x: "%10.5f" % x})


def contact_test_2d(frg_inf, bin_bnd, n_perm=1000, verbose=True, sigma=1.0):

    # find covered bins
    n_frg = frg_inf.shape[0]
    frg_bdx = np.zeros(n_frg, dtype=np.int64)
    is_fw = frg_inf[:, 4] == 1
    frg_bdx[ is_fw] = np.searchsorted(bin_bnd[:, 0], frg_inf[ is_fw, 2], side='left') - 1
    frg_bdx[~is_fw] = np.searchsorted(bin_bnd[:, 0], frg_inf[~is_fw, 3], side='left') - 1

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
    del frg_inf

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
    decay_prob = get_decay_prob(rds2bin, n_bin, sigma=sigma)

    # estimate expected distributions
    print('Estimating expected distributions:')
    all_rids = np.arange(n_read)
    exp_obj = []
    for bi in range(n_bin):
        exp_obj.append([OnlineStats() for _ in range(n_bin)])
    for ei in range(n_perm):
        if ei % 1 == 0:
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
    zscr_mat = np.zeros([n_bin, n_bin])
    for bi in range(n_bin):
        for bj in range(n_bin):
            if bi == bj:
                continue
            if exp_obj[bi][bj].std != 0:
                zscr_mat[bi, bj] = (obs_smt[bi, bj] - exp_obj[bi][bj].mean) / exp_obj[bi][bj].std

    # plot
    exp_avg = np.zeros([n_bin, n_bin])
    exp_std = np.zeros([n_bin, n_bin])
    for bi in range(n_bin):
        for bj in range(n_bin):
            exp_avg[bi, bj] = exp_obj[bi][bj].mean
            exp_std[bi, bj] = exp_obj[bi][bj].std
    from matplotlib import pyplot as plt
    plt.close('all')
    plt.figure(figsize=[15, 13])
    for pi, mat in enumerate([obs_smt, exp_avg, exp_std, zscr_mat]):
        ax = plt.subplot(2, 2, pi + 1)
        img_h = ax.imshow(mat)
        plt.colorbar(ax=ax, mappable=img_h)



def get_decay_prob(rd2bins, n_bin, sigma):

    # fill in coverage matrix
    cvg_mat = np.zeros([n_bin, n_bin])
    n_read = 0
    for read in rd2bins:
        if len(read) == 0:
            continue
        n_read += 1
        hit_lst = np.unique(flatten(read))
        for bi in hit_lst:
            for bj in hit_lst:
                cvg_mat[bi, bj] += 1
    valid_rows = np.where(np.sum(cvg_mat, axis=1) > n_read * 0.01)[0]
    cvg_mat = cvg_mat[valid_rows, :]
    n_row = cvg_mat.shape[0]

    # normalize rows to sum=1, removing mean from columns
    cvg_row1 = cvg_mat / np.sum(cvg_mat, axis=1).reshape(-1, 1)
    cvg_both1 = cvg_row1 - np.mean(cvg_row1, axis=0)

    # rolling profiles to align their "view point" on top of each other
    cvg_rol = np.zeros([n_row, n_bin])
    for ri in range(n_row):
        cvg_rol[ri, :] = np.roll(cvg_both1[ri, :], n_bin // 2 - valid_rows[ri])

    # smoothening the profiles
    rol_smt = np.zeros([n_row, n_bin])
    kernel = get_gauss_kernel(size=11, sigma=sigma, ndim=1)
    print('Smoothing ROI decay profiles by: {:s}'.format(', '.join(['{:0.4f}'.format(k) for k in kernel])))
    for ri in range(n_row):
        rol_smt[ri, :] = np.convolve(cvg_rol[ri, :], kernel, mode='same')
    smt_stk = np.vstack([rol_smt[:, n_bin // 2:], np.fliplr(rol_smt[:, 1:n_bin // 2 + 1])])
    smt_stk = np.hstack([smt_stk, np.zeros_like(smt_stk)])
    stk_med = np.median(smt_stk, axis=0)

    # from matplotlib import pyplot as plt
    # plt.close('all')
    # ax = plt.figure(figsize=(15, 3)).gca()
    # stk_std = np.std(smt_stk, axis=0)
    # ax.plot(np.mean(smt_stk, axis=0), color='orange', linewidth=1, label='mean profile', zorder=10)
    # ax.plot(stk_med, color='red', linewidth=1, label='median profile')
    # ax.fill_between(range(n_bin), stk_med - stk_std, stk_med + stk_std, color='red', linewidth=0.2, label='std. profile', alpha=0.1)
    # plt.legend()
    # plt.show()

    # corrections
    # stk_avg[stk_avg < 0] = 0
    stk_med = stk_med - stk_med.min()
    stk_med[np.argmin(stk_med):] = 0
    decay_prob = stk_med / np.sum(stk_med)

    return decay_prob


def contact_test_by_decay(frg_inf, pos_crd, bin_bnd, decay_prob, n_perm=1000, verbose=True, sigma=1.0):

    # re-index circles
    frg_inf[:, 0] = np.unique(frg_inf[:, 0], return_inverse=True)[1] + 1
    n_read = np.max(frg_inf[:, 0])

    # convert fragments to bin-coverage
    rd2bins = [list() for i in range(n_read + 1)]
    n_frg = frg_inf.shape[0]
    assert len(np.unique(frg_inf[:, 1])) == 1
    for fi in range(n_frg):
        bin_idx = np.where(hasOL(frg_inf[fi, 2:4], bin_bnd))[0]
        rd2bins[frg_inf[fi, 0]].append(list(bin_idx))

    # select positive/negative circles
    is_pos = hasOL(pos_crd, frg_inf[:, 1:4])
    frg_pos = frg_inf[ np.isin(frg_inf[:, 0], frg_inf[is_pos, 0]), :]
    frg_neg = frg_inf[~np.isin(frg_inf[:, 0], frg_inf[is_pos, 0]), :]
    rd2bins_pos = [rd2bins[i] for i in np.unique(frg_pos[:, 0])]
    rd2bins_neg = [rd2bins[i] for i in np.unique(frg_neg[:, 0])]
    n_pos = len(rd2bins_pos)
    n_neg = len(rd2bins_neg)
    if verbose:
        print('#reads in sets: pos={:,d} vs. neg={:,d}'.format(n_pos, n_neg))

    # make positive profile
    n_bin = bin_bnd.shape[0]
    prf_pos = np.zeros(n_bin)
    for pi in range(n_pos):
        hit_bins = flatten(rd2bins_pos[pi])
        prf_pos[hit_bins] += 1

    decay_prob = get_decay_prob(rd2bins, n_bin, sigma=sigma)

    # assign probability to neg fragments
    soi_bdx = int(np.mean(np.where(hasOL(pos_crd[1:], bin_bnd))[0]))
    frg_bdx = np.searchsorted(bin_bnd[:, 0], frg_inf[:, 2], side='left') - 1  # selection from all reads, should be from neg reads
    frg_prob = decay_prob[np.abs(soi_bdx - frg_bdx)]
    frg_prob = frg_prob / np.sum(frg_prob)

    # make background profile from negative set
    prf_rnd = np.zeros([n_perm, n_bin])
    neg_lst = range(n_neg)
    for ei in np.arange(n_perm):
        if verbose and (((ei + 1) % 200) == 0):
            print('\t{:d} randomized profiles are computed.'.format(ei + 1))
        np.random.shuffle(neg_lst)
        neg_fbdx = np.random.choice(frg_bdx, p=frg_prob, size=n_pos)
        for ni in range(n_pos):
            frg2bins_neg = rd2bins_neg[neg_lst[ni]]
            np.random.shuffle(frg2bins_neg)
            prf_rnd[ei, flatten(frg2bins_neg[1:])] += 1  # making sure one element is randomly ignored everytime
            prf_rnd[ei, neg_fbdx[ni]] += 1

    # smoothing, if needed
    if sigma != 0:
        if verbose:
            print('Smoothing profiles using Gaussian (sig={:0.2f}) ...'.format(sigma))
        from utilities import get_gauss_kernel
        kernel = get_gauss_kernel(size=7, sigma=sigma, ndim=1)
        prf_pos = np.convolve(prf_pos, kernel, mode='same')
        for ei in np.arange(n_perm):
            prf_rnd[ei, :] = np.convolve(prf_rnd[ei, :], kernel, mode='same')

    return prf_pos, prf_rnd, frg_pos, frg_neg, decay_prob
