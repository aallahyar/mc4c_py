
from __future__ import print_function
from os import path

import numpy as np

from utilities import showprogress

MIN_N_POS = 100


def estimate_decay_effect(rd2bins, n_bin, sigma):
    from utilities import flatten, get_gauss_kernel

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
    valid_rows = np.where(np.sum(cvg_mat, axis=1) >= n_read * 0.01)[0]
    cvg_mat = cvg_mat[valid_rows, :]
    n_row = cvg_mat.shape[0]

    # normalize rows to sum=1, removing mean from columns
    cvg_prob = cvg_mat / np.sum(cvg_mat, axis=1).reshape(-1, 1)
    cvg_no_mean = cvg_prob - np.mean(cvg_prob, axis=0)

    # rolling profiles to align their "view point" on top of each other
    cvg_rolled = np.zeros([n_row, n_bin])
    for ri in range(n_row):
        cvg_rolled[ri, :] = np.roll(cvg_no_mean[ri, :], n_bin // 2 - valid_rows[ri])

    # smoothening the profiles
    cvg_smooth = np.zeros([n_row, n_bin])
    kernel = get_gauss_kernel(size=11, sigma=sigma, ndim=1)
    print('Smoothing ROI decay profiles by: {:s}'.format(', '.join(['{:0.4f}'.format(k) for k in kernel])))
    for ri in range(n_row):
        cvg_smooth[ri, :] = np.convolve(cvg_rolled[ri, :], kernel, mode='same')
    cvg_stacked = np.vstack([cvg_smooth[:, n_bin // 2:], np.fliplr(cvg_smooth[:, 1:n_bin // 2 + 1])])
    cvg_stacked = np.hstack([cvg_stacked, np.zeros_like(cvg_stacked)])
    cvg_chance = np.median(cvg_stacked, axis=0)

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
    cvg_chance = cvg_chance - cvg_chance.min()
    cvg_chance[np.argmin(cvg_chance):] = 0
    decay_prob = cvg_chance / np.sum(cvg_chance)
    return decay_prob


def compute_mc_2d_associations_by_decay(frg_inf, bin_bnd, cmd_args):
    # import scipy.ndimage as ndimage

    from utilities import get_gauss_kernel, hasOL, OnlineStats, flatten, normalize_matrix

    # find covered bins for each fragment
    n_frg = frg_inf.shape[0]
    n_bin = bin_bnd.shape[0]
    bin_w = bin_bnd[0, 1] - bin_bnd[0, 0]
    frg_bdx = np.zeros(n_frg, dtype=np.int64)
    is_fw = frg_inf[:, 4] == 1
    frg_bdx[ is_fw] = np.searchsorted(bin_bnd[:, 0], frg_inf[ is_fw, 2], side='right') - 1
    frg_bdx[~is_fw] = np.searchsorted(bin_bnd[:, 0], frg_inf[~is_fw, 3], side='right') - 1
    frg_bdx[frg_bdx == -1] = 0  # TODO: Better correction for out of bounds fragments

    # reset read indices
    read_idxs = np.unique(frg_inf[:, 0], return_inverse=True)[1]
    n_read = np.max(read_idxs) + 1

    # convert fragments to bin-coverage
    print('Mapping fragments to bins, and bins to reads ...')
    frg2bin = []
    rds2bin = [list() for _ in range(n_read)]
    bin2rds = [list() for _ in range(n_bin)]
    assert len(np.unique(frg_inf[:, 1])) == 1  # make sure its cis-frags only, note: not cheking for cis=vp_chr
    for fi in range(n_frg):
        ov_bdxs = np.where(hasOL(frg_inf[fi, 2:4], bin_bnd))[0].tolist()
        frg2bin.append(ov_bdxs)
        rds2bin[read_idxs[fi]].append(ov_bdxs)
        for bin_idx in ov_bdxs:
            bin2rds[bin_idx].append(read_idxs[fi])
    for bi in range(n_bin):
        bin2rds[bi] = np.unique(bin2rds[bi])

    # compute observed coverage
    print('Calculating observed coverage ...')
    obs_org = np.zeros([n_bin, n_bin])
    bin_pids = [list() for _ in range(n_bin)]
    bin_npos = np.zeros(n_bin)
    for soi_bdx in np.arange(n_bin):
        is_pos = hasOL(bin_bnd[soi_bdx, :], frg_inf[:, 2:4], offset=bin_w)
        bin_pids[soi_bdx] = np.unique(read_idxs[is_pos]).tolist()
        bin_npos[soi_bdx] = len(bin_pids[soi_bdx])
        for rd_idx in bin_pids[soi_bdx]:
            obs_org[soi_bdx, flatten(rds2bin[rd_idx])] += 1

    print('Smoothing observed profiles using sig={:0.1f}'.format(cmd_args.sigma))
    # kernel_2d = get_gauss_kernel(size=11, sigma=args.sigma, ndim=2)
    # obs_smt = ndimage.convolve(obs_org, kernel_2d, mode='constant')
    kernel = get_gauss_kernel(size=11, sigma=cmd_args.sigma, ndim=1)
    obs_smt = np.zeros([n_bin, n_bin])
    for bi in range(n_bin):
        obs_smt[bi, :] = np.convolve(obs_org[bi, :], kernel, mode='same')

    if cmd_args.cvg_norm == 'none':
        print('No normalization is performed')
        def norm_func(x): return x.copy()
    else:
        print('Coverage normalization using: {:s}'.format(cmd_args.cvg_norm))
        if cmd_args.cvg_norm == 'iter':
            def norm_func(x): return normalize_matrix(x, method='iterative', scale=True)
        elif cmd_args.cvg_norm == 'KR':
            def norm_func(x): return normalize_matrix(x, method='KR', scale=True)
        elif cmd_args.cvg_norm == '1d':
            def norm_func(x): return normalize_matrix(x, method='1d', scale=True)
        elif cmd_args.cvg_norm == '1dNoScale':
            def norm_func(x): return normalize_matrix(x, method='1d', scale=False)
        elif cmd_args.cvg_norm == 'nrd':
            def norm_func(x):
                return x / bin_npos.reshape(-1, 1)
        else:
            raise ValueError('Unknown normalization method')
    obs_nrm = norm_func(obs_smt)

    # estimate decay profile
    print('Estimating decay profile ...')
    if cmd_args.correction == 'decay':
        decay_prob = estimate_decay_effect(rds2bin, n_bin, sigma=cmd_args.sigma)

    # estimate expected distributions
    print('Estimating expected distributions:')
    all_fids = np.arange(n_frg)
    all_rids = np.arange(n_read)
    exp_obj = []
    for bi in range(n_bin):
        exp_obj.append([OnlineStats() for _ in range(n_bin)])
    for ei in range(cmd_args.n_perm):
        if ei % 50 == 0:
            print('\tEpoch #{:04d}/{:04d}: '.format(ei, cmd_args.n_perm))

        # loop over each SOI
        bkg_org = np.zeros([n_bin, n_bin])
        for soi_bdx in range(n_bin):

            # select pos/neg reads
            neg_rids = all_rids[~np.isin(all_rids, bin_pids[soi_bdx])]
            n_pos = len(bin_pids[soi_bdx])
            n_neg = len(neg_rids)

            # assign probability to neg fragments
            if cmd_args.correction == 'decay':
                frg_prob = decay_prob[np.abs(soi_bdx - frg_bdx)]
                frg_prob = frg_prob / np.sum(frg_prob)
                rnd_frgs = [frg2bin[i] for i in np.random.choice(all_fids, p=frg_prob, size=n_pos)]
            else:
                rnd_frgs = [[]] * n_pos

            # make background coverage
            rnd_idxs = np.random.randint(n_neg, size=n_pos)
            for ni in range(n_pos):
                rnd_read = rds2bin[neg_rids[rnd_idxs[ni]]]
                rnd_nfrg = len(rnd_read)
                del_idx = np.random.randint(rnd_nfrg)
                for fi in range(rnd_nfrg):
                    if fi == del_idx:
                        bkg_org[soi_bdx, rnd_frgs[ni]] += 1
                    else:
                        bkg_org[soi_bdx, rnd_read[fi]] += 1

        # smoothing/normalizing background
        # bkg_smt = ndimage.convolve(bkg_org, kernel_2d, mode='constant')
        bkg_smt = np.zeros([n_bin, n_bin])
        for bi in range(n_bin):
            bkg_smt[bi, :] = np.convolve(bkg_org[bi, :], kernel, mode='same')
        bkg_nrm = norm_func(bkg_smt)

        # store the current epoch
        for bi in range(n_bin):
            for bj in range(n_bin):
                exp_obj[bi][bj].include(bkg_nrm[bi, bj])

    # collect background values
    bkg_avg = np.zeros([n_bin, n_bin])
    bkg_std = np.zeros([n_bin, n_bin])
    for bi in range(n_bin):
        for bj in range(n_bin):
            bkg_avg[bi, bj] = exp_obj[bi][bj].mean
            bkg_std[bi, bj] = exp_obj[bi][bj].std

    return obs_org, obs_smt, obs_nrm, bkg_avg, bkg_std, bin_npos


def compute_mc_associations_by_decay(frg_inf, pos_crd, bin_bnd, cmd_args, verbose=True):
    from utilities import hasOL, flatten

    # re-index circles
    frg_inf[:, 0] = np.unique(frg_inf[:, 0], return_inverse=True)[1] + 1
    n_read = np.max(frg_inf[:, 0])
    n_frg = frg_inf.shape[0]

    # convert fragments to bin-coverage
    frg2bins = [list() for _ in range(n_frg)]
    rd2bins = [list() for _ in range(n_read + 1)]
    n_frg = frg_inf.shape[0]
    assert len(np.unique(frg_inf[:, 1])) == 1
    for fi in range(n_frg):
        bin_idx = np.where(hasOL(frg_inf[fi, 2:4], bin_bnd))[0]
        frg2bins[fi] = bin_idx.tolist()
        rd2bins[frg_inf[fi, 0]].append(bin_idx.tolist())

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
    obs_org = np.zeros(n_bin)
    for ri in range(n_pos):
        hit_bins = flatten(rd2bins_pos[ri])
        obs_org[hit_bins] += 1

    # find fragment bin index
    frg_bdx = np.zeros(n_frg, dtype=np.int64)
    is_fw = frg_inf[:, 4] == 1
    frg_bdx[ is_fw] = np.searchsorted(bin_bnd[:, 0], frg_inf[ is_fw, 2], side='right') - 1
    frg_bdx[~is_fw] = np.searchsorted(bin_bnd[:, 0], frg_inf[~is_fw, 3], side='right') - 1
    frg_bdx[frg_bdx == -1] = 0  # TODO: Better correction for out of bounds fragments
    del is_fw

    # assign probability to neg fragments
    decay_prob = estimate_decay_effect(rd2bins, n_bin, sigma=cmd_args.sigma)
    soi_bdx = int(np.mean(np.where(hasOL(pos_crd[1:], bin_bnd))[0]))
    frg_prob = decay_prob[np.abs(soi_bdx - frg_bdx)]
    frg_prob = frg_prob / np.sum(frg_prob)

    # make background profile from negative set
    all_fids = np.arange(n_frg)
    bkg_rnd = np.zeros([cmd_args.n_perm, n_bin])
    neg_lst = range(n_neg)
    for ei in np.arange(cmd_args.n_perm):
        if verbose and (((ei + 1) % 200) == 0):
            print('\t{:d} randomized profiles are computed.'.format(ei + 1))
        np.random.shuffle(neg_lst)
        neg_fbdx = [frg2bins[i] for i in np.random.choice(all_fids, p=frg_prob, size=n_pos)]
        for ni in range(n_pos):
            frg2bins_neg = rd2bins_neg[neg_lst[ni]]
            np.random.shuffle(frg2bins_neg)
            bkg_rnd[ei, flatten(frg2bins_neg[1:])] += 1  # making sure one element is randomly ignored everytime
            bkg_rnd[ei, neg_fbdx[ni]] += 1

    # smoothing, if needed
    if cmd_args.sigma != 0:
        if verbose:
            print('Smoothing profiles using Gaussian (sig={:0.2f}) ...'.format(cmd_args.sigma))
        from utilities import get_gauss_kernel
        kernel = get_gauss_kernel(size=11, sigma=cmd_args.sigma, ndim=1)
        obs_smt = np.convolve(obs_org, kernel, mode='same')
        for ei in np.arange(cmd_args.n_perm):
            bkg_rnd[ei, :] = np.convolve(bkg_rnd[ei, :], kernel, mode='same')
    else:
        obs_smt = obs_org.copy()

    return obs_org, obs_smt, bkg_rnd, frg_pos, frg_neg, decay_prob


def compute_mc_associations(frg_inf, pos_crd, bin_bnd, n_perm, sigma, verbose=True):
    from utilities import hasOL, flatten

    # initialization
    n_bin = bin_bnd.shape[0]

    # re-index circles
    frg_inf[:, 0] = np.unique(frg_inf[:, 0], return_inverse=True)[1] + 1
    n_read = np.max(frg_inf[:, 0])

    # convert fragments to bin-coverage
    cfb_lst = [list() for _ in range(n_read + 1)]
    n_frg = frg_inf.shape[0]
    for fi in range(n_frg):
        bin_idx = np.where(hasOL(frg_inf[fi, 2:4], bin_bnd))[0]
        # if frg_inf[fi, 4] == 1:
        #     bin_idx = [bin_idx[0]]
        # else:
        #     bin_idx = [bin_idx[-1]]
        cfb_lst[frg_inf[fi, 0]].append(list(bin_idx))

    # select positive/negative circles
    is_pos = np.where(hasOL(pos_crd, frg_inf[:, 1:4]))[0]
    frg_pos = frg_inf[ np.isin(frg_inf[:, 0], frg_inf[is_pos, 0]), :]
    frg_neg = frg_inf[~np.isin(frg_inf[:, 0], frg_inf[is_pos, 0]), :]
    # bin_idx = np.where(hasOL(pos_crd[1:], bin_bnd))[0]
    # pos_ids = [rd_idx for rd_idx, rd_lst in enumerate(cfb_lst) if bin_idx[0] in flatten(rd_lst)]
    # frg_pos = frg_inf[ np.isin(frg_inf[:, 0], pos_ids), :]
    # frg_neg = frg_inf[~np.isin(frg_inf[:, 0], pos_ids), :]
    cfb_pos = [cfb_lst[i] for i in np.unique(frg_pos[:, 0])]
    cfb_neg = [cfb_lst[i] for i in np.unique(frg_neg[:, 0])]
    n_pos = len(cfb_pos)
    n_neg = len(cfb_neg)
    if verbose:
        print('#reads in sets: pos={:,d} vs. neg={:,d}'.format(n_pos, n_neg))

    # make positive profile
    prf_org = np.zeros(n_bin)
    for pi in range(n_pos):
        bin_lst = flatten(cfb_pos[pi])
        prf_org[bin_lst] += 1

    # make background profile from negative set
    prf_rnd = np.zeros([n_perm, n_bin])
    neg_lst = range(n_neg)
    for ei in np.arange(n_perm):
        if verbose and (((ei + 1) % 200) == 0):
            print('\t{:d} randomized profiles are computed.'.format(ei + 1))
        np.random.shuffle(neg_lst)
        for ni in range(n_pos):
            f2b_rnd = cfb_neg[neg_lst[ni]]
            np.random.shuffle(f2b_rnd)
            prf_rnd[ei, flatten(f2b_rnd[1:])] += 1  # making sure one element is randomly ignored everytime

    # smoothing if needed
    if sigma != 0:
        if verbose:
            print('Smoothing profiles using Gaussian (sig={:0.2f}) ...'.format(sigma))
        from utilities import get_gauss_kernel
        kernel = get_gauss_kernel(size=11, sigma=sigma, ndim=1)
        prf_smt = np.convolve(prf_org, kernel, mode='same')
        for ei in np.arange(n_perm):
            prf_rnd[ei, :] = np.convolve(prf_rnd[ei, :], kernel, mode='same')
    else:
        prf_smt = prf_org.copy()

    return prf_org, prf_smt, prf_rnd, frg_pos, frg_neg


def perform_vpsoi_analysis(config_lst, soi_name, min_n_frg):
    import platform
    import matplotlib
    if platform.system() == 'Linux':
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt, patches
    from matplotlib.colors import LinearSegmentedColormap

    from utilities import load_mc4c, load_annotation, hasOL, flatten

    # initialization
    run_id = ','.join([config['run_id'] for config in config_lst])
    configs = config_lst[0]
    if configs['output_file'] is None:
        configs['output_file'] = path.join(configs['output_dir'],
                                           'analysis_atVP-SOI_{:s}_{:s}_'.format(run_id, soi_name) +
                                           'sig{:0.2f}_corr-{:s}_'.format(configs['cmd_args'].sigma, configs['cmd_args'].correction) +
                                           'np{:0.2f}k_zlm{:0.1f}.pdf'.format(configs['cmd_args'].n_perm / 1e3, configs['cmd_args'].zscr_lim))
    edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
    bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
    bin_cen = np.mean(bin_bnd, axis=1, dtype=np.int64)
    bin_w = bin_bnd[0, 1] - bin_bnd[0, 0]
    x_lim = [configs['roi_start'], configs['roi_end']]
    y_lim = [0, 10]

    # load MC-4C data
    frg_dp = load_mc4c(config_lst, unique_only=True, valid_only=True, min_mq=20, reindex_reads=True, verbose=True)
    frg_np = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'Strand']].values
    del frg_dp

    # select within roi fragments
    vp_crd = [configs['vp_cnum'], configs['vp_start'], configs['vp_end']]
    roi_crd = [configs['vp_cnum'], configs['roi_start'], configs['roi_end']]
    is_vp = hasOL(vp_crd, frg_np[:, 1:4])
    is_roi = hasOL(roi_crd, frg_np[:, 1:4])
    frg_roi = frg_np[~is_vp & is_roi, :]
    del frg_np

    # filter small circles (>1 roi-frg, ex.)
    cir_size = np.bincount(frg_roi[:, 0])[frg_roi[:, 0]]
    frg_inf = frg_roi[cir_size >= min_n_frg, :]
    frg_inf[:, 0] = np.unique(frg_inf[:, 0], return_inverse=True)[1] + 1
    n_read = len(np.unique(frg_inf[:, 0]))

    # convert fragments to bin-coverage
    cfb_lst = [list() for i in range(n_read + 1)]
    n_frg = frg_inf.shape[0]
    for fi in range(n_frg):
        bin_idx = np.where(hasOL(frg_inf[fi, 2:4], bin_bnd))[0]
        cfb_lst[frg_inf[fi, 0]].append(bin_idx.tolist())

    # filter circles for (>1 bin cvg)
    valid_lst = []
    for rd_nid in range(1, n_read + 1):
        fb_lst = cfb_lst[rd_nid]
        bin_cvg = np.unique(flatten(fb_lst))
        if len(bin_cvg) > 1:
            valid_lst.append(rd_nid)
    frg_inf = frg_inf[np.isin(frg_inf[:, 0], valid_lst), :]
    frg_inf[:, 0] = np.unique(frg_inf[:, 0], return_inverse=True)[1] + 1
    n_read = np.max(frg_inf[:, 0])

    # get soi info
    ant_pd = load_annotation(configs['genome_build'], roi_crd=roi_crd).reset_index(drop=True)
    n_ant = ant_pd.shape[0]
    is_in = np.where(np.isin(ant_pd['ant_name'], soi_name))[0]
    assert len(is_in) == 1
    soi_pd = ant_pd.loc[is_in[0], :].copy()
    soi_crd = [soi_pd['ant_cnum'], soi_pd['ant_pos'] - int(bin_w * 1.5), soi_pd['ant_pos'] + int(bin_w * 1.5)]
    ant_pd['zscore'] = np.nan
    if hasOL(soi_crd, vp_crd)[0]:
        print('[w] Selected SOI coordinate overlaps with the view point. Ignoring the analysis')
        return ant_pd

    # compute positive profile and backgrounds
    print('Computing expected profile for bins:')
    if configs['cmd_args'].correction == 'decay':
        prf_org, prf_frq, prf_rnd, frg_pos, frg_neg, decay_prob = compute_mc_associations_by_decay(frg_inf, soi_crd, bin_bnd, cmd_args=configs['cmd_args'])
    else:
        prf_org, prf_frq, prf_rnd, frg_pos, frg_neg = compute_mc_associations(frg_inf, soi_crd, bin_bnd, n_perm=configs['cmd_args'].n_perm, sigma=configs['cmd_args'].sigma)
    n_pos = len(np.unique(frg_pos[:, 0]))
    prf_obs = prf_frq * 100.0 / n_pos
    print('{:,d} reads are found to cover '.format(n_pos) +
          '{:s} area ({:s}:{:d}-{:d})'.format(soi_pd['ant_name'], soi_pd['ant_chr'], soi_crd[1], soi_crd[2]))

    # check enough #pos
    if n_pos < MIN_N_POS:
        print('[w] #reads in the positive set is insufficient (n={:d}, required >{:d})'.format(n_pos, MIN_N_POS))
        print('Analysis is ignored ...')
        return ant_pd

    # compute scores
    nrm_rnd = prf_rnd * 100.0 / n_pos
    prf_exp = np.mean(nrm_rnd, axis=0)
    prf_std = np.std(nrm_rnd, axis=0, ddof=0)
    np.seterr(all='ignore')
    bin_scr = np.divide(prf_obs - prf_exp, prf_std)
    np.seterr(all=None)

    # set vp bins to nan
    vp_bnd = [configs['vp_start'], configs['vp_end']]
    is_vp = hasOL(vp_bnd, bin_bnd)
    bin_scr[is_vp] = np.nan

    # compute score for annotations
    print('Computing expected profile for annotations:')
    ant_pos = ant_pd['ant_pos'].values.reshape(-1, 1)
    ant_bnd = np.hstack([ant_pos - int(bin_w * 1.5), ant_pos + int(bin_w * 1.5)])
    ant_scr = np.full(shape=n_ant, fill_value=np.nan)
    for ai in range(n_ant):
        ov_idxs = np.where(hasOL(ant_bnd[ai, :], bin_bnd))[0]
        ov_idxs = ov_idxs[~np.isnan(bin_scr[ov_idxs])]
        ov_sim = 1 / np.abs(np.mean(ant_bnd[ai, :]) - bin_cen[ov_idxs])
        ov_sim = ov_sim / np.sum(ov_sim)
        ant_scr[ai] = np.sum(bin_scr[ov_idxs] * ov_sim)

    # set vp score to nan
    is_vp = hasOL(vp_bnd, ant_bnd)
    is_soi = hasOL(soi_crd[1:3], ant_bnd)
    ant_scr[is_vp | is_soi] = np.nan
    ant_pd['zscore'] = ant_scr

    # plotting
    fig = plt.figure(figsize=(15, 3))
    ax_prf = plt.subplot2grid((20, 40), (0, 0), rowspan=19, colspan=39)
    ax_cmp = plt.subplot2grid((20, 40), (0, 39), rowspan=10, colspan=1)
    ax_scr = plt.subplot2grid((20, 40), (19, 0), rowspan=1, colspan=39)

    # set up colorbar
    clr_lst = ['#ff1a1a', '#ff7575', '#ffcccc', '#ffffff', '#ffffff', '#ffffff', '#ccdfff', '#3d84ff', '#3900f5']
    clr_map = LinearSegmentedColormap.from_list('test', clr_lst, N=9)
    clr_map.set_bad('gray', 0.05)
    norm = matplotlib.colors.Normalize(vmin=-configs['cmd_args'].zscr_lim, vmax=configs['cmd_args'].zscr_lim)
    cbar_h = matplotlib.colorbar.ColorbarBase(ax_cmp, cmap=clr_map, norm=norm)
    # cbar_h.ax.tick_params(labelsize=12)
    cbar_h.ax.set_ylabel('z-score', rotation=90)

    # profile plot
    ax_prf.plot(bin_cen, prf_obs, color='#5757ff', linewidth=1, zorder=3)
    ax_prf.plot(bin_cen, prf_exp, color='#cccccc', linewidth=1, zorder=2)
    if configs['cmd_args'].correction == 'decay':
        soi_cen = np.mean(soi_crd[1:])
        ax_prf.plot(bin_cen + (soi_cen - bin_cen[0]), decay_prob * 100, color='#377d22', linewidth=0.5, alpha=0.5, zorder=200)
        ax_prf.plot(bin_cen - (bin_cen[-1] - soi_cen), decay_prob[::-1] * 100, color='#377d22', linewidth=0.5, alpha=0.5, zorder=200)
    ax_prf.fill_between(bin_cen, prf_exp - prf_std, prf_exp + prf_std, color='#ebebeb', linewidth=0.2, zorder=1)

    ax_prf.add_patch(patches.Rectangle([vp_bnd[0], y_lim[0]], vp_bnd[1] - vp_bnd[0], y_lim[1] - y_lim[0],
                                       edgecolor='None', facecolor='orange', zorder=100))
    ax_prf.add_patch(patches.Rectangle([soi_crd[1], y_lim[0]], soi_crd[2] - soi_crd[1], y_lim[1] - y_lim[0],
                                       edgecolor='None', facecolor='#377d22', zorder=100))
    ax_prf.set_xlim(x_lim)
    ax_prf.set_ylim(y_lim)
    ax_prf.set_xticks([])

    # add score plot
    ax_scr.imshow(bin_scr.reshape(1, -1), extent=x_lim + [-1, 1], aspect='auto',
                  cmap=clr_map, vmin=configs['zscr_lim'][0], vmax=configs['zscr_lim'][1])
    ax_scr.set_xlim(x_lim)
    ax_scr.set_yticks([])
    # ax_scr.spines['top'].set_color('none')

    # add annotations
    for ai in range(n_ant):
        cmap_normed = matplotlib.cm.ScalarMappable(norm=norm, cmap=clr_map)
        ax_prf.text(ant_pos[ai], y_lim[1], ant_pd.loc[ai, 'ant_name'],
                    horizontalalignment='center', verticalalignment='bottom', rotation=60)
        ax_prf.plot(ant_pos[[ai, ai]], y_lim, ':', color='#bfbfbf', linewidth=1, alpha=0.4)

        if not np.isnan(ant_scr[ai]):
            ax_prf.add_patch(patches.Rectangle([ant_bnd[ai, 0], y_lim[1] - 0.15], ant_bnd[ai, 1] - ant_bnd[ai, 0], 0.15,
                                               edgecolor='None', facecolor=cmap_normed.to_rgba(ant_scr[ai]), zorder=10))
            ax_prf.text(ant_pos[ai], y_lim[1] - 0.2, '{:+0.1f}'.format(ant_scr[ai]),
                        horizontalalignment='center', verticalalignment='top', fontweight='bold', fontsize=6)

    # final adjustments
    x_ticks = np.linspace(configs['roi_start'], configs['roi_end'], 15, dtype=np.int64)
    y_ticks = ax_prf.get_yticks()
    x_tick_label = ['{:0.2f}m'.format(x / 1e6) for x in x_ticks]
    y_tick_label = ['{:0.0f}%'.format(y) for y in y_ticks]
    ax_scr.set_xticks(x_ticks)
    ax_scr.set_xticklabels(x_tick_label)
    ax_prf.set_yticklabels(y_tick_label)
    ax_prf.set_ylabel('Percentage of reads')
    ax_prf.set_title('VP-SOI from {:s}, SOI={:s}\n'.format(run_id, soi_name) +
                     '#read (#roiFrg>{:d}, ex. vp)={:,d}, #pos={:d}, '.format(min_n_frg - 1, n_read, n_pos) +
                     'correction={:s}\nsigma={:0.2f}; '.format(configs['cmd_args'].correction, configs['cmd_args'].sigma) +
                     'bin-w={:0.0f}; soi-w={:0.0f}; '.format(bin_w, ant_bnd[0, 1] - ant_bnd[0, 0]) +
                     '#perm={:d}\n\n\n'.format(configs['cmd_args'].n_perm)
                     )
    plt.savefig(configs['output_file'], bbox_inches='tight')
    return ant_pd


def perform_soisoi_analysis(config_lst, min_n_frg, n_perm):
    import platform
    import matplotlib
    if platform.system() == 'Linux':
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    from utilities import load_mc4c, load_annotation, hasOL, flatten

    # initialization
    run_id = ','.join([config['run_id'] for config in config_lst])
    if config_lst[0]['output_file'] is None:
        config_lst[0]['output_file'] = path.join(config_lst[0]['output_dir'],
                                                 'analysis_atSOI-SOI_{:s}.pdf'.format(run_id))
    edge_lst = np.linspace(config_lst[0]['roi_start'], config_lst[0]['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
    bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
    bin_w = bin_bnd[0, 1] - bin_bnd[0, 0]
    del edge_lst

    # load MC-HC data
    frg_dp = load_mc4c(config_lst, unique_only=True, valid_only=True, min_mq=20, reindex_reads=True, verbose=True)
    frg_np = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd']].values
    del frg_dp

    # select within roi fragments
    vp_crd = [config_lst[0]['vp_cnum'], config_lst[0]['vp_start'], config_lst[0]['vp_end']]
    roi_crd = [config_lst[0]['vp_cnum'], config_lst[0]['roi_start'], config_lst[0]['roi_end']]
    is_vp = hasOL(vp_crd, frg_np[:, 1:4])
    is_roi = hasOL(roi_crd, frg_np[:, 1:4])
    frg_roi = frg_np[~is_vp & is_roi, :]
    del frg_np

    # filter small read (>1 roi-frg, ex.)
    cir_size = np.bincount(frg_roi[:, 0])[frg_roi[:, 0]]
    frg_inf = frg_roi[cir_size >= min_n_frg, :]
    frg_inf[:, 0] = np.unique(frg_inf[:, 0], return_inverse=True)[1] + 1
    n_read = len(np.unique(frg_inf[:, 0]))

    # convert fragments to bin-coverage
    cfb_lst = [list() for i in range(n_read + 1)]
    n_frg = frg_inf.shape[0]
    for fi in range(n_frg):
        bin_idx = np.where(hasOL(frg_inf[fi, 2:4], bin_bnd))[0]
        cfb_lst[frg_inf[fi, 0]].append(bin_idx.tolist())

    # filter reads for (>1 bin cvg)
    valid_lst = []
    for rd_nid in range(1, n_read + 1):
        fb_lst = cfb_lst[rd_nid]
        bin_cvg = np.unique(flatten(fb_lst))
        if len(bin_cvg) > 1:
            valid_lst.append(rd_nid)
    frg_inf = frg_inf[np.isin(frg_inf[:, 0], valid_lst), :]

    # Downsample and re-index
    # rnd_rid = np.random.choice(np.unique(frg_inf[:, 0]), 8618, replace=False)  ### random selection
    # frg_inf = frg_inf[np.isin(frg_inf[:, 0], rnd_rid), :]
    frg_inf[:, 0] = np.unique(frg_inf[:, 0], return_inverse=True)[1] + 1
    n_read = np.max(frg_inf[:, 0])

    # loop over each SOI
    ant_pd = load_annotation(config_lst[0]['genome_build'], roi_crd=roi_crd).reset_index(drop=True)
    n_ant = ant_pd.shape[0]
    ant_name_lst = ant_pd['ant_name'].values
    ant_scr = np.full(shape=[n_ant, n_ant], fill_value=np.nan)
    n_pos = np.zeros(n_ant, dtype=np.int)
    x_tick_lbl = []
    for ai in range(n_ant):
        soi_pd = ant_pd.loc[ai, :]
        soi_crd = [soi_pd['ant_cnum'], soi_pd['ant_pos'] - int(bin_w * 1.5), soi_pd['ant_pos'] + int(bin_w * 1.5)]
        if hasOL(vp_crd[1:], soi_crd[1:]):
            x_tick_lbl.append(ant_name_lst[ai])
            continue

        # compute score for annotations
        print('Computing expected profile for {:s}:'.format(soi_pd['ant_name']))
        ant_pos = ant_pd['ant_pos'].values.reshape(-1, 1)
        ant_bnd = np.hstack([ant_pos - int(bin_w * 1.5), ant_pos + int(bin_w * 1.5)])
        ant_org, ant_obs, soi_rnd, frg_pos = compute_mc_associations(frg_inf, soi_crd, ant_bnd, n_perm=n_perm, sigma=0)[:3]
        n_pos[ai] = len(np.unique(frg_pos[:, 0]))
        x_tick_lbl.append('{:s}\n#{:,d}'.format(ant_name_lst[ai], n_pos[ai]))
        del frg_pos

        # check number of positive reads
        if n_pos[ai] <= MIN_N_POS:
            print('[w] #reads (n={:d}) in the positive set is insufficient '.format(n_pos[ai]) + \
                  '(required >{:d}). This analysis is ignored ...'.format(MIN_N_POS))
            continue

        # calculate expected profile
        ant_exp = np.mean(soi_rnd, axis=0)
        ant_std = np.std(soi_rnd, axis=0, ddof=0)
        np.seterr(all='ignore')
        ant_scr[:, ai] = np.divide(ant_obs - ant_exp, ant_std)
        np.seterr(all=None)

        # set vp score to nan
        is_vp = hasOL(vp_crd[1:], ant_bnd)
        is_soi = hasOL(soi_crd[1:3], ant_bnd)
        ant_scr[is_vp | is_soi, ai] = np.nan

    # plotting
    plt.figure(figsize=(8, 7))
    ax_scr = plt.subplot2grid((40, 40), (0,  0), rowspan=39, colspan=39)
    ax_cmp = plt.subplot2grid((40, 40), (0, 39), rowspan=20, colspan=1)

    # set up colorbar
    c_lim = [-6, 6]
    clr_lst = ['#ff1a1a', '#ff7575', '#ffcccc', '#ffffff', '#ffffff', '#ffffff', '#ccdfff', '#3d84ff', '#3900f5']
    clr_map = LinearSegmentedColormap.from_list('test', clr_lst, N=9)
    clr_map.set_bad('gray', 0.2)
    norm = matplotlib.colors.Normalize(vmin=c_lim[0], vmax=c_lim[1])
    cbar_h = matplotlib.colorbar.ColorbarBase(ax_cmp, cmap=clr_map, norm=norm)
    # cbar_h.ax.tick_params(labelsize=12)
    cbar_h.ax.set_ylabel('z-score', rotation=90)
    cbar_edge = np.round(cbar_h.cmap(norm(c_lim)), decimals=2)

    # add score scatter matrix
    x_lim = [0, n_ant]
    img_h = ax_scr.imshow(ant_scr, extent=x_lim + x_lim, cmap=clr_map,
                  vmin=c_lim[0], vmax=c_lim[1], interpolation='nearest', origin='bottom')
    ax_scr.set_xlim(x_lim)
    ax_scr.set_ylim(x_lim)

    # add score values to each box
    for ai in range(n_ant):
        for aj in range(n_ant):
            if np.isnan(ant_scr[ai, aj]):
                continue
            ant_clr = np.round(img_h.cmap(img_h.norm(ant_scr[ai, aj])), decimals=2)
            if np.array_equal(ant_clr, cbar_edge[0]) or np.array_equal(ant_clr, cbar_edge[1]):
                txt_clr = '#ffffff'
            else:
                txt_clr = '#000000'
            ax_scr.text(aj + 0.5, ai + 0.5, '{:+0.1f}'.format(ant_scr[ai, aj]), color=txt_clr,
                        horizontalalignment='center', verticalalignment='center', fontsize=12)

    # final adjustments
    ax_scr.set_xticks(np.arange(n_ant) + 0.5)
    ax_scr.set_yticks(np.arange(n_ant) + 0.5)
    ax_scr.set_xticklabels(x_tick_lbl)
    ax_scr.set_yticklabels(ant_name_lst)
    ax_scr.set_xlabel('Selected SOIs')
    ax_scr.set_title('Association matrix from {:s}\n'.format(run_id) +
                     '#read (#roiFrg>{:d}, ex. vp)={:,d}, '.format(min_n_frg - 1, n_read) +
                     'bin-w={:d}; #perm={:d}'.format(config_lst[0]['bin_width'], n_perm)
                     )
    plt.savefig(config_lst[0]['output_file'], bbox_inches='tight')


def perform_at_across_roi(config_lst, min_n_frg):
    import platform
    import matplotlib
    if platform.system() == 'Linux':
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt, patches, cm
    from matplotlib.colors import LinearSegmentedColormap

    from utilities import load_mc4c, load_annotation, hasOL, flatten, limit_to_roi

    # initialization
    run_id = ','.join([config['run_id'] for config in config_lst])
    configs = config_lst[0]
    roi_w = configs['roi_end'] - configs['roi_start']
    if configs['output_file'] is None:
        if configs['cmd_args'].downsample:
            run_id += '_ds{:d}'.format(configs['cmd_args'].downsample)
        configs['output_file'] = path.join(configs['output_dir'],
                                           'analysis_atAcrossROI_{:s}_'.format(run_id) +
                                           'rw{:0.1f}kb_sig{:0.2f}_'.format(roi_w / 1e3, configs['cmd_args'].sigma) +
                                           'nrm-{:s}_corr-{:s}_'.format(configs['cmd_args'].cvg_norm, configs['cmd_args'].correction) +
                                           'np{:0.2f}k_zlm{:0.0f}.pdf'.format(configs['cmd_args'].n_perm / 1e3, configs['cmd_args'].zscr_lim))

    # create bin list
    edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
    bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
    bin_w = bin_bnd[0, 1] - bin_bnd[0, 0]
    n_bin = bin_bnd.shape[0]

    # make block list
    bin_cen = np.mean(bin_bnd, axis=1, dtype=np.int64).reshape(-1, 1)
    blk_crd = np.hstack([np.repeat(configs['vp_cnum'], n_bin).reshape(-1, 1), bin_cen - int(bin_w * 1.5), bin_cen + int(bin_w * 1.5) - 1])
    blk_w = blk_crd[0, 2] - blk_crd[0, 1]
    n_blk = blk_crd.shape[0]
    del edge_lst

    # define areas
    roi_cen = np.mean([np.min(configs['prm_start']), np.max(configs['prm_end'])], dtype=np.int)
    vp_crd = np.array([configs['vp_cnum'], roi_cen - int(bin_w * 1.5), roi_cen + int(bin_w * 1.5)])
    roi_crd = [configs['vp_cnum'], configs['roi_start'], configs['roi_end']]

    # load MC-4C data
    frg_dp = load_mc4c(config_lst, unique_only=True, valid_only=True, min_mq=20, reindex_reads=True, verbose=True)
    read_all = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'Strand']].values
    del frg_dp

    # select >2 roi-fragments
    read_inf = limit_to_roi(read_all, vp_crd=vp_crd, roi_crd=roi_crd, min_n_frg=min_n_frg)
    del read_all

    # re-index reads
    read_inf[:, 0] = np.unique(read_inf[:, 0], return_inverse=True)[1] + 1
    n_read = len(np.unique(read_inf[:, 0]))

    # convert fragments to bin-coverage
    print('Mapping reads to bins ...')
    cfb_lst = [list() for i in range(n_read + 1)]
    n_frg = read_inf.shape[0]
    for fi in range(n_frg):
        bin_idx = np.where(hasOL(read_inf[fi, 2:4], bin_bnd))[0]
        cfb_lst[read_inf[fi, 0]].append(bin_idx.tolist())

    # filter circles for (>1 bin cvg)
    print('Selecting only reads with >1 bins covered')
    valid_lst = []
    for rd_nid in range(1, n_read + 1):
        fb_lst = cfb_lst[rd_nid]
        bin_cvg = np.unique(flatten(fb_lst))
        if len(bin_cvg) > 1:
            valid_lst.append(rd_nid)
    read_inf = read_inf[np.isin(read_inf[:, 0], valid_lst), :]

    # subsample reads
    if configs['cmd_args'].downsample:
        n_read = len(np.unique(read_inf[:, 0]))
        print('[i] Downsampling {:,d} informative reads to {:d} reads.'.format(n_read, configs['cmd_args'].downsample))
        rnd_ids = np.random.choice(np.unique(read_inf[:, 0]), configs['cmd_args'].downsample, replace=False)
        read_inf = read_inf[np.isin(read_inf[:, 0], rnd_ids), :]

    # reindexing reads
    read_inf[:, 0] = np.unique(read_inf[:, 0], return_inverse=True)[1] + 1
    n_read = np.max(read_inf[:, 0])
    print('{:,d} reads are left after bin-coverage filter.'.format(n_read))

    # get soi info
    ant_pd = load_annotation(configs['genome_build'], roi_crd=roi_crd)
    ant_bnd = np.hstack([ant_pd[['ant_pos']].values, ant_pd[['ant_pos']].values])

    # choose the model
    print('Computing expected profile using "{:s}" model, '.format(configs['correction']), end='')
    print('over {:d} bins, required coverage: {:d} reads'.format(n_bin, MIN_N_POS))
    obs_org, obs_smt, obs_nrm, bkg_avg, bkg_std, bin_npos = compute_mc_2d_associations_by_decay(read_inf, bin_bnd, cmd_args=configs['cmd_args'])
    # if configs['cmd_args'].correction == 'decay':
    #     print('over {:d} bins, required coverage: {:d} reads'.format(n_bin, MIN_N_POS))
    #     obs_org, obs_smt, obs_nrm, bkg_avg, bkg_std, bin_npos = compute_mc_2d_associations_by_decay(read_inf, bin_bnd, cmd_args=configs['cmd_args'])
    # else:
    #     print('over {:d} blocks, required coverage: {:d} reads'.format(n_blk, MIN_N_POS))
    #
    #     # compute score for annotations
    #     obs_org = np.full([n_blk, n_blk], fill_value=np.nan)
    #     obs_smt = np.full([n_blk, n_blk], fill_value=np.nan)
    #     bkg_avg = np.full([n_blk, n_blk], fill_value=np.nan)
    #     bkg_std = np.full([n_blk, n_blk], fill_value=np.nan)
    #     bin_npos = np.zeros(n_bin)
    #     for bi in range(n_blk):
    #         showprogress(bi, n_blk, n_step=20)
    #         if hasOL(blk_crd[bi, :], vp_crd)[0]:
    #             continue
    #         obs_org[bi, :], obs_smt[bi, :], bkg_rnd, frg_pos = compute_mc_associations(read_inf, blk_crd[bi, :], bin_bnd, n_perm=configs['cmd_args'].n_perm, sigma=configs['cmd_args'].sigma, verbose=False)[:4]
    #         bin_npos[bi] = len(np.unique(frg_pos[:, 0]))
    #         bkg_avg[bi, :] = np.mean(bkg_rnd, axis=0)
    #         bkg_std[bi, :] = np.std(bkg_rnd, axis=0, ddof=0)
    #     obs_nrm = obs_smt.copy()

    # compute z-scores
    n_ignored = 0
    zscr_mat = np.full([n_bin, n_bin], fill_value=np.nan)
    for bi in range(n_bin):
        if bin_npos[bi] < 100:
            n_ignored += 1
            continue
        for bj in range(n_bin):
            if np.abs(bi - bj) <= 3:
                continue
            if bkg_std[bi, bj] != 0:
                zscr_mat[bi, bj] = (obs_nrm[bi, bj] - bkg_avg[bi, bj]) / bkg_std[bi, bj]
    print('[w] {:d}/{:d} blocks are ignored due to low coverage.'.format(n_ignored, n_blk))

    # export to excel file
    if configs['cmd_args'].to_tsv:
        import pandas as pd

        tsv_fname = configs['output_file'][:-4] + '.tsv'
        print('Exporting z-scores to .tsv file: {:s}'.format(tsv_fname))

        # adjust ticks
        y_tick_lbl = [[] for _ in range(n_bin)]
        for ai in range(ant_bnd.shape[0]):
            ov_idxs = np.where(hasOL(ant_bnd[ai, :], bin_bnd))[0]
            for ov_idx in ov_idxs:
                y_tick_lbl[ov_idx].append(ant_pd.at[ai, 'ant_name'])
        y_tick_lbl = [','.join(y_lbl) for y_lbl in y_tick_lbl]

        # storing
        zscr_pd = pd.DataFrame(zscr_mat, columns=[bin_cen.flatten(), y_tick_lbl], index=[bin_cen.flatten(), y_tick_lbl])
        zscr_pd.to_csv(tsv_fname, sep='\t', index=True, header=True)

    # plotting the scores
    plt.figure(figsize=(25, 14))
    x_lim = [configs['roi_start'], configs['roi_end']]
    img_names = ['Original', 'Smoothed', 'Normalized/Observed', 'Expected', 'Standard deviation', 'z-score']
    zclr_lst = ['#ff1a1a', '#ff7575', '#ffcccc', '#ffffff', '#ffffff', '#ffffff', '#ccdfff', '#3d84ff', '#3900f5']
    cmap_lst = [cm.get_cmap('hot_r', 20), cm.get_cmap('hot_r', 20), cm.get_cmap('hot_r', 20), cm.get_cmap('hot_r', 20),
                cm.get_cmap('summer_r', 20), LinearSegmentedColormap.from_list('bwwr', zclr_lst, 9)]
    for img_idx, img in enumerate([obs_org, obs_smt, obs_nrm, bkg_avg, bkg_std, zscr_mat]):
        ax = plt.subplot(2, 3, img_idx + 1)

        # plot the image
        cmap = cmap_lst[img_idx]
        cmap.set_bad('#b2b2a3', 0.4)
        img_h = plt.imshow(img, interpolation=None, cmap=cmap, extent=x_lim + x_lim, origin='lower')

        # adjustments
        if img_names[img_idx] in ['Original', 'Smoothed']:
            plt.clim(np.nanpercentile(obs_smt, [20, 95]))
        elif img_names[img_idx] in ['Normalized/Observed', 'Expected']:
            plt.clim([np.nanpercentile(obs_nrm, 20), np.nanpercentile(obs_nrm, 95)])
        elif img_names[img_idx] in ['Standard deviation']:
            plt.clim(np.nanpercentile(bkg_std, [60, 95]))
        elif img_names[img_idx] in ['z-score']:
            plt.clim(configs['zscr_lim'])

        # add vp patches
        ax.add_patch(patches.Rectangle([x_lim[0], vp_crd[1]], roi_w, vp_crd[2] - vp_crd[1],
                                       linewidth=0, edgecolor='None', facecolor='orange'))
        ax.add_patch(patches.Rectangle([vp_crd[1], x_lim[0]], vp_crd[2] - vp_crd[1], roi_w,
                                       linewidth=0, edgecolor='None', facecolor='orange'))

        plt.colorbar(img_h, fraction=0.046, pad=0.04)  # , extend='both'
        ax.set_xticks(ant_pd['ant_pos'])
        ax.set_yticks(ant_pd['ant_pos'])
        ax.set_xticklabels(ant_pd['ant_name'], rotation=45)
        ax.set_yticklabels(ant_pd['ant_name'])
        # ax.set_xlabel('Coverage/profile')
        ax.set_ylabel('Selected SOIs')
        # ax.tick_params(length=0)
        ax.set_xlim(x_lim)
        ax.set_ylim(x_lim)
        ax.invert_yaxis()
        ax.set_title(img_names[img_idx])

    # final adjustments
    plt.suptitle('Association matrix from {:s}\n'.format(run_id) +
                 '#read (#roiFrg>{:d}, ex. vp)={:,d}; '.format(min_n_frg - 1, n_read) +
                 'sigma={:0.2f}; cvg_norm={:s}, correction={:s}\n'.format(configs['cmd_args'].sigma, configs['cmd_args'].cvg_norm, configs['cmd_args'].correction) +
                 'bin-w={:0.0f}; block-w={:0.0f}; #perm={:d}'.format(bin_w, blk_w, configs['cmd_args'].n_perm)
                 )
    plt.subplots_adjust(wspace=0.25, hspace=0.15, top=0.91)
    plt.savefig(configs['output_file'], bbox_inches='tight')
    plt.close()
    print('ROI-ROI z-scores are plotted in {:s}'.format(configs['output_file']))
