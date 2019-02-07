import numpy as np


def perform_mc_analysis(configs, min_n_frg=2):
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt, patches
    from matplotlib.colors import LinearSegmentedColormap

    from utilities import hasOL
    from mc4c_tools import load_mc4c, load_annotation

    # initialization
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/plt_mcTest_' + configs['run_id'] + '.pdf'
    edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
    bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
    n_bin = bin_bnd.shape[0]
    n_epoch = 1000
    x_lim = [configs['roi_start'], configs['roi_end']]

    # load MC-HC data
    frg_dp = load_mc4c(configs, only_unique=True, only_valid=True, min_mq=20, reindex_reads=False)
    frg_np = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd']].values
    del frg_dp

    # select within roi fragments
    vp_crd = [configs['vp_cnum'], configs['vp_start'], configs['vp_end']]
    roi_crd = [configs['vp_cnum'], configs['roi_start'], configs['roi_end']]
    is_vp = hasOL(vp_crd, frg_np[:, 1:4])
    is_roi = hasOL(roi_crd, frg_np[:, 1:4])
    frg_roi = frg_np[~is_vp & is_roi, :]
    del frg_np

    # filter small circles
    cir_size = np.bincount(frg_roi[:, 0])[frg_roi[:, 0]]
    frg_roi = frg_roi[cir_size >= min_n_frg, :]
    n_read = len(np.unique(frg_roi[:, 0]))

    # re-index circles
    frg_roi[:, 0] = np.unique(frg_roi[:, 0], return_inverse=True)[1]

    # convert reads to bin coverage
    cvg_lst = [list() for i in range(n_read)]
    for fi in range(frg_roi.shape[0]):
        bin_idx = np.where(hasOL(frg_roi[fi, 2:4], bin_bnd))[0]
        cvg_lst[frg_roi[fi, 0]].extend(bin_idx)
    cvg_lst = [np.unique(cvg_lst[i]) for i in range(n_read)]

    # looping over bins
    print 'Performing the MC analysis using {:d} reads ...'.format(n_read)
    mat_freq = np.full([n_bin, n_bin], fill_value=np.nan)
    mat_zscr = np.full([n_bin, n_bin], fill_value=np.nan)
    for bi in range(n_bin):
        if bi % (n_bin / 10) == 0:
            print '{:0.0f}%,'.format(bi * 100.0 / n_bin),
        is_pos = hasOL(bin_bnd[bi, :], frg_roi[:, 2:4])
        frg_pos = frg_roi[np.isin(frg_roi[:, 0], frg_roi[is_pos, 0]), :]
        frg_neg = frg_roi[~np.isin(frg_roi[:, 0], frg_pos[:, 0]), :]
        ids_pos = np.unique(frg_pos[:, 0])
        ids_neg = np.unique(frg_neg[:, 0])
        n_pos = len(ids_pos)
        n_neg = len(ids_neg)
        assert n_pos <= n_neg
        if n_pos < 100:
            continue

        # calculate the background
        rnd_freq = np.zeros([n_epoch, n_bin])
        for ei in np.arange(n_epoch):
            rnd_lst = np.random.choice(ids_neg, n_pos, replace=False)
            for rd_idx in rnd_lst:
                bin_cvg = cvg_lst[rd_idx]
                n_cvg = len(bin_cvg)
                rnd_freq[ei, bin_cvg] += 1
                rnd_freq[ei, bin_cvg[np.random.randint(n_cvg)]] -= 1

        # calculate observed
        for bj in range(bi+1, n_bin):
            is_cov = hasOL(bin_bnd[bj, :], frg_pos[:, 2:4])
            mat_freq[bi, bj] = len(np.unique(frg_pos[is_cov, 0]))

            zscr_avg = np.mean(rnd_freq[:, bj])
            zscr_std = np.std(rnd_freq[:, bj])
            if zscr_std == 0:
                continue
            mat_zscr[bi, bj] = (mat_freq[bi, bj] - zscr_avg) / zscr_std
            mat_zscr[bj, bi] = mat_zscr[bi, bj]

    # set vp bins to nan
    is_vp = hasOL([configs['vp_start'], configs['vp_end']], bin_bnd)
    mat_zscr[is_vp, :] = np.nan
    mat_zscr[:, is_vp] = np.nan
    vp_bnd = [bin_bnd[is_vp, 0][0], bin_bnd[is_vp, 1][-1]]

    # plotting
    plt.figure(figsize=(17, 9))
    clr_lst = ['#ff1a1a', '#ff8a8a', '#ffffff', '#ffffff', '#ffffff', '#8ab5ff', '#3900f5']
    clr_map = LinearSegmentedColormap.from_list('test', clr_lst, N=10)
    clr_map.set_bad('gray', 0.05)
    plt.imshow(mat_zscr, extent=x_lim + x_lim, cmap=clr_map, origin='bottom', interpolation='nearest')
    plt.gca().add_patch(patches.Rectangle([vp_bnd[0], x_lim[0]], vp_bnd[1] - vp_bnd[0], x_lim[1] - x_lim[0],
                                          linewidth=0, edgecolor='None', facecolor='orange'))
    plt.gca().add_patch(patches.Rectangle([x_lim[0], vp_bnd[0]], x_lim[1] - x_lim[0], vp_bnd[1] - vp_bnd[0],
                                          linewidth=0, edgecolor='None', facecolor='orange'))
    cbar_h = plt.colorbar()
    cbar_h.ax.tick_params(labelsize=14)
    plt.clim(-6, 6)

    # add annotations
    ant_pd = load_annotation(configs['genome_build'], roi_crd=[configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
    for ai in range(ant_pd.shape[0]):
        ant_pos = ant_pd.loc[ai, 'ant_pos']
        plt.text(ant_pos, x_lim[1], ant_pd.loc[ai, 'ant_name'],
                 horizontalalignment='left', verticalalignment='bottom', rotation=60)
        plt.text(x_lim[1], ant_pos, ' ' + ant_pd.loc[ai, 'ant_name'],
                 horizontalalignment='left', verticalalignment='center')
        plt.plot([ant_pos, ant_pos], x_lim, ':', color='#bfbfbf', linewidth=1, alpha=0.4)
        plt.plot(x_lim, [ant_pos, ant_pos], ':', color='#bfbfbf', linewidth=1, alpha=0.4)

    # final adjustments
    plt.xlim(x_lim)
    plt.ylim(x_lim)
    x_ticks = np.linspace(configs['roi_start'], configs['roi_end'], 7, dtype=np.int64)
    x_tick_label = ['{:0.2f}m'.format(x / 1e6) for x in x_ticks]
    plt.xticks(x_ticks, x_tick_label, rotation=0, horizontalalignment='center')
    plt.yticks(x_ticks, x_tick_label, rotation=0)
    plt.title('Multicontact matrix, {:s}\n'.format(configs['run_id']) +
              '#read (#roiFrg>{:d}, ex. vp)={:,d}\n\n\n'.format(min_n_frg - 1, n_read)
              )
    plt.savefig(configs['output_file'], bbox_inches='tight')


def perform_vpsoi_analysis(configs, soi_name, min_n_frg=2, n_perm=10):
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt, patches
    from matplotlib.colors import LinearSegmentedColormap

    from utilities import hasOL
    from mc4c_tools import load_mc4c, load_annotation

    # initialization
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/plt_vpsoi_' + configs['run_id'] + '.pdf'
    edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
    bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
    bin_cen = np.mean(bin_bnd, axis=1, dtype=np.int64)
    x_lim = [configs['roi_start'], configs['roi_end']]
    y_lim = [0, 10]
    # y_lim = plt.ylim()

    # load MC-HC data
    frg_dp = load_mc4c(configs, only_unique=True, only_valid=True, min_mq=20, reindex_reads=False)
    frg_np = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd']].values
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
    n_read = len(np.unique(frg_inf[:, 0]))

    # get soi info
    ant_pd = load_annotation(configs['genome_build'], roi_crd=roi_crd)
    is_in = np.where(np.isin(ant_pd['ant_name'], soi_name))[0]
    assert len(is_in) == 1
    soi_pd = ant_pd.loc[is_in[0], :]
    pos_crd = [soi_pd['ant_cnum'], soi_pd['ant_pos'] - 1500, soi_pd['ant_pos'] + 1500]

    # compute positive profile and backgrounds
    prf_pos, prf_rnd, frg_pos, frg_neg = compute_mc_associations(frg_inf, pos_crd, bin_bnd, n_perm=n_perm)
    n_pos = len(np.unique(frg_pos[:, 0]))
    nrm_pos = prf_pos * 100.0 / n_pos

    # compute scores
    nrm_rnd = prf_rnd * 100.0 / n_pos
    nrm_exp = np.mean(nrm_rnd, axis=0)
    nrm_std = np.std(nrm_rnd, axis=0, ddof=0)
    bin_scr = np.divide(nrm_pos - nrm_exp, nrm_std)

    # set vp bins to nan
    is_vp = hasOL([configs['vp_start'], configs['vp_end']], bin_bnd)
    bin_scr[is_vp] = np.nan
    vp_bnd = [bin_bnd[is_vp, 0][0], bin_bnd[is_vp, 1][-1]]

    # plotting
    plt.figure(figsize=(15, 4))
    plt.plot(bin_cen, nrm_pos, color='#5757ff', linewidth=1)
    plt.plot(bin_cen, nrm_exp, color='#cccccc', linewidth=1)
    plt.fill_between(bin_cen, nrm_exp - nrm_std, nrm_exp + nrm_std, color='#ebebeb', linewidth=0.2)


    # plt.subplot(1, 2, vi + 1)
    # clr_lst = ['#ff1a1a', '#ff8a8a', '#ffffff', '#ffffff', '#ffffff', '#8ab5ff', '#3900f5']
    # clr_map = LinearSegmentedColormap.from_list('test', clr_lst, N=10)
    # clr_map.set_bad('gray', 0.05)
    # plt.imshow(mat_zscr, extent=x_lim + x_lim, cmap=clr_map, origin='bottom', interpolation='nearest')
    plt.gca().add_patch(patches.Rectangle([vp_bnd[0], y_lim[0]], vp_bnd[1] - vp_bnd[0], y_lim[1] - y_lim[0],
                                          linewidth=0, edgecolor='None', facecolor='orange', zorder=10))
    plt.gca().add_patch(patches.Rectangle([pos_crd[1], y_lim[0]], pos_crd[2] - pos_crd[1], y_lim[1] - y_lim[0],
                                          linewidth=0, edgecolor='None', facecolor='green', zorder=10))
    # cbar_h = plt.colorbar()
    # cbar_h.ax.tick_params(labelsize=14)
    # plt.clim(-6, 6)

    # add annotations
    for ai in range(ant_pd.shape[0]):
        ant_pos = ant_pd.loc[ai, 'ant_pos']
        plt.text(ant_pos, y_lim[1], ant_pd.loc[ai, 'ant_name'],
                 horizontalalignment='center', verticalalignment='bottom', rotation=60)
        plt.plot([ant_pos, ant_pos], y_lim, ':', color='#bfbfbf', linewidth=1, alpha=0.4)

    # final adjustments
    plt.xlim(x_lim)
    plt.ylim(y_lim)
    x_ticks = np.linspace(configs['roi_start'], configs['roi_end'], 7, dtype=np.int64)
    y_ticks = plt.yticks()[0]
    x_tick_label = ['{:0.2f}m'.format(x / 1e6) for x in x_ticks]
    y_tick_label = ['{:0.0f}%'.format(y) for y in y_ticks]
    plt.xticks(x_ticks, x_tick_label, rotation=0, horizontalalignment='center')
    plt.yticks(y_ticks, y_tick_label, rotation=0, horizontalalignment='right')
    plt.ylabel('Percentage of reads')
    plt.title('VP-SOI for {:s}, in {:s}\n'.format(soi_name, configs['run_id']) +
              '#read (#roiFrg>{:d}, ex. vp)={:,d}\n'.format(min_n_frg - 1, n_read) +
              '#pos = {:d}\n\n\n'.format(n_pos)
              )
    # plt.subplots_adjust(top=1.4, wspace=0.5)
    plt.savefig(configs['output_file'], bbox_inches='tight')


def compute_mc_associations(frg_inf, pos_crd, bin_bnd, n_perm=1000):
    from utilities import hasOL, flatten

    # initialization
    n_bin = bin_bnd.shape[0]

    # re-index circles
    frg_inf[:, 0] = np.unique(frg_inf[:, 0], return_inverse=True)[1] + 1
    n_read = np.max(frg_inf[:, 0])

    # convert fragments to bin-coverage
    cfb_lst = [list() for i in range(n_read + 1)]
    n_frg = frg_inf.shape[0]
    for fi in range(n_frg):
        bin_idx = np.where(hasOL(frg_inf[fi, 2:4], bin_bnd))[0]
        cfb_lst[frg_inf[fi, 0]].append(list(bin_idx))

    # filter circles for (>1 bin cvg)
    valid_lst = []
    for rd_nid in range(1, n_read):
        fb_lst = cfb_lst[rd_nid]
        bin_lst = np.unique(flatten(fb_lst))
        if len(bin_lst) > 1:
            valid_lst.append(rd_nid)
    frg_inf = frg_inf[np.isin(frg_inf[:, 0], valid_lst), :]

    # select positive/negative circles
    is_pos = np.where(hasOL(pos_crd, frg_inf[:, 1:4]))[0]
    frg_pos = frg_inf[ np.isin(frg_inf[:, 0], frg_inf[is_pos, 0]), :]
    frg_neg = frg_inf[~np.isin(frg_inf[:, 0], frg_inf[is_pos, 0]), :]
    cfb_pos = [cfb_lst[i] for i in np.unique(frg_pos[:, 0])]
    cfb_neg = [cfb_lst[i] for i in np.unique(frg_neg[:, 0])]
    n_pos = len(np.unique(frg_pos[:, 0]))
    n_neg = len(np.unique(frg_neg[:, 0]))
    assert n_pos == len(cfb_pos)
    assert n_neg == len(cfb_neg)

    # make positive profile
    prf_pos = np.zeros(n_bin)
    for pi in range(n_pos):
        bin_lst = np.unique(flatten(cfb_pos[pi]))
        for bi in bin_lst:
            prf_pos[bi] += 1

    # make background profile from negative set
    prf_rnd = np.zeros([n_perm, n_bin])
    print 'Computing expected profile:'
    for ei in np.arange(n_perm):
        if ((ei + 1) % 200) == 0:
            print '\t{:d} randomized profiles are computed.'.format(ei + 1)
        rnd_lst = np.random.permutation(n_neg)[:n_pos]
        for rd_idx in rnd_lst:
            f2b_rnd = cfb_neg[rd_idx]
            n_frg = len(f2b_rnd)
            if n_frg > 1:
                frg_red = [f2b_rnd[i] for i in np.random.permutation(n_frg)[1:]]
                for i in np.unique(flatten(frg_red)):
                    prf_rnd[ei, i] += 1

    return prf_pos, prf_rnd, frg_pos, frg_neg






