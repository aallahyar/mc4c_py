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

    # looping over bins
    print 'Performing the MC analysis using {:d} reads ...'.format(n_read)
    mat_freq = np.full([n_bin, n_bin], fill_value=np.nan)
    mat_zscr = np.full([n_bin, n_bin], fill_value=np.nan)
    vec_freq = np.zeros(n_bin)
    for bi in range(n_bin):
        is_pos = hasOL(bin_bnd[bi, :], frg_roi[:, 2:4])
        frg_pos = frg_roi[np.isin(frg_roi[:, 0], frg_roi[is_pos, 0]), :]
        frg_neg = frg_roi[~np.isin(frg_roi[:, 0], frg_pos[:, 0]), :]
        ids_pos = np.unique(frg_pos[:, 0])
        ids_neg = np.unique(frg_neg[:, 0])
        n_pos = len(ids_pos)
        n_neg = len(ids_neg)
        assert n_pos <= n_neg

        for bj in range(bi+1, n_bin):
            is_cov = hasOL(bin_bnd[bj, :], frg_pos[:, 2:4])
            mat_freq[bi, bj] = len(np.unique(frg_pos[is_cov, 0]))
            if mat_freq[bi, bj] < 5:
                continue

            for ei in range(n_epoch):
                ids_rnd = np.random.choice(ids_neg, n_pos, replace=False)
                frg_rnd = frg_neg[np.isin(frg_neg[:, 0], ids_rnd), :]
                is_cov = hasOL(bin_bnd[bj, :], frg_rnd[:, 2:4])
                vec_freq[bj] = len(np.unique(frg_rnd[is_cov, 0]))

            zscr_avg = np.mean(vec_freq)
            zscr_std = np.std(vec_freq)
            mat_zscr[bi, bj] = (mat_freq[bi, bj] - zscr_avg) / zscr_std

    # set vp bins to nan
    is_vp = hasOL([configs['vp_start'], configs['vp_end']], bin_bnd)
    mat_zscr[is_vp, :] = np.nan
    mat_zscr[:, is_vp] = np.nan
    vp_bnd = [bin_bnd[is_vp, 0][0], bin_bnd[is_vp, 1][-1]]

    # plotting
    plt.figure(figsize=(15, 5))
    clr_lst = ['#3900f5', '#8ab5ff', '#ffffff', '#ffffff', '#ffffff', '#ff8a8a', '#ff1a1a']
    clr_map = LinearSegmentedColormap.from_list('test', clr_lst, N=10)
    plt.imshow(mat_zscr, extent=x_lim + x_lim, cmap=clr_map, origin='upper', interpolation='nearest')
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
                 horizontalalignment='center', verticalalignment='bottom')
        plt.text(x_lim[1], ant_pos, ' ' + ant_pd.loc[ai, 'ant_name'],
                 horizontalalignment='left', verticalalignment='center')
        plt.plot([ant_pos, ant_pos], x_lim, ':', color='#bfbfbf', linewidth=1, alpha=0.4)
        plt.plot(x_lim, [ant_pos, ant_pos], ':', color='#bfbfbf', linewidth=1, alpha=0.4)

    # final adjustments
    plt.xlim(x_lim)
    plt.ylim(x_lim)
    x_ticks = np.linspace(configs['roi_start'], configs['roi_end'], 7, dtype=np.int64)
    x_tick_label = ['{:0.2f}m'.format(x / 1e6) for x in x_ticks]
    plt.xticks(x_ticks, x_tick_label, rotation=20)
    plt.yticks(x_ticks, x_tick_label, rotation=0)
    plt.title('Multicontact matrix, {:s}\n'.format(configs['run_id']) +
              '#read (#roiFrg>{:d}, ex. vp)={:,d}\n'.format(min_n_frg - 1, n_read)
              )
    plt.savefig(configs['output_file'], bbox_inches='tight')

