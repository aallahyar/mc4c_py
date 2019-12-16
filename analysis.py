from __future__ import print_function
import numpy as np
from utilities import showprogress

MIN_N_POS = 100


def compute_mc_associations(frg_inf, pos_crd, bin_bnd, n_perm=1000, pos_ids=None, verbose=True):
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

    # select positive/negative circles
    if pos_ids is not None:
        assert len(pos_crd) == 0
        is_pos = np.isin(frg_inf[:, 0], pos_ids)
    else:
        is_pos = np.where(hasOL(pos_crd, frg_inf[:, 1:4]))[0]
    frg_pos = frg_inf[ np.isin(frg_inf[:, 0], frg_inf[is_pos, 0]), :]
    frg_neg = frg_inf[~np.isin(frg_inf[:, 0], frg_inf[is_pos, 0]), :]
    cfb_pos = [cfb_lst[i] for i in np.unique(frg_pos[:, 0])]
    cfb_neg = [cfb_lst[i] for i in np.unique(frg_neg[:, 0])]
    n_pos = len(cfb_pos)
    n_neg = len(cfb_neg)

    # make positive profile
    prf_pos = np.zeros(n_bin)
    for pi in range(n_pos):
        bin_lst = flatten(cfb_pos[pi])
        prf_pos[bin_lst] += 1

    # make background profile from negative set
    prf_rnd = np.zeros([n_perm, n_bin])
    neg_lst = range(n_neg)
    for ei in np.arange(n_perm):
        if verbose and (((ei + 1) % 200) == 0):
            print('\t{:d} randomized profiles are computed.'.format(ei + 1))
        np.random.shuffle(neg_lst)
        for rd_idx in neg_lst[:n_pos]:
            f2b_rnd = cfb_neg[rd_idx]
            np.random.shuffle(f2b_rnd)
            prf_rnd[ei, flatten(f2b_rnd[1:])] += 1  # making sure one element is randomly removed everytime

    return prf_pos, prf_rnd, frg_pos, frg_neg


def perform_mc_analysis(configs, min_n_frg=2):
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt, patches
    from matplotlib.colors import LinearSegmentedColormap

    from utilities import load_mc4c, load_annotation, hasOL

    # initialization
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/analysis_mcTest_' + configs['run_id'] + '.pdf'
    edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
    bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
    n_bin = bin_bnd.shape[0]
    n_epoch = 1000
    x_lim = [configs['roi_start'], configs['roi_end']]

    # load MC-HC data
    frg_dp = load_mc4c(configs, unique_only=True, valid_only=True, min_mq=20, reindex_reads=False)
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
    print('Performing the MC analysis using {:d} reads ...'.format(n_read))
    mat_freq = np.full([n_bin, n_bin], fill_value=np.nan)
    mat_zscr = np.full([n_bin, n_bin], fill_value=np.nan)
    for bi in range(n_bin):
        if bi % (n_bin / 10) == 0:
            print('{:0.0f}%,'.format(bi * 100.0 / n_bin), end='')
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


def perform_vpsoi_analysis(config_lst, soi_name, min_n_frg=2, n_perm=1000):
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
        configs['output_file'] = configs['output_dir'] + '/analysis_atVP-SOI_{:s}_{:s}.pdf'.format(run_id, soi_name)
    edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
    bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
    bin_cen = np.mean(bin_bnd, axis=1, dtype=np.int64)
    bin_w = bin_bnd[0, 1] - bin_bnd[0, 0]
    x_lim = [configs['roi_start'], configs['roi_end']]
    y_lim = [0, 10]

    # load MC-HC data
    frg_dp = load_mc4c(config_lst, unique_only=True, valid_only=True, min_mq=20, reindex_reads=True, verbose=True)
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
    soi_pd = ant_pd.loc[is_in[0], :]
    soi_crd = [soi_pd['ant_cnum'], soi_pd['ant_pos'] - int(bin_w * 1.5), soi_pd['ant_pos'] + int(bin_w * 1.5)]
    if hasOL(soi_crd, vp_crd)[0]:
        print('[w] Selected SOI coordinate overlaps with the view point. Ignoring the analysis')
        return

    # compute positive profile and backgrounds
    print('Computing expected profile for bins:')
    prf_frq, prf_rnd, frg_pos, frg_neg = compute_mc_associations(frg_inf, soi_crd, bin_bnd, n_perm=n_perm)
    n_pos = len(np.unique(frg_pos[:, 0]))
    prf_obs = prf_frq * 100.0 / n_pos
    print('{:,d} reads are found to cover '.format(n_pos) + \
          '{:s} area ({:s}:{:d}-{:d})'.format(soi_pd['ant_name'], soi_pd['ant_chr'], soi_crd[1], soi_crd[2]))

    # check enough #pos
    if n_pos < MIN_N_POS:
        print('[w] #reads in the positive set is insufficient (n={:d}, required >{:d})'.format(n_pos, MIN_N_POS))
        print('Analysis is ignored ...')
        return

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
    ant_obs, soi_rnd = compute_mc_associations(frg_inf, soi_crd, ant_bnd, n_perm=n_perm)[:2]
    ant_exp = np.mean(soi_rnd, axis=0)
    ant_std = np.std(soi_rnd, axis=0, ddof=0)
    np.seterr(all='ignore')
    ant_scr = np.divide(ant_obs - ant_exp, ant_std)
    np.seterr(all=None)

    # set vp score to nan
    is_vp = hasOL(vp_bnd, ant_bnd)
    is_soi = hasOL(soi_crd[1:3], ant_bnd)
    ant_scr[is_vp | is_soi] = np.nan

    # plotting
    fig = plt.figure(figsize=(15, 3))
    ax_prf = plt.subplot2grid((20, 40), (0, 0), rowspan=19, colspan=39)
    ax_cmp = plt.subplot2grid((20, 40), (0, 39), rowspan=10, colspan=1)
    ax_scr = plt.subplot2grid((20, 40), (19, 0), rowspan=1, colspan=39)

    # set up colorbar
    c_lim = [-6, 6]
    clr_lst = ['#ff1a1a', '#ff7575', '#ffcccc', '#ffffff', '#ffffff', '#ffffff', '#ccdfff', '#3d84ff', '#3900f5']
    clr_map = LinearSegmentedColormap.from_list('test', clr_lst, N=9)
    clr_map.set_bad('gray', 0.05)
    norm = matplotlib.colors.Normalize(vmin=c_lim[0], vmax=c_lim[1])
    cbar_h = matplotlib.colorbar.ColorbarBase(ax_cmp, cmap=clr_map, norm=norm)
    # cbar_h.ax.tick_params(labelsize=12)
    cbar_h.ax.set_ylabel('z-score', rotation=90)

    # profile plot
    ax_prf.plot(bin_cen, prf_obs, color='#5757ff', linewidth=1)
    ax_prf.plot(bin_cen, prf_exp, color='#cccccc', linewidth=1)
    ax_prf.fill_between(bin_cen, prf_exp - prf_std, prf_exp + prf_std, color='#ebebeb', linewidth=0.2)

    ax_prf.add_patch(patches.Rectangle([vp_bnd[0], y_lim[0]], vp_bnd[1] - vp_bnd[0], y_lim[1] - y_lim[0],
                                       edgecolor='None', facecolor='orange', zorder=100))
    ax_prf.add_patch(patches.Rectangle([soi_crd[1], y_lim[0]], soi_crd[2] - soi_crd[1], y_lim[1] - y_lim[0],
                                       edgecolor='None', facecolor='green', zorder=100))
    ax_prf.set_xlim(x_lim)
    ax_prf.set_ylim(y_lim)
    ax_prf.set_xticks([])

    # add score plot
    ax_scr.imshow(bin_scr.reshape(1, -1), extent=x_lim + [-500, 500], cmap=clr_map,
                  vmin=c_lim[0], vmax=c_lim[1], interpolation='nearest')
    ax_scr.set_xlim(x_lim)
    ax_scr.set_yticks([])

    # add annotations
    for ai in range(n_ant):
        ax_prf.text(ant_pos[ai], y_lim[1], ant_pd.loc[ai, 'ant_name'],
                    horizontalalignment='center', verticalalignment='bottom', rotation=60)
        ax_prf.plot(ant_pos[[ai, ai]], y_lim, ':', color='#bfbfbf', linewidth=1, alpha=0.4)

        if not np.isnan(ant_scr[ai]):
            ax_prf.add_patch(patches.Rectangle([ant_bnd[ai, 0], y_lim[1]-0.15], ant_bnd[ai, 1] - ant_bnd[ai, 0], 0.15,
                                               edgecolor='None', facecolor=clr_map(ant_scr[ai]), zorder=10))
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
    ax_prf.set_title('VP-SOI from {:s}, using as SOI {:s}\n'.format(run_id, soi_name) +
                     '#read (#roiFrg>{:d}, ex. vp)={:,d}, #pos={:d}\n'.format(min_n_frg - 1, n_read, n_pos) +
                     'bin-w={:0.0f}; soi-w={:0.0f}; #perm={:d}\n\n\n'.format(bin_w, np.diff(ant_bnd[0]), n_perm)
                     )
    plt.savefig(configs['output_file'], bbox_inches='tight')


def perform_soisoi_analysis(config_lst, min_n_frg=2, n_perm=1000):
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
        config_lst[0]['output_file'] = config_lst[0]['output_dir'] + '/analysis_atSOI-SOI_{:s}.pdf'.format(run_id)
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
        ant_obs, soi_rnd, frg_pos = compute_mc_associations(frg_inf, soi_crd, ant_bnd, n_perm=n_perm)[:3]
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


def perform_at_across_roi(config_lst, min_n_frg=2, n_perm=1000, downsample=None, xls_export=False):
    import platform
    import matplotlib
    if platform.system() == 'Linux':
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt, patches
    from matplotlib.colors import LinearSegmentedColormap

    from utilities import load_mc4c, load_annotation, hasOL, flatten, limit_to_roi

    # initialization
    run_id = ','.join([config['run_id'] for config in config_lst])
    configs = config_lst[0]
    if configs['output_file'] is None:
        roi_w = configs['roi_end'] - configs['roi_start']
        if downsample:
            run_id += '_ds{:d}'.format(downsample)
        configs['output_file'] = configs['output_dir'] + \
                                 '/analysis_atAcrossROI_{:s}_'.format(run_id) + \
                                 'rw{:0.1f}kb_np{:0.1f}k_'.format(roi_w / 1e3, n_perm / 1e3) + \
                                 'zlm{:0.0f},{:0.0f}.pdf'.format(*configs['zscr_lim'])

    # create bin list
    edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
    bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
    bin_w = bin_bnd[0, 1] - bin_bnd[0, 0]
    n_bin = bin_bnd.shape[0]

    # make block list
    bin_cen = np.mean(bin_bnd, axis=1, dtype=np.int64).reshape(-1, 1)
    # blk_crd = np.hstack([np.repeat(configs['vp_cnum'], n_bin / 3).reshape(-1, 1), edge_lst[:-3:3], edge_lst[3::3] - 1])
    blk_crd = np.hstack([np.repeat(configs['vp_cnum'], n_bin).reshape(-1, 1), bin_cen - int(bin_w * 1.5), bin_cen + int(bin_w * 1.5) - 1])
    blk_w = blk_crd[0, 2] - blk_crd[0, 1]
    n_blk = blk_crd.shape[0]
    del edge_lst

    # define areas
    roi_cen = np.mean([np.min(configs['prm_start']), np.max(configs['prm_end'])], dtype=np.int)
    vp_crd = np.array([configs['vp_cnum'], roi_cen - int(bin_w * 1.5), roi_cen + int(bin_w * 1.5)])
    roi_crd = [configs['vp_cnum'], configs['roi_start'], configs['roi_end']]

    # load MC-HC data
    frg_dp = load_mc4c(config_lst, unique_only=True, valid_only=True, min_mq=20, reindex_reads=True, verbose=True)
    read_all = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd']].values
    del frg_dp

    # select >2 roi-fragments
    read_inf = limit_to_roi(read_all[:, :4], vp_crd=vp_crd, roi_crd=roi_crd, min_n_frg=min_n_frg)
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
    'Selecting only reads with >1 bins covered'
    valid_lst = []
    for rd_nid in range(1, n_read + 1):
        fb_lst = cfb_lst[rd_nid]
        bin_cvg = np.unique(flatten(fb_lst))
        if len(bin_cvg) > 1:
            valid_lst.append(rd_nid)
    read_inf = read_inf[np.isin(read_inf[:, 0], valid_lst), :]

    # subsample reads
    if downsample:
        n_read = len(np.unique(read_inf[:, 0]))
        print('[i] Downsampling {:,d} informative reads to {:d} reads.'.format(n_read, downsample))
        rnd_ids = np.random.choice(np.unique(read_inf[:, 0]), downsample, replace=False)
        read_inf = read_inf[np.isin(read_inf[:, 0], rnd_ids), :]

    # reindexing reads
    read_inf[:, 0] = np.unique(read_inf[:, 0], return_inverse=True)[1] + 1
    n_read = np.max(read_inf[:, 0])
    print('{:,d} reads are left after bin-coverage filter.'.format(n_read))

    # get soi info
    ant_pd = load_annotation(configs['genome_build'], roi_crd=roi_crd)
    ant_bnd = np.hstack([ant_pd[['ant_pos']].values, ant_pd[['ant_pos']].values])

    # compute score for annotations
    print('Computing expected profile for {:d} blocks (required coverage: {:d} reads):'.format(n_blk, MIN_N_POS))
    blk_scr = np.full([n_blk, n_blk], fill_value=np.nan)
    # x_tick_lbl = [' '] * n_blk
    y_tick_lbl = [' '] * n_blk
    n_ignored = 0
    for bi in range(n_blk):
        showprogress(bi, n_blk, n_step=20)

        # add axes labels
        ant_idx = np.where(hasOL(blk_crd[bi, 1:], ant_bnd, offset=0))[0]
        if len(ant_idx) > 0:
            ant_name = ','.join([ant_pd.loc[i, 'ant_name'] for i in ant_idx])
            # x_tick_lbl[bi] = ('{:s}, #{:0.0f}'.format(ant_name, n_pos))
            y_tick_lbl[bi] = ant_name
        # else:
            # x_tick_lbl[bi] = ('#{:0.0f}'.format(n_pos))

        # ignore if vp
        if hasOL(blk_crd[bi, :], vp_crd, offset=blk_w)[0]:
            continue

        # compute the observe and background
        blk_obs, blk_rnd, read_pos = compute_mc_associations(read_inf, blk_crd[bi, :], blk_crd[:, 1:],
                                                             n_perm=n_perm, verbose=False)[:3]
        n_pos = len(np.unique(read_pos[:, 0]))
        if n_pos < MIN_N_POS:
            n_ignored += 1
            continue

        # compute the scores
        blk_exp = np.mean(blk_rnd, axis=0)
        blk_std = np.std(blk_rnd, axis=0, ddof=0)
        np.seterr(all='ignore')
        blk_scr[:, bi] = np.divide(blk_obs - blk_exp, blk_std)
        np.seterr(all=None)

        # remove scores overlapping with positive set
        is_nei = hasOL(blk_crd[bi, 1:], blk_crd[:, 1:], offset=blk_w)
        blk_scr[is_nei, bi] = np.nan

    if n_ignored != 0:
        print('[w] {:d}/{:d} blocks are ignored due to low coverage.'.format(n_ignored, n_blk))

    # set self scores to nan
    # np.fill_diagonal(blk_scr, val=np.nan)

    # clean up tick labels

    # plotting the scores
    plt.figure(figsize=(15, 13))
    ax_scr = plt.subplot2grid((40, 40), (0, 0), rowspan=39, colspan=39)
    ax_cmp = plt.subplot2grid((40, 40), (0, 39), rowspan=20, colspan=1)

    # set up color bar
    clr_lst = ['#ff1a1a', '#ff7575', '#ffcccc', '#ffffff', '#ffffff', '#ffffff', '#ccdfff', '#3d84ff', '#3900f5']
    clr_map = LinearSegmentedColormap.from_list('test', clr_lst, N=9)
    clr_map.set_bad('gray', 0.1)
    norm = matplotlib.colors.Normalize(vmin=configs['zscr_lim'][0], vmax=configs['zscr_lim'][1])
    cbar_h = matplotlib.colorbar.ColorbarBase(ax_cmp, cmap=clr_map, norm=norm)
    # cbar_h.ax.tick_params(labelsize=12)
    cbar_h.ax.set_ylabel('z-score', rotation=90)
    cbar_edge = np.round(cbar_h.cmap(norm(configs['zscr_lim'])), decimals=2)

    # add score scatter matrix
    x_lim = [0, n_blk]
    ax_scr.imshow(blk_scr, extent=x_lim + x_lim, cmap=clr_map,
                  vmin=configs['zscr_lim'][0], vmax=configs['zscr_lim'][1], interpolation='nearest', origin='lower')
    ax_scr.set_xlim(x_lim)
    ax_scr.set_ylim(x_lim)
    ax_scr.invert_yaxis()

    # add vp patches
    vp_idx = np.where(hasOL(vp_crd, blk_crd, offset=blk_w))[0]
    ax_scr.add_patch(patches.Rectangle([0, vp_idx[0]], n_blk, vp_idx[-1] - vp_idx[0],
                                       linewidth=0, edgecolor='None', facecolor='orange'))
    ax_scr.add_patch(patches.Rectangle([vp_idx[0], 0], vp_idx[-1] - vp_idx[0], n_blk,
                                       linewidth=0, edgecolor='None', facecolor='orange'))

    # add score values to each box
    # for bi in range(n_blk):
    #     for bj in range(n_blk):
    #         if np.isnan(blk_scr[bi, bj]):
    #             continue
    #         ant_clr = np.round(img_h.cmap(img_h.norm(blk_scr[bi, bj])), decimals=2)
    #         if np.array_equal(ant_clr, cbar_edge[0]) or np.array_equal(ant_clr, cbar_edge[1]):
    #             txt_clr = '#ffffff'
    #         else:
    #             txt_clr = '#000000'
    #         ax_scr.text(bj + 0.5, bi + 0.5, '{:+0.1f}'.format(blk_scr[bi, bj]), color=txt_clr,
    #                     horizontalalignment='center', verticalalignment='center', fontsize=12)

    # adjust ticks
    for lbl in np.unique(y_tick_lbl):
        if lbl == ' ':
            continue
        idx_lst = np.where(np.isin(y_tick_lbl, lbl))[0]
        if len(idx_lst) > 1:
            kpt_idx = np.mean(idx_lst, dtype=np.int)
            for idx in idx_lst:
                y_tick_lbl[idx] = 'l'
            y_tick_lbl[kpt_idx] = lbl + ' '

    # final adjustments
    ax_scr.set_xticks(np.arange(n_blk) + 0.5)
    ax_scr.set_yticks(np.arange(n_blk) + 0.5)
    ax_scr.set_xticklabels(y_tick_lbl, rotation=90)
    ax_scr.set_yticklabels(y_tick_lbl)
    ax_scr.set_xlabel('Selected SOIs')
    ax_scr.tick_params(length=0)
    ax_scr.set_title('Association matrix from {:s}\n'.format(run_id) +
                     '#read (#roiFrg>{:d}, ex. vp)={:,d}, '.format(min_n_frg - 1, n_read) +
                     'bin-w={:0.0f}; block-w={:0.0f}; #perm={:d}'.format(bin_w, blk_w, n_perm)
                     )
    plt.savefig(configs['output_file'], bbox_inches='tight')
    plt.close()
    print('ROI-ROI z-scores are plotted in {:s}'.format(configs['output_file']))

    # export to excel file
    if xls_export:
        xls_fname = configs['output_file'][:-4] + '.xlsx'
        print('Exporting z-scores to excel sheet: {:s}'.format(xls_fname))

        import pandas as pd
        zscr_pd = pd.DataFrame(blk_scr, columns=[bin_cen.flatten(), y_tick_lbl], index=[bin_cen.flatten(), y_tick_lbl])
        zscr_pd.to_excel(xls_fname, sheet_name='z-scores')



