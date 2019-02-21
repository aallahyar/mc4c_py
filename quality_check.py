import numpy as np


def plot_readSizeDistribution(configs):
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import pysam

    # initializations
    if configs['input_file'] is None:
        configs['input_file'] = './reads/rd_' + configs['run_id'] + '.fasta.gz'
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/plt_ReadSizeDistribution_' + configs['run_id'] + '.pdf'
    MAX_SIZE = 8000
    edge_lst = np.linspace(0, MAX_SIZE, 81)
    n_bin = len(edge_lst) - 1

    # loop over reads
    nbp_total = 0
    nrd_all = 0
    nbp_lrg = 0
    nrd_lrg = 0
    size_dist = np.zeros(n_bin, dtype=np.int64)
    print 'Calculating read size distribution for: {:s}'.format(configs['input_file'])
    with pysam.FastxFile(configs['input_file']) as gz_fid:
        for rd_idx, read in enumerate(gz_fid):
            if rd_idx % 50000 == 0:
                print('{:,d} reads are processed.'.format(rd_idx))
            seq_size = len(read.sequence)

            nrd_all += 1
            nbp_total = nbp_total + seq_size
            if seq_size > 1500:
                nrd_lrg += 1
                nbp_lrg = nbp_lrg + seq_size

            if seq_size >= MAX_SIZE:
                seq_size = MAX_SIZE - 1
            bin_idx = np.digitize(seq_size, edge_lst) - 1
            size_dist[bin_idx] += 1

    # plotting
    plt.figure(figsize=(7, 5))
    plt.bar(range(n_bin), size_dist, width=0.95)
    x_ticks_idx = range(0, n_bin, 5)
    plt.xticks(x_ticks_idx, ['{:0.1f}'.format(edge_lst[i] / 1e3) for i in x_ticks_idx],
               rotation=0, fontsize=10)
    plt.xlim([-1, n_bin])
    plt.xlabel('#base pairs (kbp)')

    y_ticks = plt.yticks()[0]
    y_tick_lbl = ['{:0,.0f}k'.format(x / 1e3) for x in y_ticks]
    plt.yticks(y_ticks, y_tick_lbl)
    plt.ylabel('#reads')

    plt.title('Read size distribution, {:s}\n'.format(configs['run_id']) +
              '#read={:,d}; #read (>1.5kb)={:,d}\n'.format(nrd_all, nrd_lrg) +
              '#bases={:,d}; #bases (>1.5kb)={:,d}'.format(nbp_total, nbp_lrg)
              )

    plt.savefig(configs['output_file'], bbox_inches='tight')


def plot_frg_size_distribution(configs):
    import numpy as np
    import gzip
    from os import path
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    from utilities import load_mc4c, get_chr_info, get_re_info

    # initialization
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/plt_frgSizeDistribution_' + configs['run_id'] + '.pdf'
    MAX_SIZE = 1500
    edge_lst = np.linspace(0, MAX_SIZE, 31)
    n_bin = len(edge_lst) - 1

    # get chr info
    chr_lst = get_chr_info(genome_str=configs['genome_build'], property='chr_name')
    n_chr = len(chr_lst)

    # Read res-enz positions
    re_fname = './renzs/{:s}_{:s}.npz'.format(configs['genome_build'], '-'.join(configs['re_name']))
    print 'Loading RE positions from: {:s}'.format(re_fname)
    if not path.isfile(re_fname):
        from utilities import extract_re_positions
        extract_re_positions(genome_str=configs['genome_build'], re_name_lst=configs['re_name'],
                             output_fname=re_fname, ref_fasta=configs['reference_fasta'])
    re_pos_lst = get_re_info(re_name='-'.join(configs['re_name']), property='pos', genome_str=configs['genome_build'])

    # compute ref fragment size
    # lcl_area = [np.min(configs['prm_start']) - 5e6, np.max(configs['prm_end']) + 5e6]
    chr_idx = np.where(np.isin(chr_lst, configs['vp_chr']))[0]
    assert len(chr_idx) == 1
    re_pos = re_pos_lst[chr_idx[0]]
    # re_lcl = re_pos[(re_pos > lcl_area[0]) & (re_pos < lcl_area[1])]
    re_size = np.diff(re_pos) + 1
    n_re_ref = len(re_size)
    re_size[re_size >= MAX_SIZE] = MAX_SIZE - 1
    bin_idx = np.digitize(re_size, edge_lst) - 1
    dist_ref = np.bincount(bin_idx, minlength=n_bin)

    # Load MC-HC data
    frg_dp = load_mc4c(configs, min_mq=0, reindex_reads=False, unique_only=False, valid_only=True)
    frg_np = frg_dp[['Chr', 'MapStart', 'MapEnd', 'MQ']].values
    del frg_dp
    print 'Total of {:,d} mapped fragments are loaded:'.format(frg_np.shape[0])

    # calculate mapped fragment sizes
    frg_size = frg_np[:, 2] - frg_np[:, 1] + 1
    is_mq01 = frg_np[:, 3] >= 1
    is_mq20 = frg_np[:, 3] >= 20
    n_bp_mq01 = np.sum(frg_size[is_mq01])
    n_bp_mq20 = np.sum(frg_size[is_mq20])
    n_frg_mq01 = np.sum(is_mq01)
    n_frg_mq20 = np.sum(is_mq20)
    del frg_np

    # calculate fragment size distribution
    frg_size[frg_size >= MAX_SIZE] = MAX_SIZE - 1
    bin_idx = np.digitize(frg_size, edge_lst) - 1
    dist_mq01 = np.bincount(bin_idx, minlength=n_bin)
    bin_idx = np.digitize(frg_size[is_mq20], edge_lst) - 1
    dist_mq20 = np.bincount(bin_idx, minlength=n_bin)
    # del frg_size

    # calculate raw fragment size
    frg_fname = './fragments/frg_{:s}.fasta.gz'.format(configs['run_id'])
    print 'Scanning raw fragments in {:s}'.format(frg_fname)
    dist_mq00 = np.zeros(n_bin, dtype=np.int64)
    n_bp_mq00 = 0
    n_frg_mq00 = 0
    # seq_size_lst = []
    with gzip.open(frg_fname, 'r') as splt_fid:
        while True:
            frg_sid = splt_fid.readline()
            frg_seq = splt_fid.readline().rstrip('\n')
            if frg_sid == '':
                break
            if n_frg_mq00 % 250000 == 0:
                print('{:,d} fragments are processed.'.format(n_frg_mq00))

            seq_size = len(frg_seq)
            n_bp_mq00 += seq_size
            # seq_size_lst.append(seq_size)
            if seq_size >= MAX_SIZE:
                seq_size = MAX_SIZE - 1

            bin_idx = np.digitize(seq_size, edge_lst) - 1
            dist_mq00[bin_idx] += 1
            n_frg_mq00 += 1

    # Plotting
    plt.figure(figsize=(7, 5))
    ref_h = plt.bar(range(n_bin), dist_ref * float(n_frg_mq00) / np.sum(dist_ref), width=1.00, color='#dddddd')
    q00_h = plt.bar(range(n_bin), dist_mq00, width=0.90, color='#c9ae18', alpha=0.5)
    q01_h = plt.bar(range(n_bin), dist_mq01, width=0.70, color='#4766ff', alpha=1.0)
    q20_h = plt.bar(range(n_bin), dist_mq20, width=0.50, color='#1be600')
    plt.legend([ref_h, q00_h, q01_h, q20_h], [
        'Ref (#frg={:0,.0f}k), normed'.format(n_re_ref / 1e3),
        'Raw (#frg={:0,.0f}k)'.format(n_frg_mq00 / 1e3),
        'MQ1', 'MQ20'])
    plt.xlim([-1, n_bin + 1])
    x_ticks_idx = range(1, n_bin, 2)
    plt.xticks(x_ticks_idx, ['{:0.0f}'.format(edge_lst[i + 1]) for i in x_ticks_idx], rotation=35)
    plt.xlabel('#base pairs')

    y_ticks = plt.yticks()[0]
    y_tick_lbl = ['{:0,.0f}k'.format(y / 1e3) for y in y_ticks]
    plt.yticks(y_ticks, y_tick_lbl)
    plt.ylabel('#Fragments')

    plt.title('Fragment size distribution, {:s}\n'.format(configs['run_id']) +
              '#bp mapped ' +
              'MQ1={:0,.1f}m ({:0.1f}%); '.format(n_bp_mq01 / 1e6, n_bp_mq01 * 1e2 / n_bp_mq00) +
              'MQ20={:0,.1f}m ({:0.1f}%)\n'.format(n_bp_mq20 / 1e6, n_bp_mq20 * 1e2 / n_bp_mq00) +
              '#frg mapped ' +
              'MQ1={:0,.1f}k ({:0.1f}%); '.format(n_frg_mq01 / 1e3, n_frg_mq01 * 1e2 / n_frg_mq00) +
              'MQ20={:0,.1f}k ({:0.1f}%)'.format(n_frg_mq20 / 1e3, n_frg_mq20 * 1e2 / n_frg_mq00))
    plt.savefig(configs['output_file'], bbox_inches='tight')

    # from scipy.stats import gamma
    # shape, loc, scale = gamma.fit(seq_size_lst, floc=0)
    # print(shape, loc, scale)
    # x = np.linspace(50, MAX_SIZE, 30)
    # y = gamma.pdf(x, shape, loc, scale)
    # raw_h = plt.plot(range(n_bin), y * 2e6)
    # # ln_h[0].remove()
    # 
    # # shape, loc, scale = gamma.fit(frg_size, floc=0)
    # k = gamma.pdf(x, shape + 1, loc, scale)
    # map_h = plt.plot(range(n_bin), k * 1e6)
    # # ln_h[0].remove()
    #
    # raw_h[0].remove()
    # map_h[0].remove()



def plot_chrCvg(configs):
    import numpy as np
    import gzip
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    from utilities import load_mc4c, get_chr_info

    # initialization
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/plt_chrCoverage_' + configs['run_id'] + '.pdf'

    # get chr info
    chr_lst = get_chr_info(genome_str=configs['genome_build'], property='chr_name')
    n_chr = len(chr_lst)

    # Load MC-HC data
    frg_dp = load_mc4c(configs, min_mq=0, reindex_reads=False, unique_only=False, valid_only=True)
    frg_np = frg_dp[['Chr', 'MapStart', 'MapEnd', 'MQ']].values
    del frg_dp
    print 'Total of {:,d} mapped fragments are loaded:'.format(frg_np.shape[0])

    # calculate chromosome coverage
    is_mq01 = frg_np[:, 3] >= 1
    is_mq20 = frg_np[:, 3] >= 20
    ccvg_mq01 = np.bincount(frg_np[is_mq01, 0] - 1, minlength=n_chr)
    ccvg_mq20 = np.bincount(frg_np[is_mq20, 0] - 1, minlength=n_chr)

    frg_size = frg_np[:, 2] - frg_np[:, 1] + 1
    n_bp_mq01 = np.sum(frg_size[is_mq01])
    n_bp_mq20 = np.sum(frg_size[is_mq20])
    n_frg_mq01 = np.sum(is_mq01)
    n_frg_mq20 = np.sum(is_mq20)
    del frg_np

    # calculate raw fragment size
    frg_fname = './fragments/frg_{:s}.fasta.gz'.format(configs['run_id'])
    print 'Scanning raw fragments in {:s}'.format(frg_fname)
    n_bp_mq00 = 0
    n_frg_mq00 = 0
    with gzip.open(frg_fname, 'r') as splt_fid:
        while True:
            frg_sid = splt_fid.readline()
            frg_seq = splt_fid.readline().rstrip('\n')
            if frg_sid == '':
                break
            if n_frg_mq00 % 250000 == 0:
                print('{:,d} fragments are processed.'.format(n_frg_mq00))

            n_bp_mq00 += len(frg_seq)
            n_frg_mq00 += 1

    # Plotting
    plt.figure(figsize=(7, 5))
    q01_h = plt.bar(range(n_chr), ccvg_mq01, color='#bc85ff', width=0.80, alpha=0.7)
    q20_h = plt.bar(range(n_chr), ccvg_mq20, color='#8929ff', width=0.60, alpha=0.7)
    plt.legend([q01_h, q20_h], ['MQ>=1', 'MQ>=20'])
    plt.xticks(range(n_chr), chr_lst, rotation=35, ha='right')
    plt.xlim([-0.5, n_chr - 0.5])
    y_ticks = plt.yticks()[0]
    y_tick_lbl = ['{:0,.1f}k'.format(y / 1e3) for y in y_ticks]
    plt.yticks(y_ticks, y_tick_lbl)
    plt.ylabel('#Fragments')
    plt.title('Chromosome coverage, {:s}\n'.format(configs['run_id']) +
              '#bp mapped ' +
              'MQ1={:0,.1f}m ({:0.1f}%); '.format(n_bp_mq01 / 1e6, n_bp_mq01 * 1e2 / n_bp_mq00) +
              'MQ20={:0,.1f}m ({:0.1f}%)\n'.format(n_bp_mq20 / 1e6, n_bp_mq20 * 1e2 / n_bp_mq00) +
              '#frg mapped ' +
              'MQ1={:0,.1f}k ({:0.1f}%); '.format(n_frg_mq01 / 1e3, n_frg_mq01 * 1e2 / n_frg_mq00) +
              'MQ20={:0,.1f}k ({:0.1f}%)'.format(n_frg_mq20 / 1e3, n_frg_mq20 * 1e2 / n_frg_mq00))

    plt.savefig(configs['output_file'], bbox_inches='tight')


def plot_cirSizeDistribution(configs, roi_only=True, uniq_only=True):
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt, cm

    from utilities import accum_array, load_mc4c

    # initialization
    MAX_SIZE = 8
    edge_lst = np.linspace(1, MAX_SIZE, num=MAX_SIZE)
    n_edge = len(edge_lst)

    # Load MC-HC data
    frg_dp = load_mc4c(configs, min_mq=20, reindex_reads=True, unique_only=uniq_only)
    frg_np = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'MQ', 'ReadLength']].values
    del frg_dp

    # select requested fragments
    if uniq_only:
        filter_lst = ['uniq']
    else:
        filter_lst = []
    if roi_only:
        from utilities import hasOL
        vp_crd = np.array([configs['vp_cnum'], configs['vp_start'], configs['vp_end']])
        roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
        is_vp = hasOL(vp_crd, frg_np[:, 1:4], offset=0)
        is_roi = hasOL(roi_crd, frg_np[:, 1:4], offset=0)
        frg_np = frg_np[~is_vp & is_roi, :]
        filter_lst += ['roi', 'ex.vp']

    # group circles
    read_grp = accum_array(frg_np[:, 0] - 1, frg_np, rebuild_index=True)
    n_grp = len(read_grp)

    # Looping over circles
    size_dist = np.zeros([4, n_edge], dtype=np.int64)
    print 'Computing circle size from {:d} reads:'.format(n_grp)
    for read_idx, frg_set in enumerate(read_grp):
        if read_idx % 50000 == 0:
            print('\t{:,d}/{:,d} Reads are processed.'.format(read_idx, n_grp))
        n_frg = frg_set.shape[0]
        if n_frg == 0:
            continue
        n_bp = frg_set[0, 5]

        if n_frg > MAX_SIZE:
            n_frg = MAX_SIZE
        bin_idx = np.digitize(n_frg, edge_lst) - 1

        if n_bp < 1500:
            size_dist[0, bin_idx] += 1
        elif n_bp < 8000:
            size_dist[1, bin_idx] += 1
        else:
            size_dist[2, bin_idx] += 1
        size_dist[3, bin_idx] += 1

    # calculate measures
    n_map0 = np.sum(size_dist[3, :])
    n_map1 = np.sum(size_dist[3, 1:])
    n_map2 = np.sum(size_dist[3, 2:])

    # Plotting
    clr_map = [cm.Blues(x) for x in np.linspace(0.3, 1.0, 3)] + [(1.0, 0.5, 0.25)]
    plt.figure(figsize=(7, 5))
    plt_h = [None] * 4
    for cls_idx in range(4):
        plt_h[cls_idx] = plt.bar(edge_lst, size_dist[cls_idx, :] * 100.0 / np.sum(size_dist[cls_idx, :]),
                                 width=0.95 - cls_idx / 4.0, color=clr_map[cls_idx])[0]

    plt.xlim([0, MAX_SIZE + 1])
    plt.xticks(edge_lst)
    plt.xlabel('Read size (#fragment)')
    plt.ylabel('Frequency (%)')
    # plt.ylim([0, 70])
    title_msg = configs['run_id']
    if len(filter_lst) != 0:
        title_msg += ' ({:s})'.format(', '.join(filter_lst))
    title_msg += '\n#map>0={:,d};\n'.format(n_map0) + \
                 '#map>1={:,d} ({:0.0f}%); '.format(n_map1, n_map1 * 1e2 / n_map0) + \
                 '#map>2={:,d} ({:0.0f}%)'.format(n_map2, n_map2 * 1e2 / n_map0)
    plt.title(title_msg)
    plt.legend(plt_h, [
        'read #bp <1.5kb (n={:,d})'.format(np.sum(size_dist[0, :])),
        'read #bp <8kb (n={:,d})'.format(np.sum(size_dist[1, :])),
        'read #bp >8kb (n={:,d})'.format(np.sum(size_dist[2, :])),
        'All (n={:,d})'.format(np.sum(size_dist[3, :]))
    ])

    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/plt_CirSizeDistribution_' + configs['run_id']
        if roi_only or uniq_only:
             configs['output_file'] += '_{:s}.pdf'.format('-'.join(filter_lst))
        else:
            configs['output_file'] += '.pdf'
    plt.savefig(configs['output_file'], bbox_inches='tight')


def plot_overallProfile(configs, min_n_frg=2):
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt, patches

    from utilities import hasOL, load_mc4c, load_annotation

    # initialization
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/plt_OverallProfile_' + configs['run_id'] + '.pdf'
    edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=201, dtype=np.int64).reshape(-1, 1)
    bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
    bin_width = bin_bnd[0, 1] - bin_bnd[0, 0]
    bin_cen = np.mean(bin_bnd, axis=1)
    n_bin = bin_bnd.shape[0]
    del edge_lst
    vp_crd = np.array([configs['vp_cnum'], configs['vp_start'], configs['vp_end']])
    roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])

    # loop over datasets
    bin_frq = np.zeros([2, n_bin], dtype=np.int)
    n_read = np.zeros(2, dtype=np.int)
    for di in range(2):

        # load MC-HC data
        frg_dp = load_mc4c(configs, unique_only=di != 0, valid_only=True, min_mq=20, reindex_reads=True)
        frg_np = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd']].values

        # filter small circles
        is_vp = hasOL(vp_crd, frg_np[:, 1:4], offset=0)
        is_roi = hasOL(roi_crd, frg_np[:, 1:4], offset=0)
        frg_nvp = frg_np[~is_vp & is_roi, :]
        cir_size = np.bincount(frg_nvp[:, 0])[frg_nvp[:, 0]]
        is_inf = np.isin(frg_np[:, 0], frg_nvp[cir_size >= min_n_frg, 0])
        frg_inf = frg_np[is_inf, :]

        # select within roi fragments
        is_roi = hasOL(roi_crd, frg_inf[:, 1:4], offset=0)
        frg_roi = frg_inf[is_roi, :]
        n_read[di] = len(np.unique(frg_roi[:, 0]))

        # looping over bins
        for bi in range(n_bin):
            is_in = hasOL(bin_bnd[bi, :], frg_roi[:, 2:4])
            bin_frq[di, bi] = len(np.unique(frg_roi[is_in, 0]))  # each circle can contribute only once to a bin

    # set vp bins to nan
    # is_vp = hasOL([configs['vp_start'], configs['vp_end']], bin_bnd)
    # bin_frq[:, is_vp] = np.nan
    vpb_idx = hasOL([configs['vp_start'], configs['vp_end']], bin_bnd)
    vpd_bnd = [bin_bnd[vpb_idx][0, 0], bin_bnd[vpb_idx][-1, 1]]

    # plotting
    plt.figure(figsize=(15, 3))
    plt_h = [None] * 2
    clr_map = ['#d0d0d0', '#43ff14']
    bin_nrm = np.zeros([2, n_bin])
    for di in range(2):
        bin_nrm[di, :] = bin_frq[di, :] * 100.0 / n_read[di]
        bin_nrm[di, vpb_idx] = np.nan

        plt_h[di] = plt.bar(bin_cen, bin_nrm[di, :], width=bin_width, color=clr_map[di], alpha=0.7)

    # add vp area
    y_lim = [0, np.nanmax(bin_nrm) * 1.1]
    plt.gca().add_patch(patches.Rectangle([vpd_bnd[0], 0], vpd_bnd[1] - vpd_bnd[0], y_lim[1],
                                          linewidth=0, edgecolor='None', facecolor='orange'))

    # add annotations
    ant_pd = load_annotation(configs['genome_build'], roi_crd=roi_crd).reset_index(drop=True)
    for ai in range(ant_pd.shape[0]):
        ant_pos = ant_pd.loc[ai, 'ant_pos']
        plt.text(ant_pos, y_lim[1], ant_pd.loc[ai, 'ant_name'],
                 horizontalalignment='center', verticalalignment='bottom')
        plt.plot([ant_pos, ant_pos], y_lim, ':', color='#bfbfbf', linewidth=1, alpha=0.5)

    # final adjustments
    plt.xlim([configs['roi_start'], configs['roi_end']])
    x_ticks = np.linspace(configs['roi_start'], configs['roi_end'], 20, dtype=np.int64)
    x_tick_label = ['{:0.2f}m'.format(x / 1e6) for x in x_ticks]
    plt.xticks(x_ticks, x_tick_label, rotation=20)
    plt.ylabel('Frequency (% of reads)')
    plt.ylim(y_lim)
    plt.legend(plt_h, [
        'All reads (n={:0,.0f})'.format(n_read[0]),
        'Unique reads (n={:0,.0f})'.format(n_read[1])
    ])
    plt.title('Overall profile (#roiFrg>{:d}, ex. vp), {:s}\n'.format(min_n_frg - 1, configs['run_id']))
    plt.savefig(configs['output_file'], bbox_inches='tight')


def find_optimal_roi(configs, set_optimal_roi=False, min_cvg=2, min_mq=20):
    import pandas as pd
    from matplotlib import pyplot as plt
    from scipy.stats import spearmanr

    from utilities import load_configs, load_mc4c, limit_to_roi, get_chr_info, hasOL, get_nreads_per_bin
    from pre_process import remove_duplicates_by_umi

    # initialization
    print '%% Finding optimal ROI by analysing local coverage ...'
    if configs['output_file'] is None:
        configs['output_file'] = './plots/plt_OptimalROI_' + configs['run_id'] + '.pdf'
    vp_crd = np.array([configs['vp_cnum'], configs['vp_start'], configs['vp_end']])
    roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])

    # load dataset
    mc4c_pd = load_mc4c(configs, unique_only=False, valid_only=True, min_mq=min_mq, reindex_reads=True, verbose=True)
    header_lst = ['ReadID', 'Chr', 'ExtStart', 'ExtEnd']
    read_all = mc4c_pd[header_lst].values
    del mc4c_pd

    # select >1 #frg
    read_inf = limit_to_roi(read_all[:, :4], vp_crd=vp_crd, roi_crd=roi_crd, min_n_frg=2)

    # select reads with #traceable fragment > 1
    roi_size = configs['roi_end'] - configs['roi_start']
    lcl_crd = np.array([configs['vp_cnum'], configs['roi_start'] - roi_size, configs['roi_end'] + roi_size])
    is_lcl = hasOL(lcl_crd, read_inf[:, 1:4], offset=0)
    umi_def = read_inf[~is_lcl, :].copy()

    # remove duplicates
    unq_set, lcl_info = remove_duplicates_by_umi(umi_def)
    is_unq = np.isin(read_all[:, 0], unq_set[:, 0])
    pcr_def = read_all[is_unq, :].copy()
    n_def = len(np.unique(pcr_def[:, 0]))

    # use trans-umi only
    is_cis = lcl_crd[0] == read_inf[:, 1]
    umi_trs = read_inf[~is_cis, :].copy()
    unq_set, trs_info = remove_duplicates_by_umi(umi_trs)
    is_unq = np.isin(read_all[:, 0], unq_set[:, 0])
    pcr_trs = read_all[is_unq, :].copy()
    n_trs = len(np.unique(pcr_def[:, 0]))

    # split vp chromosome into blocks
    chr_size = get_chr_info(genome_str=configs['genome_build'], property='chr_size')
    blk_lst = np.arange(0, chr_size[configs['vp_cnum'] - 1], 30e3, dtype=np.int64).reshape(-1, 1)
    blk_lst[-1] = chr_size[configs['vp_cnum'] - 1]
    n_blk = len(blk_lst) - 1
    blk_crd = np.hstack([np.repeat(configs['vp_cnum'], n_blk).reshape(-1, 1), blk_lst[:-1], blk_lst[1:] - 1])
    blk_w = blk_crd[0, 1] - blk_crd[0, 0]
    blk_cen = np.mean(blk_crd[:, 1:3], axis=1)

    # compute coverage over chromosome
    cvg_def = get_nreads_per_bin(pcr_def[:, :4], bin_crd=blk_crd)[0]
    cvg_trs = get_nreads_per_bin(pcr_trs[:, :4], bin_crd=blk_crd)[0]
    nrm_def = pd.Series(cvg_def * 1e2 / n_def).rolling(5).mean().values
    nrm_trs = pd.Series(cvg_trs * 1e2 / n_trs).rolling(5).mean().values

    # select highly covered region
    high_cvg_idx = np.where(nrm_trs > min_cvg)[0]
    roi_crd_adj = np.array([configs['vp_cnum'], blk_crd[high_cvg_idx[0], 0], blk_crd[high_cvg_idx[-1], 1]])
    edge_lst_adj = np.linspace(roi_crd_adj[1], roi_crd_adj[2], num=201, dtype=np.int64).reshape(-1, 1)
    bin_w_adj = edge_lst_adj[1] - edge_lst_adj[0]
    roi_cen = np.mean([configs['vp_start'], configs['vp_end']])
    vp_crd_adj = np.array([configs['vp_cnum'], roi_cen - int(bin_w_adj * 1.5), roi_cen + int(bin_w_adj * 1.5)])

    # use adjusted area to compute duplicates
    adj_w = roi_crd_adj[2] - roi_crd_adj[1]
    is_adj = hasOL([roi_crd_adj[0], roi_crd_adj[1] - adj_w, roi_crd_adj[2] + adj_w], read_inf[:, 1:4], offset=0)
    umi_adj = read_inf[~is_adj, :].copy()
    adj_set, adj_info = remove_duplicates_by_umi(umi_adj)
    is_adj = np.isin(read_all[:, 0], adj_set[:, 0])
    pcr_adj = read_all[is_adj, :].copy()

    # compute coverage over chromosome using adjusted region
    cvg_adj = get_nreads_per_bin(pcr_adj[:, :4], bin_crd=blk_crd)[0]
    nrm_adj = pd.Series(cvg_adj * 1e2 / n_def).rolling(5).mean().values

    # compute roi profiles
    n_bin = 200
    edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=n_bin + 1, dtype=np.int64).reshape(-1, 1)
    bin_crd = np.hstack([np.repeat(configs['vp_cnum'], n_bin).reshape(-1, 1), edge_lst[:-1], edge_lst[1:] - 1])
    bin_w = edge_lst[1] - edge_lst[0]
    prf_def, n_def = get_nreads_per_bin(pcr_def[:, :4], bin_crd=bin_crd, min_n_frg=2)
    prf_trs, n_trs = get_nreads_per_bin(pcr_trs[:, :4], bin_crd=bin_crd, min_n_frg=2)
    prf_adj, n_adj = get_nreads_per_bin(pcr_adj[:, :4], bin_crd=bin_crd, min_n_frg=2)

    # plotting
    fig = plt.figure(figsize=(15, 5))
    ax_crr = plt.subplot2grid((1, 3), (0, 0), rowspan=1, colspan=1)
    ax_prf = plt.subplot2grid((1, 3), (0, 1), rowspan=1, colspan=2)

    # plot correlations
    ax_crr.plot(prf_def * 1e2 / n_def, prf_trs * 1e2 / n_trs, 'o', color='gray')
    ax_crr.plot(prf_adj * 1e2 / n_adj, prf_trs * 1e2 / n_trs, 'x', color='blue')
    ax_crr.xlim([0, 10])
    ax_crr.ylim([0, 10])
    ax_crr.title('ROI coverage correlations, trans UMI vs far-cis + trans UMI' +
              'spr-corr: {:0.5f}'.format(spearmanr([prf_def, prf_adj]).correlation))

    # plot roi profiles
    y_lim = [0, 10]
    plt.plot([roi_crd_adj[1], roi_crd_adj[1]], [0, 1], color='blue')
    plt.plot([roi_crd_adj[2], roi_crd_adj[2]], [0, 1], color='blue')
    plt.text(roi_crd_adj[1], y_lim[1] * 0.8, '{:,d}> '.format(roi_crd_adj[1]), horizontalalignment='right', color='red')
    plt.text(roi_crd_adj[2], y_lim[1] * 0.8, ' <{:,d}'.format(roi_crd_adj[2]), horizontalalignment='left', color='red')
    plt.xlim([roi_crd_adj[1] - adj_w * 3, roi_crd_adj[2] + adj_w * 3])
    plt.ylim(y_lim)

    plt_h = [None] * 3
    plt_h[0] = plt.plot(blk_cen, nrm_def, ':', color='gray')[0]
    plt_h[1] = plt.plot(blk_cen, nrm_trs, ':', color='gray')[0]
    plt_h[2] = plt.plot(blk_cen, nrm_adj, ':', color='gray')[0]
    plt.legend(plt_h, [
        'Far-cis (n={:0.0f})'.format(n_def),
        'Trans only (n={:0.0f})'.format(n_trs),
        'Adjusted (n={:0.0f})'.format(n_adj)
    ])
    plt.suptitle('{:s}\n'.format(configs['run_id']) +
                 '#block={:d}, block_w={:0.0f}k, bin_w={:0.0f}\n'.format(n_blk, blk_w / 1e3, bin_w))
    plt.savefig(configs['output_file'], bbox_inches='tight')

