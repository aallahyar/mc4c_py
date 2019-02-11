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
        configs['input_file'] = './fastqs/raw_' + configs['run_id'] + '.fastq.gz'
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
    dist_ref = np.zeros(n_bin, dtype=np.int64)
    n_re_ref = 0
    for chr_idx, chr_name in enumerate(chr_lst):
        re_size = np.diff(re_pos_lst[chr_idx]) + 1
        n_re_ref += len(re_size)
        re_size[re_size >= MAX_SIZE] = MAX_SIZE - 1

        bin_idx = np.digitize(re_size, edge_lst) - 1
        dist_ref += np.bincount(bin_idx, minlength=n_bin)

    # Load MC-HC data
    frg_dp = load_mc4c(configs, min_mq=0, reindex_reads=False, uniq_only=False, valid_only=True)
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
    del frg_size

    # calculate raw fragment size
    frg_fname = './fragments/frg_{:s}.fasta.gz'.format(configs['run_id'])
    print 'Scanning raw fragments in {:s}'.format(frg_fname)
    dist_mq00 = np.zeros(n_bin, dtype=np.int64)
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

            seq_size = len(frg_seq)
            n_bp_mq00 += seq_size
            if seq_size >= MAX_SIZE:
                seq_size = MAX_SIZE - 1

            bin_idx = np.digitize(seq_size, edge_lst) - 1
            dist_mq00[bin_idx] += 1
            n_frg_mq00 += 1

    # Plotting
    plt.figure(figsize=(7, 5))
    ref_h = plt.bar(range(n_bin), dist_ref * float(n_frg_mq00) / np.sum(dist_ref), width=1.00, color='#dddddd')
    q00_h = plt.bar(range(n_bin), dist_mq00, width=0.90, color='#aaaaaa')
    q01_h = plt.bar(range(n_bin), dist_mq01, width=0.70)
    q20_h = plt.bar(range(n_bin), dist_mq20, width=0.50)
    plt.legend([ref_h, q00_h, q01_h, q20_h], [
        'Ref (#frg={:0,.0f}k), normed'.format(n_re_ref / 1e3),
        'Raw (#frg={:0,.0f}k)'.format(n_frg_mq00 / 1e3),
        'MQ1', 'MQ20'])
    plt.xlim([-1, n_bin + 1])
    x_ticks_idx = range(1, n_bin, 2)
    plt.xticks(x_ticks_idx, ['{:0.0f}'.format(edge_lst[i + 1]) for i in x_ticks_idx], rotation=35)
    plt.xlabel('#base pairs')

    y_ticks = plt.yticks()[0]
    y_tick_lbl = ['{:0,.1f}k'.format(y / 1e3) for y in y_ticks]
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
    frg_dp = load_mc4c(configs, min_mq=0, reindex_reads=False, uniq_only=False, valid_only=True)
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
    frg_dp = load_mc4c(configs, min_mq=20, reindex_reads=True, uniq_only=uniq_only)
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
    read_grp = accum_array(frg_np[:, 0] - 1, frg_np)
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

    # Plotting
    clr_map = [cm.Blues(x) for x in np.linspace(0.3, 1.0, 3)] + [(1.0, 0.5, 0.25)]
    plt.figure(figsize=(10, 5))
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
    title_msg += '\n#mapped={:,d}; '.format(np.sum(size_dist[3, :])) + \
                 '#map>1={:,d}; '.format(np.sum(size_dist[3, 1:])) + \
                 '#map>2={:,d}'.format(np.sum(size_dist[3, 2:]))
    plt.title(title_msg)
    plt.legend(plt_h, [
        'read size <1.5kb (n={:,d})'.format(np.sum(size_dist[0, :])),
        'read size <8kb (n={:,d})'.format(np.sum(size_dist[1, :])),
        'read size >8kb (n={:,d})'.format(np.sum(size_dist[2, :])),
        'All (n={:,d})'.format(np.sum(size_dist[3, :]))
    ])

    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/plt_CirSizeDistribution_' + configs['run_id']
        if roi_only or uniq_only:
             configs['output_file'] += '_{:s}.pdf'.format('-'.join(filter_lst))
        else:
            configs['output_file'] += '.pdf'
    plt.savefig(configs['output_file'], bbox_inches='tight')


def plot_overallProfile(configs, uniq_only=True, MIN_N_FRG=2):
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
    n_bin = len(bin_bnd)
    del edge_lst

    # load MC-HC data
    frg_dp = load_mc4c(configs, uniq_only=uniq_only, valid_only=True, min_mq=20, reindex_reads=True)
    frg_np = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd']].values
    del frg_dp

    # select within roi fragments
    is_roi = hasOL([configs['vp_cnum'], configs['roi_start'], configs['roi_end']], frg_np[:, 1:4])
    frg_roi = frg_np[is_roi, :]
    del frg_np

    # filter small circles
    cir_size = np.bincount(frg_roi[:, 0])[frg_roi[:, 0]]
    frg_roi = frg_roi[cir_size >= MIN_N_FRG, :]
    n_read = len(np.unique(frg_roi[:, 0]))

    # looping over bins
    bin_freq = np.zeros(n_bin)
    for bi in range(n_bin):
        is_in = hasOL(bin_bnd[bi, :], frg_roi[:, 2:4])
        bin_freq[bi] = len(np.unique(frg_roi[is_in, 0]))  # each circle can contribute only once to a bin

    # set vp bins to nan
    is_vp = hasOL([configs['vp_start'], configs['vp_end']], bin_bnd)
    bin_freq[is_vp] = np.nan
    vp_bnd = [bin_bnd[is_vp, 0][0], bin_bnd[is_vp, 1][-1]]

    # plotting
    plt.figure(figsize=(15, 5))
    bin_cen = np.mean(bin_bnd, axis=1)
    bin_nrm = bin_freq * 100.0 / np.nansum(bin_freq)
    plt.bar(bin_cen, bin_nrm, width=bin_width, color='#43ff14')

    # add vp area
    y_lim = [0, np.nanmax(bin_nrm) * 1.1]
    plt.gca().add_patch(patches.Rectangle([vp_bnd[0], 0], vp_bnd[1] - vp_bnd[0], y_lim[1],
                      linewidth=0, edgecolor='None', facecolor='orange'))

    # add annotations
    ant_pd = load_annotation(configs['genome_build'], roi_crd=[configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
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
    plt.ylabel('Frequency (%)')
    plt.ylim(y_lim)
    plt.title('Overall profile, {:s}\n'.format(configs['run_id']) +
              '#read (#roiFrg>{:d}, ex. vp)={:,d}\n'.format(MIN_N_FRG - 1, n_read)
              )
    plt.savefig(configs['output_file'], bbox_inches='tight')

