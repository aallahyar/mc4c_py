import numpy as np


def plot_readSizeDistribution(configs):
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


def find_optimal_roi(config_lst, min_cvg=2, min_mq=20):
    import pandas as pd
    from matplotlib import pyplot as plt
    from scipy.stats import spearmanr

    from utilities import load_mc4c, limit_to_roi, get_chr_info, hasOL, get_nreads_per_bin, load_annotation
    from pre_process import remove_duplicates_by_umi

    # initialization
    print '%% Finding optimal ROI by analysing local coverage ...'
    run_id = ','.join([config['run_id'] for config in config_lst])
    configs = config_lst[0]
    if configs['output_file'] is None:
        configs['output_file'] = './plots/plt_FindROIbyCoverage_' + run_id + '.pdf'
    def_vp_crd = np.array([configs['vp_cnum'], configs['vp_start'], configs['vp_end']])
    def_roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
    def_roi_cen = int(np.mean([configs['vp_start'], configs['vp_end']]))
    def_roi_w = configs['roi_end'] - configs['roi_start']
    blk_w = 30000
    n_bin = 200

    # split vp chromosome into blocks
    chr_size = get_chr_info(genome_str=configs['genome_build'], property='chr_size')
    blk_lst = np.arange(0, chr_size[configs['vp_cnum'] - 1], blk_w, dtype=np.int64).reshape(-1, 1)
    blk_lst[-1] = chr_size[configs['vp_cnum'] - 1]
    n_blk = len(blk_lst) - 1
    blk_crd = np.hstack([np.repeat(configs['vp_cnum'], n_blk).reshape(-1, 1), blk_lst[:-1], blk_lst[1:] - 1])
    blk_cen = np.mean(blk_crd[:, 1:3], axis=1)

    # use raw reads and extract UMIs
    print '[i] Extracting UMIs from each dataset ...'
    read_all = np.empty([0, 4], dtype=np.int)
    def_pcr = np.empty([0, 4], dtype=np.int)
    trs_pcr = np.empty([0, 4], dtype=np.int)
    max_n_read = 10000000
    for idx, cfg in enumerate(config_lst):

        # load raw data
        print('Loading all reads from {:s} ...'.format(cfg['run_id']))
        mc4c_part = load_mc4c(cfg, unique_only=False, valid_only=True, min_mq=min_mq, reindex_reads=True, verbose=False)
        read_prt = mc4c_part[['ReadID', 'Chr', 'ExtStart', 'ExtEnd']].values

        # use default approach, use >1 roi-frg
        has_inf = limit_to_roi(read_prt[:, :4], vp_crd=def_vp_crd, roi_crd=def_roi_crd, min_n_frg=2)
        def_read_inf = read_prt[np.isin(read_prt[:, 0], has_inf[:, 0]), :].copy()
        def_umi_set = def_read_inf[def_read_inf[:, 1] != configs['vp_cnum'], :].copy()
        has_uid = remove_duplicates_by_umi(def_umi_set)[0]
        def_pcr_prt = read_prt[np.isin(read_prt[:, 0], has_uid[:, 0]), :].copy()
        del has_inf, has_uid
        print '\t{:,d} unique reads are added using [far-cis + trans] UMIs.'.format(len(np.unique(def_pcr_prt[:, 0])))

        # select >1 cis-frg
        trs_vp_crd = np.array([configs['vp_cnum'], def_roi_cen - 5000, def_roi_cen + 5000])
        is_cis = read_prt[:, 1] == configs['vp_cnum']
        is_vp = hasOL(trs_vp_crd, read_prt[:, 1:4])
        trs_read_m1c = read_prt[~is_vp & is_cis, :].copy()
        read_size = np.bincount(trs_read_m1c[:, 0], minlength=np.max(trs_read_m1c[:, 0]) + 1)[trs_read_m1c[:, 0]]
        trs_read_m1c = read_prt[np.isin(read_prt[:, 0], trs_read_m1c[read_size >= 2, 0]), :].copy()
        del read_size

        # select and process UMIs
        trs_umi_set = trs_read_m1c[trs_read_m1c[:, 1] != configs['vp_cnum'], :].copy()
        trs_uid = remove_duplicates_by_umi(trs_umi_set)[0]
        trs_pcr_prt = read_prt[np.isin(read_prt[:, 0], trs_uid[:, 0]), :].copy()
        print '\t{:,d} unique reads are added using [trans only] UMIs.'.format(len(np.unique(trs_pcr_prt[:, 0])))

        # add data specific identifiers
        assert np.max(read_prt[:, 0]) < max_n_read
        def_pcr_prt[:, 0] = def_pcr_prt[:, 0] + (idx + 1) * max_n_read
        trs_pcr_prt[:, 0] = trs_pcr_prt[:, 0] + (idx + 1) * max_n_read
        read_prt[:, 0] = read_prt[:, 0] + (idx + 1) * max_n_read

        # appending
        read_all = np.vstack([read_all, read_prt])
        def_pcr = np.vstack([def_pcr, def_pcr_prt])
        trs_pcr = np.vstack([trs_pcr, trs_pcr_prt])

    # compute coverage over chromosome
    def_cvg, n_def = get_nreads_per_bin(def_pcr[:, :4], bin_crd=blk_crd, min_n_frg=2)
    trs_cvg, n_trs = get_nreads_per_bin(trs_pcr[:, :4], bin_crd=blk_crd, min_n_frg=2)
    def_nrm = pd.Series(def_cvg * 1e2 / n_def).rolling(5, center=True).mean().values
    trs_nrm = pd.Series(trs_cvg * 1e2 / n_trs).rolling(5, center=True).mean().values

    # select highly covered region
    np.seterr(all='ignore')
    cvd_idx = np.where(trs_nrm > min_cvg)[0]
    np.seterr(all=None)
    adj_roi_crd = np.array([configs['vp_cnum'], blk_crd[cvd_idx[0], 1], blk_crd[cvd_idx[-1], 2]])
    adj_roi_w = adj_roi_crd[2] - adj_roi_crd[1]
    adj_edge_lst = np.linspace(adj_roi_crd[1], adj_roi_crd[2], num=n_bin + 1, dtype=np.int64).reshape(-1, 1)
    adj_bin_w = adj_edge_lst[1, 0] - adj_edge_lst[0, 0]
    adj_vp_crd = np.array([configs['vp_cnum'], def_roi_cen - int(adj_bin_w * 1.5), def_roi_cen + int(adj_bin_w * 1.5)])

    # select informative reads
    has_inf = limit_to_roi(read_all[:, :4], vp_crd=adj_vp_crd, roi_crd=adj_roi_crd, min_n_frg=2)
    adj_read_inf = read_all[np.isin(read_all[:, 0], has_inf[:, 0]), :].copy()
    adj_lcl_crd = [adj_roi_crd[0], adj_roi_crd[1] - adj_roi_w, adj_roi_crd[2] + adj_roi_w]
    is_lcl = hasOL(adj_lcl_crd, adj_read_inf[:, 1:4], offset=0)
    adj_umi = adj_read_inf[~is_lcl, :].copy()
    adj_uid = remove_duplicates_by_umi(adj_umi)[0]
    is_adj = np.isin(adj_read_inf[:, 0], adj_uid[:, 0])
    adj_pcr = adj_read_inf[is_adj, :].copy()

    # compute coverage over chromosome using adjusted region
    adj_cvg, n_adj = get_nreads_per_bin(adj_pcr[:, :4], bin_crd=blk_crd, min_n_frg=2)
    adj_nrm = pd.Series(adj_cvg * 1e2 / n_adj).rolling(5, center=True).mean().values

    # compute roi profiles
    def_edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=n_bin + 1, dtype=np.int64).reshape(-1, 1)
    def_bin_crd = np.hstack([np.repeat(configs['vp_cnum'], n_bin).reshape(-1, 1), def_edge_lst[:-1], def_edge_lst[1:] - 1])
    def_bin_w = def_bin_crd[0, 2] - def_bin_crd[0, 1]
    def_prf, n_def = get_nreads_per_bin(def_pcr[:, :4], bin_crd=def_bin_crd, min_n_frg=2)
    trs_prf, n_trs = get_nreads_per_bin(trs_pcr[:, :4], bin_crd=def_bin_crd, min_n_frg=2)
    adj_prf, n_adj = get_nreads_per_bin(adj_pcr[:, :4], bin_crd=def_bin_crd, min_n_frg=2)

    # plotting
    plt.figure(figsize=(25, 5))
    ax_crr = plt.subplot2grid((1, 4), (0, 0), rowspan=1, colspan=1)
    ax_prf = plt.subplot2grid((1, 4), (0, 1), rowspan=1, colspan=3)
    plt_h = [None] * 4
    x_lim = [configs['roi_start'] - def_roi_w * 2, configs['roi_end'] + def_roi_w * 2]
    y_lim = [0, 10]

    # plot correlations
    plt_h[0] = ax_crr.plot(def_prf * 1e2 / n_def, trs_prf * 1e2 / n_trs, 'x', color='#ffad14')[0]
    plt_h[1] = ax_crr.plot(def_prf * 1e2 / n_def, adj_prf * 1e2 / n_adj, 'o', color='#2e2eff')[0]
    ax_crr.set_xlim(y_lim)
    ax_crr.set_ylim(y_lim)
    ax_crr.set_xlabel('Default ROI')
    ax_crr.set_ylabel('Trans / Adjusted ROI')
    ax_crr.set_title('ROI coverage Spearman correlations\n' +
                     'def-UMI vs. trs-UMI: {:0.5f}\n'.format(spearmanr(def_prf, trs_prf).correlation) +
                     'def-UMI vs. adj-UMI: {:0.5f}'.format(spearmanr(def_prf, adj_prf).correlation))
    ax_crr.legend(plt_h[:2], ['Default vs Trans profile', 'Default vs. Adjusted profile'])

    # plot roi profiles
    ax_prf.plot([def_roi_crd[1], def_roi_crd[1]], y_lim, color='#bcbcbc')
    ax_prf.plot([def_roi_crd[2], def_roi_crd[2]], y_lim, color='#bcbcbc')
    ax_prf.plot([adj_roi_crd[1], adj_roi_crd[1]], y_lim, color='#a3d1ff')
    ax_prf.plot([adj_roi_crd[2], adj_roi_crd[2]], y_lim, color='#a3d1ff')
    ax_prf.text(def_roi_crd[1], y_lim[1] * 0.9, '{:,d}> '.format(def_roi_crd[1]), horizontalalignment='right', color='#9c9c9c')
    ax_prf.text(def_roi_crd[2], y_lim[1] * 0.9, ' <{:,d}'.format(def_roi_crd[2]), horizontalalignment='left', color='#9c9c9c')
    ax_prf.text(adj_roi_crd[1], y_lim[1] * 0.8, '{:,d}> '.format(adj_roi_crd[1]), horizontalalignment='right', color='#52a8ff')
    ax_prf.text(adj_roi_crd[2], y_lim[1] * 0.8, ' <{:,d}'.format(adj_roi_crd[2]), horizontalalignment='left', color='#52a8ff')

    plt_h[0] = ax_prf.plot(blk_cen, def_nrm, '--', color='#777777')[0]
    plt_h[1] = ax_prf.plot(blk_cen, trs_nrm, ':o', color='#ffad14', alpha=0.8, markersize=4, markeredgecolor=None)[0]
    plt_h[2] = ax_prf.plot(blk_cen, adj_nrm, '-',  color='#2e2eff')[0]

    # add annotations
    ant_pd = load_annotation(configs['genome_build'], roi_crd=[configs['vp_cnum']] + x_lim)
    for ai in range(ant_pd.shape[0]):
        ant_pos = ant_pd.loc[ai, 'ant_pos']
        if (ai == 0) or (np.abs(ant_pd.loc[ai - 1, 'ant_pos'] - ant_pos) > def_roi_w / 50.0):
            ax_prf.text(ant_pos, y_lim[0], ' ' + ant_pd.loc[ai, 'ant_name'],
                        horizontalalignment='center', verticalalignment='bottom', rotation=90, fontsize=6)
        ax_prf.plot([ant_pos, ant_pos], [y_lim[0] + y_lim[1] * 0.05, y_lim[1]],
                    ':', color='#bfbfbf', linewidth=1, alpha=0.3)
    plt_h[3] = ax_prf.plot([x_lim[0], x_lim[1]], [min_cvg, min_cvg], '--', color='#ff8f8f')[0]

    ax_prf.set_xlim(x_lim)
    ax_prf.set_ylim(y_lim)
    x_ticks = np.linspace(x_lim[0], x_lim[1], 25, dtype=np.int64)
    x_tick_label = ['{:0.3f}m'.format(x / 1e6) for x in x_ticks]
    plt.xticks(x_ticks, x_tick_label, rotation=20)
    plt.ylabel('Frequency (% of reads)')
    ax_prf.legend(plt_h, [
        'Far-cis + trans (n={:0.0f})'.format(n_def),
        'Trans only (n={:0.0f})'.format(n_trs),
        'Adjusted (n={:0.0f})'.format(n_adj),
        'Coverage threshold ({:0.1f}%)'.format(min_cvg)
    ], loc='center left')

    plt.title('{:s}\n'.format(run_id) +
              '#block={:d}, block_w={:0.0f}k\n'.format(n_blk, blk_w / 1e3) +
              'bin-w (def, adjusted): {:0,.0f}bp; {:0,.0f}bp'.format(def_bin_w, adj_bin_w) + '\n' +
              'roi-w (def, adjusted): {:0,.0f}bp; {:0,.0f}bp'.format(def_roi_w, adj_roi_w))
    plt.savefig(configs['output_file'], bbox_inches='tight')

