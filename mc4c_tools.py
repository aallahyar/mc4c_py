import numpy as np
import pandas as pd
import h5py


def load_mc4c(config_lst, target_field='frg_np', data_path='./datasets/',
              min_mq=20, only_unique=True, reindex_reads=True, max_rows=np.inf):
    MAX_N_CIR = 1000000000000
    out_pd = pd.DataFrame()
    if not isinstance(config_lst, list):
        config_lst = [config_lst]

    header_lst = []
    for vi, configs in enumerate(config_lst):
        if configs['input_file'] is None:
            configs['input_file'] = data_path + '/mc4c_{:s}_uniq.hdf5'.format(configs['run_id'])
        print('Loading mc4c dataset: {:s}'.format(configs['input_file']))

        h5_fid = h5py.File(configs['input_file'], 'r')
        if np.isinf(max_rows):
            data_np = h5_fid[target_field].value
        else:
            print 'Selecting only top [{:d}] rows in dataset'.format(max_rows)
            data_np = h5_fid[target_field][:max_rows]

        header_lst = list(h5_fid[target_field + '_header_lst'].value)
        h5_fid.close()
        part_pd = pd.DataFrame(data_np, columns=header_lst)

        # Filtering fragments
        if min_mq:
            part_pd = part_pd.loc[part_pd['MQ'] >= min_mq].copy()
        if only_unique:
            part_pd = part_pd.loc[part_pd['IsUnique'] > 0].copy()

        # Adjust Read IDs
        assert np.max(part_pd['ReadID']) < MAX_N_CIR
        part_pd['ReadID'] = part_pd['ReadID'] + (vi + 1) * MAX_N_CIR

        # Append the part
        out_pd = out_pd.append(part_pd, ignore_index=True)
    out_pd = out_pd[header_lst]

    if reindex_reads:
        header_lst.append('ReadID_original')
        out_pd[header_lst[-1]] = out_pd['ReadID'].copy()
        out_pd['ReadID'] = np.unique(out_pd['ReadID'], return_inverse=True)[1] + 1
        print 'Got [{:,d}] reads and [{:,d}] fragments after re-indexing.'.format(
            np.max(out_pd['ReadID']), out_pd.shape[0])

    return out_pd[header_lst]


def load_configs(cfg_fname):
    """ Read configurations from given file, put it into a dict

    :param cfg_fname: takes a path to a tab-separated file with one variable name and value per line, multiple values
    are seprated by ";").

    :returns: Dictionary where keys are based on the first column with values in a list.
    """
    from os import path
    from utilities import get_chr_info

    # check if config_file is a file
    if cfg_fname[-4:] != '.cfg':
        # print 'Configuration file does not end with ".cfg". Assuming the given name as a run ID.'
        cfg_fname = './configs/cfg_' + cfg_fname + '.cfg'

    # add global config file to list
    if path.isfile('./mc4c.cfg'):
        cfg_flist = ['./mc4c.cfg', cfg_fname]
    else:
        cfg_flist = [cfg_fname]

    # Load global and then given configs
    configs = dict()
    for fname in cfg_flist:
        with open(fname, 'r') as cfg_fid:
            for line in cfg_fid:
                if (line[0] == '#') or (len(line) == 1):
                    continue
                columns = line.rstrip('\n').split('\t')
                assert len(columns) == 2
                fld_lst = columns[1].split(';')
                if len(fld_lst) == 1:
                    configs[columns[0]] = fld_lst[0]
                else:
                    configs[columns[0]] = fld_lst

    # Conversions
    for cnf_name in ['vp_start', 'vp_end', 'roi_start', 'roi_end']:
        configs[cnf_name] = int(configs[cnf_name])
    for cnf_name in ['prm_start', 'prm_end']:
        configs[cnf_name] = [int(value) for value in configs[cnf_name]]
    for cnf_name in ['bwa_index_path', 'ref_genome_file']:
        configs[cnf_name] = configs[cnf_name].replace('%REF%', configs['genome_build'])
    chr_lst = get_chr_info(genome_str=configs['genome_build'], property='chr_name')
    chr_map = dict(zip(chr_lst, range(1, len(chr_lst) + 1)))
    configs['vp_cnum'] = chr_map[configs['vp_chr']]

    # Check configs that should be of equal length
    linked_configs = [
        ['prm_seq','prm_start','prm_end'],
        ['re_name','re_seq'],
    ]
    for cnf_set in linked_configs:
        assert len(set([len(configs[x]) for x in cnf_set])) == 1, \
            'Error: different lengths for linked configs:'+','.join(str(x) for x in cnf_set)

    return configs


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
        configs['output_file'] = configs['output_dir'] + '/plt_' + configs['run_id'] + '_ReadSizeDistribution.pdf'
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
    with pysam.FastxFile(configs['input_file'], persist=False) as gz_fid:
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


def plot_cvgDistribution(configs):
    import numpy as np
    import gzip
    from os import path
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt

    from utilities import get_chr_info, get_re_info

    # initialization
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/plt_' + configs['run_id'] + '_cvgDistribution.pdf'
    MAX_SIZE = 1500
    edge_lst = np.linspace(0, MAX_SIZE, 31)
    n_bin = len(edge_lst) - 1
    re_set_name = '-'.join(configs['re_name'])

    # get chr info
    chr_lst = get_chr_info(genome_str=configs['genome_build'], property='chr_name')
    n_chr = len(chr_lst)

    # Read res-enz positions
    re_fname = './renzs/{:s}_{:s}.npz'.format(configs['genome_build'], '-'.join(configs['re_name']))
    print 'Loading RE positions from: {:s}'.format(re_fname)
    if not path.isfile(re_fname):
        from utilities import extract_re_positions
        extract_re_positions(genome_str=configs['genome_build'], re_name_lst=configs['re_name'],
                             output_fname=re_fname, ref_fasta=configs['ref_genome_file'])
    re_pos_lst = get_re_info(re_name=re_set_name, property='pos', genome_str=configs['genome_build'])

    # compute ref fragment size
    ref_re_sd = np.zeros(n_bin, dtype=np.int64)
    for chr_idx, chr_name in enumerate(chr_lst):
        re_size = np.diff(re_pos_lst[chr_idx])
        re_size[re_size >= MAX_SIZE] = MAX_SIZE - 1

        bin_idx = np.digitize(re_size, edge_lst) - 1
        ref_re_sd += np.bincount(bin_idx, minlength=n_bin)

    # Load MC-HC data
    frg_dp = load_mc4c(configs, min_mq=0, reindex_reads=False)
    frg_np = frg_dp[['Chr', 'ExtStart', 'ExtEnd', 'MQ']].values
    del frg_dp

    # calculate chromosome coverage
    frg_size = frg_np[:, 2] - frg_np[:, 1] + 1
    is_mq20 = frg_np[:, 3] >= 20
    ccvg_mq01 = np.bincount(frg_np[:, 0] - 1, minlength=n_chr)
    ccvg_mq20 = np.bincount(frg_np[is_mq20, 0] - 1, minlength=n_chr)
    del frg_np
    n_bp_mq01 = np.sum(frg_size)
    n_bp_mq20 = np.sum(frg_size[is_mq20])

    # calculate fragment size distribution
    frg_size[frg_size >= MAX_SIZE] = MAX_SIZE - 1
    bin_idx = np.digitize(frg_size, edge_lst) - 1
    dist_mq01 = np.bincount(bin_idx, minlength=n_bin)
    n_mq01 = len(frg_size)

    bin_idx = np.digitize(frg_size[is_mq20], edge_lst) - 1
    dist_mq20 = np.bincount(bin_idx, minlength=n_bin)
    n_mq20 = np.sum(is_mq20)
    del frg_size

    # calculate raw fragment size
    frg_fname = './fragments/frg_{:s}.fasta.gz'.format(configs['run_id'])
    print 'Reading {:s}'.format(frg_fname)
    dist_mq00 = np.zeros(n_bin, dtype=np.int64)
    n_mq00 = 0
    with gzip.open(frg_fname, 'r') as splt_fid:
        while True:
            frg_sid = splt_fid.readline()
            frg_seq = splt_fid.readline().rstrip('\n')
            if frg_sid == '':
                break
            if n_mq00 % 100000 == 0:
                print('{:,d} fragments are read.'.format(n_mq00))

            seq_size = len(frg_seq)
            if seq_size >= MAX_SIZE:
                seq_size = MAX_SIZE - 1

            bin_idx = np.digitize(seq_size, edge_lst) - 1
            dist_mq00[bin_idx] += 1
            n_mq00 += 1

    # Plotting
    plt.figure(figsize=(15, 5))
    plt.subplot(1, 2, 1)
    q01_h = plt.bar(range(n_chr), ccvg_mq01, color='#bc85ff', width=0.80, alpha=0.7)
    q20_h = plt.bar(range(n_chr), ccvg_mq20, color='#8929ff', width=0.60, alpha=0.7)
    plt.legend([q01_h, q20_h], [
        'MQ>=1, #mapped bp: {:,d}'.format(n_bp_mq01),
        'MQ>=20, #mapped bp: {:,d}'.format(n_bp_mq20)
    ])
    plt.xticks(range(n_chr), chr_lst, rotation=35, ha='right')
    plt.xlim([-0.5, n_chr - 0.5])
    y_ticks = plt.yticks()[0]
    y_tick_lbl = ['{:0,.1f}m'.format(y / 1e6) for y in y_ticks]
    plt.yticks(y_ticks, y_tick_lbl)
    plt.ylabel('#Fragments')
    plt.title('Chromosome coverage, {:s}\n'.format(configs['run_id']) +
              '#unmapped MQ1={:0,.1f}m; MQ20={:0,.1f}m'.format((n_mq00 - n_mq01) / 1e6, (n_mq00 - n_mq20) / 1e6))

    plt.subplot(1, 2, 2)
    ref_h = plt.bar(range(n_bin), ref_re_sd * float(n_mq00) / np.sum(ref_re_sd), width=1.00, color='#dddddd')
    q00_h = plt.bar(range(n_bin), dist_mq00, width=0.90, color='#aaaaaa')
    q01_h = plt.bar(range(n_bin), dist_mq01, width=0.70)
    q20_h = plt.bar(range(n_bin), dist_mq20, width=0.50)
    plt.legend([ref_h, q00_h, q01_h, q20_h], [
        'Ref (#frg={:0,.1f}m), normed'.format(np.sum(ref_re_sd) / 1e6),
        'Raw (#frg={:0,.1f}m)'.format(n_mq00 / 1e6),
        'MQ1 (#frg={:0,.1f}m)'.format(n_mq01 / 1e6),
        'MQ20 (#frg={:0,.1f}m)'.format(n_mq20 / 1e6)])
    plt.xlim([-1, n_bin + 1])
    x_ticks_idx = range(1, n_bin, 2)
    plt.xticks(x_ticks_idx, ['{:0.0f}'.format(edge_lst[i + 1]) for i in x_ticks_idx], rotation=35)
    plt.xlabel('#base pairs')

    y_ticks = plt.yticks()[0]
    y_tick_lbl = ['{:0,.1f}m'.format(y / 1e6) for y in y_ticks]
    plt.yticks(y_ticks, y_tick_lbl)
    plt.ylabel('#Fragments')

    plt.title('Fragment size distribution, {:s}'.format(configs['run_id']) +
              '\nMap success rate MQ1={:0.1f}%; MQ20={:0.1f}%]'.format(n_mq01 * 100 / n_mq00, n_mq20 * 100 / n_mq00))
    plt.savefig(configs['output_file'], bbox_inches='tight')


def plot_cirSizeDistribution(configs):
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt, cm

    from utilities import accum_array, hasOL

    # initialization
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/plt_' + configs['run_id'] + '_CirSizeDistribution.pdf'
    MAX_SIZE = 8
    edge_lst = np.linspace(1, MAX_SIZE, num=MAX_SIZE)
    n_edge = len(edge_lst)

    # Load MC-HC data
    frg_dp = load_mc4c(configs, min_mq=20, reindex_reads=True)
    frg_np = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'MQ', 'ReadLength']].values
    del frg_dp

    # group circles
    read_grp = accum_array(frg_np[:, 0] - 1, frg_np)
    n_grp = len(read_grp)

    # Looping over circles
    size_dist = np.zeros([4, n_edge], dtype=np.int64)
    print 'Computing circle size from {:d} reads:'.format(n_grp)
    for read_idx, frg_set in enumerate(read_grp):
        if read_idx % 10000 == 0:
            print('\t{:,d}/{:,d} Reads are processed.'.format(read_idx, n_grp))
        n_frg = frg_set.shape[0]
        if n_frg == 0:
            continue

        if n_frg > MAX_SIZE:
            n_frg = MAX_SIZE
        bin_idx = np.digitize(n_frg, edge_lst) - 1

        if frg_set[0, 5] < 3000:
            read_cls = 0
        elif frg_set[0, 5] < 7000:
            read_cls = 1
        else:
            read_cls = 2

        size_dist[read_cls, bin_idx] += 1
    size_dist[3, :] = np.sum(size_dist, axis=0)

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
    plt.ylim([0, 70])
    plt.title(configs['run_id'] + '\n' +
              '#mapped={:,d}; '.format(np.sum(size_dist[3, :])) +
              '#map>1={:,d}; '.format(np.sum(size_dist[3, 1:])) +
              '#map>2={:,d}'.format(np.sum(size_dist[3, 2:]))
              )
    plt.legend(plt_h, [
        'read size <3kb (n={:,d})'.format(np.sum(size_dist[0, :])),
        'read size <7kb (n={:,d})'.format(np.sum(size_dist[1, :])),
        'read size >7kb (n={:,d})'.format(np.sum(size_dist[2, :])),
        'All (n={:,d})'.format(np.sum(size_dist[3, :]))
    ])
    plt.savefig(configs['output_file'], bbox_inches='tight')



