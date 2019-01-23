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
            configs['input_file'] = data_path + 'mc4c_{:s}_uniq.hdf5'.format(configs['run_id'])
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

    # Convert to Integer
    for cnf_name in ['vp_start', 'vp_end', 'roi_start', 'roi_end']:
        configs[cnf_name] = int(configs[cnf_name])

    for cnf_name in ['prm_start', 'prm_end']:
        configs[cnf_name] = [int(value) for value in configs[cnf_name]]

    # Check configs that should be of equal length
    linked_configs = [
        ['prm_seq','prm_start','prm_end'],
        ['re_name','re_seq'],
    ]
    for cnf_set in linked_configs:
        assert len(set([len(configs[x]) for x in cnf_set])) == 1, \
            'Error: different lengths for linked configs:'+','.join(str(x) for x in cnf_set)

    # convert chr name to chromosome number
    chr_lst = get_chr_info(genome_str=configs['genome_build'], property='chr_name')
    chr_map = dict(zip(chr_lst, range(1, len(chr_lst) + 1)))
    configs['vp_cnum'] = chr_map[configs['vp_chr']]

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
        configs['output_file'] = configs['output_dir'] + '/rep_' + configs['run_id'] + '_ReadSizeDistribution.pdf'
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


def plot_cirSizeDistribution(configs):
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt, cm

    from utilities import accum_array, hasOL

    # initialization
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/rep_' + configs['run_id'] + '_CirSizeDistribution.pdf'
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


