import numpy as np
import pandas as pd
import h5py


def load_mc4c(config_lst, target_field='frg_np', data_path='./datasets/',
              min_mq=20, only_unique=True, reindex_reads=True, max_rows=np.inf):
    MAX_N_CIR = 100000000
    out_pd = pd.DataFrame()
    if not isinstance(config_lst, list):
        config_lst = [config_lst]

    header_lst = []
    for vi, vp_info in enumerate(config_lst):
        tsv_name = data_path + 'frg_{:s}_{:s}.hdf5'.format(vp_info['run_id'], vp_info['genome'])
        print('Loading [{:d}: {:s}] data from: {:s}'.format(vp_info['row_index'], vp_info['run_id'], tsv_name))

        h5_fid = h5py.File(tsv_name, 'r')
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
            part_pd = part_pd.loc[part_pd['isUnique'] > 0].copy()

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


def plot_ReadSizeDistribution(configs):
    import platform
    if platform.system() == 'Linux':
        import matplotlib
        matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    from os import path
    import pysam

    # initializations
    if configs['input_file'] is None:
        configs['input_file'] = './fastqs/raw_' + configs['run_id'] + '.fastq.gz'
    if configs['output_file'] is None:
        configs['output_file'] = configs['output_dir'] + '/rep_' + configs['run_id'] + '_ReadSizeDistribution.pdf'
    MAX_SIZE = 25000
    edge_lst = np.linspace(0, MAX_SIZE, 51)
    n_bin = len(edge_lst) - 1

    # loop over reads
    nbp_total = 0
    nbp_inf = 0
    nrd_inf = 0
    n_read = 0
    size_dist = np.zeros(n_bin, dtype=np.int64)
    print 'Computing read size for: {:s}'.format(configs['input_file'])
    with pysam.FastxFile(configs['input_file'], persist=False) as gz_fid:
        for rd_idx, read in enumerate(gz_fid):
            if rd_idx % 50000 == 0:
                print('{:,d} reads are processed.'.format(rd_idx))
            seq_size = len(read.sequence)

            n_read += 1
            nbp_total = nbp_total + seq_size
            if seq_size > 1500:
                nrd_inf += 1
                nbp_inf = nbp_inf + seq_size

            if seq_size >= MAX_SIZE:
                seq_size = MAX_SIZE - 1
            bin_idx = np.digitize(seq_size, edge_lst) - 1
            size_dist[bin_idx] += 1

    # plotting
    plt.figure(figsize=(7, 5))
    plt.bar(range(n_bin), size_dist, width=0.95)
    x_ticks_idx = range(1, n_bin, 2)
    plt.xticks(x_ticks_idx, ['{:0.0f}'.format(edge_lst[i + 1] / 1e3) for i in x_ticks_idx],
               rotation=0, fontsize=10)
    plt.xlim([-1, n_bin])
    plt.xlabel('#base pairs (kbp)')

    y_ticks = plt.yticks()[0]
    y_tick_lbl = ['{:0,.0f}k'.format(x / 1e3) for x in y_ticks]
    plt.yticks(y_ticks, y_tick_lbl)
    plt.ylabel('#Reads')

    plt.title('Read size distribution, {:s}\n'.format(configs['run_id']) +
              '#read={:,d}; #read (>1.5kb)={:,d}\n'.format(n_read, nrd_inf) +
              '#bases={:,d}; #bases (>1.5kb)={:,d}'.format(nbp_total, nbp_inf)
              )

    plt.savefig(configs['output_file'], bbox_inches='tight')


def plot_FragSizeDistribution(configs):
    print 'Fragment size distribution'

