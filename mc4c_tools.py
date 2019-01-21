import numpy as np
import pandas as pd
import h5py


def get_mchc_data(config_lst, target_field='frg_np', data_path='./mc4c_files/',
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


def load_configs(cnf_fname):
    """ Read configurations from given file, put it into a dict

    :param cnf_fname: takes a path to a tab-separated file with one variable name and value per line, multiple values
    are seprated by ";").

    :returns: Dictionary where keys are based on the first column with values in a list.
    """
    from utilities import get_chr_info

    # Load global and then given configs
    configs = dict()
    for fname in ['./mc4c.cnf', cnf_fname]:
        with open(fname, 'r') as cnf_fid:
            for line in cnf_fid:
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
