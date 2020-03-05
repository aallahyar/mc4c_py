import numpy as np


def get_chr_info(genome_str, property='chr_name'):
    # chr size: https://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
    chr_details = dict({
        'hg19': dict({
            'chr_name': [
                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
                'chrX', 'chrY', 'chrM'
            ],
            'chr_size': [
                249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,
                135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                155270560, 59373566, 16571
            ]
        }),
        'hg38': dict({
            'chr_name': [
                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
                'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM'
            ],
            'chr_size': [
                248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
                133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569
            ]
        }),
        'mm9': dict({
            'chr_name': [
                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                'chrX', 'chrY', 'chrM'
            ],
            'chr_size': [
                197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172,
                129993255, 121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430,
                166650296, 15902555, 16299
            ]
        }),
        'mm10': dict({
            'chr_name': [
                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                'chrX', 'chrY', 'chrM'
            ],
            'chr_size': [
                195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110,
                130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566,
                171031299, 91744698, 16299,
            ]
        })
    })
    return chr_details[genome_str][property]


def get_re_info(re_name='DpnII', property='seq', genome_str=None):
    re_details = dict({
        'DpnII': dict({'seq': 'GATC'}),
        'MboI': dict({'seq': 'GATC'}),
        'Csp6I': dict({'seq': 'GTAC'}),
        'NlaIII': dict({'seq': 'CATG'}),
        'XbaI': dict({'seq': 'TCTAGA'}),
        'BamHI': dict({'seq': 'GGATCC'}),
        'SacI': dict({'seq': 'GAGCTC'}),
        'PstI': dict({'seq': 'CTGCAG'}),
        'HindIII': dict({'seq': 'AAGCTT'})
    })

    if property == 'pos':
        re_fname = './renzs/{:s}_{:s}.npz'.format(genome_str, re_name)
        chr_lst = get_chr_info(genome_str=genome_str, property='chr_name')

        re_data = np.load(re_fname)['arr_0']
        assert np.array_equal(re_data[1], chr_lst)
        assert re_data[2] == genome_str
        return re_data[0]
    else:
        return re_details[re_name][property]


def extract_re_positions(genome_str, re_name_lst, output_fname=None, ref_fasta=None):
    from os import path, makedirs
    import pysam
    import re

    # Initialization
    chr_lst = get_chr_info(genome_str=genome_str, property='chr_name')
    chr_map = dict(zip(chr_lst, np.arange(len(chr_lst))))

    if output_fname is None:
        output_fname = './renzs/{:s}_{:s}.npz'.format(genome_str, '-'.join(re_name_lst))
    if path.isfile(output_fname):
        print('[w] Restriction enzyme file exists: ' + output_fname)
        return
    if not path.isdir(path.dirname(output_fname)):
        makedirs(path.dirname(output_fname))
    if ref_fasta is None:
        ref_fasta = '../../../datasets/reference_genomes/' + genome_str + '/chrAll.fa.gz'
    print('Searching in the reference genome defined in: ' + ref_fasta)

    # get re sequences
    seq_lst = []
    for re_name in re_name_lst:
        seq_lst.append(get_re_info(genome_str=genome_str, re_name=re_name, property='seq'))
    re_regex = '|'.join(seq_lst)

    # Loop over chromosomes
    re_pos_lst = [None] * len(chr_lst)
    chr_lst_loaded = [None] * len(chr_lst)
    with pysam.FastxFile(path.expanduser(ref_fasta)) as ref_fid:
        print('Scanning chromosomes for restriction recognition sequences: {:s}'.format(', '.join(seq_lst)))
        for chr_ind, chr in enumerate(ref_fid):
            if not chr.name in chr_lst:
                print('\t{:s} is ignored,'.format(chr.name))
                continue
            print('\t{:s},'.format(chr.name))

            cut_sites = []
            for frg in re.finditer(re_regex, chr.sequence, re.IGNORECASE):
                cut_sites.append(frg.start() + 1)
            re_pos_lst[chr_map[chr.name]] = np.array(cut_sites, dtype=np.uint32)
            chr_lst_loaded[chr_map[chr.name]] = chr.name
        if not np.array_equal(chr_lst, chr_lst_loaded):
            raise Exception('[e] Inconsistent reference genome!')
        print('')

    # Save the result
    np.savez(output_fname, [re_pos_lst, chr_lst_loaded, genome_str])


def get_fasta_sequence(genome, chromosome, pos_start, pos_end):
    import urllib2
    from xml.etree import ElementTree
    # example: 'mm9', '1', 13052650, 13052691

    message = 'http://genome.ucsc.edu/cgi-bin/das/{:s}/dna?segment={:s}:{:d},{:d}'.format(
            genome, chromosome, pos_start, pos_end)
    response_xml = urllib2.urlopen(message)
    html = response_xml.read()  # I'm going to assume a safe XML here
    response_tree = ElementTree.fromstring(html)
    return response_tree[0][0].text.replace('\n', '').replace('\r', '')


def seq_complement(seq):
    from string import maketrans

    trans_tbl = maketrans('TCGAtcga', 'AGCTagct')
    return seq.translate(trans_tbl)


def seq_rev_comp(seq):
    return seq_complement(seq)[::-1]


def hasOL(que_item, ref_lst, include_ref_left=False, include_ref_right=False, offset=0):
    if isinstance(que_item, list):
        que_item = np.array(que_item)
    if isinstance(ref_lst, list):
        ref_lst = np.array(ref_lst)
    if ref_lst.ndim == 1:
        ref_lst = ref_lst.reshape(1, -1)
    que_ncol = que_item.shape[0]
    ref_nrow = ref_lst.shape[0]
    assert que_item.ndim == 1, 'Query must be only one element'
    assert que_ncol == ref_lst.shape[1], 'Inconsistency between number of columns in query and reference.'

    result = np.ones(ref_nrow, dtype=bool)
    crd_ind = 0
    if que_ncol == 4:  # Orientation
        result = que_item[3] == ref_lst[:, 3]
    if que_ncol >= 3:  # Chromosome
        result = np.logical_and(result, que_item[0] == ref_lst[:, 0])
        crd_ind = 1
    if include_ref_left:
        OvlL = ref_lst[:, crd_ind] <= que_item[crd_ind+1] + offset
    else:
        OvlL = ref_lst[:, crd_ind] <  que_item[crd_ind+1] + offset
    if include_ref_right:
        OvlR = ref_lst[:, crd_ind+1] >= que_item[crd_ind] - offset
    else:
        OvlR = ref_lst[:, crd_ind+1] >  que_item[crd_ind] - offset
    result = np.logical_and(result, np.logical_and(OvlL, OvlR))
    return result


def which(program):
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def accum_array(group_idx, arr, func=None, default_value=None, min_n_group=None, rebuild_index=False):
    """groups a by indices, and then applies func to each group in turn.
    e.g. func example: [func=lambda g: g] or [func=np.sum] or None for speed
    based on https://github.com/ml31415/numpy-groupies
    """

    if rebuild_index:
        group_idx = np.unique(group_idx.copy(), return_inverse=True)[1]
    if not min_n_group:
        min_n_group = np.max(group_idx) + 1

    order_group_idx = np.argsort(group_idx, kind='mergesort')
    counts = np.bincount(group_idx, minlength=min_n_group)

    if isinstance(arr, np.ndarray):
        groups = np.split(arr[order_group_idx], np.cumsum(counts)[:-1], axis=0)
    else:  # If arr is a Pandas DataFrame
        groups = np.split(arr.loc[order_group_idx,:], np.cumsum(counts)[:-1], axis=0)

    if func:
        ret = [default_value] * min_n_group
        for i, grp in enumerate(groups):
            if len(grp) > 0:
                ret[i] = func(grp)
        return ret
    else:
        return groups


def flatten(nested_lst):
    out_lst = []
    for item in nested_lst:
        if isinstance(item, list):
            out_lst.extend(flatten(item))
        else:
            out_lst.append(item)
    return out_lst


def showprogress(index, n_iter, n_step=10, output_format='{:1.0f}%,'):
    index = index + 1
    if ((index % (n_iter / float(n_step))) - ((index - 1) % (n_iter / float(n_step))) < 0) or (
            n_iter / float(n_step) <= 1):
        print(output_format.format(index * 100 / n_iter)),
        if index == n_iter:
            print


def get_gauss_kernel(size, sigma, ndim=1):
    if ndim == 1:
        if sigma == 0:
            kernel = np.zeros(size)
            kernel[size // 2] = 1
        else:
            kernel = np.exp(-((np.arange(-size // 2 + 1, size // 2 + 1) ** 2) / (2.0 * sigma ** 2)))
    else:
        if sigma == 0:
            kernel = np.zeros([size, size])
            kernel[size // 2, size // 2] = 1
        else:
            dx, dy = np.mgrid[-size // 2 + 1:size // 2 + 1, -size // 2 + 1:size // 2 + 1]
            kernel = np.exp(-((dx ** 2 + dy ** 2) / (2.0 * sigma ** 2)))
    return kernel / kernel.sum()


class OnlineStats(object):
    """
    Welford's algorithm computes the sample variance incrementally.
    from: https://stackoverflow.com/questions/5543651/computing-standard-deviation-in-a-stream
    example:
        from _utilities import OnlineStats
        ostat = OnlineStats(ddof=0)
        for item in item_lst:
            bkg_ostat.include(item)
        assert np.allclose(ostat.mean, np.mean(item_lst))
        assert np.allclose(ostat.std, np.std(item_lst))
    """

    def __init__(self, iterable=None, ddof=0):
        self.ddof, self.n, self.mean, self.M2 = ddof, 0, 0.0, 0.0
        if iterable is not None:
            for datum in iterable:
                self.include(datum)

    def include(self, datum):
        self.n += 1
        self.delta = datum - self.mean
        self.mean += self.delta / self.n
        self.M2 += self.delta * (datum - self.mean)

    @property
    def variance(self):
        return self.M2 / (self.n - self.ddof)

    @property
    def std(self):
        return np.sqrt(self.variance)


def normalize_matrix(mat, method):
    if method == 'KR':
        from krbalancing import kr_balancing
        from scipy import sparse
        mat_nrm = sparse.csr_matrix(mat.copy())
        kr = kr_balancing(mat_nrm.shape[0],
                          mat_nrm.shape[1],
                          mat_nrm.count_nonzero(),
                          mat_nrm.indptr.astype(np.int64, copy=False),
                          mat_nrm.indices.astype(np.int64, copy=False),
                          mat_nrm.data.astype(np.float64, copy=False))
        kr.computeKR()
        mat_nrm = kr.get_normalised_matrix(True).toarray()
        mat_nrm += np.triu(mat_nrm, 1).T
    elif method == 'iterative':
        n_epoch = 100
        # p = np.ones(n_epoch, dtype=np.float)
        # p = np.linspace(0, 1, n_epoch, dtype=np.float)
        # p = np.logspace(np.log10(1.0 / n_epoch), 0, n_epoch, dtype=np.float)
        mat_nrm = mat.copy()
        for ei in range(n_epoch):
            row_sum = np.nansum(mat_nrm, axis=0)
            mat_nrm = mat_nrm / row_sum.reshape(1, -1)
            col_sum = np.nansum(mat_nrm, axis=1)
            mat_nrm = mat_nrm / col_sum.reshape(-1, 1)
            # mat_nrm = (1 - p[ei]) * mat_nrm + p[ei] * mat_nrm / row_sum.reshape(1, -1)
            # mat_nrm = (1 - p[ei]) * mat_nrm + p[ei] * mat_nrm / col_sum.reshape(-1, 1)
    return mat_nrm

################### MC-4C related functions #########################

def load_annotation(genome_str, roi_crd=None):
    import pandas as pd

    # load annotation
    inp_fname = './annotations/ant_{:s}.tsv'.format(genome_str)
    ant_pd = pd.read_csv(inp_fname, delimiter='\t', comment='#')

    # convert map to chr_nums
    chr_lst = get_chr_info(genome_str=genome_str, property='chr_name')
    chr_map = dict(zip(chr_lst, range(1, len(chr_lst) + 1)))
    ant_pd['ant_cnum'] = ant_pd['ant_chr'].map(chr_map)

    # filter annotations outside ROI
    if roi_crd is not None:
        is_in = (ant_pd['ant_cnum'] == roi_crd[0]) & \
                (ant_pd['ant_pos'] >= roi_crd[1]) & \
                (ant_pd['ant_pos'] <= roi_crd[2])
        ant_pd = ant_pd.loc[is_in]

    return ant_pd.reset_index(drop=True)


def load_configs(input_fname, max_n_configs=None):
    """ Read configurations from given file, put it into a dict

    :param input_fname: takes a path to a tab-separated file (or a "config_id") with one variable name and value
    per line, multiple values are seprated by ";").

    :returns: Dictionary where keys are based on the first column with values in a list.
    """
    from os import path

    # check number of given configs
    cfg_file_list = input_fname.split(',')
    if max_n_configs is not None:
        assert len(cfg_file_list) <= max_n_configs, \
            'Maximum of {:d} configs are allowed to be loaded.'.format(max_n_configs)

    # loop over configs
    config_lst = []
    for cfg_fname in cfg_file_list:

        # check if config_file is a file
        if cfg_fname[-4:] != '.cfg':
            cfg_fname = './configs/cfg_' + cfg_fname + '.cfg'
        assert path.isfile(cfg_fname), 'Configuration file could not be found: {:s}'.format(cfg_fname)

        # Load global and then given configs
        configs = dict()
        for fname in ['./mc4c.cfg', cfg_fname]:
            if not path.isfile(fname):
                continue
            with open(fname, 'r') as cfg_fid:
                for line in cfg_fid:
                    if (line[0] == '#') or (len(line) == 1):
                        continue
                    columns = line.rstrip('\n').split('\t')
                    assert len(columns) == 2
                    fld_lst = columns[1].split(',')
                    if len(fld_lst) == 1:
                        configs[columns[0]] = fld_lst[0]
                    else:
                        configs[columns[0]] = fld_lst

        # conversions
        for cfg_name in ['vp_start', 'vp_end', 'roi_start', 'roi_end']:
            if cfg_name in configs.keys():
                configs[cfg_name] = int(configs[cfg_name])
        for cfg_name in ['prm_start', 'prm_end']:
            configs[cfg_name] = [int(value) for value in configs[cfg_name]]
        for cfg_name in ['bwa_index', 'reference_fasta']:
            if cfg_name in configs.keys():
                configs[cfg_name] = configs[cfg_name].replace('%REF%', configs['genome_build'])

        # get chromosome info
        chr_lst = get_chr_info(genome_str=configs['genome_build'], property='chr_name')
        chr_map = dict(zip(chr_lst, range(1, len(chr_lst) + 1)))
        configs['vp_cnum'] = chr_map[configs['vp_chr']]

        # check configs that should be of equal length
        linked_configs = [
            ['prm_seq','prm_start','prm_end'],
            ['re_name','re_seq'],
        ]
        for cnf_set in linked_configs:
            assert len(set([len(configs[x]) for x in cnf_set])) == 1, \
                'Error: different lengths for linked configs:'+','.join(str(x) for x in cnf_set)

        # set default if needed
        roi_cen = int(np.mean([np.min(configs['prm_start']), np.max(configs['prm_end'])]))
        # TODO: Shall we allow the user to choose bin size?
        if 'n_bin' not in configs.keys():
            configs['n_bin'] = 200
        if not np.all([key in configs.keys() for key in ['roi_start', 'roi_end']]):
            configs['roi_start'] = roi_cen - 1000000
            configs['roi_end'] = roi_cen + 1000000
        edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=configs['n_bin'] + 1, dtype=np.int64)
        configs['bin_width'] = edge_lst[1] - edge_lst[0]
        if not np.all([key in configs.keys() for key in ['vp_start', 'vp_end']]):
            configs['vp_start'] = roi_cen - int(configs['bin_width'] * 1.5)
            configs['vp_end'] = roi_cen + int(configs['bin_width'] * 1.5)

        assert (configs['roi_end'] - configs['roi_start'] < 2e6), '[e] ROI can not be defined to be larger than 2mb!'
        assert (configs['roi_end'] - configs['roi_start'] > 120e3), '[e] ROI can not be defined to be smaller than 120kb!'
        assert (configs['n_bin'] >= 100) and (configs['n_bin'] <= 300), \
            '[e] #bins={:d}, #bins should be in the interval of 100 <= #bin <= 300'.format(configs['n_bin'])

        # add to list of configs
        config_lst.append(configs.copy())
    return config_lst


def load_mc4c(config_lst, target_field='frg_np', data_path='./datasets/', verbose=True,
              min_mq=20, valid_only=True, unique_only=True, reindex_reads=True, max_rows=np.inf):
    import pandas as pd
    import h5py

    MAX_N_CIR = 1000000000000
    out_pd = pd.DataFrame()
    if not isinstance(config_lst, list):
        config_lst = [config_lst]

    header_lst = []
    for cfg_idx, configs in enumerate(config_lst):
        if unique_only:
            inp_fname = data_path + '/mc4c_{:s}_{:s}_uniq.hdf5'.format(configs['run_id'], configs['genome_build'])
        else:
            inp_fname = data_path + '/mc4c_{:s}_{:s}_all.hdf5'.format(configs['run_id'], configs['genome_build'])
        if verbose:
            print('Loading {:s} dataset ...'.format(inp_fname))

        h5_fid = h5py.File(inp_fname, 'r')
        if np.isinf(max_rows):
            data_np = h5_fid[target_field][()]
        else:
            print('Selecting only top [{:d}] rows in the dataset'.format(max_rows))
            data_np = h5_fid[target_field][:max_rows]

        header_lst = list(h5_fid[target_field + '_header_lst'][()])
        h5_fid.close()
        part_pd = pd.DataFrame(data_np, columns=header_lst)

        # Filtering fragments
        if min_mq > 0:
            part_pd = part_pd.loc[part_pd['MQ'] >= min_mq]
        if valid_only:
            is_val = np.bitwise_and(part_pd['Flag'], 1) == 0
            part_pd = part_pd.loc[is_val, :]

        # Adjust Read IDs
        assert np.max(part_pd['ReadID']) < MAX_N_CIR
        part_pd['ReadID'] = part_pd['ReadID'] + (cfg_idx + 1) * MAX_N_CIR
        if verbose and (len(config_lst) > 1):
            print('\tGot [{:,d}] reads and [{:,d}] fragments.'.format(len(np.unique(part_pd['ReadID'])), part_pd.shape[0]))

        # Append the part
        if cfg_idx != 0:
            assert np.array_equal(out_pd.columns, part_pd.columns), '[e] Merged datasets do not have identical headers.'
        out_pd = out_pd.append(part_pd, ignore_index=True)
    out_pd = out_pd[header_lst]

    if reindex_reads:
        if verbose:
            print('Reindexing reads ...')
        header_lst.append('ReadID_original')
        out_pd[header_lst[-1]] = out_pd['ReadID'].copy()
        out_pd['ReadID'] = np.unique(out_pd['ReadID'], return_inverse=True)[1] + 1
    if verbose:
        print('In total, [{:,d}] reads and [{:,d}] fragments are loaded.'.format(len(np.unique(out_pd['ReadID'])), out_pd.shape[0]))

    return out_pd[header_lst]


def limit_to_roi(reads, vp_crd=None, roi_crd=None, min_n_frg=2):
    # Reads format: ReadID, Chr, StartCrd, EndCrd
    assert reads.shape[1] >= 4
    n_frg = reads.shape[0]
    if n_frg == 0:
        return np.empty([0, reads.shape[1]])

    is_val = np.ones(n_frg, dtype=np.bool)
    if vp_crd is not None:
        assert len(vp_crd) == 3
        is_val = is_val & ~hasOL(vp_crd, reads[:, 1:4], offset=0)
    if roi_crd is not None:
        assert len(roi_crd) == 3
        is_val = is_val & hasOL(roi_crd, reads[:, 1:4], offset=0)
    reads_roi = reads[is_val, :].copy()

    if min_n_frg is not None:
        read_size = np.bincount(reads_roi[:, 0], minlength=np.max(reads_roi[:, 0]) + 1)[reads_roi[:, 0]]
        reads_roi = reads_roi[read_size >= min_n_frg, :]

    return reads_roi


def get_nreads_per_bin(reads, bin_crd=None, n_bin=None, boundary=None, min_n_frg=None):
    # Reads format: ReadID, Chr, StartCrd, EndCrd
    # Bin format: Chr, StartCrd, EndCrd
    assert reads.shape[1] == 4

    if boundary is None:
        boundary = [bin_crd[0, 0], bin_crd[0, 1], bin_crd[-1, 2]]
    if min_n_frg is not None:
        assert len(boundary) == 3
        reads = limit_to_roi(reads, vp_crd=None, roi_crd=boundary, min_n_frg=min_n_frg)

    if n_bin is not None:
        edge_lst = np.linspace(boundary[1], boundary[2], num=n_bin + 1, dtype=np.int64).reshape(-1, 1)
        bin_crd = np.hstack([np.repeat(boundary[0], n_bin).reshape(-1, 1), edge_lst[:-1], edge_lst[1:] - 1])
    else:
        n_bin = bin_crd.shape[0]
    assert bin_crd.shape[1] == 3
    n_read = len(np.unique(reads[:, 0]))

    # looping over bins
    bin_cvg = np.zeros(n_bin, dtype=np.int)
    for bi in range(n_bin):
        is_in = hasOL(bin_crd[bi, :], reads[:, 1:4])
        bin_cvg[bi] = len(np.unique(reads[is_in, 0]))

    return bin_cvg, n_read


