import numpy as np


def get_chr_info(genome_str, property='chr_name'):
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
        'Csp6I': dict({'seq': 'GTAC'}),
        'NlaIII': dict({'seq': 'CATG'}),
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
    import numpy as np
    from os.path import isfile
    import pysam
    import re

    # Initialization
    chr_lst = get_chr_info(genome_str=genome_str, property='chr_name')
    chr_map = dict(zip(chr_lst, np.arange(len(chr_lst))))

    if output_fname is None:
        output_fname = './renzs/{:s}_{:s}.npz'.format(genome_str, '-'.join(re_name_lst))
    if isfile(output_fname):
        print '[w] Restriction enzyme file exists: ' + output_fname
        return
    if ref_fasta is None:
        ref_fasta = '../../../datasets/reference_genomes/' + genome_str + '/chrAll.fa'
    print 'Searching in reference: ' + ref_fasta

    # get re sequences
    seq_lst = []
    for re_name in re_name_lst:
        seq_lst.append(get_re_info(genome_str=genome_str, re_name=re_name, property='seq'))
    re_regex = '|'.join(seq_lst)

    # Loop over chromosomes
    re_pos_lst = [None] * len(chr_lst)
    chr_lst_loaded = [None] * len(chr_lst)
    with pysam.FastxFile(ref_fasta) as ref_fid:
        for chr_ind, chr in enumerate(ref_fid):
            if not chr.name in chr_lst:
                print '\t[{:s}] is ignored.'.format(chr.name)
                continue
            print '\tScanning [{:s}]'.format(chr.name)

            cut_sites = []
            for frg in re.finditer(re_regex, chr.sequence, re.IGNORECASE):
                cut_sites.append(frg.start() + 1)
            re_pos_lst[chr_map[chr.name]] = np.array(cut_sites, dtype=np.uint32)
            chr_lst_loaded[chr_map[chr.name]] = chr.name
        if not np.array_equal(chr_lst, chr_lst_loaded):
            raise Exception('[e] Inconsistent reference genome!')

    # Save the result
    np.savez(output_fname, [re_pos_lst, chr_lst_loaded, genome_str])


def getFastaSequence(genome, chromosome, pos_start, pos_end):
    import urllib2
    from xml.etree import ElementTree

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
    que_dim = que_item.shape[0]
    [n_ref, ref_dim] = np.shape(ref_lst)
    result = np.ones(n_ref, dtype=bool)
    if que_dim<>ref_dim or que_item.ndim<>1:
        raise ValueError('Query or reference are inconsistent')
    crd_ind = 0

    if que_dim == 4:  # Orientation
        result = que_item[3] == ref_lst[:, 3]
    if que_dim >= 3:  # Chromosome
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

