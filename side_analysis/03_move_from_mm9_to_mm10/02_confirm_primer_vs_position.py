#! /usr/bin/env python2

import sys
from glob import glob
from os import path

sys.path.insert(0, '../../')
from utilities import load_configs

from utilities import get_fasta_sequence, seq_rev_comp as rc

# initialization
cfg_dir = '../../configs/'
# cfg_dir = path.expanduser('~/Downloads/movingfrommm9tomm10/')
glob_ptrn = cfg_dir + 'cfg_Prdm14_*.cfg'
file_lst = glob(glob_ptrn)

for cfg_fname in file_lst:
    print('Checking: {:s}'.format(cfg_fname))
    configs = load_configs(cfg_fname)[0]

    fw_got = get_fasta_sequence(configs['genome_build'], configs['vp_chr'][3:], configs['prm_start'][0], configs['prm_end'][0]).upper()
    assert configs['prm_seq'][1] == fw_got

    rv_got = get_fasta_sequence(configs['genome_build'], configs['vp_chr'][3:], configs['prm_start'][1], configs['prm_end'][1]).upper()
    assert rc(configs['prm_seq'][0]) == rv_got

    print(configs['prm_start'][1] - configs['prm_start'][0])

    