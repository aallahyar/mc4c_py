#! /usr/bin/env python2

import sys

import numpy as np

sys.path.insert(0, '../../')
from utilities import load_configs

from utilities import get_fasta_sequence, seq_rev_comp as rc

# initialization
run_id = 'Prdm14_RB_LB-DEL'
cfg_fname = '../../configs/cfg_' + run_id + '.cfg'
configs = load_configs(cfg_fname)[0]
# roi_seq_mm9 = get_fasta_sequence('mm9', configs['vp_chr'][3:], configs['roi_start'], configs['roi_end']).upper()
with open('./roi_seq_mm9.txt', 'r') as file:
	roi_seq_mm9 = file.read().replace('\n', '')
mm9_prm_fwSIdx = roi_seq_mm9.find(configs['prm_seq'][0])
mm9_vpf_endIdx = roi_seq_mm9.find(configs['re_seq'][0], mm9_prm_fwSIdx) + len(configs['re_seq'][0])
mm9_vpf_begIdx = roi_seq_mm9.rfind(configs['re_seq'][0], 0, mm9_prm_fwSIdx)
mm9_vpf_seq = roi_seq_mm9[mm9_vpf_begIdx:mm9_vpf_endIdx]
assert mm9_vpf_seq[5:-5].find(configs['re_seq'][0]) == -1
mm9_prm_fwSeq = roi_seq_mm9[mm9_prm_fwSIdx:mm9_vpf_endIdx]
print('=== mm9: {:s}:{:,d}-{:,d}'.format(configs['vp_chr'][3:], configs['roi_start'], configs['roi_end']))
print('FW:\n{:s}\n{:s}\n'.format(configs['prm_seq'][0][:30], mm9_prm_fwSeq[:30]))

mm9_prm_rvSIdx = mm9_vpf_seq.find(rc(configs['prm_seq'][1])) + len(configs['prm_seq'][1])
mm9_prm_rvSeq = mm9_vpf_seq[:mm9_prm_rvSIdx]
mm9_tmp = ' ' * (len(mm9_prm_rvSeq) - len(configs['prm_seq'][1])) + rc(configs['prm_seq'][1])
print('Rev:\n{:s}\n{:s}'.format(mm9_tmp[-30:], mm9_prm_rvSeq[-30:]))

# switch to mm10
mm10_bnd = [12840000, 13185000]
# roi_seq_mm10 = get_fasta_sequence('mm10', configs['vp_chr'][3:], mm10_bnd[0], mm10_bnd[1]).upper()
with open('./roi_seq_mm10.txt', 'r') as file:
	mm10_roi_seq = file.read().replace('\n', '')
	
mm10_prm_fwSIdx = mm10_roi_seq.find(configs['prm_seq'][0])
mm10_vpf_begIdx = mm10_roi_seq.rfind(configs['re_seq'][0], 0, mm10_prm_fwSIdx)
mm10_vpf_endIdx = mm10_roi_seq.find(configs['re_seq'][0], mm10_prm_fwSIdx) + len(configs['re_seq'][0])
mm10_vpf_seq = mm10_roi_seq[mm10_vpf_begIdx:mm10_vpf_endIdx]
assert mm10_vpf_seq[5:-5].find(configs['re_seq'][0]) == -1

mm10_prm_fwSeq = mm10_roi_seq[mm10_prm_fwSIdx:mm10_vpf_endIdx]
mm10_prm_rvSIdx = mm10_roi_seq.find(rc(configs['prm_seq'][1]))
mm10_prm_rvSeq = mm10_roi_seq[mm10_vpf_begIdx:mm10_prm_rvSIdx + len(configs['prm_seq'][1])]
assert mm10_prm_fwSeq[:len(configs['prm_seq'][0])] == configs['prm_seq'][0]
assert mm10_prm_rvSeq[-len(configs['prm_seq'][1]):] == rc(configs['prm_seq'][1])

mm10_prm_fwPosBeg = mm10_bnd[0] + mm10_prm_fwSIdx
mm10_prm_fwPosEnd = mm10_prm_fwPosBeg + len(configs['prm_seq'][0]) - 1
mm10_prm_rvPosBeg = mm10_bnd[0] + mm10_prm_rvSIdx
mm10_prm_rvPosEnd = mm10_prm_rvPosBeg + len(configs['prm_seq'][1]) - 1
assert configs['prm_seq'][0] == get_fasta_sequence('mm10', configs['vp_chr'][3:], mm10_prm_fwPosBeg, mm10_prm_fwPosEnd).upper()
assert rc(configs['prm_seq'][1]) == get_fasta_sequence('mm10', configs['vp_chr'][3:], mm10_prm_rvPosBeg, mm10_prm_rvPosEnd).upper()

print('=== mm10: {:s}:{:,d}-{:,d}'.format(configs['vp_chr'][3:], mm10_bnd[0], mm10_bnd[1]))
print('\tFw: {:d} - {:d}'.format(mm10_prm_fwPosBeg, mm10_prm_fwPosEnd))
print('\tRv: {:d} - {:d}'.format(mm10_prm_rvPosBeg, mm10_prm_rvPosEnd))

