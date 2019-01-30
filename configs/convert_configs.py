#! /usr/bin/env python

import numpy as np
import pandas as pd

np.set_printoptions(linewidth=180, threshold=5000)  # , suppress=True, formatter={'float_kind':'{:0.5f}'.format}
pd.set_option('display.width', 180)
pd.set_option('display.max_rows', 200)
pd.set_option('display.max_columns', 25)

# initialization
src_fname = '../../64_Cleaning_Up/Dataset_info.tsv'
run_lst = ['LVR-BMaj-96x', 'LVR-BMaj-NP']

# load source file
vpi_pd = pd.read_csv(src_fname, delimiter='\t')

# loop over rows
for run_id in run_lst:
    di = np.where(vpi_pd['id'] == run_id)[0]
    assert len(di) == 1

    vp_info = vpi_pd.loc[di[0]]
    out_fname = './cfg_{:s}.cfg'.format(vp_info['id'])
    print 'writing config file to: ' + out_fname
    with open(out_fname, 'w') as out_fid:
        out_fid.write('run_id\t{:s}\n'.format(vp_info['id']))
        out_fid.write('vp_chr\tchr{:0.0f}\n'.format(vp_info['Target_Chr']))
        out_fid.write('prm_start\t{:0.0f};{:0.0f}\n'.format(vp_info['PR1_Be'], vp_info['PR2_Be']))
        out_fid.write('prm_end\t{:0.0f};{:0.0f}\n'.format(vp_info['PR1_En'], vp_info['PR2_En']))
        out_fid.write('prm_seq\t{:s};{:s}\n'.format(vp_info['PR1_Seq'], vp_info['PR2_Seq']))
        out_fid.write('re_name\t{:s};{:s}\n'.format(vp_info['Digest1'], vp_info['Digest_2']))
        out_fid.write('re_seq\t{:s};{:s}\n'.format(vp_info['Digest_1_Seq'], vp_info['Digest_2_Seq']))
        out_fid.write('roi_start\t{:0.0f}\n'.format(vp_info['Win_be']))
        out_fid.write('roi_end\t{:0.0f}\n'.format(vp_info['Win_en']))
        out_fid.write('genome_build\t{:s}\n'.format(vp_info['Genome_Build']))


