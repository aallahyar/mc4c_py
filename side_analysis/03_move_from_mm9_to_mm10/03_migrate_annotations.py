#! /usr/bin/env python2

import sys

import pandas as pd

sys.path.insert(0, '../../')
from utilities import get_fasta_sequence

# initialization
genome = 'mm9'
ant_fname = '../../annotations/ant_{:s}.tsv'.format(genome)
roi_bnd = ['chr1', 12840000, 13185000]

# load and filter
ant_pd = pd.read_csv(ant_fname, sep='\t', comment='#')
is_in = (ant_pd['ant_chr'] == roi_bnd[0]) & \
		(ant_pd['ant_pos'] >= roi_bnd[1]) & \
		(ant_pd['ant_pos'] <= roi_bnd[2])
ant_pd = ant_pd.loc[is_in].reset_index(drop=True)

# loop over annotations
for ant_idx, ant_item in ant_pd.iterrows():
	print('>{:s}_{:,d}'.format(*ant_item[['ant_name', 'ant_pos']]))
	ant_seq = get_fasta_sequence(genome, ant_item['ant_chr'], ant_item['ant_pos'] - 30, ant_item['ant_pos'] + 30).upper()
	print('{:s}'.format(ant_seq))
