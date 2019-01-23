#! /usr/bin/env python

import argparse
import sys
import gzip
from os import path, makedirs
import numpy as np
import pandas as pd

import mc4c_tools

flag_DEBUG = True
np.set_printoptions(linewidth=180, threshold=5000)  # , suppress=True, formatter={'float_kind':'{:0.5f}'.format}
pd.set_option('display.width', 180)
pd.set_option('display.max_rows', 200)
pd.set_option('display.max_columns', 15)


def initialize_run(args):
    # try:
    #     configs = mc4c_tools.load_configs(args.config_file)
    # except:
    #     raise Exception('Configuration file is missing required fields or not formatted correctly.')

    # TODO:
    # test refrence genome path
    # test existance of bwa
    # test if bwa index is present
    # check and make RE position file if non existant
    # test existance of samtools
    print '[TO DO] Checking the pipeline for related files for given experiment'


def setReadIds(args):
    print '%% Assigning traceable identifiers to reads ...'

    configs = mc4c_tools.load_configs(args.config_file)

    # initialize
    if args.input_file is None:
        args.input_file = './fastqs/raw_' + configs['run_id'] + '.fastq.gz'
    if args.output_file is None:
        args.output_file = './reads/rd_' + configs['run_id'] + '.fasta.gz'
    # assert not path.isfile(args.output_file), '[e] output file already exists: {:s}'.format(args.output_file)
    if not path.isdir(path.dirname(args.output_file)):
        makedirs(path.dirname(args.output_file))
    print('Writing reads with traceable identifiers to: {:s}'.format(args.output_file))

    # loop over reads
    with gzip.open(args.output_file, 'w') as out_fid:
        inp_flst = args.input_file.split(';')
        print 'Total of [{:d}] files are given as input.'.format(len(inp_flst))
        for inp_fidx, inp_fname in enumerate(inp_flst):
            print('\t{:d}. Reading from: {:s}'.format(inp_fidx + 1, inp_fname))
            rd_idx = 0
            with gzip.open(inp_fname, 'r') as inp_fid:
                while True:
                    rd_idx = rd_idx + 1
                    rd_oid = inp_fid.readline().rstrip('\n')
                    rd_seq = inp_fid.readline().rstrip('\n')
                    rd_plus = inp_fid.readline().rstrip('\n')
                    rd_pred = inp_fid.readline().rstrip('\n')
                    if rd_oid == '':
                        break
                    if rd_oid[0] != '@' or rd_plus != '+':
                        raise Exception('[e] the input file is corrupted.\n' +
                                        'Read #{:d}:\n'.format(rd_idx) +
                                        '\tID: [{:s}],\n\tplus: [{:s}]'.format(rd_oid, rd_plus))
                    if rd_idx % 10000 == 0:
                        print('\t\tprocessed {:,d} reads.'.format(rd_idx))

                    rd_sid = 'Fl.Id:{:d};Rd.Id:{:d};Rd.Ln:{:d}'.format(inp_fidx + 1, rd_idx, len(rd_seq))
                    out_fid.write('>' + rd_sid + '\n')
                    out_fid.write(rd_seq + '\n')
    print '[i] Read identifier conversion is completed successfully.'


def splitReads(args):
    from utilities import get_re_info
    import re

    print '%% Splitting reads into fragments ...'
    configs = mc4c_tools.load_configs(args.config_file)

    if args.input_file is None:
        args.input_file = './reads/rd_' + configs['run_id'] + '.fasta.gz'
    if args.output_file is None:
        args.output_file = './frgments/frg_' + configs['run_id'] + '.fasta.gz'
    if not path.isdir(path.dirname(args.output_file)):
        makedirs(path.dirname(args.output_file))
    # assert not path.isfile(args.output_file), '[e] Output file already exists: {:s}'.format(args.output_file)
    print('Reading reads from: {:s}'.format(args.input_file))
    print('Writing fragments to: {:s}'.format(args.output_file))

    MAX_FRG_SIZE = 2000
    re_seq_lst = [get_re_info(re_name=re_name, property='seq') for re_name in configs['re_name']] + ['$']
    regex_ptr = '|'.join(re_seq_lst)

    rd_ind = 1
    frg_ind = 1
    n_reduced = 0
    with gzip.open(args.input_file, 'r') as inp_fid, \
            gzip.open(args.output_file, 'w') as out_fid:
        while True:
            rd_sid = inp_fid.readline().rstrip('\n')
            rd_seq = inp_fid.readline().rstrip('\n')
            if rd_sid == '':
                break
            if rd_ind % 10000 == 0:
                print('\tprocessed {:,d} reads and produced {:,d} fragments.'.format(rd_ind, frg_ind))

            frg_be = 0
            for res_enz in re.finditer(regex_ptr, rd_seq):
                frg_en = res_enz.end()
                if frg_en - frg_be > MAX_FRG_SIZE:
                    n_reduced += 1
                    frg_en = MAX_FRG_SIZE
                out_fid.write('{:s};Fr.Id:{:d};Fr.SBp:{:d};Fr.EBp:{:d}\n'.format(rd_sid, frg_ind, frg_be + 1, frg_en))
                out_fid.write(rd_seq[frg_be:frg_en] + '\n')
                frg_ind = frg_ind + 1
                frg_be = res_enz.start()
            rd_ind = rd_ind + 1
    if n_reduced != 0:
        print '[i] [{:,d}] fragments are reduced to {:,d}bp.'.format(n_reduced, MAX_FRG_SIZE)
    print '[i] Total of {:,d} reads and {:,d} fragments are produced successfully.'.format(rd_ind, frg_ind)


def mapFragments(args):
    print '%% Mapping fragments to genome ...'

    configs = mc4c_tools.load_configs(args.config_file)

    # Map split fragments to genome
    if args.input_file is None:
        args.input_file = './fragments/frg_' + configs['run_id'] + '.fasta.gz'
    if args.output_file is None:
        args.output_file = './bams/bam_{:s}.bam'.format(configs['run_id'])
    if not path.isdir(path.dirname(args.output_file)):
        makedirs(path.dirname(args.output_file))
    # assert not path.isfile(args.output_file)
    print('Reading fragments from: {:s}'.format(args.input_file))
    print('Writing mapped fragments to: {:s}'.format(args.output_file))

    cmd_str = \
        configs['bwa_path'] + ' bwasw -b 5 -q 2 -r 1 -z 5 -T 15 -t {:d} '.format(args.n_thread) + \
        configs['bwa_index_path'].replace('%%', configs['genome_build']) + ' ' + args.input_file + \
        ' | samtools view -q 1 -hbS - ' + \
        '> ' + args.output_file

    if args.return_command:
        print '{:s}'.format(cmd_str)
    else:
        print 'Running BWA by: {:s}'.format(cmd_str)
        import subprocess
        map_prs = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE)
        std_out, std_err = map_prs.communicate()
        if std_err is not None:
            raise Exception('[e] Error: BWA failed to run properly.')
        print '[i] Fragments are mapped to genome successfully.'


def processMappedFragments(args):
    import h5py
    from os import remove
    from pandas import read_csv
    import pysam

    from utilities import get_chr_info, hasOL

    configs = mc4c_tools.load_configs(args.config_file)

    if args.input_file is None:
        args.input_file = './bams/bam_{:s}.bam'.format(configs['run_id'])
    if args.output_file is None:
        args.output_file = './datasets/mc4c_' + configs['run_id'] + '_all.hdf5'
    if not path.isdir(path.dirname(args.output_file)):
        makedirs(path.dirname(args.output_file))
    # assert not path.isfile(args.output_file), '[e] Output file already exists: {:s}'.format(args.output_file)
    print('Reading fragments from: {:s}'.format(args.input_file))

    chr_lst = get_chr_info(configs['genome_build'], 'chr_name')
    chr_map = dict(zip(chr_lst, np.arange(len(chr_lst)) + 1))

    # extend fragment coordinates to reference genome
    re_pos_fname = './renz_files/{:s}_{:s}.npz'.format(configs['genome_build'], '-'.join(configs['re_name']))
    re_pos, re_chr_lst = np.load(re_pos_fname)['arr_0']
    assert np.array_equal(re_chr_lst, chr_lst)

    # Read file line by line
    ReadID_old = -1
    FrgID_old = -1
    n_processed = 0
    header_lst = ['ReadID', 'Chr', 'MapStart', 'MapEnd', 'Strand', 'ExtStart', 'ExtEnd', 'MQ',
                  'FileID', 'FrgID', 'SeqStart', 'SeqEnd', 'ReadLength']
    n_header = len(header_lst)
    tmp_fname = args.output_file + '.tmp'
    print('Writing processed fragments to a temprary file first: {:s}'.format(tmp_fname))
    with pysam.AlignmentFile(args.input_file, 'rb') as bam_fid, gzip.open(tmp_fname, 'wb') as gz_fid:
        gz_fid.write('\t'.join(header_lst) + '\n')
        frg_template = '\t'.join(['{:d}'] * n_header) + '\n'

        frg_set = np.empty([0, n_header], dtype=np.int64)
        for que_idx, que_line in enumerate(bam_fid):
            if que_idx % 50000 == 0:
                print('\tprocessed {:,d} fragments in {:,d} reads.'.format(que_idx, n_processed))
            if (np.bitwise_and(que_line.flag, 0x800) == 0x800) or (que_line.reference_name not in chr_lst):
                continue
            FileID, ReadID, ReadLength, FrgID, SeqStart, SeqEnd = \
                [int(x.split(':')[1]) for x in que_line.query_name.split(';')]
            MapChrNid = chr_map[que_line.reference_name]
            MapStart = que_line.reference_start
            MapEnd = que_line.reference_end
            MapStrand = 1 - (que_line.is_reverse * 2)

            # extending coordinates to nearby restriction site
            nei_left = np.searchsorted(re_pos[MapChrNid - 1], MapStart, side='left') - 1
            if np.abs(re_pos[MapChrNid - 1][nei_left + 1] - MapStart) < 10:
                nei_left = nei_left + 1
            nei_right = np.searchsorted(re_pos[MapChrNid - 1], MapEnd, side='left')
            if np.abs(MapEnd - re_pos[MapChrNid - 1][nei_right - 1]) < 10:
                nei_right = nei_right - 1

            if nei_left == nei_right:
                dist_left = np.abs(MapStart - re_pos[MapChrNid - 1][nei_left])
                dist_right = np.abs(MapEnd - re_pos[MapChrNid - 1][nei_right])
                if dist_right < dist_left:
                    nei_left = nei_left - 1
                else:
                    nei_right = nei_right + 1
            ExtStart = re_pos[MapChrNid - 1][nei_left]
            ExtEnd = re_pos[MapChrNid - 1][nei_right]
            # TODO: Needs to account for unmapped fragments

            # combine into an array
            frg_info = np.array([
                ReadID, MapChrNid, MapStart, MapEnd, MapStrand, ExtStart, ExtEnd, que_line.mapping_quality,
                FileID, FrgID, SeqStart, SeqEnd, ReadLength]).reshape([1, -1])

            # Check order of fragments
            if FrgID < FrgID_old:
                raise Exception('Order of fragments are lost.')
            FrgID_old = FrgID

            # Check if read has ended
            if ReadID == ReadID_old:
                frg_set = np.vstack([frg_set, frg_info])
                continue

            # Save the read
            fi = 0
            while fi < frg_set.shape[0] - 1:
                if hasOL(frg_set[fi, [1, 5, 6]], frg_set[fi + 1:fi + 2, [1, 5, 6]], offset=20)[0]:
                    # this still occurs if A1- and A2+ are adjacent (ignoring strand)
                    frg_set[fi, 2] = np.min(frg_set[fi:fi + 2, 2])
                    frg_set[fi, 3] = np.max(frg_set[fi:fi + 2, 3])
                    frg_set[fi, 5] = np.min(frg_set[fi:fi + 2, 5])
                    frg_set[fi, 6] = np.max(frg_set[fi:fi + 2, 6])
                    frg_set[fi, 7] = np.max(frg_set[fi:fi + 2, 7])
                    frg_set[fi, 10] = np.min(frg_set[fi:fi + 2, 10])
                    frg_set[fi, 11] = np.max(frg_set[fi:fi + 2, 11])
                    frg_set = np.delete(frg_set, fi + 1, axis=0)
                else:
                    fi = fi + 1

            for frg in frg_set:
                gz_fid.write(frg_template.format(*frg))
            frg_set = frg_info.copy()
            ReadID_old = ReadID

            n_processed += 1

        if frg_set.shape[0] != 0:  # Saving the last read after file has finished
            for frg in frg_set:
                gz_fid.write(frg_template.format(*frg))

    # Load fragments in pandas format and sort fragments within a read according to their relative positions
    print('Loading temporary file: {:s}'.format(tmp_fname))
    frg_pd = read_csv(tmp_fname, delimiter='\t', compression='gzip')
    frg_pd = frg_pd.iloc[np.lexsort([frg_pd['SeqStart'], frg_pd['ReadID']])]

    print('Writing dataset to: {:s}'.format(args.output_file))
    h5_fid = h5py.File(args.output_file, 'w')
    if frg_pd.shape[0] > 5000:
        h5_fid.create_dataset('frg_np', data=frg_pd.values, compression='gzip', compression_opts=5,
                              chunks=(5000, frg_pd.shape[1]))
    else:
        h5_fid.create_dataset('frg_np', data=frg_pd.values, compression='gzip', compression_opts=5)
    h5_fid.create_dataset('frg_np_header_lst', data=list(frg_pd.columns.values))
    h5_fid.create_dataset('chr_lst', data=chr_lst)
    h5_fid.close()

    print('Removing temporary file: {:s}'.format(tmp_fname))
    remove(tmp_fname)
    print '[i] MC4C dataset is created successfully.'


def removeDuplicates(args):
    import h5py
    from utilities import hasOL

    configs = mc4c_tools.load_configs(args.config_file)

    if args.input_file is None:
        args.input_file = './datasets/mc4c_' + configs['run_id'] + '_all.hdf5'
    if args.output_file is None:
        args.output_file = './datasets/mc4c_' + configs['run_id'] + '_uniq.hdf5'
    if not path.isdir(path.dirname(args.output_file)):
        makedirs(path.dirname(args.output_file))
    print('Reading MC4C dataset from: {:s}'.format(args.input_file))
    print('Writing unique MC4C reads to: {:s}'.format(args.output_file))

    # load mc4c data
    h5_fid = h5py.File(args.input_file, 'r')
    target_field = 'frg_np'
    data_np = h5_fid[target_field].value
    header_lst = list(h5_fid[target_field + '_header_lst'].value)
    mc4c_pd = pd.DataFrame(data_np, columns=header_lst)
    chr_lst = list(h5_fid['chr_lst'].value)
    h5_fid.close()
    MAX_ReadID = np.max(mc4c_pd['ReadID'])
    print 'There are {:d} reads in the dataset.'.format(len(np.unique(mc4c_pd['ReadID'])))

    # Filtering reads according to their MQ
    print 'Ignoring fragments with MQ < {:d}:'.format(args.min_mq)
    frg_trs = mc4c_pd[['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'MQ']].values
    frg_trs = frg_trs[frg_trs[:, 4] >= args.min_mq, :4]
    print '{:d} reads are left after MQ filtering.'.format(len(np.unique(frg_trs[:, 0])))

    # select far-cis/trans fragments
    roi_size = configs['roi_end'] - configs['roi_start']
    local_area = np.array([configs['vp_cnum'], configs['roi_start'] - roi_size, configs['roi_end'] + roi_size])
    is_lcl = hasOL(local_area, frg_trs[:, 1:4])
    frg_trs = frg_trs[np.isin(frg_trs[:, 0], frg_trs[~is_lcl, 0]), :]
    print 'Selecting for reads with #trans-fragment > 0: {:d} reads are selected.'.format(len(np.unique(frg_trs[:, 0])))

    # count #roi-frag per read
    roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
    is_roi = hasOL(roi_crd, frg_trs[:, 1:4], offset=0)
    frg_roi = frg_trs[is_roi, :]
    roi_freq = np.bincount(frg_roi[:, 0], minlength=MAX_ReadID + 1).reshape(-1, 1)

    # count #tras-frag per read
    trs_freq = np.bincount(frg_trs[:, 0], minlength=MAX_ReadID + 1).reshape(-1, 1)
    frg_trs = np.hstack([frg_trs, roi_freq[frg_trs[:, 0], :], trs_freq[frg_trs[:, 0], :]])

    # filter reads with #roi-frg > 1
    frg_trs = frg_trs[frg_trs[:, 4] > 1, :]
    print 'Selecting for reads with #roi-fragment > 1: {:d} reads are left.'.format(len(np.unique(frg_trs[:, 0])))

    # sort reads according to #trans
    trs_sid = np.lexsort([frg_trs[:, 0], frg_trs[:, -1]])[::-1]
    frg_trs = frg_trs[trs_sid, :]

    # loop over trans fragments
    trs_idx = 0
    print 'Scanning for duplicated trans-fragments:'
    while trs_idx < frg_trs.shape[0]:
        if trs_idx % 10000 == 0:
            print '\tscanned {:,d} fragments for duplicates.'.format(trs_idx)
        has_ol = hasOL(frg_trs[trs_idx, 1:4], frg_trs[:, 1:4], offset=-10)
        if np.sum(has_ol) > 1:
            # select duplicates
            frg_dup = frg_trs[np.isin(frg_trs[:, 0], frg_trs[has_ol, 0]), :]

            # keep largest read
            keep_rid = frg_dup[np.argmax(frg_dup[:, 4]), 0]
            frg_dup = frg_dup[frg_dup[:, 0] != keep_rid, :]

            # remove extra duplicates
            frg_trs = frg_trs[~ np.isin(frg_trs[:, 0], frg_dup[:, 0]), :]
        trs_idx = trs_idx + 1
    print 'Result statistics:'
    print '\t#reads: {:,d} --> {:,d}'.format(len(np.unique(mc4c_pd['ReadID'])), len(np.unique(frg_trs[:, 0])))
    print '\t#fragments: {:,d} --> {:,d}'.format(mc4c_pd['ReadID'].shape[0], frg_trs.shape[0])

    # select and save unique reads
    is_uniq = np.isin(mc4c_pd['ReadID'], frg_trs[:, 0])
    uniq_pd = mc4c_pd.loc[is_uniq, :]
    print('Writing dataset to: {:s}'.format(args.output_file))
    h5_fid = h5py.File(args.output_file, 'w')
    h5_fid.create_dataset('frg_np', data=uniq_pd.values, compression='gzip', compression_opts=5)
    h5_fid.create_dataset('frg_np_header_lst', data=list(uniq_pd.columns.values))
    h5_fid.create_dataset('chr_lst', data=chr_lst)
    h5_fid.close()
    print '[i] PCR duplicates are removed from MC4C dataset successfully.'


def getSumRep(args):
    import mc4c_tools

    configs = mc4c_tools.load_configs(args.config_file)

    if args.input_file is None:
        args.input_file = './fastqs/raw_' + configs['run_id'] + '.fastq.gz'
    if args.output_file is None:
        configs['output_dir'] = './plots/'
    else:
        configs['output_dir'] = path.dirname(args.output_file)
    if not path.isdir(configs['output_dir']):
        makedirs(configs['output_dir'])
    configs['input_file'] = args.input_file
    configs['output_file'] = args.output_file
    print('Reading MC4C dataset from: {:s}'.format(args.input_file))


    # read size distribution
    if args.report_type == 'readSizeDist':
        mc4c_tools.plot_ReadSizeDistribution(configs)
    elif args.report_type == 'fragSizeDist':
        mc4c_tools.plot_FragSizeDistribution(configs)
    else:
        raise Exception()

# Huge wall of argparse text starts here
def main():
    """ Everything in here is to interpret calls from a command line.
        Anything being run will just call a similarly named function
        above.
    """

    parser = argparse.ArgumentParser(description="MC4C pipeline for processing multi-contact 4C data")
    subparsers = parser.add_subparsers()

    # init command
    # parser_init = subparsers.add_parser('initRun',
    #                                       description='Initialize and prepare the pipeline for given configuration')
    # parser_init.add_argument('config_file', metavar='config-file',
    #                            type=str,
    #                            help='Configuration file containing experiment specific details')
    # parser_init.set_defaults(func=initialize_run)

    # Set read identifiers
    parser_readid = subparsers.add_parser('setReadIds',
                                          description='Defining identifiers for sequenced reads')
    parser_readid.add_argument('config_file', metavar='config-file',
                               type=str,
                               help='Configuration file containing experiment specific details')
    parser_readid.add_argument('--input-file',
                               default=None,
                               type=str,
                               help='Input file (in FASTQ format) containing raw sequenced reads')
    parser_readid.add_argument('--output-file',
                               default=None,
                               type=str,
                               help='Output file (in FASTA format) containing sequenced reads with traceable IDs')
    parser_readid.set_defaults(func=setReadIds)

    # Split reads into fragments
    parser_readSplt = subparsers.add_parser('splitReads',
                                            description='Splitting reads into fragments using restriction enzyme recognition sequence')
    parser_readSplt.add_argument('config_file', metavar='config-file',
                                 type=str,
                                 help='Configuration file containing experiment specific details')
    parser_readSplt.add_argument('--input-file',
                                 default=None,
                                 type=str,
                                 help='Input file (in FASTA format) containing reads with traceable IDs.')
    parser_readSplt.add_argument('--output-file',
                                 default=None,
                                 type=str,
                                 help='Output file (in FASTA format) containing fragments with traceable IDs')
    parser_readSplt.set_defaults(func=splitReads)

    # Mapping fragments to reference genome
    parser_mapFrg = subparsers.add_parser('mapFragments',
                                          description='Mapping fragments to reference genome')
    parser_mapFrg.add_argument('config_file', metavar='config-file',
                               type=str,
                               help='Configuration file containing experiment specific details')
    parser_mapFrg.add_argument('--input-file',
                               default=None,
                               type=str,
                               help='Input file (in FASTA format) containing fragments with traceable IDs')
    parser_mapFrg.add_argument('--output-file',
                               default=None,
                               type=str,
                               help='Output file (in BAM format) containing fragments with traceable IDs')
    parser_mapFrg.add_argument('--n_thread',
                               default=6,
                               type=int,
                               help='Number of threads should be used by the aligner')
    parser_mapFrg.add_argument('--return_command',
                               action="store_true",
                               help='Return only mapping command instead of running it ' +
                                    '(useful for running the pipeline in a cluster)')
    parser_mapFrg.set_defaults(func=mapFragments)

    # Process mapped fragments
    parser_mkDataset = subparsers.add_parser('makeDataset',
                                             description='Processed the mapped fragments and create a MC-4C dataset')
    parser_mkDataset.add_argument('config_file', metavar='config-file',
                                  type=str,
                                  help='Configuration file containing experiment specific details')
    parser_mkDataset.add_argument('--input-file',
                                  default=None,
                                  type=str,
                                  help='Input file (in BAM format) containing fragments with traceable IDs')
    parser_mkDataset.add_argument('--output-file',
                                  default=None,
                                  type=str,
                                  help='Output file (in HDF5 format) containing processed fragments')
    parser_mkDataset.set_defaults(func=processMappedFragments)

    # Remove PCR duplicated
    parser_remDup = subparsers.add_parser('removeDuplicates',
                                          description='Remove PCR duplicates from a given MC-4C dataset')
    parser_remDup.add_argument('config_file', metavar='config-file',
                               type=str,
                               help='Configuration file containing experiment specific details')
    parser_remDup.add_argument('--input-file',
                               default=None,
                               type=str,
                               help='Input file (in HDF5 format) containing MC4C data.')
    parser_remDup.add_argument('--output-file',
                               default=None,
                               type=str,
                               help='Output file (in HDF5 format) containing MC4C data.')
    parser_remDup.add_argument('--min-mq',
                               default=20,
                               type=int,
                               help='Minimum mapping quality (MQ) to consider a fragment as confidently mapped.')
    parser_remDup.set_defaults(func=removeDuplicates)

    # produce statistics plots
    parser_sumReport = subparsers.add_parser('getSumRep',
                                            description='Generate various summary reports about a MC-4C dataset.')
    parser_sumReport.add_argument('report_type', metavar='report-type',
                                 choices=['readSizeDist', 'fragSizeDist'],
                                 type=str,
                                 help='Type of summary report that needs to be generated')
    parser_sumReport.add_argument('config_file', metavar='config-file',
                                 type=str,
                                 help='Configuration file containing experiment specific details')
    parser_sumReport.add_argument('--input-file',
                                 default=None,
                                 type=str,
                                 help='Input file (in HDF5 format) containing MC4C data.')
    parser_sumReport.add_argument('--output-file',
                                 default=None,
                                 type=str,
                                 help='Output file (in PDF format) containing the requested summary report.')
    parser_sumReport.set_defaults(func=getSumRep)

    if flag_DEBUG:
        # sys.argv = ['./mc4c.py', 'init', './cfg_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'setReadIds', './cnf_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'splitReads', './cnf_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'mapFragments', './cnf_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'makeDataset', './cnf_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'removeDuplicates', './cnf_files/cfg_LVR-BMaj.cnf']
        sys.argv = ['./mc4c.py', 'getSumRep', 'readSizeDist', 'LVR-BMaj']
    args = parser.parse_args(sys.argv[1:])
    args.func(args)


if __name__ == '__main__':
    main()
