#! /usr/bin/env python

import argparse
import sys
import gzip
from os import path, makedirs, environ
import numpy as np
import pandas as pd

import mc4c_tools

import platform  # ###
# flag_DEBUG = platform.system() != 'Linux'

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


def processMC4C(args):
    from copy import deepcopy
    print '%% Processing MC-4C library ...'
    setattr(args, 'input_file', None)
    setattr(args, 'output_file', None)
    setattr(args, 'return_command', False)

    setReadIds(deepcopy(args))
    splitReads(deepcopy(args))
    mapFragments(deepcopy(args))
    processMappedFragments(deepcopy(args))
    removeDuplicates(deepcopy(args))
    print '[i] Processing MC-4C experiment is completed successfully.'

def setReadIds(args):
    print '%% Assigning traceable identifiers to reads ...'

    configs = mc4c_tools.load_configs(args.config_file, max_n_configs=1)[0]

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
                    if rd_idx % 5000 == 0:
                        print('\t\tprocessed {:,d} reads.'.format(rd_idx))

                    rd_sid = 'Fl.Id:{:d};Rd.Id:{:d};Rd.Ln:{:d}'.format(inp_fidx + 1, rd_idx, len(rd_seq))
                    out_fid.write('>' + rd_sid + '\n')
                    out_fid.write(rd_seq + '\n')
    print '[i] Read identifier conversion is completed successfully.'


def splitReads(args):
    from utilities import get_re_info
    import re

    print '%% Splitting reads into fragments ...'
    configs = mc4c_tools.load_configs(args.config_file, max_n_configs=1)[0]

    if args.input_file is None:
        args.input_file = './reads/rd_' + configs['run_id'] + '.fasta.gz'
    if args.output_file is None:
        args.output_file = './fragments/frg_' + configs['run_id'] + '.fasta.gz'
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
            if rd_ind % 50000 == 0:
                print('\tprocessed {:,d} reads and produced {:,d} fragments.'.format(rd_ind, frg_ind))

            frg_be = 0
            for res_enz in re.finditer(regex_ptr, rd_seq):
                frg_en = res_enz.end()
                if frg_en - frg_be > MAX_FRG_SIZE:
                    n_reduced += 1
                    frg_en = frg_be + MAX_FRG_SIZE
                out_fid.write('{:s};Fr.Id:{:d};Fr.SBp:{:d};Fr.EBp:{:d}\n'.format(rd_sid, frg_ind, frg_be + 1, frg_en))
                out_fid.write(rd_seq[frg_be:frg_en] + '\n')
                frg_ind = frg_ind + 1
                frg_be = res_enz.start()
            rd_ind = rd_ind + 1
    if n_reduced != 0:
        print '[i] [{:,d}] fragments are reduced to {:,d}bp.'.format(n_reduced, MAX_FRG_SIZE)
    print '[i] Total of {:,d} reads and {:,d} fragments are produced successfully.'.format(rd_ind - 1, frg_ind - 1)


def mapFragments(args):
    if not args.return_command:
        print '%% Mapping fragments to genome ...'

    configs = mc4c_tools.load_configs(args.config_file, max_n_configs=1)[0]

    # Map split fragments to genome
    if args.input_file is None:
        args.input_file = './fragments/frg_' + configs['run_id'] + '.fasta.gz'
    if args.output_file is None:
        args.output_file = './bams/bam_{:s}.bam'.format(configs['run_id'])
    if not path.isdir(path.dirname(args.output_file)):
        makedirs(path.dirname(args.output_file))
    # assert not path.isfile(args.output_file)
    if not args.return_command:
        print('Reading fragments from: {:s}'.format(args.input_file))
        print('Writing mapped fragments to: {:s}'.format(args.output_file))

    cmd_str = \
        configs['bwa'] + ' bwasw -b 5 -q 2 -r 1 -z 5 -T 15 -t {:d} '.format(args.n_thread) + \
        configs['bwa_index'] + ' ' + args.input_file + \
        ' | samtools view -q 1 -hbS - ' + \
        '> ' + args.output_file

    if args.return_command:
        print '{:s}'.format(cmd_str)
    else:
        print 'Running: {:s}'.format(cmd_str)
        import subprocess
        map_prs = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        std_out, std_err = map_prs.communicate()
        # TODO: A better error handling here would be nice
        assert std_err.split('\n')[-2][:18] == '[main] Real time: ', \
            'bwa failed to run properly, see below:\n{:s}'.format(std_err)
        print '[i] Fragments are mapped to genome successfully.'


def processMappedFragments(args):
    import h5py
    from os import remove, path
    from pandas import read_csv
    import pysam

    from utilities import get_chr_info, hasOL

    print '%% Creating a MC-4C dataset from mapped fragments ...'
    configs = mc4c_tools.load_configs(args.config_file, max_n_configs=1)[0]

    if args.input_file is None:
        args.input_file = './bams/bam_{:s}.bam'.format(configs['run_id'])
    if args.output_file is None:
        args.output_file = './datasets/mc4c_' + configs['run_id'] + '_all.hdf5'
    if not path.isdir(path.dirname(args.output_file)):
        makedirs(path.dirname(args.output_file))
    # assert not path.isfile(args.output_file), '[e] Output file already exists: {:s}'.format(args.output_file)
    print('Reading fragments from: {:s}'.format(args.input_file))

    # get chromosome information
    chr_lst = get_chr_info(configs['genome_build'], 'chr_name')
    chr_size = get_chr_info(configs['genome_build'], 'chr_size')
    chr_map = dict(zip(chr_lst, np.arange(len(chr_lst)) + 1))

    # loading corresponding restriction fragment positions
    re_pos_fname = './renzs/{:s}_{:s}.npz'.format(configs['genome_build'], '-'.join(configs['re_name']))
    if not path.isfile(re_pos_fname):
        from utilities import extract_re_positions
        print 'Database of restriction enzyme cut sites is not found. ' + \
              'Scanning the reference genome to create this database ...'
        extract_re_positions(genome_str=configs['genome_build'], re_name_lst=configs['re_name'],
                             ref_fasta=configs['reference_fasta'])
    re_pos, re_chr_lst, re_genome_str = np.load(re_pos_fname)['arr_0']
    assert np.array_equal(re_chr_lst, chr_lst)
    assert configs['genome_build'] == re_genome_str
    for ri in range(len(re_pos)):  # extending positions to length of chromosomes
        re_pos[ri] = np.hstack([0, re_pos[ri], chr_size[ri]])

    # extend fragment coordinates to reference genome
    n_processed = 0
    header_lst = ['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'Strand', 'MapStart', 'MapEnd', 'MQ',
                  'FileID', 'FrgID', 'SeqStart', 'SeqEnd', 'ReadLength', 'IsValid']
    n_header = len(header_lst)
    tmp_fname = args.output_file + '.tmp'
    print('Writing processed fragments to a temporary file first: {:s}'.format(tmp_fname))
    with pysam.AlignmentFile(args.input_file, 'rb') as bam_fid, gzip.open(tmp_fname, 'wb') as gz_fid:
        gz_fid.write('\t'.join(header_lst) + '\n')
        frg_template = '\t'.join(['{:d}'] * n_header) + '\n'

        frg_set = np.empty([0, n_header], dtype=np.int64)
        for que_idx, que_line in enumerate(bam_fid):
            if que_idx % 100000 == 0:
                print('\tprocessed {:,d} fragments in {:,d} reads.'.format(que_idx, n_processed))
            if (np.bitwise_and(que_line.flag, 0x800) == 0x800) or (que_line.reference_name not in chr_lst):
                continue
            FileID, ReadID, ReadLength, FrgID, SeqStart, SeqEnd = \
                [int(x.split(':')[1]) for x in que_line.query_name.split(';')]
            MapChrNum = chr_map[que_line.reference_name]
            MapStart = que_line.reference_start
            MapEnd = que_line.reference_end
            MapStrand = 1 - (que_line.is_reverse * 2)
            if que_idx == 0:
                FrgID_old = FrgID
                ReadID_old = ReadID
            # TODO: Unmapped fragments are ignored here

            # extending coordinates to nearby restriction site
            n_re = len(re_pos[MapChrNum - 1])
            nei_left = np.searchsorted(re_pos[MapChrNum - 1], MapStart, side='left') - 1
            if (nei_left < n_re - 1) and (np.abs(re_pos[MapChrNum - 1][nei_left + 1] - MapStart) < 10):
                nei_left = nei_left + 1
            nei_right = np.searchsorted(re_pos[MapChrNum - 1], MapEnd, side='left')
            if (nei_right > 0) and (np.abs(MapEnd - re_pos[MapChrNum - 1][nei_right - 1]) < 10):
                nei_right = nei_right - 1

            if nei_left == nei_right:
                dist_left = np.abs(MapStart - re_pos[MapChrNum - 1][nei_left])
                dist_right = np.abs(MapEnd - re_pos[MapChrNum - 1][nei_right])
                if dist_right < dist_left:
                    nei_left = nei_left - 1
                else:
                    nei_right = nei_right + 1

            ExtStart = re_pos[MapChrNum - 1][nei_left]
            try:
                ExtEnd = re_pos[MapChrNum - 1][nei_right]
            except:
                if nei_right == len(re_pos[MapChrNum - 1]):
                    ExtEnd = re_pos[MapChrNum - 1][-1] + 100
                else:
                    raise Exception('Error in: {:s}'.format(que_line.query_name))
            # TODO: Why the mapped coordinates are after end of chromosome
            # happened in: frg_BRN-BMaj-96x.fasta.gz
            # frg_id: Fl.Id:1;Rd.Id:1819914;Rd.Ln:2908;Fr.Id:12027752;Fr.SBp:1291;Fr.EBp:1700
            # frg_seq: GATCATATGGGCAGAAACATCAACATAATATGATTAAAATCAATAAATCATAAATACTCCACAAGTAAAATTTTACTTGTAAAATATGATTAAAGCATGTCTACCAAAATAATCCATCTAATCTTGATACTATCTTACACTACTCTATCATAAGATTTATTGAATATGGCATTTCAGAAAACATCACTAGTGTTCTGTGCATATCAGAACGCCAGTTTACACATATATCAACTATGGAAACAAATCAAGGGATACCAGCATATGAATTGATGCATAAATTGCTACGCTTTAACGTTTCTAAGGTATCAGGTGCAAGACACTTGTGTTACTCATCTGAAGCACCACTTACAATGCAACAAATAATACATTGTAGCCTCCCTCAGCCTTGAAGCAAGTGGTATCTGATGATC

            # adjust fragment seq according to skipped bases
            cigar_tup = que_line.cigartuples
            if MapStrand == -1:
                cigar_tup.reverse()
            if (cigar_tup[0][0] == 4) or (cigar_tup[0][0] == 5):
                SeqStart = SeqStart + cigar_tup[0][1]
            if (cigar_tup[-1][0] == 4) or (cigar_tup[-1][0] == 5):
                SeqEnd = SeqEnd - cigar_tup[-1][1]

            # combine into an array
            frg_info = np.array([
                ReadID, MapChrNum, ExtStart, ExtEnd, MapStrand, MapStart, MapEnd, que_line.mapping_quality,
                FileID, FrgID, SeqStart, SeqEnd, ReadLength, 1]).reshape([1, -1])

            # Check order of fragments
            if FrgID < FrgID_old:
                raise Exception('Order of fragments are lost.')
            FrgID_old = FrgID

            # Check if read has ended
            if ReadID == ReadID_old:
                frg_set = np.vstack([frg_set, frg_info])
                continue

            # sort the read according to seqs
            frg_set = frg_set[np.argsort(frg_set[:, 10]), :]

            # merge adjacent fragments
            fi = 0
            while fi < frg_set.shape[0] - 1:
                if hasOL(frg_set[fi, 1:5], frg_set[fi + 1:fi + 2, 1:5], offset=20)[0]:
                    frg_set[fi, 2] = np.min(frg_set[fi:fi + 2, 2])
                    frg_set[fi, 3] = np.max(frg_set[fi:fi + 2, 3])
                    assert frg_set[fi, 4] == frg_set[fi+1, 4]

                    frg_set[fi, 5] = np.min(frg_set[fi:fi + 2, 5])
                    frg_set[fi, 6] = np.max(frg_set[fi:fi + 2, 6])

                    frg_set[fi, 7] = np.max(frg_set[fi:fi + 2, 7])

                    frg_set[fi, 10] = np.min(frg_set[fi:fi + 2, 10])
                    frg_set[fi, 11] = np.max(frg_set[fi:fi + 2, 11])
                    frg_set = np.delete(frg_set, fi + 1, axis=0)
                else:
                    fi = fi + 1

            # make sure fragments are not overlapping (ignoring the strand)
            n_frg = frg_set.shape[0]
            for fi in range(n_frg):
                if frg_set[fi, 13] == 0:
                    continue
                for fj in range(fi+1, n_frg):
                    if frg_set[fj, 13] == 0:
                        continue

                    # TODO: This should also check to see if primers are found in the middle of reads (read-fusion)
                    if hasOL(frg_set[fi, 1:4], frg_set[fj:fj + 1, 1:4], offset=-20)[0]:
                        if frg_set[fi, 7] >= frg_set[fj, 7]:
                            frg_set[fj, 13] = 0
                        else:
                            frg_set[fi, 13] = 0

            # save the read
            for frg in frg_set:
                gz_fid.write(frg_template.format(*frg))
            frg_set = frg_info.copy()
            ReadID_old = ReadID

            n_processed += 1

        if frg_set.shape[0] != 0:  # Saving the last read after file has finished
            for frg in frg_set:
                gz_fid.write(frg_template.format(*frg))

    # Load fragments in pandas format and sort fragments within a read according to their relative positions
    print('Loading the temporary file: {:s}'.format(tmp_fname))
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

    print '%% Removing pcr duplicates from a MC-4C dataset ...'
    configs = mc4c_tools.load_configs(args.config_file, max_n_configs=1)[0]

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

    # filtering reads according to their MQ
    read_all = mc4c_pd[['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'MQ', 'IsValid']].values
    is_mapped = read_all[:, 4] >= args.min_mq
    is_valid = read_all[:, 5] == 1
    read_all = read_all[is_mapped & is_valid, :4]
    print 'Selected non-overlapping fragments with MQ >= {:d}: {:d} reads are left.'.format(
        args.min_mq, len(np.unique(read_all[:, 0])))
    del is_mapped, is_valid

    # select informative reads (#frg > 1), ignoring VP fragments
    vp_crd = np.array([configs['vp_cnum'], configs['vp_start'], configs['vp_end']])
    roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
    is_vp = hasOL(vp_crd, read_all[:, 1:4], offset=0)
    is_roi = hasOL(roi_crd, read_all[:, 1:4], offset=0)
    frg_roi = read_all[~is_vp & is_roi, :]
    read_n_roi = np.bincount(frg_roi[:, 0], minlength=MAX_ReadID + 1)
    is_inf = np.isin(read_all[:, 0], frg_roi[read_n_roi[frg_roi[:, 0]] > 1, 0])
    read_inf = np.hstack([read_all[is_inf, :], read_n_roi[read_all[is_inf, 0]].reshape(-1, 1)])
    print 'Selected reads #cis fragment > 1: {:d} reads are selected.'.format(len(np.unique(read_inf[:, 0])))

    # select reads with #traceable fragment > 1
    roi_size = configs['roi_end'] - configs['roi_start']
    lcl_crd = np.array([configs['vp_cnum'], configs['roi_start'] - roi_size, configs['roi_end'] + roi_size])
    is_lcl = hasOL(lcl_crd, read_inf[:, 1:4], offset=0)
    frg_trs = read_inf[~is_lcl, :]
    print 'Selected reads with #trans fragment > 0: {:d} reads are selected.'.format(len(np.unique(frg_trs[:, 0])))

    # make duplicate list of fragments
    frg_uid, frg_idx, frg_cnt = np.unique(frg_trs[:, 1:4], axis=0, return_index=True, return_counts=True)
    frg_idx = np.argsort(frg_idx)
    dup_info = np.hstack([frg_uid[frg_idx, :], frg_cnt[frg_idx].reshape(-1, 1)])

    # sort trans-fragments according to #duplicates
    dup_info = dup_info[np.lexsort([-dup_info[:, -1]]), :]

    # loop over trans fragments
    dup_idx = 0
    print 'Scanning for duplicated trans-fragments:'
    while dup_idx < dup_info.shape[0]:
        if dup_idx % 1000 == 0:
            print '\tscanned {:,d} trans-fragments, '.format(dup_idx) + \
                  '{:,d} reads are still unique.'.format(len(np.unique(frg_trs[:, 0])))
        has_ol = hasOL(dup_info[dup_idx, :3], frg_trs[:, 1:4], offset=0)
        if np.sum(has_ol) > 1:
            # select duplicates
            dup_set = frg_trs[has_ol, :]

            # keep largest read according to #roi fragments
            keep_rid = dup_set[np.argmax(dup_set[:, -1]), 0]
            dup_set = dup_set[dup_set[:, 0] != keep_rid, :]

            # remove extra duplicates
            frg_trs = frg_trs[~ np.isin(frg_trs[:, 0], dup_set[:, 0]), :]
        dup_idx = dup_idx + 1
    print 'Result statistics (before --> after filtering):'
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

    # load config files
    configs = mc4c_tools.load_configs(args.config_file, max_n_configs=1)[0]

    if args.output_file is None:
        configs['output_dir'] = './plots/'
    else:
        configs['output_dir'] = path.dirname(args.output_file)
    if not path.isdir(configs['output_dir']):
        makedirs(configs['output_dir'])
    configs['input_file'] = args.input_file
    configs['output_file'] = args.output_file

    # call the requested function
    if args.report_type == 'readSizeDist':
        mc4c_tools.plot_readSizeDistribution(configs)
    elif args.report_type =='cvgDist':
        mc4c_tools.plot_cvgDistribution(configs)
    elif args.report_type == 'cirSizeDist':
        mc4c_tools.plot_cirSizeDistribution(configs, roi_only=args.roi_only)
    elif args.report_type == 'overallProfile':
        mc4c_tools.plot_overallProfile(configs, MIN_N_FRG=2)
    else:
        raise Exception()


def perform_analysis(args):
    import mc4c_analysis

    configs = mc4c_tools.load_configs(args.config_file, max_n_configs=1)[0]
    if args.output_file is None:
        configs['output_dir'] = './plots/'
    else:
        configs['output_dir'] = path.dirname(args.output_file)
    if not path.isdir(configs['output_dir']):
        makedirs(configs['output_dir'])
    configs['input_file'] = args.input_file
    configs['output_file'] = args.output_file

    # call the requested function
    if args.analysis_type == 'mcTest':
        mc4c_analysis.perform_mc_analysis(configs)
    if args.analysis_type == 'vpSoi':
        mc4c_analysis.perform_vpsoi_analysis(configs, soi_name=args.ant_name, n_perm=args.n_perm)
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

    # Runs all pre-processing steps
    parser_process = subparsers.add_parser('process', description='Performs the entire pre-processing ' +
                                                                  'steps required for a MC-4C experiment.')
    parser_process.add_argument('config_file', metavar='config-file', type=str,
                               help='Configuration file containing experiment specific details')
    parser_process.add_argument('--n_thread', default=6,type=int,
                               help='Number of threads should be used by the aligner')
    parser_process.add_argument('--min-mq', default=20, type=int,
                               help='Minimum mapping quality (MQ) to consider a fragment as confidently mapped.')
    parser_process.set_defaults(func=processMC4C)

    # Set read identifiers
    parser_readid = subparsers.add_parser('setReadIds', description='Defining identifiers for sequenced reads')
    parser_readid.add_argument('config_file', metavar='config-file', type=str,
                               help='Configuration file containing experiment specific details')
    parser_readid.add_argument('--input-file', default=None, type=str,
                               help='Input file (in FASTQ format) containing raw sequenced reads')
    parser_readid.add_argument('--output-file', default=None, type=str,
                               help='Output file (in FASTA format) containing sequenced reads with traceable IDs')
    parser_readid.set_defaults(func=setReadIds)

    # Split reads into fragments
    parser_readSplt = subparsers.add_parser('splitReads', description='Splitting reads into fragments using ' +
                                                                      'restriction enzyme recognition sequence')
    parser_readSplt.add_argument('config_file', metavar='config-file', type=str,
                                 help='Configuration file containing experiment specific details')
    parser_readSplt.add_argument('--input-file', default=None, type=str,
                                 help='Input file (in FASTA format) containing reads with traceable IDs.')
    parser_readSplt.add_argument('--output-file', default=None, type=str,
                                 help='Output file (in FASTA format) containing fragments with traceable IDs')
    parser_readSplt.set_defaults(func=splitReads)

    # Mapping fragments to reference genome
    parser_mapFrg = subparsers.add_parser('mapFragments', description='Mapping fragments to reference genome')
    parser_mapFrg.add_argument('config_file', metavar='config-file', type=str,
                               help='Configuration file containing experiment specific details')
    parser_mapFrg.add_argument('--input-file', default=None, type=str,
                               help='Input file (in FASTA format) containing fragments with traceable IDs')
    parser_mapFrg.add_argument('--output-file', default=None, type=str,
                               help='Output file (in BAM format) containing fragments with traceable IDs')
    parser_mapFrg.add_argument('--n_thread', default=6, type=int,
                               help='Number of threads should be used by the aligner')
    parser_mapFrg.add_argument('--return_command', action="store_true",
                               help='Return only mapping command instead of running it ' +
                                    '(useful for running the pipeline in a cluster)')
    parser_mapFrg.set_defaults(func=mapFragments)

    # Process mapped fragments
    parser_mkDataset = subparsers.add_parser('makeDataset',
                                             description='Processed the mapped fragments and create a MC-4C dataset')
    parser_mkDataset.add_argument('config_file', metavar='config-file', type=str,
                                  help='Configuration file containing experiment specific details')
    parser_mkDataset.add_argument('--input-file', default=None, type=str,
                                  help='Input file (in BAM format) containing fragments with traceable IDs')
    parser_mkDataset.add_argument('--output-file', default=None, type=str,
                                  help='Output file (in HDF5 format) containing processed fragments')
    parser_mkDataset.set_defaults(func=processMappedFragments)

    # Remove PCR duplicated
    parser_remDup = subparsers.add_parser('removeDuplicates',
                                          description='Remove PCR duplicates from a given MC-4C dataset')
    parser_remDup.add_argument('config_file', metavar='config-file', type=str,
                               help='Configuration file containing experiment specific details')
    parser_remDup.add_argument('--input-file', default=None, type=str,
                               help='Input file (in HDF5 format) containing MC4C data.')
    parser_remDup.add_argument('--output-file', default=None, type=str,
                               help='Output file (in HDF5 format) containing MC4C data.')
    parser_remDup.add_argument('--min-mq', default=20, type=int,
                               help='Minimum mapping quality (MQ) to consider a fragment as confidently mapped.')
    parser_remDup.set_defaults(func=removeDuplicates)

    # produce statistics plots
    parser_sumReport = subparsers.add_parser('getSumRep',
                                            description='Generate various summary reports about a MC-4C dataset.')
    parser_sumReport.add_argument('report_type', type=str,
                                 choices=['cvgDist', 'readSizeDist', 'cirSizeDist', 'overallProfile'],
                                 help='Type of summary report that needs to be generated')
    parser_sumReport.add_argument('config_file', metavar='config-file', type=str,
                                 help='Configuration file containing experiment specific details')
    parser_sumReport.add_argument('--input-file', default=None, type=str,
                                 help='Input file (in HDF5 format) containing MC4C data.')
    parser_sumReport.add_argument('--output-file', default=None, type=str,
                                 help='Output file (in PDF format) containing the requested summary report.')
    parser_sumReport.add_argument('--roi-only', action="store_true",
                                  help='Limits the requested summary report to be generated from roi-fragments only.')
    parser_sumReport.set_defaults(func=getSumRep)

    # perform basic analysis
    parser_analysis = subparsers.add_parser('analysis',
                                             description='Perform analysis on a MC-4C dataset.')
    parser_analysis.add_argument('analysis_type', choices=['mcTest', 'vpSoi'], type=str,
                                  help='Type of analysis that needs to be performed')
    parser_analysis.add_argument('config_file', metavar='config-file', type=str,
                                  help='Configuration file containing experiment specific details')
    parser_analysis.add_argument('ant_name', metavar='ant-name', type=str,
                                 help='Name of annotation for which VP-SOI plot needs to be computed. ' +
                                      'Only used for VP-SOI (i.e. "vpSoi") analysis')
    parser_analysis.add_argument('--input-file', default=None, type=str,
                                  help='Input file (in HDF5 format) containing MC4C data.')
    parser_analysis.add_argument('--output-file', default=None, type=str,
                                  help='Output file (in PDF format) containing the result of the requested analysis.')
    parser_analysis.add_argument('--roi-only', action="store_true",
                                  help='Limits the requested analysis to be generated from roi-fragments only. ' +
                                       'By default this flag is set to TRUE.')
    parser_analysis.add_argument('--n-perm', default=1000, type=int,
                                 help='Number of profiles that needs to be drawn from negative reads (i.e. reads ' +
                                      'that contain no fragment from site of interest) to produce the expected profile.')
    parser_analysis.set_defaults(func=perform_analysis)

    if 'FROM_PYCHARM' in environ:
        Warning('We are in the PyCharm!')
        # sys.argv = ['./mc4c.py', 'process', 'LVR-BMaj-PB']
        # sys.argv = ['./mc4c.py', 'init', './cfg_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'setReadIds', './cnf_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'splitReads', 'LVR-BMaj']
        # sys.argv = ['./mc4c.py', 'mapFragments', 'LVR-BMaj']
        # sys.argv = ['./mc4c.py', 'makeDataset', 'LVR-BMaj-96x']
        # sys.argv = ['./mc4c.py', 'removeDuplicates', 'LVR-BMaj']
        # sys.argv = ['./mc4c.py', 'getSumRep', 'readSizeDist', 'K562-GATA1']
        # sys.argv = ['./mc4c.py', 'getSumRep', 'cvgDist', 'K562-GATA1']
        # sys.argv = ['./mc4c.py', 'getSumRep', 'cirSizeDist', 'K562-WplD-10x', '--roi-only']
        # sys.argv = ['./mc4c.py', 'getSumRep', 'overallProfile', 'K562-WplD-10x']
        # sys.argv = ['./mc4c.py', 'analysis', 'mcTest', 'K562-WplD-10x']
        sys.argv = ['./mc4c.py', 'analysis', 'vpSoi', '--n-perm=10', 'LVR-BMaj-96x', 'HS2']

    args = parser.parse_args(sys.argv[1:])
    args.func(args)


if __name__ == '__main__':
    main()

# cluster run example:
# qsub -P hub_laat -N mc4c -l h_rt=05:00:00 -l h_vmem=50G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 ./mc4c.py setReadIds LVR-BMaj"
