import numpy as np
from os import path, makedirs


def remove_duplicates_by_umi(umi_set, verbose=False):
    from utilities import hasOL

    # make duplicate list of fragments
    frg_uid, frg_idx, frg_cnt = np.unique(umi_set[:, 1:4], axis=0, return_index=True, return_counts=True)
    frg_idx = np.argsort(frg_idx)
    frg_umi = np.hstack([frg_uid[frg_idx, :], frg_cnt[frg_idx].reshape(-1, 1)])

    # sort trans-fragments according to #duplicates
    frg_umi = frg_umi[np.lexsort([-frg_umi[:, -1]]), :]
    n_umi = frg_umi.shape[0]
    # np.savetxt('/Users/aallahyar/Downloads/python.txt', frg_umi, delimiter='\t', fmt='%0.0f')

    # loop over trans fragments
    umi_idx = 0
    duplicate_info = []
    if verbose:
        print 'Scanning {:,d} UMIs for duplicates:'.format(n_umi)
    while umi_idx < n_umi:
        if verbose and (umi_idx % 1000 == 0):
            print '\tscanned {:,d} trans-fragments, '.format(umi_idx) + \
                  '{:,d} reads are still unique.'.format(len(np.unique(umi_set[:, 0])))
        has_ol = hasOL(frg_umi[umi_idx, :3], umi_set[:, 1:4], offset=0)
        n_ovl = len(np.unique(umi_set[has_ol, 0]))
        if n_ovl > 1:
            # select duplicates
            dup_set = umi_set[has_ol, :].copy()

            # keep largest read according to #roi fragments
            keep_rid = dup_set[np.argmax(dup_set[:, -1]), 0]
            dup_set = dup_set[dup_set[:, 0] != keep_rid, :]

            # remove extra duplicates
            umi_set = umi_set[~ np.isin(umi_set[:, 0], dup_set[:, 0]), :]
        elif n_ovl == 1:
            keep_rid = umi_set[has_ol, 0][0]
            dup_set = umi_set[has_ol, :].copy()
        else:
            keep_rid = -1
            dup_set = np.empty([0, 5])

        # save information
        duplicate_info.append((keep_rid, dup_set[:, 0].copy(), frg_umi[umi_idx, :].copy()))
        umi_idx = umi_idx + 1

    return umi_set, duplicate_info


##########################################################
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
    import gzip

    from utilities import load_configs
    print '%% Assigning traceable identifiers to reads ...'

    configs = load_configs(args.config_file, max_n_configs=1)[0]

    # initialize
    if args.input_file is None:
        args.input_file = './fastqs/fq_' + configs['run_id'] + '.fastq.gz'
    if args.output_file is None:
        args.output_file = './reads/rd_' + configs['run_id'] + '.fasta.gz'
    # assert not path.isfile(args.output_file), '[e] output file already exists: {:s}'.format(args.output_file)
    if not path.isdir(path.dirname(args.output_file)):
        makedirs(path.dirname(args.output_file))
    print('Writing reads with traceable identifiers to: {:s}'.format(args.output_file))

    # loop over reads
    with gzip.open(args.output_file, 'w') as out_fid:
        inp_flst = args.input_file.split(',')
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
    # TODO: check to see if #lines in inp_file equals #lines read


def splitReads(args):
    import re
    import gzip

    from utilities import load_configs, get_re_info

    print '%% Splitting reads into fragments ...'
    configs = load_configs(args.config_file, max_n_configs=1)[0]

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

    from utilities import load_configs
    if not args.return_command:
        print '%% Mapping fragments to genome ...'

    configs = load_configs(args.config_file, max_n_configs=1)[0]

    # Map split fragments to genome
    if args.input_file is None:
        args.input_file = './fragments/frg_' + configs['run_id'] + '.fasta.gz'
    if args.output_file is None:
        args.output_file = './bams/bam_{:s}.bam'.format(configs['run_id'])
    if not path.isdir(path.dirname(args.output_file)):
        makedirs(path.dirname(args.output_file))
    if not args.return_command:
        print('Reading fragments from: {:s}'.format(args.input_file))
        print('Writing mapped fragments to: {:s}'.format(args.output_file))

    # prepare the command
    cmd_str = \
        configs['bwa'] + ' bwasw -b 5 -q 2 -r 1 -z 5 -T 15 -t {:d} '.format(args.n_thread) + \
        configs['bwa_index'] + ' ' + args.input_file + \
        ' | samtools view -q 1 -hbS - ' + \
        '> ' + args.output_file
    if args.return_command:
        print '{:s}'.format(cmd_str)
    else:
        assert path.isfile(args.input_file), '[e] Input file could not be found: {:s}'.format(args.input_file)
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
    from os import remove
    from pandas import read_csv
    import pysam
    import gzip

    from utilities import load_configs, get_chr_info, hasOL

    print '%% Creating an MC-4C dataset from mapped fragments ...'
    configs = load_configs(args.config_file, max_n_configs=1)[0]

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

    # identify vp_fragment coordinates
    idx_lft = np.searchsorted(re_pos[configs['vp_cnum'] - 1], np.min(configs['prm_start']), side='right') - 1
    idx_rgt = np.searchsorted(re_pos[configs['vp_cnum'] - 1], np.max(configs['prm_end']) - len(configs['re_seq'][0]), side='left')
    vp_frg = [configs['vp_cnum'], re_pos[configs['vp_cnum'] - 1][idx_lft], re_pos[configs['vp_cnum'] - 1][idx_rgt]]
    if idx_lft + 1 == idx_rgt:
        print '[w] Can not map primer positions on a single fragment.'

    # define fragment headers
    header_lst = ['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'Strand', 'MapStart', 'MapEnd', 'MQ',
                  'FileID', 'FrgID', 'SeqStart', 'SeqEnd', 'ReadLength', 'Flag']
    n_header = len(header_lst)

    # Processing fragments:
    # 1. extending fragment coordinates to nearest restriction site in the reference genome
    # 2. checking overlap within reads
    # 3. checking for read-fusion events
    tmp_fname = args.output_file + '.tmp'
    print('Writing processed fragments to a temporary file first: {:s}'.format(tmp_fname))
    n_fusion = 0
    n_processed = 0
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
            # TODO: Note that unmapped fragments are ignored here

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
                ExtEnd = re_pos[MapChrNum - 1][nei_right] - 1  # Adjacent fragments in ref should not overlap with 1bp
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
                FileID, FrgID, SeqStart, SeqEnd, ReadLength, 0]).reshape([1, -1])

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
                if np.bitwise_and(frg_set[fi, 13], 1) == 1:
                    continue
                for fj in range(fi+1, n_frg):
                    if np.bitwise_and(frg_set[fj, 13], 1) == 1:
                        continue

                    if hasOL(frg_set[fi, 1:4], frg_set[fj:fj + 1, 1:4], offset=-20)[0]:
                        if frg_set[fi, 7] >= frg_set[fj, 7]:
                            frg_set[fj, 13] = np.bitwise_or(frg_set[fj, 13], 1)
                        else:
                            frg_set[fi, 13] = np.bitwise_or(frg_set[fi, 13], 1)

            # check for read-fusion events
            for fi in range(1, frg_set.shape[0] - 1):
                is_prv_vp = hasOL(vp_frg, frg_set[fi - 1, [1, 5, 6]], offset=-10)[0]
                is_cur_vp = hasOL(vp_frg, frg_set[fi    , [1, 5, 6]], offset=-10)[0]
                is_nxt_vp = hasOL(vp_frg, frg_set[fi + 1, [1, 5, 6]], offset=-10)[0]
                if ~is_prv_vp & is_cur_vp & ~is_nxt_vp:
                    frg_set[:, 13] = np.bitwise_or(frg_set[:, 13], 2)
                    n_fusion += 1
                    break

            # save the read
            for frg in frg_set:
                gz_fid.write(frg_template.format(*frg))
            frg_set = frg_info.copy()
            ReadID_old = ReadID

            n_processed += 1

        if frg_set.shape[0] != 0:  # Saving the last read after file has finished
            for frg in frg_set:
                gz_fid.write(frg_template.format(*frg))
            n_processed += 1

    if n_fusion != 0:
        print '[w] {:,d} fused reads ({:0.1f}% of total) are identified and flagged.'.format(
            n_fusion, n_fusion * 100.0 / n_processed)

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

#
# def find_optimal_roi(args):
#     import pandas as pd
#     from matplotlib import pyplot as plt
#     from scipy.stats import spearmanr
#
#     from utilities import load_mc4c, limit_to_roi, get_chr_info, hasOL, get_nreads_per_bin, load_annotation, load_configs
#
#     print '%% Finding optimal ROI by analysing local coverage ...'
#     config_lst = load_configs(args.config_file)
#
#     # initialization
#     run_id = ','.join([config['run_id'] for config in config_lst])
#     configs = config_lst[0]
#     if args.output_file is None:
#         args.output_file = './plots/plt_selectROI-ByCoverage_' + run_id + '.pdf'
#     if not path.isdir(path.dirname(args.output_file)):
#         makedirs(path.dirname(args.output_file))
#
#     def_vp_crd = np.array([configs['vp_cnum'], configs['vp_start'], configs['vp_end']])
#     def_roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
#     def_roi_cen = int(np.mean([configs['vp_start'], configs['vp_end']]))
#     def_roi_w = configs['roi_end'] - configs['roi_start']
#     blk_w = 30000
#     n_bin = 200
#
#     # split vp chromosome into blocks
#     chr_size = get_chr_info(genome_str=configs['genome_build'], property='chr_size')
#     blk_lst = np.arange(0, chr_size[configs['vp_cnum'] - 1], blk_w, dtype=np.int64).reshape(-1, 1)
#     blk_lst[-1] = chr_size[configs['vp_cnum'] - 1]
#     n_blk = len(blk_lst) - 1
#     blk_crd = np.hstack([np.repeat(configs['vp_cnum'], n_blk).reshape(-1, 1), blk_lst[:-1], blk_lst[1:] - 1])
#     blk_cen = np.mean(blk_crd[:, 1:3], axis=1)
#
#     # use raw reads and extract UMIs
#     print 'Extracting UMIs from given dataset(s) ...'
#     read_all = np.empty([0, 4], dtype=np.int)
#     def_pcr = np.empty([0, 4], dtype=np.int)
#     trs_pcr = np.empty([0, 4], dtype=np.int)
#     max_n_read = 10000000
#     for idx, cfg in enumerate(config_lst):
#
#         # load raw data
#         print('Loading all reads from {:s} ...'.format(cfg['run_id']))
#         mc4c_part = load_mc4c(cfg, unique_only=False, valid_only=True, min_mq=args.min_mq, reindex_reads=True, verbose=False)
#         read_prt = mc4c_part[['ReadID', 'Chr', 'ExtStart', 'ExtEnd']].values
#
#         # use default approach, use >1 roi-frg
#         has_inf = limit_to_roi(read_prt[:, :4], vp_crd=def_vp_crd, roi_crd=def_roi_crd, min_n_frg=2)
#         def_read_inf = read_prt[np.isin(read_prt[:, 0], has_inf[:, 0]), :].copy()
#         def_umi_set = def_read_inf[def_read_inf[:, 1] != configs['vp_cnum'], :].copy()
#         has_uid = remove_duplicates_by_umi(def_umi_set)[0]
#         def_pcr_prt = read_prt[np.isin(read_prt[:, 0], has_uid[:, 0]), :].copy()
#         del has_inf, has_uid
#         print '\t{:,d} unique reads are added using [far-cis + trans] UMIs.'.format(len(np.unique(def_pcr_prt[:, 0])))
#
#         # select >1 cis-frg
#         trs_vp_crd = np.array([configs['vp_cnum'], def_roi_cen - 5000, def_roi_cen + 5000])
#         is_cis = read_prt[:, 1] == configs['vp_cnum']
#         is_vp = hasOL(trs_vp_crd, read_prt[:, 1:4])
#         trs_read_m1c = read_prt[~is_vp & is_cis, :].copy()
#         read_size = np.bincount(trs_read_m1c[:, 0], minlength=np.max(trs_read_m1c[:, 0]) + 1)[trs_read_m1c[:, 0]]
#         trs_read_m1c = read_prt[np.isin(read_prt[:, 0], trs_read_m1c[read_size >= 2, 0]), :].copy()
#         del read_size
#
#         # select and process UMIs
#         trs_umi_set = trs_read_m1c[trs_read_m1c[:, 1] != configs['vp_cnum'], :].copy()
#         trs_uid = remove_duplicates_by_umi(trs_umi_set)[0]
#         trs_pcr_prt = read_prt[np.isin(read_prt[:, 0], trs_uid[:, 0]), :].copy()
#         print '\t{:,d} unique reads are added using [trans only] UMIs.'.format(len(np.unique(trs_pcr_prt[:, 0])))
#
#         # add data specific identifiers
#         assert np.max(read_prt[:, 0]) < max_n_read
#         def_pcr_prt[:, 0] = def_pcr_prt[:, 0] + (idx + 1) * max_n_read
#         trs_pcr_prt[:, 0] = trs_pcr_prt[:, 0] + (idx + 1) * max_n_read
#         read_prt[:, 0] = read_prt[:, 0] + (idx + 1) * max_n_read
#
#         # appending
#         read_all = np.vstack([read_all, read_prt])
#         def_pcr = np.vstack([def_pcr, def_pcr_prt])
#         trs_pcr = np.vstack([trs_pcr, trs_pcr_prt])
#
#     # compute coverage over chromosome
#     def_cvg, n_def = get_nreads_per_bin(def_pcr[:, :4], bin_crd=blk_crd, min_n_frg=0)
#     trs_cvg, n_trs = get_nreads_per_bin(trs_pcr[:, :4], bin_crd=blk_crd, min_n_frg=0)
#     def_nrm = pd.Series(def_cvg * 1e2 / n_def).rolling(5, center=True).mean().values
#     trs_nrm = pd.Series(trs_cvg * 1e2 / n_trs).rolling(5, center=True).mean().values
#
#     # select highly covered region
#     np.seterr(all='ignore')
#     cvd_idx = np.where(trs_nrm > args.min_cvg)[0]
#     np.seterr(all=None)
#     adj_roi_crd = np.array([configs['vp_cnum'], blk_crd[cvd_idx[0], 1], blk_crd[cvd_idx[-1], 2]])
#     adj_roi_w = adj_roi_crd[2] - adj_roi_crd[1]
#     adj_edge_lst = np.linspace(adj_roi_crd[1], adj_roi_crd[2], num=n_bin + 1, dtype=np.int64).reshape(-1, 1)
#     adj_bin_w = adj_edge_lst[1, 0] - adj_edge_lst[0, 0]
#     adj_vp_crd = np.array([configs['vp_cnum'], def_roi_cen - int(adj_bin_w * 1.5), def_roi_cen + int(adj_bin_w * 1.5)])
#     print 'Sufficient (>{:0.1f}%) coverage is found between -> '.format(args.min_cvg) + \
#           '{:s}:{:,d}-{:,d}'.format(configs['vp_chr'], adj_roi_crd[1], adj_roi_crd[2])
#
#     # select informative reads
#     has_inf = limit_to_roi(read_all[:, :4], vp_crd=adj_vp_crd, roi_crd=adj_roi_crd, min_n_frg=2)
#     adj_read_inf = read_all[np.isin(read_all[:, 0], has_inf[:, 0]), :].copy()
#     adj_lcl_crd = [adj_roi_crd[0], adj_roi_crd[1] - adj_roi_w, adj_roi_crd[2] + adj_roi_w]
#     is_lcl = hasOL(adj_lcl_crd, adj_read_inf[:, 1:4], offset=0)
#     adj_umi = adj_read_inf[~is_lcl, :].copy()
#     adj_uid = remove_duplicates_by_umi(adj_umi)[0]
#     is_adj = np.isin(adj_read_inf[:, 0], adj_uid[:, 0])
#     adj_pcr = adj_read_inf[is_adj, :].copy()
#
#     # compute coverage over chromosome using adjusted region
#     adj_cvg, n_adj = get_nreads_per_bin(adj_pcr[:, :4], bin_crd=blk_crd, min_n_frg=2)
#     adj_nrm = pd.Series(adj_cvg * 1e2 / n_adj).rolling(5, center=True).mean().values
#
#     # compute roi profiles
#     def_edge_lst = np.linspace(configs['roi_start'], configs['roi_end'], num=n_bin + 1, dtype=np.int64).reshape(-1, 1)
#     def_bin_crd = np.hstack([np.repeat(configs['vp_cnum'], n_bin).reshape(-1, 1), def_edge_lst[:-1], def_edge_lst[1:] - 1])
#     def_bin_w = def_bin_crd[0, 2] - def_bin_crd[0, 1]
#     def_prf, n_def = get_nreads_per_bin(def_pcr[:, :4], bin_crd=def_bin_crd, min_n_frg=2)
#     trs_prf, n_trs = get_nreads_per_bin(trs_pcr[:, :4], bin_crd=def_bin_crd, min_n_frg=2)
#     adj_prf, n_adj = get_nreads_per_bin(adj_pcr[:, :4], bin_crd=def_bin_crd, min_n_frg=2)
#
#     # plotting
#     plt.figure(figsize=(25, 5))
#     ax_crr = plt.subplot2grid((1, 4), (0, 0), rowspan=1, colspan=1)
#     ax_prf = plt.subplot2grid((1, 4), (0, 1), rowspan=1, colspan=3)
#     plt_h = [None] * 4
#     x_lim = [configs['roi_start'] - def_roi_w * 2, configs['roi_end'] + def_roi_w * 2]
#     y_lim = [0, 10]
#
#     # plot correlations
#     plt_h[0] = ax_crr.plot(def_prf * 1e2 / n_def, trs_prf * 1e2 / n_trs, 'x', color='#ffad14', alpha=0.9)[0]
#     plt_h[1] = ax_crr.plot(def_prf * 1e2 / n_def, adj_prf * 1e2 / n_adj, 'o', color='#2e2eff', alpha=0.5, markeredgecolor='None')[0]
#     ax_crr.set_xlim(y_lim)
#     ax_crr.set_ylim(y_lim)
#     ax_crr.set_xlabel('Default ROI')
#     ax_crr.set_ylabel('Trans / Adjusted ROI')
#     ax_crr.set_title('ROI coverage Spearman correlations\n' +
#                      'def-UMI vs. trs-UMI: {:0.5f}\n'.format(spearmanr(def_prf, trs_prf).correlation) +
#                      'def-UMI vs. adj-UMI: {:0.5f}'.format(spearmanr(def_prf, adj_prf).correlation))
#     ax_crr.legend(plt_h[:2], ['Default vs Trans profile', 'Default vs. Adjusted profile'])
#
#     # plot roi profiles
#     ax_prf.plot([def_roi_crd[1], def_roi_crd[1]], y_lim, color='#bcbcbc')
#     ax_prf.plot([def_roi_crd[2], def_roi_crd[2]], y_lim, color='#bcbcbc')
#     ax_prf.plot([adj_roi_crd[1], adj_roi_crd[1]], y_lim, color='#a3d1ff')
#     ax_prf.plot([adj_roi_crd[2], adj_roi_crd[2]], y_lim, color='#a3d1ff')
#     ax_prf.text(def_roi_crd[1], y_lim[1] * 0.9, '{:,d}> '.format(def_roi_crd[1]), horizontalalignment='right', color='#9c9c9c')
#     ax_prf.text(def_roi_crd[2], y_lim[1] * 0.9, ' <{:,d}'.format(def_roi_crd[2]), horizontalalignment='left', color='#9c9c9c')
#     ax_prf.text(adj_roi_crd[1], y_lim[1] * 0.8, '{:,d}> '.format(adj_roi_crd[1]), horizontalalignment='right', color='#52a8ff')
#     ax_prf.text(adj_roi_crd[2], y_lim[1] * 0.8, ' <{:,d}'.format(adj_roi_crd[2]), horizontalalignment='left', color='#52a8ff')
#
#     plt_h[0] = ax_prf.plot(blk_cen, def_nrm, '--', color='#777777')[0]
#     plt_h[1] = ax_prf.plot(blk_cen, trs_nrm, ':o', color='#ffad14', alpha=0.8, markersize=4, markeredgecolor=None)[0]
#     plt_h[2] = ax_prf.plot(blk_cen, adj_nrm, '-',  color='#2e2eff')[0]
#
#     # add annotations
#     ant_pd = load_annotation(configs['genome_build'], roi_crd=[configs['vp_cnum']] + x_lim)
#     for ai in range(ant_pd.shape[0]):
#         ant_pos = ant_pd.loc[ai, 'ant_pos']
#         if (ai == 0) or (np.abs(ant_pd.loc[ai - 1, 'ant_pos'] - ant_pos) > def_roi_w / 50.0):
#             ax_prf.text(ant_pos, y_lim[0], ' ' + ant_pd.loc[ai, 'ant_name'],
#                         horizontalalignment='center', verticalalignment='bottom', rotation=90, fontsize=6)
#         ax_prf.plot([ant_pos, ant_pos], [y_lim[0] + y_lim[1] * 0.05, y_lim[1]],
#                     ':', color='#bfbfbf', linewidth=1, alpha=0.3)
#     plt_h[3] = ax_prf.plot([x_lim[0], x_lim[1]], [args.min_cvg, args.min_cvg], '--', color='#ff8f8f')[0]
#
#     ax_prf.set_xlim(x_lim)
#     ax_prf.set_ylim(y_lim)
#     x_ticks = np.linspace(x_lim[0], x_lim[1], 25, dtype=np.int64)
#     x_tick_label = ['{:0.3f}m'.format(x / 1e6) for x in x_ticks]
#     plt.xticks(x_ticks, x_tick_label, rotation=20)
#     plt.ylabel('Frequency (% of reads)')
#     ax_prf.legend(plt_h, [
#         'Far-cis + trans (n={:0.0f})'.format(n_def),
#         'Trans only (n={:0.0f})'.format(n_trs),
#         'Adjusted (n={:0.0f})'.format(n_adj),
#         'Coverage threshold ({:0.1f}%)'.format(args.min_cvg)
#     ], loc='center left')
#
#     plt.title('{:s}\n'.format(run_id) +
#               '#block={:d}, block_w={:0.0f}k\n'.format(n_blk, blk_w / 1e3) +
#               'bin-w (def, adjusted): {:0,.0f}bp; {:0,.0f}bp'.format(def_bin_w, adj_bin_w) + '\n' +
#               'roi-w (def, adjusted): {:0,.0f}bp; {:0,.0f}bp'.format(def_roi_w, adj_roi_w))
#     plt.savefig(args.output_file, bbox_inches='tight')
#     print 'Coverage profiles are plotted in: {:s}'.format(args.output_file)


def removeDuplicates(args):
    import h5py
    from os import path, makedirs
    import pandas as pd

    from utilities import load_configs, hasOL

    print '%% Removing pcr duplicates from an MC-4C dataset ...'
    configs = load_configs(args.config_file, max_n_configs=1)[0]

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
    data_np = h5_fid[target_field][()]
    header_lst = list(h5_fid[target_field + '_header_lst'][()])
    mc4c_pd = pd.DataFrame(data_np, columns=header_lst)
    chr_lst = list(h5_fid['chr_lst'][()])
    h5_fid.close()
    MAX_ReadID = np.max(mc4c_pd['ReadID'])
    print 'There are {:d} reads in the dataset.'.format(len(np.unique(mc4c_pd['ReadID'])))

    # filtering reads according to their MQ
    read_all = mc4c_pd[['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'MQ', 'Flag']].values
    is_mapped = read_all[:, 4] >= args.min_mq
    is_valid = np.bitwise_and(read_all[:, 5], 1) == 0
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
    umi_set = read_inf[~is_lcl, :]
    print 'Selected reads with #trans fragment > 0: {:d} reads are selected.'.format(len(np.unique(umi_set[:, 0])))

    # remove duplicates
    unq_set, duplicate_info = remove_duplicates_by_umi(umi_set, verbose=True)
    print 'Result statistics (before --> after filtering):'
    print '\t#reads: {:,d} --> {:,d}'.format(len(np.unique(mc4c_pd['ReadID'])), len(np.unique(unq_set[:, 0])))
    print '\t#fragments: {:,d} --> {:,d}'.format(mc4c_pd['ReadID'].shape[0], unq_set.shape[0])

    # select and save unique reads
    is_uniq = np.isin(mc4c_pd['ReadID'], unq_set[:, 0])
    uniq_pd = mc4c_pd.loc[is_uniq, :]
    print('Writing dataset to: {:s}'.format(args.output_file))
    h5_fid = h5py.File(args.output_file, 'w')
    h5_fid.create_dataset('frg_np', data=uniq_pd.values, compression='gzip', compression_opts=5)
    h5_fid.create_dataset('frg_np_header_lst', data=list(uniq_pd.columns.values))
    h5_fid.create_dataset('chr_lst', data=chr_lst)

    dup_dt = np.dtype([('kept_id', np.int64),
                       ('dup_ids', h5py.special_dtype(vlen=np.int64)),
                       ('dup_frg', h5py.special_dtype(vlen=np.int64))])
    h5_fid.create_dataset('duplicate_info', data=np.array(duplicate_info, dtype=dup_dt), dtype=dup_dt)

    h5_fid.close()
    print '[i] PCR duplicates are removed from MC4C dataset successfully.'

