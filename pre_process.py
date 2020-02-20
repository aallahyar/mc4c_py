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
        print('Scanning {:,d} UMIs for duplicates:'.format(n_umi))
    while umi_idx < n_umi:
        if verbose and (umi_idx % 1000 == 0):
            print('\tscanned {:,d} trans-fragments, '.format(umi_idx) + \
                  '{:,d} reads are still unique.'.format(len(np.unique(umi_set[:, 0]))))
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


def processMC4C(args):
    from copy import deepcopy

    print('%% Processing MC-4C library ...')
    setattr(args, 'input_file', None)
    setattr(args, 'output_file', None)
    setattr(args, 'return_command', False)

    setReadIds(deepcopy(args))
    splitReads(deepcopy(args))
    mapFragments(deepcopy(args))
    processMappedFragments(deepcopy(args))
    removeDuplicates(deepcopy(args))
    print('[i] Processing MC-4C experiment is completed successfully.')


def setReadIds(args):
    import gzip

    from utilities import load_configs
    print('%% Assigning traceable identifiers to reads ...')

    configs = load_configs(args.config_file, max_n_configs=1)[0]

    # initialize
    if args.input_file is None:
        args.input_file = './fastqs/' + configs['run_id'] + '.fastq.gz'
    if args.output_file is None:
        args.output_file = './reads/rd_' + configs['run_id'] + '.fasta.gz'
    # assert not path.isfile(args.output_file), '[e] output file already exists: {:s}'.format(args.output_file)
    if not path.isdir(path.dirname(args.output_file)):
        makedirs(path.dirname(args.output_file))
    print('Writing reads with traceable identifiers to: {:s}'.format(args.output_file))

    # loop over reads
    with gzip.open(args.output_file, 'w') as out_fid:
        inp_flst = args.input_file.split(',')
        print('Total of [{:d}] files are given as input.'.format(len(inp_flst)))
        for inp_fidx, inp_fname in enumerate(inp_flst):
            print('\t{:d}. Reading from: {:s}'.format(inp_fidx + 1, inp_fname))
            rd_nid = 1
            with gzip.open(inp_fname, 'r') as inp_fid:
                while True:
                    rd_oid = inp_fid.readline()
                    if rd_oid == '':
                        break
                    rd_oid = rd_oid.rstrip('\n')

                    rd_seq = inp_fid.readline().rstrip('\n')
                    rd_plus = inp_fid.readline().rstrip('\n')
                    rd_pred = inp_fid.readline().rstrip('\n')
                    n_nt = len(rd_seq)
                    if (rd_oid[0] != '@') or (rd_plus != '+') or (n_nt == 0):
                        raise Exception('[e] the input file is corrupted.\n' +
                                        'Read #{:d}:\n'.format(rd_nid) +
                                        '\tID: [{:s}],\n\tplus: [{:s}]'.format(rd_oid, rd_plus))
                    if rd_nid % 50000 == 0:
                        print('\t\tprocessed {:,d} reads.'.format(rd_nid))

                    rd_sid = 'Fl.Id:{:d};Rd.Id:{:d};Rd.Ln:{:d}'.format(inp_fidx + 1, rd_nid, n_nt)
                    out_fid.write('>' + rd_sid + '\n')
                    out_fid.write(rd_seq + '\n')
                    rd_nid = rd_nid + 1
    print('[i] Read identifier conversion is completed successfully.')


def splitReads(args):
    import re
    import gzip

    from utilities import load_configs, get_re_info

    print('%% Splitting reads into fragments ...')
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
        print('[i] {:,d} fragments are reduced to {:,d}bp.'.format(n_reduced, MAX_FRG_SIZE))
    print('[i] Total of {:,d} reads and {:,d} fragments are produced successfully.'.format(rd_ind - 1, frg_ind - 1))


def mapFragments(args):

    from utilities import load_configs
    if not args.return_command:
        print('%% Mapping fragments to genome ...')

    configs = load_configs(args.config_file, max_n_configs=1)[0]
    map_argument = '-b 5 -q 2 -r 1 -z 5 -T 15'

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
        configs['bwa'] + ' bwasw {:s} -t {:d} '.format(map_argument, args.n_thread) + \
        configs['bwa_index'] + ' ' + args.input_file + \
        ' | samtools sort -n | samtools view -q 1 -hbS - ' + \
        '> ' + args.output_file
    if args.return_command:
        print('{:s}'.format(cmd_str))
    else:
        assert path.isfile(args.input_file), '[e] Input file could not be found: {:s}'.format(args.input_file)
        print('Running: {:s}'.format(cmd_str))
        import subprocess
        map_prs = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        std_out, std_err = map_prs.communicate()
        # TODO: A better error handling here would be nice
        assert std_err.split('\n')[-2][:18] == '[main] Real time: ', \
            'bwa failed to run properly, see below:\n{:s}'.format(std_err)
        print('[i] Fragments are mapped to genome successfully.')


def processMappedFragments(args):
    import h5py
    from os import remove
    from pandas import read_csv
    import pysam
    import gzip

    from utilities import load_configs, get_chr_info, hasOL

    print('%% Creating an MC-4C dataset from mapped fragments ...')
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
        print('Database of restriction enzyme cut sites is not found. ' + \
              'Scanning the reference genome to create this database ...')
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
    if idx_lft + 1 != idx_rgt:
        print('[w] Can not map primer positions on a single fragment.')

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
            if que_idx % 500000 == 0:
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
        print('[w] {:,d} fused reads ({:0.1f}% of total) are identified and flagged.'.format(n_fusion, n_fusion * 100.0 / n_processed))

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
    print('[i] MC4C dataset is created successfully.')


def removeDuplicates(args):
    import h5py
    from os import path, makedirs
    import pandas as pd

    from utilities import load_configs, hasOL

    print('%% Removing pcr duplicates from an MC-4C dataset ...')
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
    print('There are {:d} reads in the dataset.'.format(len(np.unique(mc4c_pd['ReadID']))))

    # filtering reads according to their MQ
    read_all = mc4c_pd[['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'MQ', 'Flag']].values
    is_mapped = read_all[:, 4] >= args.min_mq
    is_valid = np.bitwise_and(read_all[:, 5], 1) == 0
    read_all = read_all[is_mapped & is_valid, :4]
    print('Selected non-overlapping fragments with MQ >= {:d}: {:d} reads are left.'.format(
        args.min_mq, len(np.unique(read_all[:, 0]))))
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
    print('Selected reads #cis fragment > 1: {:d} reads are selected.'.format(len(np.unique(read_inf[:, 0]))))

    # select reads with #traceable fragment > 1
    roi_size = configs['roi_end'] - configs['roi_start']
    lcl_crd = np.array([configs['vp_cnum'], configs['roi_start'] - roi_size, configs['roi_end'] + roi_size])
    is_lcl = hasOL(lcl_crd, read_inf[:, 1:4], offset=0)
    umi_set = read_inf[~is_lcl, :]
    print('Selected reads with #trans fragment > 0: {:d} reads are selected.'.format(len(np.unique(umi_set[:, 0]))))

    # remove duplicates
    unq_set, duplicate_info = remove_duplicates_by_umi(umi_set, verbose=True)
    print('Result statistics (before --> after filtering):')
    print('\t#reads: {:,d} --> {:,d}'.format(len(np.unique(mc4c_pd['ReadID'])), len(np.unique(unq_set[:, 0]))))
    print('\t#fragments: {:,d} --> {:,d}'.format(mc4c_pd['ReadID'].shape[0], unq_set.shape[0]))

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
    print('[i] PCR duplicates are removed from MC4C dataset successfully.')

