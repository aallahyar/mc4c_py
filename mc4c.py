
import argparse
import sys
import gzip
from os import path, makedirs
import numpy as np
import pandas as pd

import mc4c_tools
import loger

flag_DEBUG = True

def makePrimerFasta(args):
    """ Turn primer sequences into a fasta file.
    """
    configs = mc4c_tools.load_configs(args.cnfFile)
    primerSeqs = mc4c_tools.getPrimerFragment(configs)
    mc4c_tools.writePrimerFasta(primerSeqs, args.outfile)


def cleaveReads(args):
    """ Cleave the reads by primer sequences. Requires BowTie2 information. """
    settings = mc4c_tools.load_configs(args.cnfFile)
    primerLens = [len(x) for x in settings['prm_seq']]
    primers = ['']
    primers.extend(settings['prm_seq'])
    #print primers
    prmCuts = mc4c_tools.combinePrimers(args.bamfile, primerLens)
    #print prmCuts[:10]
    mc4c_tools.applyCuts(args.fastqfile,args.outfile,prmCuts,primers)


# def splitReads(args):
#     """ Split the reads by restriction site information based on the reference genome. """
#     settings = mc4c_tools.load_config(args.cnfFile)
#     restSeqs = settings['re_seq']
#     # TODO: Substitute reference genome with reads (?)
#     mc4c_tools.findRestrictionSeqs(args.fastqfile,args.outfile,restSeqs)


def findRefRestSites(args):
    """ Determine the location of restriction sites on the reference genome. Takes a fasta file
        and stores results as a list per chromosome in a dictionary, which is saved as an npz.
    """
    from utilities import extract_re_positions

    configs = mc4c_tools.load_configs(args.cnfFile)
    extract_re_positions(genome_str=configs['genome_build'], re_name_lst=configs['re_name'])

def getRefResPositions(args):
    """ Extract a subset of restriction site positions from the reference genome. """
    settings = mc4c_tools.load_configs(args.cnfFile)
    print [settings['vp_chr']],[settings['vp_start'], settings['vp_end']]
    print 'Loading restrsites, this takes a while...'
    restrefs=np.load(args.restfile)['restrsites'].item()
    print 'Finished loading, moving on'
    result = mc4c_tools.mapToRefSite(restrefs[settings['vp_chr'][0]],[settings['vp_start'][0], settings['vp_end'][0]])

    refPosList = []

    for i in range(result[0],result[1]+1):
        #print i,restrefs[settings['vp_chr'][0]][i]
        refPosList.append(restrefs[settings['vp_chr'][0]][i])

    pdFrame = pd.DataFrame(refPosList, index=range(result[0],result[1]+1), columns=['start','stop'])

    np.savez_compressed(args.outfile,
                        pdframe=pdFrame,
                        pdcolumns=pdFrame.columns,
                        pdindex=pdFrame.index)

def exportToPlot(args):
    """ Originally written to easily import the data into interactive plotting tools.
        Converts the mapped data to a pandas dataframe and adds restriction site information.
        Additionally it creates 2 files that link between restrition sites and read ids for
        interaction down the line.
    """
    settings = mc4c_tools.load_configs(args.cnfFile)
    print 'Loading restrsites, this takes a while...'
    restrefs=np.load(args.restfile)['restrsites'].item()
    print 'Finished loading, moving on'
    byRegion,byRead,pdFrame = mc4c_tools.exportToPlot(settings,restrefs,args.bamfile)

    #dupSet = mc4c_tools.findDuplicates(settings,byRead,byRegion)
    #pdFrame['Duplicate'] = np.where(pdFrame['CircleId'].isin(dupSet), True, False)

    #print pdFrame
    np.savez_compressed(args.plotfile,
                        pdframe=pdFrame,
                        pdcolumns=pdFrame.columns,
                        pdindex=pdFrame.index)

    np.savez_compressed(args.plotfile+'_extra',
                        byregion=byRegion,
                        byread=byRead)


def markDuplicates(args):
    """ This function aims to identify reads that are most likely PCR duplicates.
        Identification is based on having overlap with eachother that is not in the viewport.
        It takes a pandas dataframe and adds a new column to the end of it.
    """
    settings = mc4c_tools.load_configs(args.cnfFile)
    exFile = np.load(args.extra)

    try:
        byRead = exFile['byread'].item()
    except KeyError:
        byRead = exFile['byreads'].item()
    byRegion = exFile['byregion'].item()

    pdFile = np.load(args.pdframe)
    pdFrame = pd.DataFrame(pdFile['pdframe'],columns=pdFile['pdcolumns'],index=pdFile['pdindex'])
    dupSet = mc4c_tools.findDuplicates(settings,byRead,byRegion)

    #df['dup']=np.where(pd.Series(df.index).isin([1,5]),True,False)
    #pdFrame['Duplicate'] = np.where(pdFrame['CircleId'].isin(dupSet), True, False)

    pdFrame['Duplicate'] = np.where(pd.Series(pdFrame.index).isin(dupSet), True, False)

    np.savez_compressed(args.outfile,
                        pdframe=pdFrame,
                        pdcolumns=pdFrame.columns,
                        pdindex=pdFrame.index)


def flattenFragments(args):
    """ This function aims to identify parts of reads that are repeats of themselves, overlapping
        the same regions multiple times.
        It takes a pandas dataframe and adds a new column to the end of it.
    """

    pdFile = np.load(args.pdframe)
    pdFrame = pd.DataFrame(pdFile['pdframe'],columns=pdFile['pdcolumns'],index=pdFile['pdindex'])
    mc4c_tools.findRepeats(pdFrame)

    print pdFrame.iloc[:10].T
    np.savez_compressed(args.outfile,
                        pdframe=pdFrame,
                        pdcolumns=pdFrame.columns,
                        pdindex=pdFrame.index)
#############################################


def initialize_run(args):
    # TODO:
    # test refrence genome path
    # test existance of bwa
    # test if bwa index is present
    # check and make RE position file if non existant
    # test existance of samtools
    print '[TO DO] Checking the pipeline for related files for given experiment'


def setReadIds(args):
    configs = mc4c_tools.load_configs(args.cnfFile)

    # initialize
    if args.input_file is None:
        args.input_file = './raw_files/raw_' + configs['run_id'] + '.fastq.gz'
    if args.output_file is None:
        args.output_file = './read_files/rd_' + configs['run_id'] + '.fasta.gz'
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

    configs = mc4c_tools.load_configs(args.cnfFile)

    if args.input_file is None:
        args.input_file = './read_files/rd_' + configs['run_id'] + '.fasta.gz'
    if args.output_file is None:
        args.output_file = './frg_files/frg_' + configs['run_id'] + '.fasta.gz'
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


def mapFragments(args):
    configs = mc4c_tools.load_configs(args.cnfFile)

    # Map split fragments to genome
    if args.input_file is None:
        args.input_file = './frg_files/frg_' + configs['run_id'] + '.fasta.gz'
    if args.output_file is None:
        args.output_file = './bam_files/bam_{:s}_{:s}.bam'.format(configs['run_id'], configs['genome_build'])
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
        import subprocess
        map_prs = subprocess.Popen(cmd_str, shell=True, stdout=subprocess.PIPE)
        std_out, std_err = map_prs.communicate()
        assert std_err == ''


def process_mapped_fragments(args):
    import h5py
    from os import remove
    from pandas import read_csv
    import pysam

    from utilities import get_chr_info, get_re_info

    configs = mc4c_tools.load_configs(args.cnfFile)

    if args.input_file is None:
        args.input_file = './bam_files/bam_{:s}_{:s}.bam'.format(configs['run_id'], configs['genome'])
    if args.output_file is None:
        args.output_file = './mc4c_files/mc4c_' + configs['run_id'] + '.hdf5'
    assert not path.isfile(args.output_file), 'Output file already exists: {:s}'.format(args.output_file)

    chr_lst = get_chr_info(configs['genome'], 'chr_name')
    chr_map = dict(zip(chr_lst, np.arange(len(chr_lst)) + 1))

    # extend fragment coordinates to reference genome
    re_pos_fname = './renz_files/{:s}_{:s}.npz'.format(configs['genome'], configs['genome'].replace(';', '-'))
    re_pos = np.load(re_pos_fname)

    # Read file line by line
    ReadID_old = -1
    FrgID_old = -1
    n_frg_info = 12
    n_processed = 0
    frg_template = '\t'.join(['{:d}'] * n_frg_info) + '\n'
    frg_set = np.empty([0, n_frg_info], dtype=np.int64)
    tmp_fname = args.output_file + '.tmp'
    print('Writing processed fragments to a temprary file first: {:s}'.format(tmp_fname))
    with pysam.AlignmentFile(args.input_file, 'rb') as bam_fid, gzip.open(tmp_fname, 'wb') as gz_fid:
        gz_fid.write(
            '\t'.join(['ReadID', 'Chr', 'RefStart', 'RefEnd', 'Strand', 'ExtStart', 'ExtEnd', 'MQ',
                       'FileID', 'FrgID', 'SeqStart', 'SeqEnd', 'ReadLength', 'TrueHop']) + '\n'
        )
        for que_idx, que_line in enumerate(bam_fid):
            if que_idx % 100000 == 0:
                print('Processed {:,d} fragments in {:,d} reads.'.format(que_idx, n_processed))

            if (np.bitwise_and(que_line.flag, 0x800) == 0x800) or (que_line.reference_name not in chr_lst):
                continue
            FileID, ReadID, ReadLength, FrgID, SeqStart, SeqEnd = \
                [int(x.split(':')[1]) for x in que_line.query_name.split(';')]
            RefChrNid = chr_map[que_line.reference_name]
            RefStart = que_line.reference_start
            RefEnd = que_line.reference_end
            RefStrand = 1 - (que_line.is_reverse * 2)

            # extending coordinates to nearby restriction site
            nei_left = np.searchsorted(re_pos[RefChrNid - 1], RefStart, side='right')
            nei_right = np.searchsorted(re_pos[RefChrNid - 1], RefEnd, side='left')
            if nei_left == nei_right:
                dist_left = RefStart - re_pos[RefChrNid - 1][nei_left]
                dist_right = RefEnd - re_pos[RefChrNid - 1][nei_right]
                if dist_right > dist_left:
                    nei_left = nei_left - 1
                else:
                    nei_right = nei_right - 1
            ExtStart = re_pos[RefChrNid - 1][nei_left]
            ExtEnd = re_pos[RefChrNid - 1][nei_right]
            # TODO: Needs to account for unmapped fragments

            # combine into an array
            frg_info = np.array([
                ReadID, RefChrNid, RefStart, RefEnd, RefStrand, ExtStart, ExtEnd, que_line.mapping_quality,
                FileID, FrgID, SeqStart, SeqEnd, ReadLength, 1]).reshape([1, -1])

            # Check order of fragments
            if FrgID < FrgID_old:
                raise Exception('Order of fragments are lost.')
            FrgID_old = FrgID

            # Check if read has ended
            if ReadID == ReadID_old:
                frg_set = np.vstack([frg_set, frg_info])
                continue

            # Save the read
            for frg in frg_set:
                gz_fid.write(frg_template.format(*frg))
            frg_set = frg_info.copy()
            ReadID_old = ReadID

            n_processed += 1

        if frg_set.shape[0] != 0:  # Saving the last read after file has finished
            for frg in frg_set:
                gz_fid.write(frg_template.format(*frg))

    # Load fragments in pandas format and sort fragments within a read according to their relative positions
    print('Reading temporary file: {:s}'.format(tmp_fname))
    frg_pd = read_csv(tmp_fname, delimiter='\t', compression='gzip')
    frg_pd = frg_pd.iloc[np.lexsort([frg_pd['SeqStart'], frg_pd['ReadID']])]

    print('Writing: {:s}'.format(args.output_file))
    h5_fid = h5py.File(args.output_file, 'w')
    if frg_pd.shape[1] > 5000:
        h5_fid.create_dataset('frg_np', data=frg_pd.values, compression='gzip', compression_opts=5,
                              chunks=(5000, frg_pd.shape[1]))
    else:
        h5_fid.create_dataset('frg_np', data=frg_pd.values, compression='gzip', compression_opts=5)
    h5_fid.create_dataset('frg_np_header_lst', data=list(frg_pd.columns.values))
    h5_fid.create_dataset('chr_lst', data=chr_lst)
    h5_fid.close()

    print('Removing temporary file: {:s}'.format(tmp_fname))
    remove(tmp_fname)


# Huge wall of argparse text starts here
def main():
    """ Everything in here is to interpret calls from a command line.
        Anything being run will just call a similarly named function
        above.
    """

    parser = argparse.ArgumentParser(description="MC4C pipeline for processing multi-contact 4C data")
    subparsers = parser.add_subparsers()

    # init command
    parser_init = subparsers.add_parser('initRun',
                                          description='Initialize and prepare the pipeline for given configuration')
    parser_init.add_argument('cnfFile',
                               type=str,
                               help='Configuration file containing experiment specific details')
    parser_init.set_defaults(func=initialize_run)

    # Set read identifiers
    parser_readid = subparsers.add_parser('setReadIds',
                                         description='Defining identifiers for sequenced reads')
    parser_readid.add_argument('cnfFile',
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
    parser_readid = subparsers.add_parser('splitReads',
                                         description='Splitting reads into fragments using restriction enzyme recognition sequence')
    parser_readid.add_argument('cnfFile',
                              type=str,
                              help='Configuration file containing experiment specific details')
    parser_readid.add_argument('--input-file',
                               default=None,
                               type=str,
                               help='Input file (in FASTA format) containing reads with traceable IDs.')
    parser_readid.add_argument('--output-file',
                               default=None,
                               type=str,
                               help='Output file (in FASTA format) containing fragments with traceable IDs')
    parser_readid.set_defaults(func=splitReads)

    # Mapping fragments to reference genome
    parser_readid = subparsers.add_parser('mapFragments',
                                          description='Mapping fragments to reference genome')
    parser_readid.add_argument('cnfFile',
                               type=str,
                               help='Configuration file containing experiment specific details')
    parser_readid.add_argument('--input-file',
                               default=None,
                               type=str,
                               help='Input file (in FASTA format) containing fragments with traceable IDs')
    parser_readid.add_argument('--output-file',
                               default=None,
                               type=str,
                               help='Output file (in BAM format) containing fragments with traceable IDs')
    parser_readid.add_argument('--n_thread',
                               default=1,
                               type=int,
                               help='Number of threads should be used by the aligner')
    parser_readid.add_argument('--return_command',
                               action="store_true",
                               help='Return only mapping command instead of running it ' +
                                    '(useful for running the pipeline in a cluster)')
    parser_readid.set_defaults(func=mapFragments)

    # Process mapped fragments
    parser_readid = subparsers.add_parser('makeDataset',
                                          description='Processed the mapped fragments and create a MC-4C dataset')
    parser_readid.add_argument('cnfFile',
                               type=str,
                               help='Configuration file containing experiment specific details')
    parser_readid.add_argument('--input-file',
                               default=None,
                               type=str,
                               help='Input file (in BAM format) containing fragments with traceable IDs')
    parser_readid.add_argument('--output-file',
                               default=None,
                               type=str,
                               help='Output file (in HDF5 format) containing processed fragments')
    parser_readid.set_defaults(func=process_mapped_fragments)

    if flag_DEBUG:
        # sys.argv = ['./mc4c.py', 'init', './cnf_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'setReadIds', './cnf_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'splitReads', './cnf_files/cfg_LVR-BMaj.cnf']
        sys.argv = ['./mc4c.py', 'mapFragments', './cnf_files/cfg_LVR-BMaj.cnf']
    args = parser.parse_args(sys.argv[1:])
    # loger.printArgs(args)
    args.func(args)


if __name__ == '__main__':
    main()
