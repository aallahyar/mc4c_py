#! /usr/bin/env python2

# cluster run example:
# qsub -P hub_laat -N mc4c -l h_rt=05:00:00 -l h_vmem=50G -pe threaded 1 ~/bulk/bin/run_script.sh "python2 ./mc4c.py setReadIds LVR-BMaj"
# downloader: aria2c -x 4 -s 4 --file-allocation=none URL

# compress raw: tar -czf ESCs-V122-PM.tar.gz ./deLaatUDNHiC5Pro3
# check integrity: tar -tzf my_tar.tar.gz > /dev/null
# combine fasq: cat *.fastq | gzip > ../../ESCs-V121-PM.fastq.gz

# fragment flag
# desc: Bits in fragment flag -> 1:overlapping frags, 2:fused-reads, 3:

# TODO: Use frequent k-mers to check primer sequences
# TODO: Check if a dup_set of reads share a second fragment in ROI

import argparse
import sys
from os import path, makedirs, environ
import numpy as np
import pandas as pd

import pre_process

# import platform  # ###
# flag_DEBUG = platform.system() != 'Linux'

np.set_printoptions(linewidth=180, threshold=5000)  # , suppress=True, formatter={'float_kind':'{:0.5f}'.format}
pd.set_option('display.width', 180)
pd.set_option('display.max_rows', 200)
pd.set_option('display.max_columns', 15)


def perform_qc(args):
    import quality_check
    from utilities import load_configs

    # load config files
    config_lst = load_configs(args.config_file)

    if args.output_file is None:
        config_lst[0]['output_dir'] = './plots/'
    else:
        config_lst[0]['output_dir'] = path.dirname(args.output_file)
    if not path.isdir(config_lst[0]['output_dir']):
        makedirs(config_lst[0]['output_dir'])
    config_lst[0]['input_file'] = args.input_file
    config_lst[0]['output_file'] = args.output_file
    config_lst[0]['sigma'] = args.sigma

    # call the requested function
    if args.report_type == 'readSizeDist':
        quality_check.plot_readSizeDistribution(config_lst)
    elif args.report_type == 'frgSizeDist':
        quality_check.plot_frg_size_distribution(config_lst)
    elif args.report_type == 'chrCvg':
        quality_check.plot_chrCvg(config_lst)
    elif args.report_type == 'cirSizeDist':
        quality_check.plot_cirSizeDistribution(config_lst, roi_only=args.roi_only, uniq_only=args.uniq_only)
    elif args.report_type == 'overallProfile':
        quality_check.plot_overallProfile(config_lst, min_n_frg=2)
    elif args.report_type == 'seqSaturation':
        assert len(config_lst) == 1, '[error] Only one dataset can be tested for saturation!'
        quality_check.plot_sequencing_saturation(config_lst[0])
    elif args.report_type == 'categorizeReads':
        quality_check.plot_reads_per_category(config_lst)
    else:
        raise Exception()
    print('[i] {:s} plot is produced successfully.'.format(args.report_type))


def perform_analysis(args):
    from copy import deepcopy
    import pandas as pd

    import analysis
    from utilities import load_configs

    config_lst = load_configs(args.config_file)
    if args.output_file is None:
        config_lst[0]['output_dir'] = './plots/'
    else:
        config_lst[0]['output_dir'] = path.dirname(args.output_file)
    if not path.isdir(config_lst[0]['output_dir']):
        makedirs(config_lst[0]['output_dir'])
    config_lst[0]['input_file'] = args.input_file
    config_lst[0]['output_file'] = args.output_file
    config_lst[0]['zscr_lim'] = [-args.zscr_lim, args.zscr_lim]
    config_lst[0]['correction'] = args.correction
    config_lst[0]['cmd_args'] = args

    # call the requested function
    if args.analysis_type == 'atVpSoi':
        if args.ant_name is None:
            from utilities import load_annotation
            roi_crd = [config_lst[0]['vp_cnum'], config_lst[0]['roi_start'], config_lst[0]['roi_end']]
            ant_pd = load_annotation(config_lst[0]['genome_build'], roi_crd=roi_crd)
            ant_name_lst = ant_pd['ant_name'].values
        else:
            ant_name_lst = args.ant_name.split(',')

        zscr_pd = pd.DataFrame()
        for ant_name in ant_name_lst:
            print('Preparing VP-SOI for [{:s}]'.format(ant_name))
            ant_pd = analysis.perform_vpsoi_analysis(deepcopy(list(config_lst)), soi_name=ant_name, min_n_frg=2)
            if len(zscr_pd) == 0:
                zscr_pd = ant_pd[['ant_name', 'ant_chr', 'ant_pos']].copy()
            assert zscr_pd[['ant_name', 'ant_chr', 'ant_pos']].equals(ant_pd[['ant_name', 'ant_chr', 'ant_pos']])
            zscr_pd['from_' + ant_name] = ant_pd['zscore']

        if args.to_tsv:
            run_id = ','.join([config['run_id'] for config in config_lst])
            tsv_fname = path.join(config_lst[0]['output_dir'],
                                  'analysis_atVP-SOI_{:s}_{:s}_'.format(run_id, ','.join(ant_name_lst)) +
                                  'sig{:0.2f}_corr-{:s}_'.format(args.sigma, args.correction) +
                                  'np{:0.1f}k.tsv'.format(args.n_perm / 1e3))
            zscr_pd.to_csv(tsv_fname, sep='\t', index=False, header=True)

    elif args.analysis_type == 'atSOISOI':
        analysis.perform_soisoi_analysis(list(config_lst), min_n_frg=2, n_perm=args.n_perm)
    elif args.analysis_type == 'atAcrossROI':
        analysis.perform_at_across_roi(list(config_lst), min_n_frg=2)
    else:
        raise Exception()
    print('[i] {:s} analysis is performed successfully.'.format(args.analysis_type))


def main():
    """ Everything in here is to interpret calls from a command line.
        Anything being run will just call a similarly named function
        above.
    """

    parser = argparse.ArgumentParser(description="MC4C pipeline for processing multi-contact 4C data")
    subparsers = parser.add_subparsers()

    # Set read identifiers
    parser_readid = subparsers.add_parser('setReadIds', description='Defining identifiers for sequenced reads')
    parser_readid.add_argument('config_file', metavar='config-file', type=str,
                               help='Configuration file containing experiment specific details')
    parser_readid.add_argument('--input-file', default=None, type=str,
                               help='Input file (in FASTQ format) containing raw sequenced reads')
    parser_readid.add_argument('--output-file', default=None, type=str,
                               help='Output file (in FASTA format) containing sequenced reads with traceable IDs')
    parser_readid.set_defaults(func=pre_process.setReadIds)

    # Split reads into fragments
    parser_readSplt = subparsers.add_parser('splitReads', description='Splitting reads into fragments using ' +
                                                                      'restriction enzyme recognition sequence')
    parser_readSplt.add_argument('config_file', metavar='config-file', type=str,
                                 help='Configuration file containing experiment specific details')
    parser_readSplt.add_argument('--input-file', default=None, type=str,
                                 help='Input file (in FASTA format) containing reads with traceable IDs.')
    parser_readSplt.add_argument('--output-file', default=None, type=str,
                                 help='Output file (in FASTA format) containing fragments with traceable IDs')
    parser_readSplt.set_defaults(func=pre_process.splitReads)

    # Mapping fragments to reference genome
    parser_mapFrg = subparsers.add_parser('mapFragments', description='Mapping fragments to reference genome')
    parser_mapFrg.add_argument('config_file', metavar='config-file', type=str,
                               help='Configuration file containing experiment specific details')
    parser_mapFrg.add_argument('--input-file', default=None, type=str,
                               help='Input file (in FASTA format) containing fragments with traceable IDs')
    parser_mapFrg.add_argument('--output-file', default=None, type=str,
                               help='Output file (in BAM format) containing mapped fragments')
    parser_mapFrg.add_argument('--n_thread', default=1, type=int,
                               help='Number of threads should be used by the aligner')
    parser_mapFrg.add_argument('--return_command', action="store_true",
                               help='Return only mapping command instead of running it (useful for running the pipeline in a cluster)')
    parser_mapFrg.set_defaults(func=pre_process.mapFragments)

    # Process mapped fragments
    parser_mkDataset = subparsers.add_parser('makeDataset',
                                             description='Processed the mapped fragments and create an MC-4C dataset')
    parser_mkDataset.add_argument('config_file', metavar='config-file', type=str,
                                  help='Configuration file containing experiment specific details')
    parser_mkDataset.add_argument('--input-file', default=None, type=str,
                                  help='Input file (in BAM format) containing fragments with traceable IDs')
    parser_mkDataset.add_argument('--output-file', default=None, type=str,
                                  help='Output file (in HDF5 format) containing processed fragments')
    parser_mkDataset.set_defaults(func=pre_process.processMappedFragments)

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
    parser_remDup.set_defaults(func=pre_process.removeDuplicates)

    # produce statistics plots
    parser_qc = subparsers.add_parser('QC',
                                            description='Generate various summary reports about an MC-4C dataset.')
    parser_qc.add_argument('report_type', type=str,
                                 choices=['readSizeDist', 'frgSizeDist', 'chrCvg', 'cirSizeDist', 'overallProfile',
                                          'seqSaturation', 'categorizeReads'],
                                 help='Type of summary report that needs to be generated')
    parser_qc.add_argument('config_file', metavar='config-file', type=str,
                                 help='Configuration file containing experiment specific details')
    parser_qc.add_argument('--input-file', default=None, type=str,
                                 help='Input file (in HDF5 format) containing MC4C data.')
    parser_qc.add_argument('--output-file', default=None, type=str,
                                 help='Output file (in PDF format) containing the requested summary report.')
    parser_qc.add_argument('--roi-only', action="store_true", default=False,
                                  help='Limits the requested summary report to be generated from roi-fragments only.')
    parser_qc.add_argument('--uniq-only', action="store_true", default=False,
                                  help='Limits the requested summary report to unique (duplicate removed) reads.')
    parser_qc.add_argument('--min-cvg', default=2, type=int,
                           help='Minimum chromosome coverage (%) used for automatic detection of ' +
                                'Region Of Interest (ROI). Default is 2.')
    parser_qc.add_argument('--sigma', default=0.0, type=float,
                           help='Sigma for Gaussian smoothing (default=0, i.e. no smoothing)')
    parser_qc.set_defaults(func=perform_qc)

    # perform basic analysis
    parser_analysis = subparsers.add_parser('analysis', description='Performs analysis on an MC-4C dataset.')
    parser_analysis.add_argument('analysis_type', choices=['mcTest', 'atVpSoi', 'atSOISOI', 'atAcrossROI'], type=str,
                                  help='Type of analysis that needs to be performed')
    parser_analysis.add_argument('config_file', metavar='config-file', type=str,
                                  help='Configuration file containing experiment specific details')
    parser_analysis.add_argument('--ant-name', default=None, type=str,
                                 help='Name of annotation for which VP-SOI plot needs to be computed. ' +
                                      'Only used for VP-SOI (i.e. "vpSoi") analysis. Every annotation within the '
                                      'Region Of Interest (ROI) will be used if this argument is not provided.')
    parser_analysis.add_argument('--input-file', default=None, type=str,
                                  help='Input file (in HDF5 format) containing MC4C data.')
    parser_analysis.add_argument('--output-file', default=None, type=str,
                                  help='Output file (in PDF format) containing the result of the requested analysis.')
    parser_analysis.add_argument('--n-perm', default=1000, type=int,
                                 help='Number of profiles that needs to be drawn from negative reads (i.e. reads ' +
                                      'that contain no fragment from site of interest) to produce the expected profile.')
    parser_analysis.add_argument('--sigma', default=0.0, type=float,
                                 help='Sigma for Gaussian smoothing (default=0, i.e. no smoothing)')
    parser_analysis.add_argument('--downsample', default=None, type=int, help='Downsample dataset before the analysis')
    parser_analysis.add_argument('--to_tsv', action="store_true", help='Store the z-scores to a tab-separated file')
    parser_analysis.add_argument('--zscr_lim', default=6, type=float)
    parser_analysis.add_argument('--correction', default='none', type=str, choices=['none', 'decay'])
    parser_analysis.add_argument('--cvg_norm', default='none', type=str, choices=['none', 'iter', 'KR', '1d', '1dNoScale'])
    parser_analysis.set_defaults(func=perform_analysis)

    #if hasattr(sys.stderr, "isatty") and sys.stderr.isatty():
    if 'PYCHARM_HOSTED' in environ:
        Warning('[w] This script is being run interactively. Running in debug mode!')
        # sys.argv = ['./mc4c.py', 'process', 'LVR-BMaj-PB']
        # sys.argv = ['./mc4c.py', 'init', './cfg_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'setReadIds', './cnf_files/cfg_LVR-BMaj.cnf']
        # sys.argv = ['./mc4c.py', 'splitReads', 'LVR-BMaj']
        # sys.argv = ['./mc4c.py', 'mapFragments', 'BMaj-test']
        # sys.argv = ['./mc4c.py', 'makeDataset', 'BMaj-test']
        # sys.argv = ['./mc4c.py', 'selectROI', 'BMaj-test']
        # sys.argv = ['./mc4c.py', 'removeDuplicates', 'BMaj-test']
        # sys.argv = ['./mc4c.py', 'QC', 'readSizeDist', 'BMaj-test,BMaj-test']
        # sys.argv = ['./mc4c.py', 'QC', 'frgSizeDist', 'BMaj-test,BMaj-test']
        # sys.argv = ['./mc4c.py', 'QC', 'chrCvg', 'BMaj-test,BMaj-test']
        # sys.argv = ['./mc4c.py', 'QC', 'cirSizeDist', 'LVR-BMaj-96x'] # , '--roi-only', '--uniq-only'
        # sys.argv = ['./mc4c.py', 'QC', 'overallProfile', 'BMaj-test,BMaj-test']

        # sys.argv = ['./mc4c.py', 'analysis', 'mcTest', 'K562-WplD-10x']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--n-perm=1000', 'BMaj-test,BMaj-test', '--ant-name', 'HS2']
        # sys.argv = ['./mc4c.py', 'analysis', 'atSOISOI', '--n-perm=1000', 'LVR-BMaj-96x-Adj']
        # sys.argv = ['./mc4c.py', 'analysis', 'atAcrossROI', '--n-perm=10', 'BMaj-test']

        # sys.argv = ['./mc4c.py', 'QC', 'seqSaturation', 'BMaj-test']

        # sys.argv = ['./mc4c.py', 'QC', 'categorizeReads', 'BMaj-test']

        # sys.argv = ['./mc4c.py', 'QC', 'overallProfile', '--sigma=1.0', 'Prdm14_Slc_WT,Prdm14_Slc_WT2,Prdm14_Slc_WT3']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--n-perm=1000', 'asMC4C_WT_B'] # , '--ant-name', 'HS2'
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', 'Prdm14_RB_LB-DEL,Prdm14_RB_LB-DEL2']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', 'K562_C11-Enh-3769_WT']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=2.0', '--n-perm=10', '--ant-name=CLOCK', 'K562_C4-Enh-1627_WT']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=default', '--ant-name=RB', 'Prdm14_Slc_WT,Prdm14_Slc_WT2,Prdm14_Slc_WT3']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=HS2', 'LVR-BMaj-96x,LVR-BMaj-NP']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=HS2', 'BRN-BMaj-96x,BRN-BMaj-96x2']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=default', '--ant-name=HS4', 'LVR-BMaj-96x,LVR-BMaj-NP']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', 'Prdm14_Slc_WT,Prdm14_Slc_WT2,Prdm14_Slc_WT3']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', 'WPL-WTD,WPL-WTD2']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', 'WPL-WTC']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=K', 'WPL-KOD,WPL-KOD2']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', 'WPL-KOC']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=H', 'WPL-KOD,WPL-KOD2']
        # sys.argv = ['./mc4c.py', 'analysis', 'atAcrossROI', '--n-perm=10', '--downsample=10000', '--to_tsv', 'LVR-BMaj-96x,LVR-BMaj-NP'] # BRN-BMaj-96x,BRN-BMaj-96x2
        # sys.argv = ['./mc4c.py', 'analysis', 'atAcrossROI', '--n-perm=10', 'Prdm14_Slc_LB-DEL']
        # sys.argv = ['./mc4c.py', 'analysis', 'atAcrossROI', '--n-perm=10', '--sigma=2.0', 'K562_C4-Enh-1627_WT']
        # sys.argv = ['./mc4c.py', 'analysis', 'atAcrossROI', '--n-perm=10', '--cvg_norm=1d', '--correction=decay', '--sigma=1.0', 'Prdm14_Slc_WT,Prdm14_Slc_WT2,Prdm14_Slc_WT3']
        # sys.argv = ['./mc4c.py', 'analysis', 'atAcrossROI', '--n-perm=10', '--correction=decay', '--zscr_lim=6', '--sigma=1.0', 'Prdm14_Slc_WT,Prdm14_Slc_WT2,Prdm14_Slc_WT3']
        # sys.argv = ['./mc4c.py', 'analysis', 'atAcrossROI', '--n-perm=10', '--sigma=1.0', '--cvg_norm=iter', '--correction=none', 'BMaj-test,BMaj-test']
        # sys.argv = ['./mc4c.py', 'analysis', 'atAcrossROI', '--n-perm=100', '--sigma=0.0', '--correction=decay', 'BRN-BMaj-96x,BRN-BMaj-96x2']
        # sys.argv = ['./mc4c.py', 'analysis', 'atAcrossROI', '--n-perm=2', '--cvg_norm=iter', '--sigma=1.0', '--correction=decay', 'LVR-BMaj-96x,LVR-BMaj-NP']
        # sys.argv = ['./mc4c.py', 'analysis', 'atAcrossROI', '--n-perm=100', '--sigma=1.0', '--correction=decay', 'LVR-HS2-96x,LVR-HS2-NP,LVR-HS2-NP2']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--to_tsv', '--ant-name=HS2', 'LVR-BMaj-96x,LVR-BMaj-NP']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=HS2', 'BRN-BMaj-96x,BRN-BMaj-96x2']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=LB', 'Prdm14_RB_WT']
        sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=tmp', 'Prdm14_Slc_WT,Prdm14_Slc_WT2,Prdm14_Slc_WT3']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=Prdm', 'Prdm14_RB_WT']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=Prdm', 'Prdm14_LB_WT']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=LB', 'Prdm14_RB_WT']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=tmp', 'Prdm14_RB_LB-DEL,Prdm14_RB_LB-DEL2']
        # sys.argv = ['./mc4c.py', 'analysis', 'atVpSoi', '--sigma=1.0', '--correction=decay', '--ant-name=tmp', 'Prdm14_Slc_LB-INV']


    args = parser.parse_args(sys.argv[1:])
    args.func(args)


if __name__ == '__main__':
    main()
