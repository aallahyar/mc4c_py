import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import spearmanr, pearsonr

from utilities import load_mc4c, limit_to_roi, hasOL, get_nreads_per_bin, load_configs, flatten

# initialization
np.set_printoptions(linewidth=180, threshold=5000)
# cnf_name = 'BMaj-test'
# cnf_name = 'LVR-BMaj-96x'
# cnf_name = '3T3-BMaj-96x,3T3-HS5-96x,BRN-BMaj-96x,BRN-BMaj-96x2,BRN-HS2-96x,BRN-HS2-NP,BRN-HS5-96x,LVR-BMaj-96x,LVR-BMaj-NP,LVR-HS2-96x,LVR-HS2-NP,LVR-HS2-NP2,LVR-HS3-96x,LVR-HS3-NP,LVR-HS5-96x,LVR-HS5-NP'
cnf_name = 'WPL-WTD,WPL-WTD2,WPL-KOD,WPL-KOD2,WPL-WTC,WPL-KOC'
n_bin_lst = np.linspace(100, 1000, 7, dtype=np.int)
n_option = len(n_bin_lst)
cls_name = ['Out left', 'In Left', 'Center', 'In right', 'Out right']
cls_clr = ['#ed0202', '#fe9f9f', '#02b66b', '#62f8fd', '#2202f2']
n_cls = len(cls_name)
show_plot = False
max_n_read = 10000000

# get configs
config_lst = load_configs(cnf_name)
configs = config_lst[0]
roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
run_id = ','.join([config['run_id'] for config in config_lst])

# load dataset
print 'Loading datasets:'
header_lst = ['ReadID', 'Chr', 'ExtStart', 'ExtEnd']
read_all = np.empty([0, 4], dtype=np.int)
for cfg_idx, cfg in enumerate(config_lst):
    if cfg_idx > 0:
        assert cfg['vp_cnum'] == config_lst[0]['vp_cnum']
        assert cfg['roi_start'] == config_lst[0]['roi_start']
        assert cfg['roi_end'] == config_lst[0]['roi_end']

    print '\tLoading {:s}:'.format(cfg['run_id']),
    mc4c_prt = load_mc4c(cfg, unique_only=True, valid_only=True, min_mq=20, reindex_reads=True, verbose=False)
    read_prt = mc4c_prt[header_lst].values
    del mc4c_prt

    # limit fragments to roi
    roi_cen = np.mean([cfg['vp_start'], cfg['vp_end']], dtype=np.int)
    vp_crd = np.array([cfg['vp_cnum'], roi_cen - 5000, roi_cen + 5000])
    read_roi = limit_to_roi(read_prt[:, :4], vp_crd=vp_crd, roi_crd=roi_crd, min_n_frg=2)

    # re-index reads
    read_roi[:, 0] = np.unique(read_roi[:, 0], return_inverse=True)[1] + 1
    n_read_prt = len(np.unique(read_roi[:, 0]))

    # add run specific read ids
    assert n_read_prt < max_n_read
    read_roi[:, 0] = read_roi[:, 0] + (cfg_idx + 1) * max_n_read
    print '[{:,d}] unique reads are added to the collection.'.format(len(np.unique(read_roi[:, 0])))

    # appending to total reads
    read_all = np.vstack([read_all, read_roi.copy()])
    del read_prt
print 'Got [{:,d}] reads in total.'.format(len(np.unique(read_all[:, 0])))

# loop over bin sets
bin_w = np.zeros(n_option, dtype=np.int)
size_score = np.zeros([n_option, 2])
size_n_test = np.zeros(n_option, dtype=np.int)
crr_exp = np.zeros(n_option)
crr_std = np.zeros(n_option)
for oi, n_bin in enumerate(n_bin_lst):

    # get untouched collection of reads
    reads = read_all.copy()
    reads[:, 0] = np.unique(reads[:, 0], return_inverse=True)[1] + 1
    n_read = len(np.unique(reads[:, 0]))

    # define bins
    edge_lst = np.linspace(roi_crd[1], roi_crd[2], num=n_bin + 1, dtype=np.int64).reshape(-1, 1)
    bin_w[oi] = edge_lst[1, 0] - edge_lst[0, 0]
    bin_crd = np.hstack([np.repeat(configs['vp_cnum'], n_bin).reshape(-1, 1), edge_lst[:-1], edge_lst[1:] - 1])
    bin_cen = np.mean(bin_crd[:, 1:], axis=1, dtype=np.int)
    n_roll = int(np.ceil(1000.0 / bin_w[oi])) + 1
    print 'Computing scores for #bin={:d}, bin-w={:0.0f}bp, #nrm={:d} ...'.format(n_bin, bin_w[oi], n_roll)

    # iterate over each bin
    bin_scr = np.full([2, n_bin], fill_value=np.nan)
    prf_roi = np.empty([0, n_bin])
    for bi in range(2, n_bin - 2):

        # select reads
        cls_crd = bin_crd[bi - 2:bi + 3, :]
        read_set = [None] * n_cls
        for ci in range(n_cls):
            is_in = hasOL(cls_crd[ci, :], reads[:, 1:4])
            read_set[ci] = reads[np.isin(reads[:, 0], reads[is_in, 0]), :].copy()

        # compute coverage profile
        cls_prf = np.zeros([n_cls, n_bin], dtype=np.int)
        cls_nrm = np.zeros([n_cls, n_bin], dtype=np.float)
        n_sel = np.zeros(n_cls, dtype=np.int)
        with np.errstate(divide='ignore', invalid='ignore'):
            for ci, read_sel in enumerate(read_set):
                cls_prf[ci, :], n_sel[ci] = get_nreads_per_bin(read_sel[:, :4], bin_crd=bin_crd, min_n_frg=2)
                cls_nrm[ci, :] = pd.Series(cls_prf[ci, :] * 1e2 / n_sel[ci]).rolling(n_roll, center=True).mean().values

        # check if coverage is enough
        if any(n_sel < 100):
            continue
        size_n_test[oi] += 1

        # removing bin selection effect
        cls_feat = cls_nrm.copy()
        is_unk = hasOL([np.min(cls_crd[:, 1:]) - 5000, np.max(cls_crd[:, 1:]) + 5000], bin_crd[:, 1:])
        cls_feat[:, is_unk] = np.nan
        prf_roi = np.vstack([prf_roi, cls_feat[2, :]])
        # cls_feat[np.isnan(cls_feat)] = 0

        # compute spr-correations
        crr_mat = spearmanr(cls_feat.T, nan_policy='omit').correlation
        # crr_mat = pd.DataFrame(cls_feat.T).corr().values
        bin_scr[0, bi] = np.mean(crr_mat[2, [1, 3]])
        bin_scr[1, bi] = np.mean(crr_mat[2, [0, 4]])

        # plotting features
        if show_plot:
            plt.figure(figsize=(20, 4))
            plt_h = [None] * n_cls
            for ci in range(n_cls):
                # plt_h[ci] = plt.plot(bin_cen, cls_nrm[ci, :], color=cls_clr[ci])[0]
                plt_h[ci] = plt.plot(bin_cen, cls_feat[ci, :], color=cls_clr[ci])[0]
            plt.legend(plt_h, ['{:s} (n={:0.0f})'.format(cls_name[i], n_sel[i]) for i in range(n_cls)])
            plt.xlim(roi_crd[1:])
            plt.ylim([0, 10])
            plt.show(block=True)
    if n_bin_lst[oi] != size_n_test[oi]:
        print '\t{:d} bins are ignored due to lack of coverage (i.e. #read < 100).'.format(n_bin_lst[oi] - size_n_test[oi])

    # compute expected correlation
    exp_mat = spearmanr(prf_roi.T, nan_policy='omit').correlation
    for di in range(n_roll):
        np.fill_diagonal(exp_mat[di:, :], np.nan)
        np.fill_diagonal(exp_mat[:, di:], np.nan)
    crr_exp[oi] = np.nanmean(exp_mat)
    crr_std[oi] = np.nanstd(exp_mat)

    # merge scores
    size_score[oi, :] = np.nanmedian(bin_scr.T, axis=0)

# plotting the scores
x_lim = [0.00, 1.00]
plt.figure(figsize=(8, 7))
plt_h = [None] * n_option
for oi in range(n_option):
    plt_h[oi] = plt.plot(size_score[oi, 1], size_score[oi, 0], '*')[0]
    plt.text(size_score[oi, 1], size_score[oi, 0], ' #{:d}'.format(n_bin_lst[oi]),
             horizontalalignment='left', verticalalignment='center')
plt.plot([x_lim[0], x_lim[1]], [x_lim[0], x_lim[1]], ':', color='#8c8c8c', alpha=0.5)
plt.xlabel('Inter SOI correlation')
plt.ylabel('Intra SOI correlation')
plt.xlim(x_lim)
plt.ylim(x_lim)
plt.legend(plt_h, ['#bin={:d}, bin-w={:d}, #test={:d}'.format(n_bin_lst[oi], bin_w[oi], size_n_test[oi])
                   for oi in range(n_option)], loc='lower right')
plt.title('Coverage correlation Inter vs. Intera SOIs\n{:s}'.format(run_id))
plt.savefig('./plots/plt_OptimalBinW_' + run_id + '.pdf', bbox_inches='tight')

