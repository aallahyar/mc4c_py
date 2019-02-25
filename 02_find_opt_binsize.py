import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import spearmanr
from utilities import showprogress

from utilities import load_mc4c, limit_to_roi, get_chr_info, hasOL, get_nreads_per_bin, load_configs, flatten
from analysis import compute_mc_associations

# initialization
cnf_name = 'BMaj-test'
n_bin_lst = np.arange(50, 100, 40, dtype=np.int)
n_param = len(n_bin_lst)
n_perm = 1000
pos_ratio = 0.1
n_ant = 10

# get configs
config_lst = load_configs(cnf_name)
configs = config_lst[0]
roi_crd = np.array([configs['vp_cnum'], configs['roi_start'], configs['roi_end']])
roi_cen = np.mean([configs['vp_start'], configs['vp_end']])
run_id = ','.join([config['run_id'] for config in config_lst])

# get chromosome info
chr_size = get_chr_info(genome_str=configs['genome_build'], property='chr_size')

# load dataset
mc4c_pd = load_mc4c(config_lst, unique_only=True, valid_only=True, min_mq=20, reindex_reads=True, verbose=True)
header_lst = ['ReadID', 'Chr', 'ExtStart', 'ExtEnd']
read_all = mc4c_pd[header_lst].values
del mc4c_pd

# loop over bin sets
fp_freq = np.zeros([n_perm, n_param], dtype=np.int)
for size_idx, n_bin in enumerate(n_bin_lst):

    # define bins and vp area
    edge_lst = np.linspace(roi_crd[1], roi_crd[2], num=n_bin + 1, dtype=np.int64).reshape(-1, 1)
    bin_w = edge_lst[1, 0] - edge_lst[0, 0]
    bin_bnd = np.hstack([edge_lst[:-1], edge_lst[1:] - 1])
    vp_crd = np.array([configs['vp_cnum'], roi_cen - int(bin_w * 1.5), roi_cen + int(bin_w * 1.5)])
    print 'Checking {:d} bins of {:d}bp size'.format(n_bin, bin_w)

    # get informative reads
    has_inf = limit_to_roi(read_all[:, :4], vp_crd=vp_crd, roi_crd=roi_crd, min_n_frg=2)
    reads = read_all[np.isin(read_all[:, 0], has_inf[:, 0]), :].copy()
    reads[:, 0] = np.unique(reads[:, 0], return_inverse=True)[1] + 1
    n_read = len(np.unique(reads[:, 0]))

    # convert fragments to bin-coverage
    cfb_lst = [list() for i in range(n_read + 1)]
    n_frg = reads.shape[0]
    for fi in range(n_frg):
        bin_idx = np.where(hasOL(reads[fi, 2:4], bin_bnd))[0]
        cfb_lst[reads[fi, 0]].append(bin_idx.tolist())

    # filter circles for (>1 bin cvg)
    valid_lst = []
    for rd_nid in range(1, n_read + 1):
        fb_lst = cfb_lst[rd_nid]
        bin_cvg = np.unique(flatten(fb_lst))
        if len(bin_cvg) > 1:
            valid_lst.append(rd_nid)
    reads = reads[np.isin(reads[:, 0], valid_lst), :]

    # downsample and re-index
    # rnd_rid = np.random.choice(np.unique(frg_inf[:, 0]), 8618, replace=False)  ### random selection
    # frg_inf = frg_inf[np.isin(frg_inf[:, 0], rnd_rid), :]
    reads[:, 0] = np.unique(reads[:, 0], return_inverse=True)[1] + 1
    reads_ids = np.unique(reads[:, 0])
    n_read = len(reads_ids)

    for pi in range(n_perm):
        showprogress(pi, n_perm)

        # random selection
        n_rnd = int(n_read * pos_ratio * 2)
        rnd_ids = np.random.choice(reads_ids, n_rnd, replace=False)
        pos_ids = rnd_ids[:n_rnd]
        read_pos = read_all[np.isin(read_all[:, 0], pos_ids), :]

        # define random annotations
        rnd_cen = np.random.randint(roi_crd[1], roi_crd[2], size=[1000, 1])
        ant_crd = np.hstack([np.repeat(configs['vp_cnum'], 1000).reshape(-1, 1),
                             rnd_cen - int(bin_w * 1.5), rnd_cen + int(bin_w * 1.5)])
        is_vp_adj = hasOL(vp_crd, ant_crd, offset=20e3)
        ant_crd = ant_crd[~is_vp_adj, :]
        assert ant_crd.shape > n_ant, 'Can not find non-overlaping SOIs'
        ant_crd = ant_crd[:n_ant]

        # compute scores
        ant_obs, soi_rnd, frg_pos = compute_mc_associations(read_pos, [], ant_crd[:, 1:3], n_perm=10, pos_ids=pos_ids)[:3]
        ant_exp = np.mean(soi_rnd, axis=0)
        ant_std = np.std(soi_rnd, axis=0, ddof=0)
        np.seterr(all='ignore')
        ant_scr = np.divide(ant_obs - ant_exp, ant_std)
        np.seterr(all=None)

        fp_freq[pi, size_idx] = np.sum(np.abs(ant_scr) > 4)

# plotting
plt.figure(figsize=(7, 8))
ax_fpr = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)
plt_h = [None] * 1

# plot correlations
fp_avg = np.mean(fp_freq, axis=0)
fp_std = np.std(fp_freq, axis=0)
x_tick_idx = range(n_param)
plt_h[0] = ax_fpr.plot(x_tick_idx, fp_avg, color='gray')[0]
ax_fpr.fill_between(x_tick_idx, fp_avg - fp_std, fp_avg + fp_std, color='#ebebeb', linewidth=0.2)

ax_fpr.set_xlim([0, n_param])
# ax_fpr.set_ylim([0, n_param])
ax_fpr.set_xlabel('Default ROI')
ax_fpr.set_ylabel('Trans / Adjusted ROI')
ax_fpr.set_title('ROI coverage Spearman correlations')
ax_fpr.legend(plt_h[:3], ['Default vs Trans profile', 'Default vs. Adjusted profile', '5MB'])

x_tick_label = ['{:d}bp'.format(item) for item in n_bin_lst]
plt.xticks(x_tick_idx, x_tick_label, rotation=20)
plt.ylabel('Frequency of False Positives (% of tests)')
# ax_fpr.legend(plt_h, [
#     '>1c reads (n=)'.format(),
#     '5al reads (n=)'.format(),
# ])

plt.title('{:s}\n'.format(run_id))
plt.savefig('./plots/plt_OptimalBinW_' + run_id + '.pdf', bbox_inches='tight')

