import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import spearmanr, pearsonr

from utilities import load_mc4c, load_configs

# initialization
if len(sys.argv) > 1:
    cfg_name = sys.argv[1]
else:
    cfg_name = 'LVR-BMaj-96x'
    # cfg_name = 'FL-PCDHaHS51-PB'
    # cfg_name = 'NPC-BMaj-PB'

# load configs
config_lst = load_configs(cfg_name)

# load MC-HC data
frg_dp = load_mc4c(config_lst, unique_only=True, valid_only=True, min_mq=20, reindex_reads=True, verbose=True)
rd_all = frg_dp[['ReadID', 'Chr', 'ExtStart', 'ExtEnd', 'ReadLength']].values
del frg_dp
n_read = len(np.unique(rd_all[:, 0]))

# calculate read sizes
read_size = np.bincount(rd_all[:, 0], minlength=np.max(rd_all[:, 0]) + 1)[rd_all[:, 0]]
rd_all = np.hstack([rd_all, read_size.reshape(-1, 1)])

# get read info
rd_info = np.unique(rd_all[:, [0, 4, 5]], axis=0)
del rd_all

# plot
plt.figure(figsize=[12, 5])
plt.cla()
nf_lst = range(1, 11)
for nf in nf_lst:
    is_in = rd_info[:, 2] == nf
    plt.boxplot(rd_info[is_in, 1], positions=[nf], widths=0.6)

plt.xlim([0, nf_lst[-1] + 1])
plt.ylim([0, 7000])
plt.xticks(nf_lst, nf_lst)
plt.title(cfg_name + '\n' +
          '#read: {:d}\n'.format(n_read) +
          'Spr-corr: {:0.1f}%'.format(spearmanr(rd_info[:, 1], rd_info[:, 2]).correlation * 100))
plt.savefig('./plots/plt_Reads_SizeVsNFrg_' + cfg_name + '.pdf', bbox_inches='tight')