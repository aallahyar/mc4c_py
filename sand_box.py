
import numpy as np
import scipy.ndimage as ndimage
from copy import copy

from utilities import get_gauss_kernel
from utilities import OnlineStats
from analysis import estimate_decay_effect

np.set_printoptions(linewidth=230, threshold=300, edgeitems=30, formatter={'float_kind': lambda x: "%8.3f" % x})



from matplotlib import pyplot as plt

plt.close('all')
plt.figure(figsize=[16, 12])
for img_idx, img in enumerate([obs_smt, bkg_smt,
                               normalize_matrix(obs_smt / bin_npos.reshape(-1, 1) / bin_npos.reshape(1, -1) / bin_npos.reshape(-1, 1) / bin_npos.reshape(1, -1), method='iter', scale=False),
                               normalize_matrix(bkg_smt / bin_npos.reshape(-1, 1) / bin_npos.reshape(1, -1) / bin_npos.reshape(-1, 1) / bin_npos.reshape(1, -1), method='iter', scale=False)]):
    plt.subplot(2, 2, img_idx + 1)
    img_h = plt.imshow(img)
    plt.colorbar(img_h)
    plt.clim(np.nanpercentile(img, [20, 95]))


# obs_org, obs_smt, bkg_rnd, n_pos, decay_prob
plt.close('all')
plt.figure(figsize=[16, 5])
plt.plot(obs_org, label='original', color='gray', alpha=0.5)
plt.plot(obs_smt, label='smoooth', color='green')
# plt.plot(np.arange(n_bin) + 122, decay_prob * n_pos, label='smoooth', color='yellow')
plt.plot(np.mean(bkg_rnd, axis=0), label='bkg_avg', color='red')
plt.ylim([0, n_pos * 0.1])
plt.legend()


plt.close('all')
plt.figure(figsize=[16, 5])
plt.plot(obs_org / bin_nrd, label='original', color='gray', alpha=0.5)
plt.plot(obs_smt / bin_nrd, label='smoooth', color='green')
# plt.plot(np.arange(n_bin) + 122, decay_prob * n_pos, label='smoooth', color='yellow')
plt.plot(np.mean(bkg_rnd, axis=0) / bin_nrd, label='bkg_avg', color='red')
# plt.ylim([0, n_pos * 0.1])
plt.legend()

