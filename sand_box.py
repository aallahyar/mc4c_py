
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