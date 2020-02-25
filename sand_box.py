
import numpy as np
import scipy.ndimage as ndimage
from copy import copy

from utilities import get_gauss_kernel
from utilities import OnlineStats
from analysis import estimate_decay_effect

np.set_printoptions(linewidth=230, threshold=300, edgeitems=30, formatter={'float_kind': lambda x: "%8.3f" % x})





