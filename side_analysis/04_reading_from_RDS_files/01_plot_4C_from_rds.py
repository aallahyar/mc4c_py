
# import rpy2.robjects as robjects
#
# readRDS = robjects.r['readRDS']
# df = readRDS('./CTCF11_10E_WT_129_VP_A.rds')
# print()

import pyreadr

result = pyreadr.read_r('./CTCF11_10E_WT_129_VP_A.rds')

write.table(x, file = "foo.csv", sep = ",", col.names = NA,
        qmethod = "double")