import numpy as np
from scipy.io import loadmat  # this is the SciPy module that loads mat-files
import pandas as pd

def load(infile, dataset=None, idKey='id'):
    mat = loadmat(infile)  # load mat-file
    #print mat.keys()
    if dataset == None:
        for key in mat.keys():
            if key[0] != '_':
                dataset = key
                break
        print 'Assuming you mean to load',dataset,'from',infile
    mdata = mat[dataset]  # variable in mat file
    mdtype = mdata.dtype  # dtypes of structures are "unsized objects"
    ndata = {n: mdata[n][0, 0] for n in mdtype.names}

    # Reconstruct the columns of the data table from just the time series
    # Use the number of intervals to test if a field is a column or metadata
    columns = [n for n, v in ndata.iteritems()]

    # There were supposed to be useful functions instead of this but those 'complained'
    # going the rude but working way now
    dataDict = dict()
    for c in columns:
        curCol = []

        for thing in mdata:
            value = thing[c][0][0]
            if type(value) == np.unicode_:
                value = str(value)
            if c == 'vp_chr':
                value = value[0]
            curCol.append(value)
        dataDict[c] = curCol

    df = pd.DataFrame(dataDict,
                      index=dataDict[idKey],
                      columns=columns)

    return df

def listIndexes(args):
    df = load(args.infile)
    for item in df.index:
        print str(item)


def main():
    import sys, os
    infile = sys.argv[1]
    df = load(infile)

    for index in df.index:
        print index

        thisType = df.loc[index]
        thisDict = thisType.to_dict()

        thisDict['prm_start'] = [thisDict['pr1_pos'][0],thisDict['pr2_pos'][0]]
        thisDict['prm_end'] = [thisDict['pr1_pos'][1],thisDict['pr2_pos'][1]]
        thisDict.pop('pr1_pos')
        thisDict.pop('pr2_pos')

        thisDict['vp_start'] = [thisDict['vp_pos'][0]]
        thisDict['vp_end'] = [thisDict['vp_pos'][1]]
        thisDict.pop('vp_pos')

        thisDict['win_start'] = [thisDict['win_bnd'][0]]
        thisDict['win_end'] = [thisDict['win_bnd'][1]]
        thisDict.pop('win_bnd')

        thisDict['pr_seq'] = [thisDict.pop('pr1_seq'),thisDict.pop('pr2_seq')]

        thisDict['re_seq'] = [thisDict.pop('re1_seq'),thisDict.pop('re2_seq')]
        thisDict['re_name'] = [thisDict.pop('re1_name'),thisDict.pop('re2_name')]

        sortedKeys = thisDict.keys()
        sortedKeys.sort()

        with open(sys.argv[2]+index+'.tsv', 'w') as f:
            delim='\t'
            for key in sortedKeys:
                value=str(thisDict[key])
                if type(thisDict[key]) == list:
                    value=delim.join([str(x) for x in thisDict[key]])
                f.write(str(key) + delim + value + os.linesep)

if __name__ == '__main__':
    main()
