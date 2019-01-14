import sys
import re
import numpy as np
from mc4c_tools import seq_rev_comp


def findSites(refFile, restSeqs, lineLen=50):
    compRestSeqs = [seq_rev_comp(x) for x in restSeqs]
    restSeqs.extend(compRestSeqs)
    reSeqs='|'.join(restSeqs)
    print reSeqs
    restSitesDict = dict()

    with open(refFile,'r') as reference:
        offset = -lineLen
        matches = []
        readSeq = ''
        for line in reference:
            if line[0] == '>':
                matches.extend([[x.start()+offset, x.end()+offset] for x in (re.finditer(reSeqs, readSeq)) if x.start()])
                readSeq = 'N'*lineLen*2
                offset = -lineLen
                matches = []
                restSitesDict[line[1:].rsplit()[0]] = matches
            else:
                readSeq = readSeq[lineLen:]+line.rsplit()[0].upper()
                matches.extend([[x.start()+offset, x.end()+offset] for x in (re.finditer(reSeqs, readSeq)) if x.start() < lineLen])
                offset += lineLen
        matches.extend([[x.start()+offset, x.end()+offset] for x in (re.finditer(reSeqs, readSeq)) if x.start()])

    return restSitesDict

np.savez_compressed(sys.argv[2], restrsites=findSites(sys.argv[1], ['GACC']))
