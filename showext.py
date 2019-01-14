import sys
import numpy as np
import collections

import matplotlib.pyplot as pp
import bisect as bs

selection = collections.defaultdict(list)

restout = np.load(sys.argv[1])
byread = restout['byread'].item()
byregion = restout['byregion'].item()

def coordsToRegionIndex(regionDict,chromosome,start,end):
    regionList = regionDict[chromosome]
    # Use bisect implementation to quickly find a matching position
    left = bs.bisect([x[0][1] for x in regionList],start)
    right = bs.bisect([x[0][0] for x in regionList],end)
    print left,right
    refLen = len(regionList)-1

    # Don't bother beyond the last position in the list
    right = min(right, refLen)
    return left, right


def getReadsByRegions(regionDict,regions):
    targetReads = set()
    for regionIndex in regions:
        region = regionDict[regionIndex[0]][regionIndex[1]]
        for read in region[1]:
            targetReads.add(read)
    return targetReads

def getRegionsByReads(readDict,reads):
    selection = collections.defaultdict(list)
    for read in reads:
        for location in readDict[read]:
            if location[0] not in selection:
                selection[location[0]] = collections.defaultdict(list)
            selection[location[0]][location[1]].append(read) #byread[read]
    return selection

def filterByRegions(regionDict,regions):
    mergeSet = set()
    for region in regions:
        # print region
        chrom = region[0]
        location = region[1]
        mergeSet = mergeSet.union(regionDict[chrom][location])
    # for chrom in regionDict:
    # for region in regionDict[chrom]:

    return mergeSet

start,end = coordsToRegionIndex(byregion,'chr7',111007928,111008141)
for index in range(start,end):
    print index,byregion['chr7'][index]
selReads = getReadsByRegions(byregion,[['chr7',index] for index in range(start,end)])
print 'selReads',selReads
touchingRegions = getRegionsByReads(byread,selReads)
print 'touchingRegions',touchingRegions
filteredReads = filterByRegions(touchingRegions,[('chr7',x) for x in [484,3960,4002,4013,4054]])
print 'filteredReads',filteredReads
filtTouchRegions = getRegionsByReads(byread,filteredReads)
print 'filtTouchRegions',filtTouchRegions

for key in byregion:
    for region in byregion[key]:
        if len(region[1])>100:
            print key,region[0]
        # for read in region[1]:
        #	 targetReads.add(read)
        # targetRegions.add((key,region[0]))
exit()

for read in targetReads:
    for location in byread[read]:
        if location[0] not in selection:
            selection[location[0]] = collections.defaultdict(list)
        selection[location[0]][location[1]].append(read) #byread[read]

for reg in targetRegions:
    print reg
print '\n'
# print selection

values = []
labels = []

sortedKeys = selection.keys()
sortedKeys.sort()
for chromosome in sortedKeys:
    sortedRegions = selection[chromosome].keys()
    sortedRegions.sort()
    for region in sortedRegions:
        reads = selection[chromosome][region]
        if len(reads) > 50:
            values.append(len(reads))
            labels.append(chromosome+':'+str(region))#str(byregion[chromosome][region][0][0])+'-'+str(byregion[chromosome][region][0][1]))
            print region,chromosome,byregion[chromosome][region][0],len(reads)#,(chromosome,byregion[chromosome][region][0]) in targetRegions
    #if len(region)>3:
    #	print region,len(selection[region])


pp.bar(range(len(values)),values,label=labels,log=True)
pp.xticks(np.arange(len(values)),labels,rotation=90)

pp.show()
