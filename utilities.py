import csv
import os.path as path

import urllib2
import xml.etree.ElementTree as ET


# I'm going to assume a safe XML here
def getFastaSequence(genome, chromosome, pos_start, pos_end):
    message = 'http://genome.ucsc.edu/cgi-bin/das/{:s}/dna?segment={:s}:{:d},{:d}'.format(
            genome, chromosome, pos_start, pos_end)
    response_xml = urllib2.urlopen(message)
    html = response_xml.read()
    response_tree = ET.fromstring(html)
    return response_tree[0][0].text.replace('\n', '').replace('\r', '')


def rowSetDataType(row, typefunc, indexes):
    for index in indexes:
        if row[index] == '':
            row[index] = None
        else:
            row[index] = typefunc(row[index])


def dataInformation(infile):
    with open(infile,'rU') as tsvFile:
        tsvIn = csv.reader(tsvFile, delimiter='\t')
        header = tsvIn.next()

        itemList = []

        for row in tsvIn:
            if row[0].startswith('#'):
                continue
            assert len(row) == 23, 'Malformed data information file, #columns != 23'
            rowSetDataType(row,int,[3,4,5,11,12,14,15,20,21])
            # for index in [3,4,5,11,12,14,15,20,21]:
            # 	if row[index] == '':
            # 		row[index] = None
            # 	else:
            # 		row[index] = int(row[index])

            itemList.append({
                'id':           row[0],
                'cell_type':    row[1],
                'vp_name':      row[2],
                'vp_start':     row[3],
                'vp_end':       row[4],
                'vp_chr':       row[5], # Is this meant to be a number?
                're1_name':     row[6],
                're1_seq':      row[7],
                're2_name':     row[8],
                're2_seq':      row[9],
                'pr1_seq':      row[10],
                'pr1_start':    row[11],
                'pr1_end':      row[12],
                'pr2_seq':      row[13],
                'pr2_start':    row[14],
                'pr2_end':      row[15],
                'loci_name':    row[16],
                'src_id':       row[17],
                'exp_date':     row[18],
                'genome_build': row[19],
                'win_start':    row[20],
                'win_end':      row[21],
                'seq_plt':      row[22]
            })

    return itemList


# Here we could use pybedtools instead but dependencies thereof make installation difficult
def prepSOfInterest(infiles):
    bedDict = dict()
    for infile in infiles:
        assert infile[-4:] == '.bed'
        thisBed = []
        with open(infile,'rU') as bedFile:
            bedIn = csv.reader(bedFile, delimiter='\t')
            header = bedIn.next()
            for row in bedIn:
                #assert len(row) >= 4
                if len(row) < 3:
                    print '\tBedfile has malformed or empty line:',infile, 'contains:', row
                    continue
                rowSetDataType(row,int,[1,2])
                thisBed.append(row) #[:4]
            # {
            # 	'ss_chr':   row[0],
            # 	'ss_start': int(row[1]),
            # 	'ss_end':   int(row[2]),
            # 	'ss_name':  row[3];
            # }

        bedDict[path.basename(infile[:-4])] = thisBed

    return bedDict


def prepAnnotation(dataInf,antBeds):
    for index, name in enumerate(set([line['loci_name'] for line in dataInf])):
        print index,name
        thisLociData = [line for line in dataInf if line['loci_name'] == name]

        assert len(set([line['genome_build'] for line in thisLociData])) == 1, 'Varying genome builds detected'
        # May need to ask the real meaning of this assertion:
        assert len(set([(line['win_start'],line['win_end']) for line in thisLociData
                        if (line['win_start'],line['win_end']) != (None,None)])) == 1, 'Varying boundaries detected'

        for cellType in [line['cell_type'] for line in thisLociData]:
            for bedKey in antBeds.keys():
                if bedKey.startswith(cellType):
                    curBed = antBeds[bedKey]
                    print cellType,bedKey

# This is where the rest of S03 should be




def prepareMeta(args):
    dataInf = dataInformation(args.infile)
    soi =  prepSOfInterest(args.soibeds)
    ant =  prepSOfInterest(args.antbeds)

    prepAnnotation(dataInf,ant)
