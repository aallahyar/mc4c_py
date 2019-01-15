import csv
import os.path as path
import numpy as np


def get_chr_info(genome_str, property='chr_name'):
    chr_details = dict({
        'hg19': dict({
            'chr_name': [
                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
                'chrX', 'chrY', 'chrM'
            ],
            'chr_size': [
                249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747,
                135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                155270560, 59373566, 16571
            ]
        }),
        'mm9': dict({
            'chr_name': [
                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                'chrX', 'chrY', 'chrM'
            ],
            'chr_size': [
                197195432, 181748087, 159599783, 155630120, 152537259, 149517037, 152524553, 131738871, 124076172,
                129993255, 121843856, 121257530, 120284312, 125194864, 103494974, 98319150, 95272651, 90772031, 61342430,
                166650296, 15902555, 16299
            ]
        }),
        'mm10': dict({
            'chr_name': [
                'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
                'chrX', 'chrY', 'chrM'
            ],
            'chr_size': [
                195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 124595110,
                130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 94987271, 90702639, 61431566,
                171031299, 91744698, 16299,
            ]
        })
    })
    return chr_details[genome_str][property]


def get_re_info(re_name='DpnII', property='seq', genome_str=None):
    re_details = dict({
        'DpnII': dict({'seq': 'GATC'}),
        'Csp6I': dict({'seq': 'GTAC'}),
        'NlaIII': dict({'seq': 'CATG'}),
        'HindIII': dict({'seq': 'AAGCTT'})
    })

    if property == 'pos':
        re_fname = './renz_files/{:s}_{:s}.npz'.format(genome_str, re_name)
        return np.load(re_fname)['arr_0'][0]
    else:
        return re_details[re_name][property]


# I'm going to assume a safe XML here
def getFastaSequence(genome, chromosome, pos_start, pos_end):
    import urllib2
    from xml.etree import ElementTree

    message = 'http://genome.ucsc.edu/cgi-bin/das/{:s}/dna?segment={:s}:{:d},{:d}'.format(
            genome, chromosome, pos_start, pos_end)
    response_xml = urllib2.urlopen(message)
    html = response_xml.read()
    response_tree = ElementTree.fromstring(html)
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
