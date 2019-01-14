import sys



# Read data from file, stuff it into a dict
settings=dict()
with open(sys.argv[1],'r') as iniFile:
    for line in iniFile:
        splitLine = line.split()
        settings[splitLine[0]] = [x for x in splitLine[1:] if x != '']

# Listed integers
for key in ['prm_start','prm_end','vp_start','vp_end','win_start','win_end']:
    settings[key]=[int(x) for x in settings[key]]

# Check lists that should be of equal length
linked=[
    ['pr_seq','prm_start','prm_end'],
    ['re_name','re_seq'],
    ['win_start','win_end'],
    ['vp_name','vp_chr','vp_start','vp_end']
]
for listed in linked:
    assert len(set([len(settings[x]) for x in listed])) == 1, 'Error: different lengths for linked data:'+','.join(str(x) for x in listed)

# Debug it
#for key in settings:
#	print key,'\t',','.join([str(x) for x in settings[key]])
