# mc4c_python
A python based approach to processing MC4C data


## Requirements

### External Tools:
To run the whole pipeline several tools from third parties are required. The following tools are what we suggest to use, including the version numbers we tested our pipeline with.
- bwa (v0.7.17-r1188)
- bowtie2 (v2.3.4.3)

### External Data:
- A reference genome


## Preparing

### Base calling
Please refer to the [wiki](https://github.com/UMCUGenetics/pymc4c/wiki/Converting-raw-signals-(i.e.-Squiggle)-to-FAST5) to convert raw signals to FAST5 reads.

### Config file
You need to have a config file that contains experiment specific details (e.g. primer sequence). You can find an example of such a file in "config_dir" folder.

### Index reference
Ensure the reference genome, `reference.fa`, is prepared for the mapper you apply later on. In the examples given here, that means it is indexed for BWA.
```
bwa index reference.fa
```

### Create primer fasta (4C Data)
Sometimes molecules attach to eachother, creating a single molecule that originates from multiple, unrelated, circles. To split these, reads are split where primer sequences map on the read. To enable mapping primers to other data, the primers need to be in fasta format, as created in this step.

```
python mc4c.py makeprimerfa \
	settings.ini \
	primers.fa
```

### Find restriction sites
The export function, described later, uses the restriction sites identified in the reference genome to define regions. This enables determining the amount of circles that overlapped such regions.
This step obtains the genomic positions where restriction sites are found and stores them for later use.

```
python mc4c.py refrestr \
	settings.ini \
	reference.fa \
	refstr.npz
```


## Running

### Rename reads / Combine and split data
Renames reads and splits data for HPC processing.

First define the data files from the samples to work with. 

```
export FILE_FASTQ=`ls path/to/*.fq`  
```

Next define where the output files are to be put. The specified path is extended with `_#.block.fa` and `_#.block.fq`, for output fasta and fastq formats respectively, where # is replaced by the index of the datablock.

```
export FILE_OUT="path/to/basename"  
```

The last variable specifies the amount of reads per file. If you desire a single file, fill in a huge number (> number of reads in total). This is what causes the datablock numbering in the file name for the previous variable definition.

```
export LINESPERFILE="20000"  
```

Now the 3 variables used are specified, run the bash script containing the actual awk code:

```
./sge_scripts/splitfq.sh  
```

### Map primers to reads (4C Data)
Sometimes circles appear to attach to eachother, creating a longer read with multiple circles. 
Therefore, reads should be cut at the restriction sites they were most likely cut at originally. 

Map the primers (made in the preparation using makeprimerfa) to the reads.

First index the reads from the samples as if they are the reference data.

```
bowtie2-build \
	-f sample_#.block.fa \
	sample_#.block
```

Next map the primers to this 'reference'.

```
bowtie2 \
	--local \
	-D 20 \
	-R 3 \
	-N 0 \
	-L 15 \
	-i S,1,0.50 \
	--rdg 2,1 \
	--rfg 2,1 \
	--mp 3,2 \
	--ma 2 -a -p 2 -f \
	-x sample_#.block \
	-U primers.fa \
	-S sample_#.block.sam
```

### Split reads by primers (4C Data)
As the primers are mapped to the reads the reads can be cut into what were likely the original circles.

```
python mc4c.py cleavereads \
	settings.ini \
	sample_#.block.sam \
	sample_#.block.fq \
	sample_#.splitpr.fq
```

### Split reads by restriction sites
This step cuts the reads into pieces that together tell what regions were close at the time the experiment was run.

> Note: The data used for input here depends on whether or not reads were previously split by mapped primers.

```
SOURCE=block # Data not split by primers
```
```
SOURCE=splitpr # Data split by primers (4C)
```

Either by defining the variable `SOURCE` as one of the above examples or by replacing the command here, split the reads by restriction sites using the `splitreads` tool:

```
python mc4c.py splitreads \
	settings.ini \
	sample_#.$SOURCE.fq \
	sample_#.splitre.fq
```

### Merge data
In case data was split previously, combine the data into a single gzipped fq file.
> Note: While not necessary if using a single file, you may want to gzip your data anyway.

```
cat *.splitre.fq | gzip > sample.splitre.fq.gz
```

### Map data to reference genome
Now most pieces of sequence that may have been separate before forming a circle together have been split into separate sequences, these sequences can be mapped to the reference genome. While any mapper should work, BWA works well for this purpose.

```
bwa bwasw \
	-b 5 \
	-q 2 \
	-r 1 \
	-z 10 \
	-T 15 \
	-t 12 \
 	-f sample.bwa.sam \
	reference.fa \
	sample.splitre.fq.gz
```

### Export results for plot tools
The mapped data now contains the region any sub-sequence mapped to, while the circles it originated from is described in the read name. Some additional filtering is done and the data is prepared for plotting using the `export` function.

```
python mc4c.py export \
	sample.bwa.sam \
	refstr.npz \
	sample.np
```
