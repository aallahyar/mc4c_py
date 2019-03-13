# MC-4C processing pipeline
A python based approach to process MC4C data. Please refer to the 
[corresponding paper](https://www.nature.com/articles/s41588-018-0161-5) for details on the method.

##Abstract:
Chromatin folding contributes to the regulation of genomic processes such as gene activity. Of note, detailed 
topological studies and genetic evidence have further indicated that individual enhancers can contact and control 
the expression of multiple genes. Conversely, single genes are often influenced by multiple enhancers.
Chromosome conformation capture (3C) methodologies aim to capture these interactions within the nucleus of the cells. 
However, due to size limitation in Illumina sequencing, conventional 3C methods split 3C (i.e. proximity ligation) 
products into smaller DNA concatemers. Therefore, sequencing such a concatemers often result in capturing pairwise 
contacts and captured multi-contact relationships in 3C products are lost. 
In contrast, MC-4C is designed to preserve and collect large concatemers of proximity ligated fragments 
for long molecule sequencing on Oxford Nanopore or Pacific Biosciences platforms, thus allowing study of 
multi-way chromatin interactions. 

The current pipeline delivers a full processing functionality for MC-4C data and includes a tailored 
statistical analysis toolbox to prob the captured interactions. The observed interactions in this statistical
toolbox are compared with a data-intrinsic background model which discerns whether contacts between more than 
two regulatory sequences are mutually exclusive or, conversely, simultaneously happening at a single allele 
resolution.


## Pipeline requirements:
The MC-4C pipeline requires the following tools:
- A Unix like shell (e.g. Bash v3.2+)
- Samtools v1.9+
- Bwa v0.7.17+
- Python v2.7+ and the following Python packages:
    - h5py v2.7.1+
    - numpy v1.13.3+
    - pandas v0.23.4+
    - pysam v0.15.1+
    - matplotlib v2.1.2+ (only for producing summary statistics)


## General remarks:

### Nomenclature: 
We will be using the following conventions. Files and folders are _italicized_ and place holders 
(i.e. variables, arguments) are enclosed by angle brackets (i.e. <config_file> represents a place holder named 
“config_file” that needs to be replaced by appropriate input from the user). Unix shell commands (such as bash, zsh or sh) 
are indicated  with “$” sign and script names (and their command line arguments) are printed in Courier New font 
(e.g. mc4c.py). Taken together, the following command:

```
$ mc4c.py
```

Runs MC-4C pipeline in a shell.


## Configuration files:
Experiment specific parameters for each MC-4C dataset are organized in a “configuration” file. Each configuration file 
is simply a tab-delimited text file with two columns, as follows:

```
param_name1 <tab> param_value1
param_name2 <tab> param_value2
```

Multiple values for a parameter are separated by a “,” sign. For example, the following configuration file specifies 
that the viewpoint lies on chromosome 7, fragments should be mapped to the mm9 mouse reference genome and 
finally GATC and AAGCTT are used as first and second restriction enzymes to prepare the library:

```
vp_chr <tab> chr7
genome_build <tab> mm9
re_seq <tab> GATC,AAGCTT
```

List of parameters that are recognized in the MC-4C pipeline and can be stored in a configuration file are shown 
in **Table.1**.

**Table.1.** Description of parameters that can be defined in a configuration file to be used for processing an MC-4C 
experiment. Required parameters are denoted by (*).

| Name | Value example | Description
| :--- | :--- | :--- |
| genome_build* | hg19 | Reference genome of interest.
| vp_chr* | chr7 | Viewpoint chromosome; the chromosome for which the viewpoint primers are designed.
| prm_start* | 110977147,110977000 | Start coordinate of primers used. Coordinates are separated by “,”.
| prm_end* | 110977171,110977000 | End coordinate of primers used. Coordinates are separated by “,”.
| prm_seq* | GATTTGTGAGCTCAGGGTTTAC,GCAGTAGTGATTCTATTCAATTTTTGGG | Sequence of primers used. Separated by “,”.
| re_name | DpnII,HindIII | Restriction enzyme name used to prepare the MC-4C library.
| re_seq* | GATC,AAGCTT | Restriction enzyme sequence used to prepare the MC-4C library.
| bwa* | /bin/bwa | Path to BWA aligner.
| bwa_index* | /bwa_indices/mm9 | Path to corresponding bwa index of reference genome.
| reference_fasta* | /genomes/mm9.fa | Path to the corresponding reference genome (in fasta format).
| roi_start | 110933500 | Start position of Region Of Interest (ROI). ROI for example will be used to define ‘far-cis’ fragments in PCR duplicate filter. This parameter will be set to 1Mb up stream of the smallest primer coordinate if not given.
| roi_end | 111066500 | End position of Region Of Interest (ROI). This parameter will be set to 1Mb downstream of the largest primer coordinate if not given.

Notes:
If a line in the configuration file starts by “#”, that line will be considered as a comment and ignored.

# Global configuration file:
Parameters that are constant across experiments (e.g. “bwa_path”) can be defined in a “global” configuration (./mc4c.cfg). Once a module is called, the MC-4C pipeline initially loads the parameters defined in this global configuration file, and then proceeds to load parameters in the experiment specific (local) configuration file. Global parameters are ignored when also entered in the local configuration file.
Extras:
If the MC-4C pipeline is used across multiple reference genomes (e.g. mm10 and hg19), the “bwa_index_path” and “ref_genome_file” parameters may include a “%REF%” placeholder, which will be replaced by the appropriate genome (i.e. genome_build) in the (global or local) configuration file. For example, the following row in the configuration file:

```
reference_fasta <tab> ~/genome/%REF%/chrAll.fa
```

will be translated to:

```
reference_fasta <tab> ~/genome/mm9/chrAll.fa
```

if the “genome_build” parameter is set to “mm9”. The user can utilize this functionality to define global configuration paths for running many different MC-4C experiments. An example setting of these global configurations can be found in the ‘./mc4c.cfg’ file.


# Modules:
The entire process of MC-4C pipeline is partitioned into modules. Each module is responsible to perform a specific task such as mapping reads to reference genome. 
Generally, each module in MC-4C receives one or more inputs (e.g. a configuration file and a FASTQ file as input), then performs the corresponding operation (e.g. mapping to reference genome) and finally produces an output file (e.g. a BAM file containing mapped fragments). These modules can be called by their name in MC-4C pipeline. For example:

```
$ mc4c.py mapFragments
```
Calls a module named “mapFragments” in the MC-4C pipeline. In this protocol, we denote the modules in boldface letters. The implemented modules in MC-4C are mentioned in Table.2.

Table.2. Modules defined in MC-4C.

| Module name | Function
| --- | ---
| setReadIds | Defines an identifier for each sequenced read.
| splitReads | Splits each read into fragments according to restriction enzyme recognition sequence.
| mapFragments | Maps the fragments to reference genome.
| makeDataset | Creates a dataset (in HDF5 format) containing mapped fragments.
| removeDuplicates | Removes duplicate reads from an MC-4C dataset.
| QC | Generates various summary report plots for an MC-4C dataset.
| analysis | Performs multi-contact analysis including VP-SOI (a single SOI vs. all SOIs), SOI-SOI (every pair of SOIs) or Cross-ROI (taking every 3-bins as a SOI) analysis.


The corresponding configuration file for an experiment can be provided as follows:

```
$ mc4c.py mapFragments <config_file>
```
 The <config_file> is simply a path to configuration file. For example:
```
$ mc4c.py mapFragments ./expr1.cfg
```
As mentioned before, each module also receives input and output file names. They can be defined as follows:
```
$ mc4c.py mapFragments <config_file> --input_file <input file> --output_file <output file>
```

For example:
```
$ mc4c.py mapFragments ./expr1.cfg --input_file ./inp.fastq.gz --output_file ./out.bam
```

maps the fragments found in “./inp.fastq.gz” file to reference genome defined in “expr1.cfg” configuration file, and then saves the results in “./out.bam” file.
Extras:
The `setReadIds` module supports multiple input files. This is useful if a single library is sequenced multiple times. To this end, separate file names by “,”. E.g. ./inp1.fastq.gz,./inp2.fastq.gz. 

The `splitReads` module supports regular expressions for restriction enzyme recognition sequence (i.e. the “re_seq” parameter in configuration file). This feature is useful if particular restriction enzymes are used to prepare an MC-4C library. For example, if ApoI restriction enzyme is used (which cuts by R^AATTY), the restriction enzyme sequence can be set to [GA]AATT[CT] to properly cut reads.

# Default directory and files:
To further reduce verbosity, MC-4C pipeline supports default paths and file names. These default paths and file names will be used if the corresponding argument is not (fully) given at the time of running a particular module. For example, in the previous example, the user can choose to provide name of the experiment instead of the full path to its configuration file when calling “mapFragments” module. This can be done as follows:
```
$ mc4c.py mapFragments <name>
```
In this case, MC-4C pipeline will look for a configuration file named “cfg_<name>.cfg” in “./configs/” folder. The input file is set by default to “./fragments/frg_<name>.fasta.gz” and the output file is set by default to “./bams/bam_<name>.bam”. List of default paths and files for each module is denoted in Table.3. Accordingly, calling MC-4C pipeline by:
```
$ mc4c.py mapFragments BMaj
```
is equivalent to:
```
$ mc4c.py mapFragments ./configs/cfg_BMaj.cfg --input_file ./fragments/frg_BMaj.fasta.gz --output_file ./bams/bam_BMaj.bam
```
Note: If the given config file name ends with “.cfg”, MC-4C pipeline assumes that the user is referring to a configuration file, not a run name.


Table.3. Default folder and file names for each module.

| Module name | Input folder and file | Output folder and file
| --- | --- | ---
| setReadIds | ./fastqs/raw_<name>.fastq.gz | ./read_files/rd_<name>.fasta.gz
| splitReads | ./reads/rd_<name>.fasta.gz | ./fragments/frg_<name>.fasta.gz
| mapFragments | ./fragments/frg_<name>.fasta.gz | ./bams/bam_<name>.bam
| makeDataset | ./bams/bam_<name>.bam | ./datasets/mc4c_<name>_all.hdf5
| removeDuplicates | ./datasets/mc4c_<name>_all.hdf5 | ./datasets/mc4c_<name>_uniq.hdf5



## Requirements

### External Tools:
To run the whole pipeline several tools from third parties are required. The following tools are what we suggest to use, including the version numbers we tested our pipeline with.
- bwa (v0.7.17-r1188)
- bowtie2 (v2.3.4.3)

### External Data:
- A reference genome


## Preparing

### Base calling (depricated)
The base calling of raw (Squiggle) data is nowadays is automatically done by Nanopore sequencing software.
In any case, if you have received Squiggle data, please refer to the 
[wiki](https://github.com/UMCUGenetics/pymc4c/wiki/Converting-raw-signals-(i.e.-Squiggle)-to-FAST5) to convert the raw 
signals to FAST5 files which contain the reads.

### Config file
Each MC-4C run has certain properties (e.g. view point position, primers sequence used, preferred reference genome, etc.) that need to be 
provided to the pipeline to properly process the data. These (run specific) "configurations" should be placed in a 
".cnf" file. Each row in this file represents a property of the run. The property name is separated from its value by 
a tab (i.e. tab-delimited format). For example, the view point position can be defined in this file in three rows as:

```
vp_chr  7
vp_start    100001000
vp_end  100002000
```


Note: If multiple values need to be given for a property, semicolon (";") can be used between these values. e.g.
```
prm_start value1;value2
```

### Index reference
Ensure the reference genome, `reference.fa`, is prepared for the mapper you apply later on. 
In the examples given here, that means it is indexed for BWA.
```
bwa index reference.fa
```

### Create primer fasta (4C Data)
Sometimes molecules attach to each other, creating a single molecule that originates from multiple, unrelated, circles. 
To split these, reads are split where primer sequences map on the read. 
To enable mapping primers to other data, the primers need to be in fasta format, as created in this step.

```
python mc4c.py makeprimerfa \
	settings.ini \
	primers.fa
```

### Find restriction sites
The export function, described later, uses the restriction sites identified in the reference genome to define regions. 
This enables determining the amount of circles that overlapped such regions.
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
