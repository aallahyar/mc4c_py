# MC-4C processing pipeline
A python based approach to process MC4C data. Please refer to the 
[corresponding paper](https://www.nature.com/articles/s41588-018-0161-5) for details on the method.

## Abstract:
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

#### Nomenclature: 
We will be using the following conventions. Files and folders are _italicized_ and place holders 
(i.e. variables, arguments) are enclosed by angle brackets (i.e. <config_file> represents a place holder named 
“config_file” that needs to be replaced by appropriate input from the user). Unix shell commands (such as bash, zsh or sh) 
are indicated  with “$” sign and script names (and their command line arguments) are printed in Courier New font 
(e.g. mc4c.py). Taken together, the following command:

```
$ mc4c.py
```

Runs MC-4C pipeline in a shell.


#### Configuration files:
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
| prm_seq* | GATTT...,GCAGT... | Sequence of primers used. Separated by “,”.
| re_name | DpnII,HindIII | Restriction enzyme name used to prepare the MC-4C library.
| re_seq* | GATC,AAGCTT | Restriction enzyme sequence used to prepare the MC-4C library.
| bwa* | /bin/bwa | Path to BWA aligner.
| bwa_index* | /bwa_indices/mm9 | Path to corresponding bwa index of reference genome.
| reference_fasta* | /genomes/mm9.fa | Path to the corresponding reference genome (in fasta format).
| roi_start | 110933500 | Start position of Region Of Interest (ROI). ROI for example will be used to define ‘far-cis’ fragments in PCR duplicate filter. This parameter will be set to 1Mb up stream of the smallest primer coordinate if not given.
| roi_end | 111066500 | End position of Region Of Interest (ROI). This parameter will be set to 1Mb downstream of the largest primer coordinate if not given.

>Notes:
> - If a line in the configuration file starts by “#”, that line will be considered as a comment and ignored.
> - The “bwa_index” and “reference_fasta” parameters may include a “%REF%” placeholder, 
 which will be replaced by the appropriate genome (i.e. according to genome_build parameter) in the configuration file. 
 This is useful for example if the MC-4C pipeline is used to process experiments using multiple reference genomes 
 (e.g. mm10 and hg19). Due this functionality, assuming that the genome is defined to be mm9 (i.e. genome_build = mm9), 
 the following row in the configuration file:
>    ```
>    reference_fasta <tab> ~/genome/%REF%/chrAll.fa
>    ```
    
>    will be translated to:
    
>    ```
>    reference_fasta <tab> ~/genome/mm9/chrAll.fa
>    ```

#### Global configuration file:
Parameters that are constant across experiments (e.g. “bwa_path”) can be defined in a “global” configuration 
(./mc4c.cfg). Once a module is called, the MC-4C pipeline initially loads the parameters defined in this global 
configuration file, and then proceeds to load parameters in the experiment specific (i.e. local) configuration file. 

> Notes:
> - Global parameters are ignored when also entered in the local configuration file.

#### Parameter definition
The most important parameter in the MC-4C pipeline which requires user attention is ROI coordiantes 
(i.e. `roi_start` and `roi_start`). Selection of a Region of Interest (ROI) requires careful consideration because 
ROI width has an important role in several parts of the MC-4C data processing pipeline. Initially, 
the ROI discerns far-cis fragments that can be used as UMIs for duplication removal procedure. 
Additionally, it identifies “informative” reads (i.e reads with more than one fragments in the ROI) 
that can be used in multi-contact analysis. Finally, the multi-contact analysis is performed by splitting the 
ROI into 200 non-overlapping and equally-spaced bins. 

The user is strongly advised to consider the following suggestions when determining the coordinates for ROI:
 - The viewpoint should always be included in the ROI. This is important as exclusion of viewpoint can change the 
 expected distribution of fragments in the ROI which in turn results in unreliable z-scores calculated in the 
 multi-contact analysis.
 - If an experiment is expected to investigate interactions within a given TAD, the ROI should encompass the 
 TAD boundaries and a few dozen kilobases of flanking sequences, as done in the [MC-4C paper](https://www.nature.com/articles/s41588-018-0161-5). 
 - The ROI size should not be smaller than 120Kb when MC-4C library is prepared using DpnII restriction enzyme 
 as the first cutter. Smaller ROIs will result in a bin width less than 600bp (i.e. 120Kb / 200) 
 which is close to expected minimum resolution of such MC-4C experiment (for human or mouse genomes).
 - The ROI can be larger if interactions between TAD boundaries are investigated as done in the 
 [MC-4C paper](https://www.nature.com/articles/s41588-018-0161-5). However, considering the expected coverage of a 
 (4bp cutter) MC-4C experiment, the ROI should not be larger than 2Mb as sufficient coverage is required for a 
 reliable multi-contact analysis and further interpretation of results (further discussed below).

#### Annotations:
Several quality control or analysis in the MC-4C pipeline use or require annotations. For example, 
in the SOI-SOI analysis, each SOI should have coordinates assigned to it to be used in the association 
test and plotting. For each analysis, only annotations within the ROI are considered. These annotations can be 
stored in a tab-delimited file which is named according to the genome build of interest (e.g. mm9 or h19). 
For example, SOIs in the mm9 genome should be stored in an annotation file named:

```bash
./annotations/ant_mm9.tsv
```

which follows the format below for coordinates:
```text
Hbb-bh1 <tab> chr7 <tab> 110990974
HS1 <tab> chr7 <tab> 111007686
HS2 <tab> chr7 <tab> 111009600
```

## Modules:
The MC-4C pipeline is partitioned into modules. Each module performs a specific task such as mapping fragments 
to the reference genome or removing PCR duplicates from a given MC-4C dataset. List of available modules in the MC-4C pipeline is given in **Table.2**.

**Table.2.** Modules defined in MC-4C.

| Module name | Function
| --- | ---
| **_setReadIds_** | Defines an identifier for each sequenced read.
| **_splitReads_** | Splits each read into fragments according to restriction enzyme recognition sequence.
| **_mapFragments_** | Maps the fragments to the reference genome.
| **_makeDataset_** | Creates a dataset (in HDF5 format) containing mapped fragments.
| **_removeDuplicates_** | Removes duplicate reads from an MC-4C dataset.
| **_QC_** | Generates various summary report plots for an MC-4C dataset.
| **_analysis_** | Performs multi-contact analysis including VP-SOI (a single SOI vs. all SOIs), SOI-SOI (every pair of SOIs) or Cross-ROI (taking every 3-bins as a SOI) analysis.


Generally, each module in the MC-4C pipeline receives one or more inputs (e.g. a configuration file and a FASTQ file), 
then performs the corresponding operation (e.g. mapping to reference genome) and finally produces an output 
file (e.g. a BAM file containing mapped fragments). These modules can be called by their name in MC-4C pipeline. 
For example:

<pre>
$ mc4c.py <b>mapFragments</b>
</pre>

Calls **_mapFragments_** module in the MC-4C pipeline.

The corresponding configuration file for an experiment can be provided as follows:

<pre>
$ mc4c.py <b>mapFragments</b> &lt;config_file&gt;
</pre>

The <config_file> is simply a path to configuration file. For example:
<pre>
$ mc4c.py <b>mapFragments</b> ./expr1.cfg
</pre>

As mentioned before, each module also receives input and output file names. They can be defined as follows:
<pre>
$ mc4c.py <b>mapFragments</b> &lt;config_file&gt; --input_file &lt;input file&gt; --output_file &lt;output file&gt;
</pre>

For example:
<pre>
$ mc4c.py <b>mapFragments</b> ./expr1.cfg --input_file ./inp.fastq.gz --output_file ./out.bam
</pre>
maps the fragments found in “./inp.fastq.gz” file to reference genome defined in “expr1.cfg” configuration file, 
and then saves the results in “./out.bam” file.

### setReadIds module:
Keeping track of reads and their captured fragments is essential for proper data processing and analysis in MC-4C pipeline. 
To facilitate this, once reads are sequenced, we assign unique identifiers to each read (and later to their originated fragments). 
This identifier is preserved across the pipeline and therefore, at each stage of data processing, the user can trace how 
a given read (or fragment) is processed. 

Read identifiers in MC-4C are a `;` separated list of characteristics including:
 - **File index (Fl.Id):** An increasing number starting from 1 which refers to sequencing runs combined. 
 This identifier is useful to track reads across multiple sequencing runs that are being processed together. 
 This will be the case if a template is sequenced multiple times and therefore corresponding reads need to be 
 combined and processed together (e.g. to remove the PCR duplicates). 
 - **Read index (Rd.Id):** An increasing number starting from 1, indicating each sequenced read.
 - **Read length (Rd.Ln):** This part refers to length of the read in terms of base pairs.

For example, a read id can be:
<pre>
Fl.Id:2;Rd.Id:2431;Rd.Ln:4463
</pre>
which refers to a read originated from second sequencing run, and is the 2431st read in the dataset with length of 
4463 base pairs. 

To execute this module, the user is required to provide a configuration file. For example:
<pre>
$ python ./mc4c.py <b>setReadIds</b> BMaj-test
</pre>

Accepted arguments in setReadIds module are as follows:
 - `--input-file`: Input file (in FASTQ format) containing raw sequenced reads (Default: determined based on the config file name).
 - `--output-file`: Output file (in FASTA format) containing sequenced reads with traceable IDs (Default: determined based on the config file name).
 
 > Notes:
> - The **_setReadIds_** module supports multiple input files. This is useful if a single library is sequenced 
 multiple times and the data is needed to be processed at the same time (to remove PCR duplicates for example).
 To achive this, the user can separate file names by “,”. E.g. ./inp1.fastq.gz,./inp2.fastq.gz. 

### splitReads module:
This module is responsible for splitting reads into fragments according to recognition sequence of restriction enzyme 
(define in the given config file) used during the MC-4C library preparation. Similar to reads, fragments produced 
in this step are tagged by a unique identifier to facilitate their tracking across MC-4C pipeline. 

Fragment identifiers are set of `;` separated fields defined as follows:
- **Read ID**: Each fragment carries identifier of the read it is originated from.
- **Fragment index (Fr.Id):** An increasing number starting from 1, indicating each fragment produced during the splitting procedure. Note that this number does not restart for each read.
- **Start Base pair (Fr.SBp):** Specifies the starting base pair number in read where the fragment is originated from
- **Start Base pair (Fr.EBp):** Specifies the end base pair number in read where the fragment is originated from

For example a fragment identifier can be:
<pre>
Fl.Id:2;Rd.Id:2431;Rd.Ln:4463;Fr.Id:15123;Fr.SBp:10;Fr.EBp:30
</pre>
which indicates that this fragment is originated from from 10th to 30th base pairs in read 
`Fl.Id:2;Rd.Id:2431;Rd.Ln:4463`.

To invoke splitting procedure, the user can call the following command in Terminal:
<pre>
$ python ./mc4c.py <b>splitReads</b> BMaj-test
</pre>

Accepted arguments in splitReads module are as follows:
- `--input-file`: Input file (in FASTA format) containing reads with IDs assigned in the setReadIds module (Default: determined based on the config file name).
- `--output-file`: Output file (in FASTA format) containing fragments with traceable IDs (Default: determined based on the config file name).

> Notes:
> - The split procedure is performed according to “exact string search” across each read. 
Due to the relatively high error rate expected in long-read sequencing technologies, splits might be introduced or missed. 
MC-4C pipeline is designed to remedy such a shortcomings as much as possible (by merging closely mapped 
fragments or an using BWA-SW aligner capable of split-mapping). According to our inspection, 
these considerations are sufficient to remedy such a false recognitions to a large extent.
> - The **_splitReads_** module supports regular expressions for restriction enzyme recognition sequence (i.e. the “re_seq” parameter in configuration file). This feature is useful if particular restriction enzymes are used to prepare an MC-4C library. For example, if ApoI restriction enzyme is used (which cuts by R^AATTY), the restriction enzyme sequence can be set to [GA]AATT[CT] to properly cut reads.

### mapFragments module:
This module maps fragments into defined reference genome (defined in the given config file). 
By default, MC-4C pipeline use BWA-SW to map the given fragments. Its mapping arguments are defined as 
follows by default in MC-4C pipeline:
- `-b 5`: Penalty of a Mismatch occurred during alignment of a fragment to the reference genome.
- `-q 2`: Penalty of opening a gap during mapping.
- `-r 1`: Gap extension penalty.
- `-z 5`: A number that is used in “z-best” heuristics implemented in bwa-sw. Briefly, BWA-SW selects the top “z” best candid alignments in the built alignment graph and further extends their alignments to identify the best alignment.
- `-T 15`: Any identified alignment with score smaller than this number is not reported by BWA-SW.

The user can further modify the following arguments in mapFragments module:
- `--input-file`: Input file (in FASTA format) containing fragments with traceable IDs (Default: determined based on the config file name).
- `--output-file`: Output file (in BAM format) containing mapped fragments (Default: determined based on the config file name).
- `--map-argument`: Mapping arguments given to BWA-SW to map fragments (Default=-b 5 -q 2 -r 1 -z 5 -T 15)
- `--n_thread`: Number of threads should be used by the aligner (Default=1).
- `--return_command`: Return only mapping command instead of running it (useful for running the pipeline in a cluster).

### makeDataset module:
This module is responsible for:
- Merging adjacently mapped fragments that are also neighbor in a read
- Marking reads that contain a view point fragment in their middle as invalid. 
- Collecting mapped fragments and storing them in a HDF5 container.

The user can call for this module using:
<pre>
$ python ./mc4c.py <b>makeDataset</b> BMaj-test
</pre>

Accepted arguments in makeDataset module are as follows:
- `--input-file`: Input file (in BAM format) containing fragments with traceable IDs (Default: determined based on the config file name).
- `--output-file`: Output file (in HDF5 format) containing processed fragments (Default: determined based on the config file name).

### removeDuplicates module:
This module is responsible for removing duplicated reads that are produced during the PCR amplification step. 
Briefly, this is done by grouping reads according to fragments mapped in any chromosome other than the viewpoint 
chromosome and removing all but one read in each group. Please refer to MC-4C manuscript for further details on this procedure. 
The user can call for this module using:
<pre>
$ python ./mc4c.py <b>removeDuplicates</b> BMaj-test
</pre>

Accepted arguments in removeDuplicates module are as follows:
- `--input-file`: Input file (in HDF5 format) containing MC4C data (Default: determined based on the config file name).
- `--output-file`: Output file (in HDF5 format) containing MC4C data (Default: determined based on the config file name).
- `--min-mq`: Minimum mapping quality (MQ) to consider a fragment as confidently mapped (Default: 20).

### Default directory and files:
To further reduce verbosity, MC-4C pipeline supports default paths and file names. These default paths and file names 
will be used if the corresponding argument is not (fully) given at the time of running a particular module. Instead,
the user can choose to provide name of the experiment instead of the full path to its configuration file as follows:
<pre>
$ mc4c.py <b>mapFragments</b> &lt;name&gt;
</pre>

Using the command above, MC-4C pipeline will look for a configuration file named “cfg_<name>.cfg” in “./configs/” folder. 
The input file is set by default to “./fragments/frg_<name>.fasta.gz” and the output file is set by default 
to “./bams/bam_<name>.bam”. If the user call the following command:

<pre>
$ mc4c.py <b>mapFragments</b> BMaj
</pre>

is equivalent to:

<pre>
$ mc4c.py <b>mapFragments</b> ./configs/cfg_BMaj.cfg --input_file ./fragments/frg_BMaj.fasta.gz --output_file ./bams/bam_BMaj.bam
</pre> 

List of default paths and input/output files for each module is denoted in **Table.3**.

**Table.3**. Default folder and file names for each module.

| Module name | Input folder and file | Output folder and file
| --- | --- | ---
| **_setReadIds_** | ./fastqs/fq_&lt;name&gt;.fastq.gz | ./read_files/rd_&lt;name&gt;.fasta.gz
| **_splitReads_** | ./reads/rd_&lt;name&gt;.fasta.gz | ./fragments/frg_&lt;name&gt;.fasta.gz
| **_mapFragments_** | ./fragments/frg_&lt;name&gt;.fasta.gz | ./bams/bam_&lt;name&gt;.bam
| **_makeDataset_** | ./bams/bam_&lt;name&gt;.bam | ./datasets/mc4c_&lt;name&gt;_all.hdf5
| **_removeDuplicates_** | ./datasets/mc4c_&lt;name&gt;_all.hdf5 | ./datasets/mc4c_&lt;name&gt;_uniq.hdf5


> Note: 
> - If the given config file name ends with “.cfg”, MC-4C pipeline assumes that the user is referring to a 
configuration file, not a run name.

## Walkthrough of MC-4C pipeline (estimated run time: 15 minutes):
    
   In this section, we provide a step by step description of how the MC-4C pipeline can be used to process sequenced reads in an MC-4C experiment. In this walkthrough, we assume that the user is using Linux operating system. Minor modifications might be required in the commands used to follow this walkthrough using Mac OSX. 
    For demonstration and testing purposes, we prepared a [test MC-4C dataset](http://dx.doi.org/10.17632/pngmgm9yr3.2). This dataset consists of a small sequencing file (i.e. r
    aw_BMaj-test.fastq.gz) and a corresponding configuration file (cfg_BMaj-test.cfg) that holds experiment specific details of this experiment. These two files can be downloaded from [here](http://dx.doi.org/10.17632/pngmgm9yr3.2).

   #### Setting up the pipeline
1) confirm that the required Python packages for MC-4C pipeline are installed by:
    ```
    $ pip install h5py numpy pandas pysam matplotlib
    ```

2) Download the latest version of MC-4C pipeline from its git repository using:
    ```
    $ wget https://github.com/aallahyar/mc4c_py/archive/master.zip
    $ unzip master.zip
    $ cd ./mc4c_py-master
    ```

3) Make an index for the reference genome of interest (mm9 in this example) using bwa by:
    ```
    $ bwa index ~/references/mm9.fa
    ```

   #### Prepare input data:
4) After downloading the aforementioned example data, create a folder named “fastqs” and move the obtained sequencing file (i.e. fq_BMaj-test.fastq.gz) to this folder:
    ```
    $ mkdir -p ./fastqs
    $ mv ~/Downloads/fq_BMaj-test.fastq.gz ./fastqs/
    ```

5) Create a folder named “configs” and move the obtained configuration file (i.e. cfg_BMaj-test.cfg) to this folder:
    ```
    $ mkdir -p ./configs
    $ mv ~/Downloads/cfg_BMaj-test.cfg ./configs/
    ```

6) Define the paths for the reference genome, bwa and the reference index (made in step 3) of the walkthrough) in the corresponding configuration file (i.e. cfg_BMaj-test.cfg). Once done, confirm that these paths are correctly defined by:
    ```
    $ grep ’bwa|bwa_index|reference_fasta’ ./configs/cfg_BMaj-test.cfg
    ```

    For example, in our system these paths are defined as follows:
    ```
    bwa /bin/bwa
    bwa_index /bwa_index/mm9.fa
    reference_fasta /reference_genomes/mm9.fa
    ```

    #### Assign read identifiers to sequenced reads (setReadIds):
    In the MC-4C pipeline, aiming to facilitate tracking of reads and fragments that originate from these reads, a unique identifier is assigned to each sequenced read. This is the first step in processing reads produced in an MC-4C library.

7) Assign a unique identifier to each sequenced read using setReadIds module by:
    <pre>
    $ python ./mc4c.py <b>setReadIds</b> BMaj-test
    </pre>

    #### Split reads into fragments (splitReads)
    Each sequenced read in an MC-4C experiment consists of multiple interacting fragments. To assist the aligner in recognizing and mapping these fragments to the reference genome, we opted to pre-split sequenced reads into fragments based on the restriction enzyme recognition sequence. 

9) Apply the splitReads module by executing the following command:
    <pre>
    $ python ./mc4c.py <b>splitReads</b> BMaj-test
    </pre>

    #### Map fragments to the reference genome (mapFragments)
10) map the produced fragments to the reference genome using the mapFragments module:
    <pre>
    $ python ./mc4c.py <b>mapFragments</b> BMaj-test
    </pre>

    #### Making an MC-4C dataset from mapped fragments (makeDataset):
    Mapped fragments are computationally extended to the nearest restriction site in the reference genome and adjacent fragments within a read that also map adjacently (<20bp) in the genome are merged.

11) apply the **_makeDataset_** module:
    <pre>
    $ python ./mc4c.py <b>makeDataset</b> BMaj-test
    </pre>

    #### Remove read duplicates (removeDuplicates)
    Multi-contact data uniquely allows the removal of PCR duplicates based on the collection of ligated fragments. However due to the (currently) limited sequencing quality of third generation sequencers, fragments that are identified in one PCR duplicate may not be identified in another copy.  Therefore we opted to use trans-ligations (and far-cis fragments) as genomic UMI’s, arguing that the chance of a viewpoint-concatemer randomly ligating to the same trans (or far-cis) fragment multiple times independently is extremely small. 

12. To remove PCR duplicate reads in an MC-4C dataset, the removeDuplicates module can be called using:
    <pre>
    $ python ./mc4c.py <b>removeDuplicates</b> BMaj-test
    </pre>

This is the final pre-processing stage in the MC-4C pipeline. 


## Citation:
Allahyar, A., Vermeulen, C., et al. (2018). Enhancer hubs and loop collisions identified from single-allele topologies. Nature genetics, 50(8), 1151. 

## Contact
For any inquiry please contact Amin Allahyar at a.allahyar@hubrecht.eu