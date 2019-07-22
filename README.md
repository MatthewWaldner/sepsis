# SepSIS

## About SepSIS

SepSIS (Separator of Strain Unique Subsequences) is an add-on tool for SPAdes used to extract subsequences unique to bacterial strains from paired, short read, bacterial datasets that have been assembled using SPAdes. SepSIS was designed to function on sequenced bacterial datasets that have originated from non-clonal samples of a single species, as well as clonal strains that are mixed by the user post-sequencing. Additionally, the user may manually mix sets of sequenced reads containing different strains. If these are used as input, SepSIS can function as a strain subsequence difference engine by using the metadata on which reads map to each assembled subsequence. 

The SepSIS pipeline takes as input one pair of .fastq short read files, and will output 3 .fasta files containing strain-specific subsequences. SepSIS attaches subsequences that are not strain-specific to one or both ends of strain-specific subsequence, to allow for the location of the strain-specific subsequences within a genome. Each of the 3 .fasta files contains subsequences that have non-strain specific subsequences on the front end of the sequence, the back end of the sequence, or both ends of the sequence, as designeated by the file name.


## Criteria for Strain-Specifc Subsequences

SepSIS has 3 runmodes that designate the type of analysis performed to designated a subsequence as strain unique. These are designated ORGANIC_P, ORGANIC_Z, and SYNTH. 

The ORGANIC_P, and ORGANIC_Z runmodes use similar criteria to designate a subsequence as strain-specific. In theory, a read set containing only 1 strain with ideal read coverage has the same level of coverage across the entire genome. Therefore, an ideal read set containing 2 or more strains will have lower levels of coverage across strain-specifc subsequences when compared to the rest of the genome. For example: if reads of 2 strains are equally mixed, the subsequences common to both strains will have a coverage of 100% and the subsequences specifc to a strain will have a coverage of 50%. 

Though coverage of short read sequenced strains will generally not be ideal, the above theory has been applied through the use of Z-Scores and percentiles in the ORGANIC_Z and ORGANIC_P runmodes respectively. ORGANIC_Z takes the coverage of all subsequences in a SPAdes assembly_graph.py, and calculates the Z-Score of a single subsequence using the coverage mean and standard deviation. Similarly, ORGANIC_P calculates a subsequence's coverage percentile from the coverage median. The subsequence Z-Score or percentile is then evaulated against a cutoff value. If the coverage is under the thresdhold, the subsequence is evaluated as strain-specific. Ideally, for runmodes ORGANIC_P, and ORGANIC_Z these cutoffs need to be set by the user based on the coverage distribution of the read sets. However, a rough estimation of valid thresholds are 20 to 30 for percentiles and -0.5244 to -0.8416 for Z-Scores. More detailed analysis of ideal thresholds will be uploaded soon.


The SYNTH runmode avoids the coverage-based heuristic of the other two modes, in favour of directly assigning strain-specifc reads to subsequences within an assembly graph. This is performed by tagging the raw (or quality controlled) reads with a strain ID. 

The cutoff for the strain-specifc 


Faults of Synth - relies upon coverage, if only one strain mapping to a common subsequence. .BAM file consequence. no spades read mapping.


## Analysis of Subgraph Types

Additionally, the user must specify the subgraph area to analyze within the assembly_graph.fastg file. These submodes are designated CYCLIC, ISOLATED, and BOTH. The recommended submode is CYCLIC. This mode analyzes only the strongly connected subgraphs of size 2 or more within the assembly graph. These subgraphs contain subsequences that have much higher average coverage than the rest of the assembly graph, and are more likely to contain correctly assembled subsequences. 


The ISOLATED submode runs the analysis on the subgraphs containing strongly connected components of size 1, meaning that none of the subgraphs in this are cyclic. However, there subgraphs will be larger than size 1 to be analyzed. The submode BOTH analyzes the entire graph.

extracted subsequences are output in the form 


## Requirements and Acknowledgements

[Python 3.0](https://www.python.org/downloads/) or greater.

[SPAdes 3.11](https://github.com/ablab/spades) or greater.

[samtools](http://www.htslib.org/)

[minimap2](https://github.com/lh3/minimap2), or another reference assembly program.

SepSIS is loosely based on the SPAdes [Recycler](https://github.com/Shamir-Lab/Recycler) utility. Two of the scripts present in SepSIS contain functions from Recycler. "make_fasta_from_fastg.py" and "recycler_utils.py" contain basic utility functions from Recycler for working with .fastg files. "recycler_utils.py" also contains a function cited by the authors of Recycler as [lh3's fast fastX reader](https://github.com/lh3/readfq/blob/master/readfq.py).


## Installation

Either download from the github broswer or enter "git clone https://github.com/MatthewWaldner/sepsis" in your terminal. All scripts are usable upon download.


## SepSIS Run Guide

### When creating manually mixed short read datasets on your computer to use SepSIS as a strain subsequence difference engine using the SYNTH runmode:

##### 1. Recomended: Trim the reads using a read trimming software for quality control.

##### 2. Run AddSampleNameToReads.py on each .fastq read file to add a sample name to the reads.

This step adds a strain ID to each read within a .fastq file for use within SepSIS.

Ex: "python AddSampleNameToReads.py path_to_read_folder/Strain1_R1.fastq Strain1 path_to_read_folder/Strain1_withName_R1.fastq" Repeat with all other read files. Do not use underscores or spaces in the sample name.

##### 3. Concatenate the separate R1 reads together and the R2 reads together in the terminal.

The read sets are combined to simulate non-clonal strain mixing. Note that if more strains than are added to the mix, the assembly starts to degrade due to the limitations of SPAdes. SepSIS works best with 2 strains, is functional with 3, and becomes more unreliable as more strains are combined.
  
Ex: "cat Strain1_withName_R1.fastq" Strain2_withName_R1.fastq" > Strain1and2_R1.fastq"
  
##### 4. Assemble the bacterial short reads using SPAdes:

Ex: "python spades.py -k 21,33,55,77,99,121 --careful --pe1-1 path_to_read_folder/Strain1and2_R1.fastq --pe1-2 path_to_read_folder/Strain1and2_R2.fastq -o path_to_output_folder/Strain1and2"

##### 5. Set the paths to minimap2 and samtools within CreateBAMFilesForContigs.py.

The script CreateBAMFilesForContigs.py 

##### 6. Run CreateBAMFilesForContigs.py using the assembly graph from SPAdes and the two read files.

Ex: "python CreateBAMFilesForContigs.py Strain1and2 path_to_output_folder/Strain1and2/assembly_graph.fastg path_to_read_folder/Strain1and2_R1.fastq path_to_read_folder/Strain1and2_R2.fastq path_to_output_folder/Strain1and2"

##### 7. Run SepSIS on the assembly_graph.fastg file using RUNMODE SYNTH

Ex: "python SepSIS.py --RUNMODE SYNTH --SUBMODE BOTH --fastgFileIn path_to_output_folder/Strain1and2/assembly_graph.fastg --ScoreValue 20 --outDirectory path_to_output_folder/Strain1and2 --bamFileIn path_to_output_folder/Strain1and2/Strain1and2.reads_minimap2_pe.sort.bam"


### When using SepSIS to extract strain-specifc subsequences from a sequenced short read dataset originating from non-clonal samples using the ORGANIC_P or ORGANIC_Z RUNMODE:

##### 1. Recomended: Trim the reads using a read trimming software for quality control.

##### 2. Assemble the bacterial short reads using SPAdes:
  
Ex: "python spades.py -k 21,33,55,77,99,121 --careful --pe1-1 path_to_read_folder/Sample1_R1.fastq --pe1-2 path_to_read_folder/Sample1_R2.fastq -o path_to_output_folder/Sample1"

##### 3. Run SepSIS on the assembly_graph.fastg file using RUNMODES ORGANIC_Z or ORGANIC_P.
  
Ex: "python SepSIS.py --RUNMODE ORGANIC_P --SUBMODE CYCLIC --fastgFileIn path_to_output_folder/Sample1/assembly_graph.fastg --ScoreValue 20 --outDirectory path_to_output_folder/Sample1"


## SepSIS Arguments

--RUNMODE, -RM : Sets the manner in which a subsequence is evaluated as strain-unique. Options: 'ORGANIC_Z', 'ORGANIC_P', or 'SYNTH'. 'SYNTH' is recommended, but is only availible to read sets mixed by the user.

--SUBMODE, -SM : Sets the graph subset to be evaluated. CYCLIC is recommended in most cases. Options: 'CYCLIC', 'ISOLATED', 'BOTH'

--fastgFileIn, -F : The path to the SPAdes output .fastg file.

--ScoreValue, -S : The thresholding value used to evaluate if a subsequence is strain-specific.

--outDirectory, -O : The path to the output directory for the output .fasta sequence files.

Optional Arguments: SUBMODE

--bamFileIn, -B : The path to the .BAM file used in the SYNTH RUNMODE.

--outSuffix, -S : An optional string to add to the output files.


## Output

SepSIS will output 3 .fasta files per run. The file names are prodiced in the format:

outDirectory + / + RUNMODE + _ + SUBMODE + _ + ENDS_TYPE + _ + outSuffix + .fasta

for example if:

outDirectory = /Users/user/Desktop

RUNMODE = SYNTH

SUBMODE = CYCLIC

outSuffix = Strain1and2

The output files will be named:

/Users/user/Desktop/SYNTH_CYCLIC_BothEnds_Strain1and2.fasta
/Users/user/Desktop/SYNTH_CYCLIC_FrontEnds_Strain1and2.fasta
/Users/user/Desktop/SYNTH_CYCLIC_BackEnds_Strain1and2.fasta

As discussed in the second paragraph of About SepSIS, each of the three files indicates which ends have non-strain-specific subsequences using the XXXXEnds.

The output files contain fasta sequences with header based upon SPAdes the .fastg headers.

The .fastg headers consist of the outgoing edge name from a node, the length of the fastg/fasta sequence contained within that .fastg node, and the coverage of that node. In SepSIS output, these nodes are concatenated using the symbols _._ with an example of 3 sequences concatenated together below. The total sequence length is the sum of the lengths in the headers. The node on the end or ends designated in the file name will not be a strain-specific subsequence. If the below sequence was in a FrontEnds file, the node EDGE_290302_length_4247_cov_59.894571 would be a non-strain-unique subsequence and would consist of the first 4247 characters of the nucleotide sequence.

\>EDGE_290302_length_4247_cov_59.894571_._EDGE_21842_length_122_cov_72.000000_._EDGE_50404_length_140_cov_176.578947_._EDGE_61376_length_174_cov_105.981132

Additionally, if the RUNMODE is SYNTH, the .fasta header will also have the strain of the strain-specific sequence in the header, separated from the node names with ___

\>Strain1___EDGE_290302_length_4247_cov_59.894571_._EDGE_21842_length_122_cov_72.000000_._EDGE_50404_length_140_cov_176.578947_._EDGE_61376_length_174_cov_105.981132
