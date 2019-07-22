# SepSIS

## About SepSIS

SepSIS (Separator of Strain Unique Subsequences) is an add-on tool for SPAdes used to extract subsequences unique to bacterial strains from paired, short read, bacterial datasets. SepSIS was designed to function on sequenced bacterial datasets that have originated from non-clonal samples of a single species, as well as clonal strains that are mixed by the user post-sequencing. By mixing sets of sequenced reads, SepSIS can function as a de novo assembled strain difference engine. This avoids any bias that may come from reference assembly. 

The SepSIS pipeline takes in a single set of paired short reads, and will output 3 .fasta files containing strain-specific subsequences. SepSIS attaches subsequences that are not strain-specific to one or both ends of strain-specific subsequence, to allow for the location of the strain-specific subsequences within a genome. Each of the 3 .fasta files contains subsequences that have non-strain specific subsequences on the front end of the sequence, the back end of the sequence, or both ends of the sequence, as designeated by the file name.

## Criteria for Strain-Specifc Subsequence

SepSIS has 3 runmodes that designate the type of analysis performed to designated a subsequence as strain unique. These are designated ORGANIC_P, ORGANIC_Z, and SYNTH. 

The ORGANIC_P, and ORGANIC_Z runmodes use similar criteria to designate a subsequence as strain-specific. In theory, a read set containing only 1 strain with ideal read coverage has the same level of coverage across the entire genome. Therefore, an ideal read set containing 2 or more strains will have lower levels of coverage across strain-specifc sequences than the 

The SYNTH runmode avoids the coverage-based heuristic of the other two modes, in favour of directly assigning strain-specifc reads to subsequences within an assembly graph. This is performed by tagging the raw (or quality controlled) reads with a strain ID 

The cutoff for the strain-specifc 

Ideally, for runmodes ORGANIC_P, and ORGANIC_Z these cutoffs need to be set by the user based on the coverage distribution of the read sets.

## Analysis of Subgraph Types

Additionally, the user must specify the subgraph area to analyze within the assembly_graph.fastg file. These submodes are designated CYCLIC, ISOLATED, and BOTH. The recommended submode is CYCLIC. This mode analyzes only the strongly connected subgraphs of size 2 or more within the assembly graph. These subgraphs contain subsequences that have much higher average coverage than the rest of the assembly graph, and are more likely to contain correctly assembled subsequences. 


The ISOLATED submode runs the analysis on the subgraphs containing SCCS of size 1, meaning that none of the subgraphs in this are cyclic. However, there subgraphs will be larger than size 1 to be analyzed. The submode BOTH analyzes the entire graph.

extracted subsequences are output in the form 

## Requirements and Acknowledgements

[Python 3.0](https://www.python.org/downloads/) or greater.

[SPAdes 3.11](https://github.com/ablab/spades) or greater.

[samtools](http://www.htslib.org/)

[minimap2](https://github.com/lh3/minimap2), or another reference assembly program.

SepSIS is loosely based on the SPAdes [Recycler](https://github.com/Shamir-Lab/Recycler) utility. Two of the scripts present in SepSIS contain functions from Recycler. "make_fasta_from_fastg.py" and "recycler_utils.py" contain basic utility functions from Recycler for working with .fastg files. "recycler_utils.py" also contains a function cited by the authors of Recycler as [lh3's fast fastX reader](https://github.com/lh3/readfq/blob/master/readfq.py).

## Installation

Either download from the github broswer or enter "git clone https://github.com/MatthewWaldner/sepsis" in your terminal. All scripts are usable upon download.

## Quick Start

##### Creating mixed short read datasets on your computer to use SepSIS as a strain subsequence difference engine:

1. Recomended: Trim the reads using a read trimming software for quality control.
2. Run 
3. Concatenate the separate R1 reads together and the R2 reads together in the terminal.
  Ex: cat Strain1_R1.fastq Strain2_R1.fastq > Strain1and2_R1.fastq
4. 

##### Using SepSIS to extract strain-specifc subsequences from a sequenced short read dataset originating from non-clonal samples:

1. Recomended: Trim the reads using a read trimming software for quality control.

2. Assemble the bacterial short reads using SPAdes:
  
  Ex: spades.py -k 21,33,55,77,99,121 --careful --pe1-1 PATH_TO_READ_FOLDER/Sample1_R1.fastq --pe1-2 PATH_TO_READ_FOLDER/Sample1_R2.fastq -o PATH_TO_OUTPUT_FOLDER/Sample1

3. Run SepSIS on the assembly_graph.fastg file using RUNMODES ORGANIC_Z or ORGANIC_P.
  
  Ex: PATH_TO_SepSIS_FOLDER/SepSIS.py --RUNMODE ORGANIC_P --SUBMODE CYCLIC --fastgFileIn PATH_TO_OUTPUT_FOLDER/Sample1/assembly_graph.fastg --ScoreValue 20 --outDirectory PATH_TO_OUTPUT_FOLDER/Sample1

## Data Preprocessing




## Preparing the .fastg File

## Preparing the .BAM File

If a reference assebmler other than minimap2 is used to generate the .BAM file, the user will have to manually use samttols to sort an index the produced .BAM file.

## Running SepSIS


Required Arguments:

--RUNMODE, -RM : Sets the manner in which a subsequence is evaluated as strain-unique. Options: 'ORGANIC_Z', 'ORGANIC_P', or 'SYNTH'. 'SYNTH' is recommended, but is only availible to read sets mixed by the user.

--SUBMODE, -SM : Sets the graph subset to be evaluated. CYCLIC is recommended in most cases. Options: 'CYCLIC', 'ISOLATED', 'BOTH'

--fastgFileIn, -F : The path to the SPAdes output .fastg file.

--ScoreValue, -S : The thresholding value used to evaluate if a subsequence is strain-specific.

--outDirectory, -O : The path to the output directory for the output .fasta sequence files.

Optional Arguments:

--bamFileIn, -B : The path to the .BAM file used in the SYNTH RUNMODE.

--outSuffix, -S : An optional string to add to the output files.


## Output

SepSIS will output 3 .fasta files per run.

