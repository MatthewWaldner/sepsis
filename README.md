# SepSIS

## About SepSIS

SepSIS (Separator of Strain Unique Subsequences) is an add-on tool for SPAdes used to extract subsequences unique to bacterial strains from paired, short read, bacterial datasets.

## Requirements and Acknowledgements

[Python 3.0](https://www.python.org/downloads/) or greater.

[SPAdes 3.11](https://github.com/ablab/spades) or greater.

[samtools](http://www.htslib.org/)

[minimap2](https://github.com/lh3/minimap2), or another reference assembly program.

SepSIS is loosely based on the SPAdes [Recycler](https://github.com/Shamir-Lab/Recycler) utility. Two of the scripts present in SepSIS contain functions from Recycler. "make_fasta_from_fastg.py" and "recycler_utils.py" contain basic utility functions from Recycler for working with .fastg files. "recycler_utils.py" also contains a function cited by the authors of Recycler as [lh3's fast fastX reader](https://github.com/lh3/readfq/blob/master/readfq.py).

## Installation

Either create a clone in the github broswer or enter "git clone https://github.com/MatthewWaldner/sepsis" in your terminal. All scripts are usable upon download.

## Quick Start

Creating mixed short read datasets on your computer to use SepSIS as a strain subsequence difference engine:

Using SepSIS to extract strain unique subsequences from a sequenced short read dataset originating from non-clonal samples:
1. Recomended: Trim the reads using a read trimming software.
2. Assemble the bacterial short reads using SPAdes:
  Ex:
3. Run SepSIS on the assembly_graph.fastg file using RUNMODES ORGANIC_Z or ORGANIC_P.
  Ex:

## Data Preprocessing



## Preparing the .fastg File

## Preparing the .BAM File

If a reference assebmler other than minimap2 is used to generate the .BAM file, the user will have to manually use samttols to sort an index the produced .BAM file.

## Running SepSIS


Required Arguments:
--RUNMODE, -RM : Sets the manner in which a subsequence is evaluated as strain-unique. Options: 'ORGANIC_Z', 'ORGANIC_P', or 'SYNTH'

--SUBMODE, -SM : Sets the graph subset to be evaluated. CYCLIC is recommended in most cases. Options: 'CYCLIC', 'ISOLATED', 'BOTH'

--fastgFileIn, -F : The path to the SPAdes output .fastg file.

--ScoreValue, -S : The thresholding value used to evaluate if a subsequence is strain-specific.

--outDirectory, -O : The path to the output directory for the output .fasta sequence files.

Optional Arguments:
--bamFileIn, -B : The path to the .BAM file used in the SYNTH RUNMODE.

--outSuffix, -S : An optional string to add to the output files.


## Output

SepSIS outputs 3 .fasta files.

