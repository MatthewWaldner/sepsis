# SepSIS

## About SepSIS

SepSIS (Separator of Strain Unique Subsequences) is an add-on tool for SPAdes used to extract subsequences unique to bacterial strains from paired, short read, bacterial datasets.

## Requirements and Acknowledgements

Python 3.0 or greater.

SPAdes 3.11 or greater.

samtools

minimap2 , or another reference assembly software.

SepSIS is loosely based on the SPAdes [Recycler](https://github.com/Shamir-Lab/Recycler) utility. Two of the scripts present in SepSIS contain functions from Recycler. "make_fasta_from_fastg.py" and "recycler_utils.py" contain basic utility functions from Recycler for working with .fastg files. "recycler_utils.py" also contains a function cited by the authors of Recycler as [lh3's fast fastX reader](https://github.com/lh3/readfq/blob/master/readfq.py).

## Installation

Either create a clone in the github broswer or enter "git clone https://github.com/MatthewWaldner/sepsis" in your terminal. All scripts are usable upon download.

## Quick Start

Creating mixed short read datasets on your computer to use SepSIS as a strain subsequence difference engine:

Using SepSIS to extract strain unique subsequecnes from a sequenced short read dataset originating from non-clonal samples:

## Data Preprocessing

## Preparing the .fastg File

## Preparing the .BAM File

## Running SepSIS

## Output
