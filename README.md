## Table of contents
* General info
* Technologies and requirements
* Setup
* How to use - options
* How to use - examples
* Contact

## General info
This is a guide on how to use splicing_detector. 
splicing_detector is a program which aim is to find potential alretnative splicing sites
based on provide annotation in General Feature Format.

The idea to create such a tool originated in Institute of Evolutionary Biology and Faculty of Biology at the University of Warsaw, Poland. 
This program is a part of Anna Korczyńska's bachelor degree project.


## Technologies and requirements
This project is created in Python 3.10.8 and uses python packages listed below:
splicing_detector requires instalation of a few python packages. These are:
	* pandas
	* gffpandas
	* numpy
	* time
	* argparse
	* intertools
Imports of these libraries are written in the script, but earlier instalation on your device is required for them to work properly.


## Setup
Download splicing_detector.py and make it executable.

On Linux 
chmod +x splicing_detector.py
On Windows using git
git add --chmod=+x splicing_detector.py


## How to use - options
There are obligatory and optional arguments. You can run the program with options:
	* -h or --help -> displays description of other options
	* -a or --annotation -> after this option path to a .gff annotation file is expected. This option is necessary to run the program. If not given a message to provide this argument is displayed.
	* -g or --groups -> after this option path to a .lst file is expected. This option is necessary to run the program. If not given a message to provide this argument is displayed.
	* -m or --minimal_number -> after this option a number of minimal common intron between genes for them to be compared is expected. This argument is optional, default is 1.
	* -t or --program_type -> after this option a string of program type is expected. The only two options are 'identical' or 'overlap'. This argumnt is optional, default is 'identical'.


## How to use - examples
There are a few examples of how to use the program and provide arguments.

Example 1 - want to find alternative splicing on an annotation with short genes or want to detect as one gene repetitive short transcripts.
Run the program on default settings providing only annotation and groups files
python splicing_detector.py -a path/to/gff/file -g path/to/groups/file
(minimal_number is 1 as default and program_type is 'identical')
the same result as
python splicing_detector.py -a path/to/gff/file -g path/to/groups/file -m 1 -t identical

Example 2 - want to detect (as one gene) only these with many common introns and find alternative splicing. Usefull if you want to filter only long genes from tha annotation file.
Run the program with changed minimal number
python splicing_detector.py -a path/to/gff/file -g path/to/groups/file -m 4
(program_type is 'identical' as default)
the same result as
python splicing_detector.py -a path/to/gff/file -g path/to/groups/file -m 4 identical

Example 3 - want to find alernative splicing in the most vide range as possible. The introns don't have to be identical, only overlapping and only one has to be in common os as in short genes with one intron alternative splicing can also be detected.
Run the program in overlap type without changing default minimal number
python splicing_detector.py -a path/to/gff/file -g path/to/groups/file -t overlap
(minimal_number is 1 as default)
the same result as
python splicing_detector.py -a path/to/gff/file -g path/to/groups/file -m 1 -t overlap

Example 4 - want to find alternative splicing also if alternative splicing occurs in every intron, but want to avoid potentially false-positive short genes.
run the program with
python splicing_detector.py -a path/to/gff/file -g path/to/groups/file -m 5 -t overlap
(minimal_number 5 for example, not to big, but big enough to filter short genes)
here there is no option to ommit any arguments and receive the same results


## Contact
If you found an error in the program or in readme file, or have an idea to make the program better, please contact aw.korczynsk@student.uw.edu.pl








