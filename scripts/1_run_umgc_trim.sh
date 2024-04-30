#!/bin/bash
set -e
set -u
set -o pipefail


for file in *.fastq.gz
do
	input=$file
	gbstrim.pl --enzyme1 psti --enzyme2 mspi --fastqfile $input --read R1 --outputfile "trimmed_ ${file}" --verbose --threads 24 --minlength 50

done
