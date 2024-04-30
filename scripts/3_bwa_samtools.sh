#!/bin/bash
set -e
set -u
set -o pipefail


for file in *.fq.gz
do
	name="$(basename $file .fq.gz)"
	bwa mem /home/wolfs064/stacks/trus/genome/GCA_019232035.1/genome.fna $file -t 48 | \  #align with BWA mem
		samtools view -b -h | \                                                           #compress alignment
		samtools sort -o /home/wolfs064/stacks/trus/alignments/first_round/${name}.bam    #sort and save in bam format
done
