# full_genetics_code
All the code used for bioinformatics and subsequent popgen analyses. I ran everything from within a conda environment to make it easier to establish the correct dependencies and version numbers.

Steps:

1) Trim padding and/or cut site residue from beginning of reads, sequencing adaptor from ends of reads, and trim all reads to length. Run '1_run_umgc_trim.sh', which calls 'gbstrim.pl' and uses python 3.7.12 and cutadapt v1.18. 

2) Run '2_process_radtags.txt' in STACKS from the command line (code in process_radtags.txt). This could be done from a bash script but it's just a one-liner, so saved as a text files instead of a shell script.

3) Download Trumpeter Swan reference genome from NCBI and create an index with the Burrows-Wheeler Alignment tool with 'bwa index genome_file_name.fna'.

4) Run '3_bwa_samtools.sh' to loop over all samples, aligning samples to the reference genome index using bwa, then compress the alignment with into bam format with a header present with 'samtools view', and then sort alignments and save to .bam output file with 'samtools sort'.

5) Run '4_refmap.txt' to call genotypes and then create an unfiltered vcf file to use in further analyses.

