Sequencing information and filtering log.

This file is to keep track of which genetic samples were removed from downstream analyses and why they were filtered out.

The names used to keep track of samples will be the unique ID give to each sample before sequencing at the UMGC

There were 2 rounds of sequencing runs:
---
First round of sequencing

Random info about the UMGC first round process:
We initially submitted a full plate of samples (96) on May 25, 2022 to UMGC. On June 30, they responded that only 42% of samples had >200ng of DNA, 22% were marginal, (100-200ng), and 35% were failing, (<100ng of DNA). The UMGC said they discovered 'a processing issue in their QC', and that their Qubit and Picogreen numbers didn't match well, so it was kinda weird. Because there was lots of uncertainty in the quantification process, we ended up giving UMGC more than the minimum requirements. On July 20 I submitted either 100pL or the rest of the remaining extraction for any samples that needed to be resubmitted. On Aug 5, they gave their quantification of the resubmission. One sample, 98, was only 10ng/pL, but the rest of the samples were over 100ng/pL, so Cody Hoffman diluted the samples by 5x (to a volume of 140 pL), and took an aliqout for downstream processing. After inquiring on the status, Cody wrote on Oct 5 that there were delays because "We use Picogreen concentrations of the indexed samples for rebalancing, but we had some concerns arise about our instrument’s accuracy and held off processing samples until that could be verified.". We received the sequence data Oct 27, 2022. On Nov 1, we got the invoice and the email included the statement "Due to the timing of the sequencing conversion, we did not perform the Nano QC". We also got the Illumina Data Release on Nov 1.

Notes from release of sequencing data (from Cody):
1) Created 96 GBS libraries using PstI + MspI.
2) Combined all libraries into a single pool and sequenced on a NextSeq P2 1x100-bp run.
3) Generated ≈ 440 M pass filter reads for the pool.
4) All but 5 of the expected barcodes and samples are detected.  Mean yield is ≈ 4.2M reads per library, though there is a lot of variety in read depth.  Most notably, one column of samples were duplicated in the pool as our indexing QC of that column came up with significantly less material than the other samples in the pool, that column was reindexed, and although normalized, the difference between the reaction batches resulted in a greater concentration (~2X) than the other samples.  Further, both replicates of those samples worked well despite the QC results.
5) Mean quality scores ≥Q30 for all libraries.

Data access
You will find the project release directory on the MSI filesystem at: 
/home/sreddy/data_delivery/umgc/2022-q4/221019_VH00601_82_AACCL5HM5/Reddy2_Project_001

Samples from the first round:
Round 1 had 96 samples. The sample ids ranges from 1 to 105, although 9 numbers between those were not represented (no samples for #25-32 and 80) because some swans had multiple DNA extractions done for the same blood card but only a single extraction (the better one) was submitted for sequencing. There were 103 samples that ended up being sequenced because Cody added some duplicate samples (10_2, 11_2, 12_2, 13_2, 14_2, 15_2, 16_2).


Samples filtered, and why:

1) After inspecting FASTQC sequencing info, 5 samples removed:
- 8    (FASTQC, low number of reads, 40,866)
- 16   (FASTQC, low number of reads, 34,104)
- 24     (FASTQC, low number of reads, 44,611)
- 40     (FASTQC, low number of reads, 42,852)
- 98     (FASTQC, low number of reads, 46,848; also possible contamination based on 10% human DNA; and FASTQ quality low)

91 samples passed on through process_radtags and used as input for ref_map

---
Second round of sequencing

Info about the 2nd round of sequencing:
I dropped off the initial submission on Dec 27, 2022. There were 138 samples across 2 plates. The first plate had 80 samples # 106-185 (Iowa, Nebraska, MN and MI feathers, Burke). The second plate had 58 samples # 186-239 (USGS plus 5 redos from the first round of sequencing: 8_2nd, 16_2nd, 24_2nd, 40_2nd, 98_2nd). On Jan 25, 2023, Cody emailed their quantification results, with 125 (91%) passing Picogreen QC with >200ng of DNA, 4 (2.9%) marginal with 100-200ng, and 8 (5.8%) failing with <100ng. Samples # 187, 205, 207, 210, and 215 all failed badly, but samples # 141, 143, and 214 were retained for sequencing. (Looks like 187 was actually left in though). The Illumina data release came on April 12, 2023.

Notes from release of sequencing data (from Cody):
1) Created 138 dual-indexed GBS libraries using enzyme combination PstI + MspI.
2) Combined all libraries into a single pool and sequenced on a 1 lane of a NextSeq P2 1x100-bp run.
3) Generated ≥ 340 M pass filter reads for the pool. 
4) All expected barcodes and samples detected and mean read depth is ≈ 2.5M.  Although, balance between barcodes is rather variable and about a dozen samples have noticeably fewer reads than the target.  Those low samples were low in the sequencing QC and while they are insufficient in the final run, I prefer to see that than to see those samples overrepresented, which has been the issue we’ve been seeing in balancing.
5) Mean quality scores ≥Q30 for all libraries.
6) Configuration of data reads and a script for trimming can be found at the following link: https://bitbucket.org/jgarbe/gbstrim.

Samples filtered, and why:

After inspecting FASTQC sequencing info, 5 samples removed:
- 107 (mean base quality score under 30 (28.8) and low number of reads, 189,950)
- 187 (mean base quality score under 30 (29.8) and low number of reads, 140,161)
- 188 (low library diversity)
- 189 (low library diversity)
- 190 (low library diversity)

133 samples passed on through process_radtags and used as input for ref_map

---

Cumulative tally: We filtered 5 samples from the first round and 5 samples from the second round, and passed 224 samples on to ref_map: 44 from the PCP (18 from Alaska, 26 from Washington), 31 from the RMP (1 Idaho, 13 Alberta, 17 Wyoming), and 148 from the IP (60 Minnesota, 20 Ohio, 7 Wisconsin, 14 Michigan, 30 Iowa, 3 Arkansas-Ontario, 15 Nebraska).

ref_map results:
Read 609478060 BAM records:
  kept 524428702 primary alignments (86.2%), of which 0 reverse reads
  skipped 31305222 primary alignments with insufficient mapping qualities (5.1%)
  skipped 8017807 excessively soft-clipped primary alignments (1.3%)
  skipped 44802960 unmapped reads (7.4%)
  skipped some suboptimal (secondary/supplementary) alignment records

  Per-sample stats (details in 'gstacks.log.distribs'):
    read 2720884.2 records/sample (281670-7855947)
    kept 53.8%-93.1% of these

Built 1089355 loci comprising 524428702 forward reads and 0 matching paired-end reads; mean insert length was 0.0 (sd: nan).

Genotyped 1089355 loci:
  effective per-sample coverage: mean=10.8x, stdev=4.6x, min=2.1x, max=28.8x
  mean number of sites per locus: 90.5
  a consistent phasing was found for 1563849 of out 1732031 (90.3%) diploid loci needing phasing

---



