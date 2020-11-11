# SearcHPV
An HPV integration point detection tool for targeted capture sequencing data

## Introdution
* SearcHPV detects HPV fusion sites on both human genome and HPV genome
* SearcHPV is able to provide locally assembled contigs for each integration events. It will report at least one and at most two contigs for each integration sites. The two contigs will provide information captured for left and right sides of the event.

## Getting started
1. Required resources
* Unix like environment
* Third-party tools:
```
Python/3.7.3 https://www.python.org/downloads/release/python-373/
samtools/1.5 https://github.com/samtools/samtools/releases/tag/1.5
BWA/0.7.15-r1140 https://github.com/lh3/bwa/releases/tag/v0.7.15
java/1.8.0_252 https://www.oracle.com/java/technologies/javase/8all-relnotes.html
Picard Tools/2.23.8 https://github.com/broadinstitute/picard/releases/tag/2.23.8
PEAR/0.9.2 https://github.com/tseemann/PEAR
CAP3/02/10/15 http://seq.cs.iastate.edu/cap3.html

```
2. Download and install
Firstly, download and install the required resources.
Then, tap these commands in your terminal:
```
pip install searcHPV

```

3. Usage
SearcHPV have four main steps. You could either run it start-to-finish or run it step-by-step.

* Usage:
```
searcHPV <options> ...
```
* Standard options:
```
 -fastq1 <str>  sequencing data: fastq/fq.gz file
 -fastq2 <str>  sequencing data: fastq/fq.gz file
 -humRef <str>  human reference genome: fasta file
 -virRef <str>  HPV reference genome: fasta file
```
* Optional options:
```
-window <int>   the length of region searching for informative reads, default=300
-output <str>   output directory, default "./"
-align  run the alignment step, step1
-genomeFusion   call the genome fusion points, step2
-assemble local assemble for each integration event, step3
-hpvFusion call the HPV fusion points, step4

```
* Examples:
1) Run it start-to-finish:
```
searcHPV -fastq1 Sample_81279.R1.fastq.gz -fastq2 Sample_81279.R2.fastq.gz -humRef hs37d5.fa -virRef HPV.fa -output /home/scratch/HPV_fusion/Sample_81279

```
2) Run it step-by-step:
```
searchHPV -align -fastq1 Sample_81279.R1.fastq.gz -fastq2 Sample_81279.R2.fastq.gz -humRef hs37d5.fa -virRef HPV.fa -output /home/scratch/HPV_fusion/Sample_81279
searchHPV -genomeFusion -fastq1 Sample_81279.R1.fastq.gz -fastq2 Sample_81279.R2.fastq.gz -humRef hs37d5.fa -virRef HPV.fa -output /home/scratch/HPV_fusion/Sample_81279
searchHPV -assemble -fastq1 Sample_81279.R1.fastq.gz -fastq2 Sample_81279.R2.fastq.gz -humRef hs37d5.fa -virRef HPV.fa -output /home/scratch/HPV_fusion/Sample_81279
searchHPV -hpvFusion -fastq1 Sample_81279.R1.fastq.gz -fastq2 Sample_81279.R2.fastq.gz -humRef hs37d5.fa -virRef HPV.fa -output /home/scratch/HPV_fusion/Sample_81279

```
## Output
1. Alignment: the marked dupliaction alignment bam file and customized reference genome.\\
2. Genome Fusion Point Calling: orignal callset, filtered callset, filtered clustered callset.\\
3. Assemble: supportive reads, contigs for each integration events (unfiltered).\\
4. HPV fusion Point Calling: alignment bam file for contigs againt human and HPV genome.\\
Final output: summary of all the integration events, contig sequences for all the integration events.

## Citation

## Contact
wenjingu@umich.edu




