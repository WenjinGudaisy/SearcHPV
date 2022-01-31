[![Documentation Status](https://readthedocs.org/projects/searchpv/badge/?version=stable)](https://searchpv.readthedocs.io/en/stable/?badge=stable)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/WenjinGudaisy/SearcHPV/blob/main/LICENSE)
[![PyPI version](https://badge.fury.io/py/searcHPV.svg)](https://badge.fury.io/py/searcHPV)
</br> 

|Host | Downloads |
|:----|:---------:|
|PyPI | [![Downloads](https://pepy.tech/badge/searchpv)](https://pepy.tech/project/searchpv)

# SearcHPV
An HPV integration point detection tool for targeted capture sequencing data

## Introdution
* SearcHPV detects HPV fusion sites on both human genome and HPV genome
* SearcHPV is able to provide locally assembled contigs for each integration events. It will report at least one and at most two contigs for each integration sites. The two contigs will provide information captured for left and right sides of the event.

## Getting started
1. Required resources
* Unix like environment


2. Download and install
Firstly, download and install the required resources.
    1) Download Anaconda >=4.11.0: https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent

    2) Download the "environment.yaml" file under this repository

    3) Creat conda environment for SearcHPV:
        ```
        conda env create -f [your_path]/environment.yaml

        ```
        This command will automatically set up all the third-party tools and packages required for SearcHPV and install latest version of SearcHPV. The name of the environment is "searcHPV".

        You can check the packages and tools in this environment by:

        ```
        conda list -n searcHPV

        ```

        You can update the environment by:
        ```
        conda env update -f [your_path]/environment.yaml

        ```



3. Usage

SearcHPV have four main steps. You could either run it start-to-finish or run it step-by-step.

* Before running SearcHPV, active the conda environment:

```
conda activate searcHPV

```

If you are running commands in a bash script, start with:

```
#!/bin/bash
source ~/anaconda3/etc/profile.d/conda.sh;
conda activate searcHPV; 
#[searcHPV commands...]
```
Note: Please check your path of "conda.sh" if you did not install Anaconda in the home directory.

* Usage of searcHPV:

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
-h, --help      show this help message and exit
-window <int>   the length of region searching for informative reads, default=300
-output <str>   output directory, default "./"
-alignment      run the alignment step, step1
-genomeFusion   call the genome fusion points, step2
-assemble local assemble for each integration event, step3
-hpvFusion call the HPV fusion points, step4
-clusterWindow <int> the length of window of clustering integration sites,default=100
-gz             if fastq files are in gz format
-poly(dn) N     poly(n), n*d(A/T/C/G), will report low confidence if contig contains poly(n), default=20
-index          index the original human and virus reference files, default=False
```

Note: If you've already indexed the virus and human reference files for BWA, Samtools, Picard, you do not need to add the "-index" option, especailly when you are running for a batch of samples that share the same virus and human reference files and you do not want to spend time on indexing references every time running a sample. The commands for indexing the virus and human reference files:

```
#activate SearcHPV conda environment first to make sure using the correct versions of tools
ref = '[path_of_your_reference_file]'
bwa index {ref}
samtools faidx {ref}
picard CreateSequenceDictionary R={ref} O={ref.replace('.fa','.dict')
```


4. Examples:

    1) Run it start-to-finish and submit a SBATCH job:
        ```
        #!/bin/bash
        #SBATCH --job-name=searcHPV
        #SBATCH --mail-user=wenjingu@umich.edu
        #SBATCH --mail-type=BEGIN,END
        #SBATCH --cpus-per-task=1
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=8
        #SBATCH --mem=40gb
        #SBATCH --time=100:00:00
        #SBATCH --account=XXXXX
        #SBATCH --partition=standard
        #SBATCH --output=searcHPV.log
        #SBATCH --error=searcHPV.err
        source ~/anaconda3/etc/profile.d/conda.sh;
        conda activate searcHPV;      
        searcHPV -fastq1 Sample_81279.R1.fastq.gz -fastq2 Sample_81279.R2.fastq.gz -humRef hs37d5.fa -virRef HPV.fa -output /home/scratch/HPV_fusion/Sample_81279 -gz -index;
        ```


    2) Run it step-by-step:


        ```
        searchHPV -alignment -fastq1 Sample_81279.R1.fastq.gz -fastq2 Sample_81279.R2.fastq.gz -humRef hs37d5.fa -virRef HPV.fa -output /home/scratch/HPV_fusion/Sample_81279 -gz -index
        searchHPV -genomeFusion -fastq1 Sample_81279.R1.fastq.gz -fastq2 Sample_81279.R2.fastq.gz -humRef hs37d5.fa -virRef HPV.fa -output /home/scratch/HPV_fusion/Sample_81279 -gz
        searchHPV -assemble -fastq1 Sample_81279.R1.fastq.gz -fastq2 Sample_81279.R2.fastq.gz -humRef hs37d5.fa -virRef HPV.fa -output /home/scratch/HPV_fusion/Sample_81279 -gz
        searchHPV -hpvFusion -fastq1 Sample_81279.R1.fastq.gz -fastq2 Sample_81279.R2.fastq.gz -humRef hs37d5.fa -virRef HPV.fa -output /home/scratch/HPV_fusion/Sample_81279 -gz

        ```
        Note: if run it step-by-step, please make sure the output directories for all steps are the same.

## Output
1. Alignment: the marked dupliaction alignment bam file and customized reference genome.\\
2. Genome Fusion Point Calling: orignal callset, filtered callset, filtered clustered callset.\\
3. Assemble: supportive reads, contigs for each integration events (unfiltered).\\
4. HPV fusion Point Calling: alignment bam file for contigs againt human and HPV genome.\\
Final outputs are under the folder "call_fusion_virus": 
summary of all the integration events : "HPVfusionPointContig.txt"
contig sequences for all the integration events: "ContigsSequence.fa"

## Citation
SearcHPV: a novel approach to identify and assemble human papillomavirus-host genomic integration events in cancer --- Accepted by Cancer

## Contact
wenjingu@umich.edu




