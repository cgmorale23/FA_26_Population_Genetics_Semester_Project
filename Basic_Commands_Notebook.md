# Copy of Lab notebook for Digital Biology, Spring 2026

## Name
Camille Morales

## Overview
Each new command and bioinformatics tool covered in class should be documented and described here.
This will serve as a useful reference for looking up commands, their sytax, what they are used for, and any other notes.
Plain text formatting is required for the lab notebook.



## shell commands

### Establish Connection to Server

**Description**: ssh
This command is used to establish a secure shell connection to a server, in our case the Big Red 200 HPCC.

**Syntax**: 
You run it like this:
 ```
ssh cgmorale@bigred200.uits.iu.edu
 ```
Password is my ordinary IU password and the characters are invisible as I type them.

### Request a Compute node

**Description**: srun
REQUEST A COMPUTE NODE - ALWAYS DO THIS

type exit to get off of a compute node

**Syntax**: 
 ```
srun -A r02453 --pty bash
 ```
 ```
srun -A c01988 --time=04:00:00 --pty bash
 ```
 ```
srun -A c01988 --mem 8gb --pty bash
 ```
If requesting more memory, --mem 8gb has to be before pty bash

**Description**: srun: the basic command to request an interactive compute node.
-A c01988: tells it to use a particular account ID, c01988 is the account we will use this semester.
--pty bash: tells it we want a "pseudo terminal" that let's us interact with the compute node.

### Navigate to my files
**Syntax**: 
```
cd /N/scratch/cgmorale/Digital_biology/
```
```
cd /N/slate/cgmorale/
```
### Copy to your computer 
**Description**: scp
YOU HAVE TO OPEN UP ANOTHER TERMINAL
NAVIGATE TO BIOINFORMATICS FOLDER WITHIN YOUR OWN COMPUTER (NOT BIGRED)
**Syntax**: 
 ```
scp cgmorale@bigred200.uits.iu.edu:/N/scratch/cgmorale/folder_for_SRR/SRR062634_2_fastqc.html .

scp cgmorale@bigred200.uits.iu.edu:/N/scratch/cgmorale/folder_for_SRR/SRR062634_1_trimmed_fastqc.html .

 ```



## Manipulating Files

### Move a File 
**Description**: moves a file to the current directory
mv can also be used to rename files

```
mv smaller_file bigger_file
```
### Copy a File

**Description**: cp

similar to moving them, except that the original is not removed as part of the operation.
cp will not copy an entire directory and all of its contents unless you add the -r flag to the command to indicate the operation should be recursive.

**Syntax**: 
```
cp <source_path> <destination_path>.
```


### Change Directory 
**Description**: moves into a directory
cd <directory_name> 

### Orient Yourself 
**Description**: pwd 
present working directory, orients you to where you are

### Make a Directory
**Description**: mkdir 
can be used to make a new directory

### List the Current Files
**Description**: ls
lists the files within that directory

### Move Backwards
**Description**: ..
navigates back a directory

### Remove a file
**Description**: rm
remove a file 

### Make a file
**Description**: type nano then the name of the new file

### Remove a directory
**Description**: rmdir
remove a directory (note: the directory must be empty)

### Full File Paths must start with /
**Description**: start of a full file path (it's a root directory). If it doesn't start with a / then it's a relative file path

### Move a folder back to the parent directory
**Description**: mv <filename> ..

### Move back to the root directory
**Description**: moves back to the root directory

**Syntax**: 
```
cd ~
```

### Open a file
**Description**:  nano <filename>
YOU HAVE TO DO THIS FIRST TO LOAD NANO
**Syntax**: 
```
module load nano
```

### Find and Replace within Nano 
**Description**: CTRL + \

### Edit text easier within nano
**Description**: use ctrl shift 6 to highlight text

**Description**: use ctrl k to delete entire lines

### Export file
**Syntax**: 
```
scp cgmorale@bigred200.uits.iu.edu:/N/u/cgmorale/exercise_1_output.txt .
```
**Description**:  on personal computer host terminal - transfers to my local files

### Denotes the current directory
**Description**: a single period .

### Moving files backwards 
**Description**: ../../ 
**Syntax**: 
```
mv cytosine_.txt ../../ (then can hit tab to see what's in there, figure out how many you need to go back)
```
## FLAGS
### Create any intermediate directories along the path if they don't exist already
**Description**: -p
NOTE: /N/scratch/ is for storing large, temporary files.

**Syntax**: 
```
mkdir -p /N/scratch/<username>/Digital_biology/Assignment_3
```

### Check how much storage you're personally allocated on scratch
**Description**: quota

### Check the size of the file 
**Description**: The -o flag shows information associated with individual files, in this case shows the size of the file 

**Syntax**: 
```
ls -o
```

### Shows the size of the file in a more user-friendly way
**Description**: 
In addition to the extra output from -o, this command makes the file sizes easier to interpret (-h), sorts by time stamp (-t), and then reverses the order of the file list (-r).

**Syntax**: 
```
ls -thor
```

### Count Lines of a File 
**Description**: Outputs three values: the total number of lines, the number of words, and the total characters. TIP: press tab to bring up the file name after typing wc

**Syntax**: 
```
wc filename
```

### Check the first lines of a file
**Description**: head

**Syntax**: 
```
head -n <x> <file>
```

### Check the end of a file
**Description**: tail prints the last x lines of a file

**Syntax**: 
```
tail -n 10 *fna
```
### KILL COMMAND 
**Description**: CTRL + c 

### learn the flags associated with wc 
**Description**: 
**Syntax**: 
```
man wc
```

### Write the standard output the number of bytes (characters) in each input file.
**Description**: 
**Syntax**: 
```
wc -c
```

### Searching withing files
**Description**: Grep searches text for lines that match a pattern and prints those lines. Can be used to search for FASTA headers, because those always start with ">"
**Syntax**: 
```
grep "string" filename
```

### Count a specific item 
**Description**: Grep can be used with the pipe command to say "search for this and then do this" For example, below example searches for FASTA headers and then counts them. 
**Syntax**: 
```
grep ">" filename | wc
```

### Exclude items in a search 
**Description**: grep and the -v flag excludes lines containing the specified string
**Syntax**:
```
grep -v "string_to_exclude" file
```

### Identify line numbers associated with a search match
**Description**: grep -n prints the line numbers associated with each match
**Syntax**: 
```
grep -n "string" file
```

### grep -i
**Description**: ignores the case (capital or not) of a search entry


## Conda

wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh  # download

bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/Miniconda  # install

source ~/Miniconda/etc/profile.d/conda.sh  # tell bash where to look for conda

conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main

conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

conda create -n "DigitalBio" --yes  # creates a new conda environment

conda activate DigitalBio  # loads the new conda env

conda install mamba -c conda-forge --yes  # installs mamba, another installer program

rm Miniconda3-latest-Linux-x86_64.sh  # clean up

## Needle

**Description**: Needle is used for alignment of sequences. you list two sequences (below is one for mouse and rat). You can specify the gap and gap extend penalty scores. Then name the outfile to be something relevant to the alignment. 

**Syntax**: 

```
needle -asequence Mus_mus.faa -bsequence Rat_nor.faa -gapopen 10 -gapextend 0.5 -outfile no_gaps_MX1_outfile
```

## Blast

### Initialize BLAST database

**Description**: A required first step for running BLAST efficiently is to "index" the target (longer) sequence or collection of sequences. This has to be done in the same directory where the file to be blasted is. Replace the file name with the correct one. 

**Syntax**: 
```
makeblastdb -in GCF_000001405.40_GRCh38.p14_genomic.fna \
            -dbtype nucl
```

### Run a a BLAST alignment 

**Syntax**:

```
blastn -query Ara_tha_transcripts.fasta\
       -db Ara_tha.fasta \
       -outfmt 7 \
       -max_target_seqs 1 \
       -qcov_hsp_perc 99 \
       -out ara_output.txt
```
**Description**: Where -query is the path to the query sequence (the smaller one).
-db is the prefix for the BLAST database. We created the database in the previous step.
-outfmt (recommended) can accept an integer that specifies the output format. For example,  "7" shows a table of summary statistics with header information, and "6" shows the same table without header information, and the default "0" visualizes each pairwise alignment (some of the other formats may show multiple-sequence alignments which we haven't covered yet).
-out (optional) specifies a path to write the output to.
-help (optional) lists possible flags (similar -help or --help flags exist for many other bioinformatics programs)
less -S can be particularly useful for visualizing the BLAST output because it lets you scroll horizontally.


### Check how many hits correspond to a chromosome
```
grep -v "^#" ara_output.txt | grep -w "NC_003070" | wc -l
```

**Description**:

# FASTQ

## Download an SRA from NCBI
USE PREFETCH
First, download this into my personal Conda for class project
```
mamba install sra-tools fastqc fastp multiqc -c bioconda
```
**Syntax**: 
```
prefetch SRR062634
```
**Description**:Note: "Run" accessions start with "SRR". NCBI also has various non-Run accession types, such as "Experiment" accession numbers beginning with "SRX".)

## Unpack the downloaded SRA
**Syntax**: 
```
fasterq-dump SRR062634 --threads 1
```
**Description**: After downloading we use a second SRA Toolkit functionality to decompress the data. It's important to specify --threads 1 to avoid using too many resources and getting in trouble with the HPC managers. This says to use a single processor instead of multiple. Later in the course we will ask for more resources—and wait our turn in the queue—in order to use additional threads.

This should produce two FASTQs with "1" and "2" in the name, which are mates from paired-end sequencing.

## Evaluate read quality 

**Description**: fastqc

**Syntax**: 
```
fastqc SRR062634_1.fastq 

fastqc SRR062634_2.fastq
```

## Look at a FASTQC output

**Description**: After you run a fastqc, our put wil be in a html format, ex: SRR7167879_fastqc.html you have to scp this onto your computer in order to view

## Trim Sequences

**Description**: fastp program does a lot for us: it (i) trims low-quality bases from the ends of reads, and filters out reads with (ii) low mean quality, (iii) adapter sequence present, or (iv) no mate (for paired-end)

**Syntax**: 
```
fastp --in1 SRR062634_1.fastq --in2 SRR062634_2.fastq --out1 SRR062634_1_trimmed.fastq --out2 SRR062634_2_trimmed.fastq --thread 1
```
OR, two at once
```
fastqc SRR062634_1.fastq SRR062634_2.fastq
```
## Trim the Front of a Sequence

**Description**: use the flag --tri_front1 then the # of bases to trim

```
fastp --in1 $acc.fastq --out1 "$acc"_trimmed.fastq --trim_front1 <num_bases> --thread 1
```
## Summarize read quality

**Description**: multiqc

**Syntax**: puts trimmed FASTQC outputs into a single report. This program reads every FASTQC output in the specified directory, so you will first need to move or copy (using cp command; similar syntax to mv) the trimmed FASTQC zip files (excluding any pre-trimming FASTQCs) into a new directory. 
```
multiqc <Trimmed_fastqc_directory/>
```

## Job Scripts

### Run the job script

**Description**: sbatch

**Syntax**: 
```
sbatch job_template.sh
```

### Check on the status
**Description**: squeue.
R -> means it's running
PD -> means waiting in lin

**Syntax**: 
```
squeue --me
```
### What to Include at the top of every job script

**Syntax**: 
```
#!/bin/bash
#SBATCH --account=c01988
#SBATCH --time=01:00:00

set -euo pipefail

module load conda
source activate /N/project/Db/Week_6/Assignment_6/Conda/envs/Db/
```
**Description**: 

#!/bin/bash - tells it to run the script using the bash shell
#SBATCH --account=c01988 - specifies the account for the class 
#SBATCH --time=01:00:00 - tells it the time to allow for the job. Format HH:MM:SS
set -euo pipefail - safety feature, tells it to
   -e → stop if any command fails
   -u → stop if you use an undefined variable
   -o pipefail → stop if any part of a pipeline fails
Then you have to load conda within the script if using jobs like sam that require it
also have to tell it the source and to activate

## Reads

### Randomly sample a read

**Description**: seqtk sample

**Syntax**: 
```
seqtk sample -s12345 accession.fastq number_of_reads > output.fastq
```
### Align Reads

**Description**: bwa mem
The -t 1 uses a single thread, -o specifies the output file.

**Syntax**: 
```
bwa mem -t 1 $reference "$acc"_1_subset.fastq "$acc"_2_subset.fastq -o "$accession".sam
```
### Index a Reference

**Description**: Indexing reorganizes the reference genome into a compressed, searchable format. After indexing, BWA can find matches in seconds and align millions of reads efficiently.

**Syntax**: 
```
bwa index reference.fasta
```
### Count unmapped reads and read depth

**Description**: use awk

Where XX is the column number and YY is the value. Quotes on YY are used for checking exact matches between strings (not for comparing numbers).

**Syntax**: 
```
awk '$XX == "YY"' <file_name>
```

### Filter reads using a bitmask

**Description**: samtools view

**Syntax**: 
```
 echo "    samtools view" 
    samtools view --exclude-flags <BITMASK> -q 20 -Sh $acc.sam -o $acc.filtered.sam
```

### Examine a .coverage file

**Description**: the columns don't align, this helps make them more neat
**Syntax**: 
```
column -t example_1.coverage | less -S
```
### Extact subsets from a sam file

**Description**: sam is the output of a bwa mem alignment. Can tell it to look at just those chromosomes, then get the mean depth coverage (column 7) and pring those values
**Syntax**: 
```
awk '$1=="NC_000023.11" || $1=="NC_000024.10" {print $1, $7}' example_1.coverage
```

### Explore alignments

**Description**: use tview to open up bam or sam(?) files
**Syntax**: 
```
samtools tview -p <chromosome>:<position> <INPUT_BAM> <REFERENCE>
```
For example, to visualize position 143,264,000 on Chr1 using a particular alignment file do this:

**Syntax**: 
```
bam=/N/scratch/chriscs/Public/Week_6/Assignment_6/Hom_sap/Extra/example.markdup.bam

samtools tview -p NC_000001.11:143264000 $bam $reference
```
note - you had to set bam first to the path to the right bam file

### Filter and process a RAD alignment
```
    #echo "    bwa mem"
    #bwa mem -t <num_threads> $reference "$acc".fastq -o "$acc".sam

    echo "    samtools view" 
    samtools view --exclude-flags <BITMASK> -q 20 -Sh $acc.sam -o $acc.filtered.sam

    echo "    sort"
    samtools sort $acc.filtered.sam -O BAM -o $acc.sorted.bam

    echo "    samtools index"
    time samtools index $acc.sorted.bam

    echo "    calculating coverage"
    samtools coverage $acc.sorted.bam > $acc.coverage
```

# Variant Calling

## BCFtools for Genome-Wide Variant Calling 

**Syntax**:
````
reference=/N/scratch/cgmorale/Digital_biology/Assignment_7/GCF_000001405.40_GRCh38.p14_genomic.fna
bamlist=bamlist.txt
bitmask=2316
outfile=variants.vcf

bcftools mpileup --fasta-ref "$reference" \
                 --bam-list "$bamlist" \
                 --skip-any-set "$bitmask" \
                 --threads 4 \
                 --output-type u \
| bcftools call --multiallelic-caller \
                --variants-only \
                --output-type v \
                --output "$outfile"
````
## Calculate number of variants in a VCF output file
**Syntax**:
````
grep -v '^#' variants.vcf | wc -l
````
## Examine Variants
**Syntax**:
````
less -S filename.vcf
````
## Remove Low Confidence Variants
**Description**: The proposed variants each have some probability of being a false-positive call: observe this in the QUAL field of the VCF.

This is a command to help you filter out sites on Chr2 with "QUAL" below 30. Could change QUAL score.
**Syntax**:
````
bcftools view -i 'QUAL>=30' <input_VCF> -O v -o <output_VCF>
````

## Remove Indels
**Description**: Some of the inferred variants are indels: observe these in the ALT column of the VCF.

In many analyses, indels are excluded. Compared to substitutions, indels are more difficult to model because they frequently involve alternate alleles of varying lengths, their mutational processes are more complex, and they are harder to align and call with high confidence. The presence of an indel also increases the likelihood of local alignment errors or other artifacts, which can reduce overall variant reliability.

It is often cleaner and much more convenient to retain only single nucleotide polymorphisms (SNPs): a SNP is a variant site where a single base differs between two sequences.

**Syntax**:
````
bcftools view -v snps <input_VCF> -O v -o <output_VCF>
````
**Description**: summarize the VCF and confirm no indels remain
````
bcftools stats <input_VCF> | less
````

# Filter sites with multiple alternate alleles
**Syntax**:
````
bcftools view -M2 <input_VCF> -O v -o <output_VCF>
````
## Assess Missingness
**Description**: Individuals missing data for a particular site have "./." for a genotype. Observe this in the first few lines of the VCF

Note: the "0.0" missing filter might be too extreme for some datasets, particularly if you analyze many specimens. So this filter will need to be tuned to each dataset

**Syntax**:
````
bcftools view -i 'F_MISSING==0.0' <input_VCF> -O v -o <output_VCF>
````

## Filter by MAF (Minor Allele Frequency)

**Description**: In some analyses (though not all) it makes sense to remove rare variants. For example, rare variants can cause structure analyses to be noisy and difficult to interpret, because the frequencies of rare alleles are estimated with high variance, and because they aren't informative about gene flow if not many individuals have the allele. Therefore, many studies define a minimum threshold for minor-allele frequency (frequency of the less common allele; MAF) to retain a SNP, e.g., 0.05 Not to mention, singletons—alleles with a single copy in the dataset—are extra likely to originate from sequencing error. 

**Syntax**:
````
bcftools view -i "MAF>=$x" $input -O v -o $output
````

**Calculating MAF**: 
n=6 (in our case, 6 individuals)
Alleles = 12 (if diploid, alleles = nx2)
Allelic frequency Singleton = 1/12 
1/12 =0.083 = MAF
Round up to 0.084

# Filter out Sex Chromosomes
**Syntax**:
````
vcfin=<input_VCF>  # use the MAF-filtered VCF

bgzip -c $vcfin > $vcfin.gz  # compression

bcftools index $vcfin.gz  # index

bcftools view -R regions.txt $vcfin.gz -o <new_VCF>
````
# Visualize inter-SNP correlations

**Description**: In class we talked about genetic linkage, where loci that are physically close together on a chromosome tend to be inherited together, while associations between distant loci get broken up by recombination. Since we are preparing to analyze population structure, it's important to filter out SNPs that are tightly linked. Imagine many SNPs clustered on a chromosome; this region would then carry extra statistical weight, biasing analyses. When applying statistics we want to aim for independent observations, not redundant observations. Through pruning SNPs that are tightly linked, it's possible to achieve approximately-independent SNPs for downstream analysis.

Go ahead ahead and run the below script that calculates the squared correlations (r^2) between pairs of SNPs that are varying genomic distances apart. Run it on the autosome-only VCF (and no singletons).

SCP the output (a PDF) and visualize. Interpret the plot and think about what this means for your SNPs (you can find the background r^2 that is expected even at large distances)

**Syntax**:
````
wget https://raw.githubusercontent.com/the-smith-lab/Digital_biology/refs/heads/main/LD/ld.py

python ld.py <input_VCF>
````
# Convert VCF to be ready for analysis with PLINK
**Description**:  First, we want to rename the the chromosome identifier codes to integers 1 through 22. Copy the below file (or give the path):

**Syntax**:
````
/N/scratch/chriscs/Public/Week_8/VCF/chrom_map.txt
````
rename the chromosomes in your VCF
**Syntax**:
````
bcftools annotate --rename-chrs chrom_map.txt <input_VCF> -o <output_VCF>
````
add names for each individual SNP, i.e., SNP-identifiers:
**Syntax**:
````
bcftools annotate --set-id +'%CHROM:%POS' <input_VCF> -o <output_VCF>
````
finally, convert PLINK—this creates multiple output files, so you only specify the output prefix:
**Syntax**:
````
plink --vcf <input_VCF> --make-bed --out <output_prefix>
````
# Prune tightly linked SNPS
**Description**: 
Based on your observations from Step 4, choose an r^2 cutoff ($r2) that is ~0.1 above the background level, intended to filter out tightly linked SNPs. Then run this command, which outputs two files: (a) a list of SNPs to keep ("<prefix>.prune.in"), and (b) a list of SNPs to exclude ("<prefix>.prune.out").

**Syntax**:
````
plink --bfile <input_prefix> \
      --indep-pairwise 10 1 $r2 \
      --out <output_prefix>
````
**Description**: Extract keeper SNPs to a new file
**Syntax**:
````
plink --vcf <input_VCF> --extract <prune.in_file> --recode vcf-iid --out <output_prefix>
````
## Calculate Divergence with Fst
Make two files: the first is called pop1.txt and contains the FIN identifiers in the VCF (bam files)
**Syntax**:
````
ERR013101.markdup.bam
ERR013116.markdup.bam
ERR015515.markdup.bam
````
and pop2.txt containing the GBD identifiers:
````
ERR016183.markdup.bam
ERR016184.markdup.bam
SRR062634.markdup.bam
````

Calcualate Fst per locus using
**Syntax**:
````
vcftools --vcf <input_VCF> \
         --weir-fst-pop pop1.txt \
         --weir-fst-pop pop2.txt \
         --out fst_results
````
Last, average the outputs to obtain an estimate of genome-wide FST
**Syntax**:
````
tail -n +2 fst_results.weir.fst | awk '{ total += $3 } END { print total/NR }'
````
**Description**: Here, the tail command here takes all lines except the first line (the header); and the awk command takes the average of the third column. 
