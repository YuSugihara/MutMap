# MutMap User Guide
#### version 2.3.9

## Table of contents
- [What is MutMap?](#what-is-mutmap)
- [Installation](#installation)
  + [Dependencies](#dependencies)
  + [Installation via bioconda](#installation-via-bioconda)
  + [Manual installation](#manual-installation)
- [Usage](#usage)
  + [Example 1 : run MutMap from FASTQ without trimming](#example-1--run-mutmap-from-fastq-without-trimming)
  + [Example 2 : run MutMap from FASTQ with trimming](#example-2--run-mutmap-from-fastq-with-trimming)
  + [Example 3 : run MutMap from BAM](#example-3--run-mutmap-from-bam)
  + [Example 4 : run MutMap from multiple FASTQs and BAMs](#example-4--run-mutmap-from-multiple-fastqs-and-bams)
  + [Example 5 : run MutPlot from VCF](#example-5--run-mutplot-from-vcf)
- [Outputs](#outputs)
- [Additional Resources](#additional-resources)
  + [MutMap commands breakdown](#mutmap-commands-breakdown)
  + [Build a custom SnpEff database](#build-a-custom-snpEff-database)
- [Frequently Asked Questions](#frequently-asked-questions)
  + [Choosing the reference genome](#choosing-the-reference-genome)
  + [When there are many contigs](#when-there-are-many-contigs)
  + [How to filter for higher-confidence SNPs](#how-to-filter-for-higher-confidence-snps)
  + [About multiple testing correction](#about-multiple-testing-correction)

## What is MutMap?
<img src="https://github.com/YuSugihara/MutMap/blob/master/images/1_logo.png" width=200>

Bulked segregant analysis, as implemented in MutMap ([Abe et al., 2012](https://www.nature.com/articles/nbt.2095)), is a powerful and efficient method to identify agronomically important loci in crop plants. MutMap requires whole-genome resequencing of a single individual from the original cultivar and the pooled sequences of F2 progeny from a cross between the original cultivar and mutant. MutMap uses the sequence of the original cultivar to polarize the site frequencies of neighbouring markers and identifies loci with an unexpected site frequency, simulating the genotype of F2 progeny. **The updated pipeline is approximately 5-8 times faster than the previous pipeline, is easier for novice users to use, and can be easily installed through bioconda with all dependencies.**

#### Citation
- Yu Sugihara, Lester Young, Hiroki Yaegashi, Satoshi Natsume, Daniel J. Shea, Hiroki Takagi, Helen Booker, Hideki Innan, Ryohei Terauchi, Akira Abe (2022). [High performance pipeline for MutMap and QTL-seq](https://doi.org/10.7717/peerj.13170). PeerJ, 10:e13170.

- Akira Abe, Shunichi Kosugi, Kentaro Yoshida, Satoshi Natsume, Hiroki Takagi, Hiroyuki Kanzaki, Hideo Matsumura, Kakoto Yoshida, Chikako Mitsuoka, Muluneh Tamiru, Hideki Innan, Liliana Cano, Sophien Kamoun & Ryohei Terauchi (2012). [Genome sequencing reveals agronomically important loci in rice using MutMap](https://www.nature.com/articles/nbt.2095). Nature Biotechnol. 30:174-179.

## Installation
### Dependencies
#### Software
- [BWA](http://bio-bwa.sourceforge.net/)
- [SAMtools](http://samtools.sourceforge.net/)
- [BCFtools](http://samtools.github.io/bcftools/)
- [SnpEff](http://snpeff.sourceforge.net/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

#### Python (>=3.5) libraries
- matplotlib
- numpy
- pandas
- seaborn (optional)

### Installation via bioconda
You can install MutMap via [bioconda](https://bioconda.github.io/index.html).
```
conda create -c bioconda -n mutmap mutmap
conda activate mutmap
```

### Manual installation
If you encounter an error during installation, you can install MutMap manually.
```
git clone https://github.com/YuSugihara/MutMap.git
cd MutMap
pip install -e .
```
You will then need to install other dependencies manually. We highly recommend installing SnpEff and Trimmomatic via bioconda.
```
conda install -c bioconda snpeff
conda install -c bioconda trimmomatic
```
After installation, please check whether SnpEff and Trimmomatic are working by using the commands below.
```
snpEff --help
trimmomatic --help
```

## Usage

If your reference genome contains more than 50 contigs, only the 50 longest contigs will be plotted.

```
mutmap -h

usage: mutmap -r <FASTA> -c <BAM|FASTQ> -b <BAM|FASTQ>
              -n <INT> -o <OUT_DIR> [-T] [-e <DATABASE>]

MutMap version 2.3.9

options:
  -h, --help         show this help message and exit
  -r , --ref         Reference FASTA file.
  -c , --cultivar    FASTQ or BAM file of cultivar. If specifying
                     FASTQ, separate paired-end files with a comma,
                     e.g., -c fastq1,fastq2. This option can be
                     used multiple times.
  -b , --bulk        FASTQ or BAM file of mutant bulk. If specifying
                     FASTQ, separate paired-end files with a comma,
                     e.g., -b fastq1,fastq2. This option can be
                     used multiple times.
  -n , --N-bulk      Number of individuals in the mutant bulk.
  -o , --out         Output directory. The specified directory must not
                     already exist.
  -t , --threads     Number of threads. If a value less than 1 is specified,
                     MutMap will use the maximum available threads. [2]
  -w , --window      Window size in kilobases (kb). [2000]
  -s , --step        Step size in kilobases (kb). [100]
  -D , --max-depth   Maximum depth of variants to be used. This cutoff
                     applies to both the cultivar and the bulk. [250]
  -d , --min-depth   Minimum depth of variants to be used. This cutoff
                     applies to both the cultivar and the bulk. [8]
  -N , --N-rep       Number of replicates for simulations to generate
                     null distribution. [5000]
  -T, --trim         Trim FASTQ files using Trimmomatic.
  -a , --adapter     FASTA file containing adapter sequences. This option
                     is used when "-T" is specified for trimming.
  --trim-params      Parameters for Trimmomatic. Input parameters
                     must be comma-separated in the following order:
                     Phred score, ILLUMINACLIP, LEADING, TRAILING,
                      SLIDINGWINDOW, MINLEN. To remove Illumina adapters,
                     specify the adapter FASTA file with "--adapter".
                     If not specified, adapter trimming will be skipped.
                     [33,<ADAPTER_FASTA>:2:30:10,20,20,4:15,75]
  -e , --snpEff      Predict causal variants using SnpEff. Check
                     available databases in SnpEff.
  --line-colors      Colors for threshold lines in plots. Specify a
                     comma-separated list in the order of SNP-index,
                     p95, and p99. ["#FE5F55,#6FD08C,#E3B505"]
  --dot-color        Color of the dots in plots. ["#11001C"]
  --mem              Maximum memory per thread when sorting BAM files;
                     suffixes K/M/G are recognized. [1G]
  -q , --min-MQ      Minimum mapping quality for mpileup. [40]
  -Q , --min-BQ      Minimum base quality for mpileup. [18]
  -C , --adjust-MQ   Adjust the mapping quality for mpileup. The default
                     setting is optimized for BWA. [50]
  -v, --version      show program's version number and exit
```

MutMap can be run from FASTQ (without or with trimming) and BAM. If you want to run MutMap from VCF, please use MutPlot (example 5). Once you run MutMap, MutMap automatically completes the subprocesses.

+ [Example 1 : run MutMap from FASTQ without trimming](#example-1--run-mutmap-from-fastq-without-trimming)
+ [Example 2 : run MutMap from FASTQ with trimming](#example-2--run-mutmap-from-fastq-with-trimming)
+ [Example 3 : run MutMap from BAM](#example-3--run-mutmap-from-bam)
+ [Example 4 : run MutMap from multiple FASTQs and BAMs](#example-4--run-mutmap-from-multiple-fastqs-and-bams)
+ [Example 5 : run MutPlot from VCF](#example-5--run-mutplot-from-vcf)



### Example 1 : run MutMap from FASTQ without trimming
```
mutmap -r reference.fasta \
       -c cultivar.1.fastq.gz,cultivar.2.fastq.gz \
       -b bulk.1.fastq.gz,bulk.2.fastq.gz \
       -n 20 \
       -o example_dir
```

`-r` : reference fasta

`-c` : FASTQs of cultivar. Please input paired-end reads separated by a comma. FASTQ files can be gzipped.

`-b` : FASTQs of bulk. Please input paired-end reads separated by a comma. FASTQ files can be gzipped.

`-n` : number of individuals in mutant bulk.

`-o` : name of output directory. The specified directory should not already exist.

### Example 2 : run MutMap from FASTQ with trimming
```
mutmap -r reference.fasta \
       -c cultivar.1.fastq.gz,cultivar.2.fastq.gz \
       -b bulk.1.fastq.gz,bulk.2.fastq.gz \
       -n 20 \
       -o example_dir \
       -T \
       -a adapter.fasta
```

`-r` : reference fasta

`-c` : FASTQs of cultivar. Please input paired-end reads separated by a comma. FASTQ files can be gzipped.

`-b` : FASTQs of mutant bulk. Please input paired-end reads separated by a comma. FASTQ files can be gzipped.

`-n` : number of individuals in mutant bulk.

`-o` : name of output directory. The specified directory should not already exist.

`-T` : trim your reads by trimmomatic.

`-a` : FASTA of adapter sequences for trimmomatic.

If you are using TrueSeq3, you can find the adapter sequences in the [Github page of Trimmomatic](https://github.com/usadellab/Trimmomatic/blob/main/adapters/TruSeq3-PE-2.fa). [This thread](https://github.com/usadellab/Trimmomatic/issues/20) is also helpful to preprare the adapter file.

### Example 3 : run MutMap from BAM
```
mutmap -r reference.fasta \
       -c cultivar.bam \
       -b bulk.bam \
       -n 20 \
       -o example_dir
```

`-r` : reference fasta

`-c` : BAM of cultivar.

`-b` : BAM of mutant bulk.

`-n` : number of individuals in mutant bulk.

`-o` : name of output directory. The specified directory should not already exist.

### Example 4 : run MutMap from multiple FASTQs and BAMs
```
mutmap -r reference.fasta \
       -c cultivar_1.1.fastq.gz,cultivar_1.2.fastq.gz \
       -c cultivar_1.bam \
       -b bulk_1.1.fastq.gz,bulk_1.2.fastq.gz \
       -b bulk_2.bam \
       -b bulk_3.bam \
       -n 20 \
       -o example_dir
```

MutMap automatically merges multiple FASTQ and BAM files. Of course, you can merge FASTQ or BAM files using `cat` or `samtools merge` before inputting them into MutMap. If you specify `-c` multiple times, please make sure that those files include only one individual. On the other hand, `-b` can include more than one individual because they are bulked samples. MutMap automatically classifies FASTQ and BAM files based on whether comma exists or not.

### Example 5 : run MutPlot from VCF
```
usage: mutplot -v <VCF> -o <OUT_DIR> -n <INT> [-w <INT>] [-s <INT>]
               [-D <INT>] [-d <INT>] [-N <INT>] [-m <FLOAT>]
               [-S <INT>] [-e <DATABASE>] [--igv] [--indel]

MutPlot version 2.3.9

options:
  -h, --help            show this help message and exit
  -v , --vcf            VCF file which contains cultivar and mutant bulk.
                        in this order. This VCF file must have AD field.
  -o , --out            Output directory. The specified directory can already
                        exist.
  -n , --N-bulk         Number of individuals in the mutant bulk.
  -w , --window         Window size in kilobases (kb). [2000]
  -s , --step           Step size in kilobases (kb). [100]
  -D , --max-depth      Maximum depth of variants to be used. This cutoff
                        applies to both the cultivar and the bulk. [250]
  -d , --min-depth      Minimum depth of variants to be used. This cutoff
                        applies to both the cultivar and the bulk. [8]
  -N , --N-rep          Number of replicates for simulations to generate
                        null distribution. [5000]
  -m , --min-SNPindex   Cutoff of minimum SNP-index for clear results. [0.3]
  -S , --strand-bias    Filter out spurious homozygous genotypes in the cultivar
                        based on strand bias. If ADF (or ADR) is higher than
                        this cutoff when ADR (or ADF) is 0, that SNP will be
                        filtered out. If you want to disable this filtering,
                        set this cutoff to 0. [7]
  -e , --snpEff         Predict causal variants using SnpEff. Check
                        available databases in SnpEff.
  --igv                 Output IGV format file to check results on IGV.
  --indel               Plot SNP-index with INDEL.
  --line-colors         Colors for threshold lines in plots. Specify a
                        comma-separated list in the order of SNP-index,
                        p95, and p99. ["#FE5F55,#6FD08C,#E3B505"]
  --dot-color           Color of the dots in plots. ["#11001C"]
  --fig-width           Width allocated in chromosome figure. [7.5]
  --fig-height          Height allocated in chromosome figure. [4.0]
  --white-space         White space between figures. (This option
                        only affects vertical direction.) [0.6]
  -f , --format         Specify the format of an output image.
                        eps/jpeg/jpg/pdf/pgf/png/rgba/svg/svgz/tif/tiff
  --version             show program's version number and exit
```
MutPlot is included in MutMap. MutMap runs MutPlot after making the VCF. Then, MutPlot will work with default parameters. If you want to change some parameters, you can use VCF inside of `(OUT_DIR/30_vcf/mutmap.vcf.gz)` to retry plotting process like below.

```
mutplot -v OUT_DIR/30_vcf/mutmap.vcf.gz \
        -o ANOTHER_DIR_NAME \
        -n 20 \
        -w 2000 \
        -s 100
```

#### Use MutPlot for a VCF which was made by yourself
In this case:
1. Ensure that your VCF includes the AD format.
2. Ensure that your VCF includes two columns of cultivar and mutant bulk in this order.

If you encounter an error, please try running MutMap from FASTQ or BAM before reporting it in the issues.

## Outputs
Inside of `OUT_DIR` is like below.
```
|-- 10_ref
|   |-- reference.fasta
|   |-- reference.fasta.amb
|   |-- reference.fasta.ann
|   |-- reference.fasta.bwt
|   |-- reference.fasta.fai
|   |-- reference.fasta.pac
|   `-- reference.fasta.sa
|-- 20_bam
|   |-- bulk.bam
|   |-- bulk.bam.bai
|   |-- cultivar.bam
|   `-- cultivar.bam.bai
|-- 30_vcf
|   |-- mutmap.vcf.gz
|   `-- mutmap.vcf.gz.tbi
|-- 40_mutmap
|   |-- snp_index.tsv
│   ├── snp_index.p95.tsv
│   ├── snp_index.p99.tsv
|   |-- sliding_window.tsv
│   ├── sliding_window.p95.tsv
│   ├── sliding_window.p99.tsv
|   `-- mutmap_plot.png
`-- log
    |-- alignment.log
    |-- bcftools.log
    |-- bwa.log
    |-- samtools.log
    `-- tabix.log
```
- If you run MutMap with trimming, you will get the directory of `00_fastq` which includes FASTQs after trimming.
- You can check the results in `40_mutmap`.
  + `snp_index.tsv` : columns in this order.
    - **CHROM** : chromosome name
    - **POSI** : position in chromosome
    - **VARIANT** : SNP or INDEL
    - **DEPTH** : depth of bulk
    - **p99** : 99% confidence interval of simulated SNP-index
    - **p95** : 95% confidence interval of simulated SNP-index
    - **SNP-index** : real SNP-index
  + `sliding_window.tsv` : columns in this order.
    - **CHROM** : chromosome name
    - **POSI** : central position of window
    - **MEAN p99** : mean of p99
    - **MEAN p95** : mean of p95
    - **MEAN SNP-index** : mean SNP-index
  + `snp_index.p95.tsv` and `snp_index.p99.tsv` contain only the SNPs that exceed the respective thresholds (95% or 99%). Similarly, `sliding_window.p95.tsv` and `sliding_window.p99.tsv` contain only the windows that exceed the respective thresholds.
  + `mutmap_plot.png` : resulting plot (like below)
    - **<span style="color: blue; ">BLUE dot</span>** : variant
    - **<span style="color: red; ">RED line</span>** : mean SNP-index
    - **<span style="color: orange; ">ORANGE line</span>** : mean p99
    - **<span style="color: green; ">GREEN line</span>** : mean p95
  + If you run MutMap with SnpEff, the following additional outputs will be generated:
    - **mutmap.snpEff.vcf**: The updated VCF file after annotation by SnpEff, located in the `40_mutmap` directory.
    - **snp_index.tsv**: This file will contain a new column, **impact**, which includes the mutation impact information predicted by SnpEff.
    - When plotting the results, variants classified as **MODERATE** by SnpEff are marked with a `+` symbol, while variants classified as **HIGH** are marked with an `x` symbol in the plot.



<img src="https://github.com/YuSugihara/MutMap/blob/master/images/2_result.png" width=600>

## Additional Resources

### MutMap commands breakdown

For a detailed breakdown of the commands used in MutMap, including explanations of each step, parameters, and best practices for troubleshooting, please refer to the [MutMap Commands Breakdown](docs/01_MutMap_Commands_Breakdown.md) document.

### Build a custom SnpEff database
If you are working with a non-model organism or your own reference genome, you may need to build a custom SnpEff database. For detailed instructions on how to build a custom SnpEff database, please refer to the [Build a Custom SnpEff Database](docs/02_Guide_to_Build_SnpEff_DB.md) document.

## Frequently Asked Questions

### Choosing the reference genome  
You can use a line that was not involved in the cross as the reference genome. In the version of MutMap published by [Abe et al., 2012](https://www.nature.com/articles/nbt.2095), the reference fasta was rebuilt using the wild-type cultivar’s reads to create a consensus sequence. However, starting from version 2, that step has been omitted. Instead, the VCF file is used to determine which mutations are derived from the wild-type cultivar and which are from the mutant.

### When there are many contigs  
The current setting has been updated to pick the top 50 contigs based on length, so only these contigs will be displayed in the plot. However, the table contains SNP-index information for all contigs, allowing you to confirm significant SNPs even for contigs not shown in the plot. You can also generate plots for these SNPs independently if needed. Since contigs smaller than the sliding window size often produce less reliable results, excluding them from the analysis should not be an issue.

### How to filter for higher-confidence SNPs  
If the phenotype is clear and the sequence data is clean, the default settings should already show some results. If you want to focus on higher-confidence SNPs, you can sequence both parents and retain only the SNPs that are clearly 0/0 in one parent and 1/1 in the other in the VCF file, which should produce cleaner results. However, in MutMap, keep in mind that you cannot input the mutant parent's sequence—only the sequence of the original cultivar can be used as input. [The linked page](https://github.com/YuSugihara/MutMap/blob/master/docs/01_MutMap_Commands_Breakdown.md) explains default MutMap commands, which involve only one parent in the workflow, but it may still be helpful as a reference when creating VCFs that include both parents.

### About multiple testing correction  
This function has been deprecated since v2.3.5. We highly recommend running MutMap without this function. However, if you would like to use this function, you can use it with versions of MutMap older than v2.3.5.
