# MutMap User Guide
#### version 2.3.5

## Table of contents
- [What is MutMap?](#what-is-mutmap)
- [Installation](#installation)
  + [Dependencies](#dependencies)
  + [Installation using bioconda](#installation-using-bioconda)
  + [Manual installation](#manual-installation)
- [Usage](#usage)
  + [Example 1 : run MutMap from FASTQ without trimming](#example-1--run-mutmap-from-fastq-without-trimming)
  + [Example 2 : run MutMap from FASTQ with trimming](#example-2--run-mutmap-from-fastq-with-trimming)
  + [Example 3 : run MutMap from BAM](#example-3--run-mutmap-from-bam)
  + [Example 4 : run MutMap from multiple FASTQs and BAMs](#example-4--run-mutmap-from-multiple-fastqs-and-bams)
  + [Example 5 : run MutPlot from VCF](#example-5--run-mutplot-from-vcf)
- [Outputs](#outputs)
- [About multiple testing correction](#about-multiple-testing-correction)
- [Built and use your own database for snpEff](#built-and-use-your-own-database-for-snpeff)

## What is MutMap?
<img src="https://github.com/YuSugihara/MutMap/blob/master/images/1_logo.png" width=200>

Bulked segregant analysis, as implemented in MutMap ([Abe et al., 2012](https://www.nature.com/articles/nbt.2095)), is a powerful and efficient method to identify agronomically important loci in crop plants. MutMap requires whole-genome resequencing of a single individual from the original cultivar and the pooled sequences of F2 progeny from a cross between the original cultivar and mutant. MutMap uses the sequence of the original cultivar to polarize the site frequencies of neighbouring markers and identifies loci with an unexpected site frequency, simulating the genotype of F2 progeny. **The updated pipeline is approximately 5-8 times faster than the previous pipeline, are easier for novice users to use and can be easily installed through bioconda with all dependencies.**

#### Citation
- Yu Sugihara, Lester Young, Hiroki Yaegashi, Satoshi Natsume, Daniel J. Shea, Hiroki Takagi, Helen Booker, Hideki Innan, Ryohei Terauchi, Akira Abe (2022). [High performance pipeline for MutMap and QTL-seq](https://doi.org/10.7717/peerj.13170). PeerJ, 10:e13170.

- Akira Abe, Shunichi Kosugi, Kentaro Yoshida, Satoshi Natsume, Hiroki Takagi, Hiroyuki Kanzaki, Hideo Matsumura, Kakoto Yoshida, Chikako Mitsuoka, Muluneh Tamiru, Hideki Innan, Liliana Cano, Sophien Kamoun & Ryohei Terauchi (2012). [Genome sequencing reveals agronomically important loci in rice using MutMap](https://www.nature.com/articles/nbt.2095). Nature Biotechnol. 30:174-179.

## Installation
### Dependencies
#### Softwares
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

### Installation using bioconda
You can install MutMap using [bioconda](https://bioconda.github.io/index.html).
```
conda create -c bioconda -n mutmap mutmap
conda activate mutmap
```

### Manual installation
If you got a error during installation, you can install MutMap, manually.
```
git clone https://github.com/YuSugihara/MutMap.git
cd MutMap
pip install -e .
```
Then you have to install other dependencies by yourself. We highly recommend you to install SnpEff and Trimmomatic using bioconda.
```
conda install -c bioconda snpeff
conda install -c bioconda trimmomatic
```
After installation, please check whether SnpEff and Trimmomatic work through the commands below.
```
snpEff --help
trimmomatic --help
```

## Usage

If your reference genome has more than 50 contigs (or chromosomes), only significant contigs will be plotted.

```
mutmap -h

usage: mutmap -r <FASTA> -c <BAM|FASTQ> -b <BAM|FASTQ>
              -n <INT> -o <OUT_DIR> [-T] [-e <DATABASE>]

MutMap version 2.3.5

options:
  -h, --help         show this help message and exit
  -r , --ref         Reference fasta.
  -c , --cultivar    fastq or bam of cultivar. If you specify
                     fastq, please separate pairs by comma,
                     e.g. -c fastq1,fastq2. You can use this
                     optiion multiple times
  -b , --bulk        fastq or bam of mutnat bulk. If you specify
                     fastq, please separate pairs by comma,
                     e.g. -b fastq1,fastq2. You can use this
                     optiion multiple times
  -n , --N-bulk      Number of individuals in mutant bulk.
  -o , --out         Output directory. Specified name must not
                     exist.
  -t , --threads     Number of threads. If you specify the number
                     below one, then MutMap will use the threads
                     as many as possible. [2]
  -w , --window      Window size (kb). [2000]
  -s , --step        Step size (kb). [100]
  -D , --max-depth   Maximum depth of variants which will be used.
                     This cutoff will be applied in both of cultivar
                     and bulk. [250]
  -d , --min-depth   Minimum depth of variants which will be used.
                     This cutoff will be applied in both of cultivar
                     and bulk. [8]
  -N , --N-rep       Number of replicates for simulation to make 
                     null distribution. [5000]
  -T, --trim         Trim fastq using trimmomatic.
  -a , --adapter     FASTA of adapter sequences. This will be used
                     when you specify "-T" for trimming.
  --trim-params      Parameters for trimmomatic. Input parameters
                     must be separated by comma with following
                     order: phred, ILLUMINACLIP, LEADING, TRAILING,
                     SLIDINGWINDOW, MINLEN. If you want to remove
                     adapters of illumina, please specify FASTA of
                     the adapter sequences with "--adapter". Specified
                     name will be inserted into <ADAPTER_FASTA>. If you
                     don't specify it, adapter trimming will be skipped.
                     [33,<ADAPTER_FASTA>:2:30:10,20,20,4:15,75]
  -e , --snpEff      Predict causal variant using SnpEff. Please
                     check available databases in SnpEff.
  --mem              Maximum memory per thread when bam sorted;
                     suffix K/M/G recognized. [1G]
  -q , --min-MQ      Minimum mapping quality in mpileup. [40]
  -Q , --min-BQ      Minimum base quality in mpileup. [18]
  -C , --adjust-MQ   "adjust-MQ" in mpileup. Default parameter
                     is suited for BWA. [50]
  -v, --version      show program's version number and exit
```

MutMap can run from FASTQ (without or with trimming) and BAM. If you want to run MutMap from VCF, please use MutPlot (example 5). Once you run MutMap, MutMap automatically complete the subprocesses.

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

`-c` : FASTQs of cultivar. Please input pair-end reads separated by comma. FASTQs can be gzipped.

`-b` : FASTQs of bulk. Please input pair-end reads separated by comma. FASTQs can be gzipped.

`-n` : number of individuals in mutant bulk.

`-o` : name of output directory. Specified name should not exist.

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

`-c` : FASTQs of cultivar. Please input pair-end reads separated by comma. FASTQs can be gzipped.

`-b` : FASTQs of bulk. Please input pair-end reads separated by comma. FASTQs can be gzipped.

`-n` : number of individuals in mutant bulk.

`-o` : name of output directory. Specified name should not exist.

`-T` : trim your reads by trimmomatic.

`-a` : FASTA of adapter sequences.

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

`-b` : BAM of bulk.

`-n` : number of individuals in mutant bulk.

`-o` : name of output directory. Specified name should not exist.

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

MutMap can automatically merge multiple FASTQs and BAMs. Of course, you can merge FASTQs or BAMs using `cat` or `samtools merge` before input them to MutMap. If you specify `-c` multiple times, please make sure that those files include only 1 individual. On the other hand, `-b` can include more than 1 individuals because those are bulked samples. MutMap can automatically classify FASTQs and BAMs from whether comma exits or not.

### Example 5 : run MutPlot from VCF
```
mutplot -h

usage: mutplot -v <VCF> -o <OUT_DIR> -n <INT> [-w <INT>] [-s <INT>]
               [-D <INT>] [-d <INT>] [-N <INT>] [-m <FLOAT>]
               [-S <INT>] [-e <DATABASE>] [--igv] [--indel]

MutPlot version 2.3.5

options:
  -h, --help            show this help message and exit
  -v , --vcf            VCF file which contains cultivar and mutant bulk.
                        in this order. This VCF file must have AD field.
  -o , --out            Output directory. Specified name can exist.
  -n , --N-bulk         Number of individuals in mutant bulk.
  -w , --window         Window size (kb). [2000]
  -s , --step           Step size (kb). [100]
  -D , --max-depth      Maximum depth of variants which will be used.
                        This cutoff will be applied in both of cultivar
                        and bulk. [250]
  -d , --min-depth      Minimum depth of variants which will be used.
                        This cutoff will be applied in both of cultivar
                        and bulk. [8]
  -N , --N-rep          Number of replicates for simulation to make 
                        null distribution. [5000]
  -m , --min-SNPindex   Cutoff of minimum SNP-index for clear results. [0.3]
  -S , --strand-bias    Filter spurious homo genotypes in cultivar using
                        strand bias. If ADF (or ADR) is higher than this
                        cutoff when ADR (or ADF) is 0, that SNP will be
                        filtered out. If you want to supress this filtering,
                        please set this cutoff to 0. [7]
  -e , --snpEff         Predict causal variant using SnpEff. Please
                        check available databases in SnpEff.
  --igv                 Output IGV format file to check results on IGV.
  --indel               Plot SNP-index with INDEL.
  --fig-width           Width allocated in chromosome figure. [7.5]
  --fig-height          Height allocated in chromosome figure. [4.0]
  --white-space         White space between figures. (This option
                        only affect vertical direction.) [0.6]
  -f , --format         Specifiy the format of an output image.
                        eps/jpeg/jpg/pdf/pgf/png/rgba/svg/svgz/tif/tiff
  --version             show program's version number and exit
```
MutPlot is included in MutMap. MutMap run MutPlot after making VCF. Then, MutPlot will work with default parameters. If you want to change some parameters, you can use VCF inside of `(OUT_DIR/30_vcf/mutmap.vcf.gz)` to retry plotting process like below.

```
mutplot -v OUT_DIR/30_vcf/mutmap.vcf.gz \
        -o ANOTHER_DIR_NAME \
        -n 20 \
        -w 2000 \
        -s 100
```

#### Use MutPlot for VCF which was made by yourself
In this case, please make sure that:
1. Your VCF include AD format.
2. Your VCF include two columns of cultivar and mutant bulk in this order.

If you got a error, please try to run MutMap from FASTQ or BAM before asking in issues.

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
|   |-- bulk.filt.bam
|   |-- bulk.filt.bam.bai
|   |-- cultivar.filt.bam
|   `-- cultivar.filt.bam.bai
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
    |-- bcftools.log
    |-- bgzip.log
    |-- bwa.log
    |-- mutplot.log
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
  + `mutmap_plot.png` : resulting plot (like below)
    - **<span style="color: blue; ">BLUE dot</span>** : variant
    - **<span style="color: red; ">RED line</span>** : mean SNP-index
    - **<span style="color: orange; ">ORANGE line</span>** : mean p99
    - **<span style="color: green; ">GREEN line</span>** : mean p95

<img src="https://github.com/YuSugihara/MutMap/blob/master/images/2_result.png" width=600>

## About multiple testing correction
Now this function is deprecated since v2.3.5.
We highly recommend running MutMap without this function.
However, if you would like to use this function, you can use it with versions of MutMap older than v2.3.5.

## Built and use your own database for snpEff
If you want to use your own database for snpEff, you need additional steps.
Here we assume that you installed MutMap via anaconda distribution, creating new environment with `conda create`.

1. Find the directory of snpEff that includes snpEff script, configuration file and database. You can find it in `/home/anaconda3/envs/{your_env_name_installed_mutmap}/share/snpeff-5.0-0/`. `anaconda3` may be `miniconda3`. Also, the version of snpeff may be different.

2. Go to this directory and follow the snpEff manual to build the database.
Don't forget to add your database info to the snpEff configuration file.
https://pcingola.github.io/SnpEff/se_buildingdb/#add-a-genome-to-the-configuration-file

3. Run MutMap with option `-e {your_database_name}`
