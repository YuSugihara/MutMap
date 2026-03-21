import os
import sys
import re
import argparse
from multiprocessing import Pool
import multiprocessing as multi
from mutmap.utils import time_stamp
from mutmap.__init__ import __version__


class Params(object):

    def __init__(self, program_name):
        self.program_name = program_name

    def set_options(self):
        if self.program_name == 'mutmap':
            parser = self.mutmap_options()
        elif self.program_name == 'mutplot':
            parser = self.mutplot_options()

        if len(sys.argv) == 1:
            args = parser.parse_args(['-h'])
        else:
            args = parser.parse_args()
        return args

    def mutmap_options(self):
        parser = argparse.ArgumentParser(description='MutMap version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('mutmap -r <FASTA> -c <BAM|FASTQ> -b <BAM|FASTQ>\n'
                        '              -n <INT> -o <OUT_DIR> [-T] [-e <DATABASE>]')

        # set options
        parser.add_argument('-r',
                            '--ref',
                            action='store',
                            required=True,
                            type=str,
                            help='Reference FASTA file.',
                            metavar='')

        parser.add_argument('-c',
                            '--cultivar',
                            action='append',
                            required=True,
                            type=str,
                            help=('FASTQ or BAM file of cultivar. If specifying\n'
                                  'FASTQ, separate paired-end files with a comma,\n'
                                  'e.g., -c fastq1,fastq2. This option can be\n'
                                  'used multiple times.'),
                            metavar='')

        parser.add_argument('-b',
                            '--bulk',
                            action='append',
                            required=True,
                            type=str,
                            help=('FASTQ or BAM file of mutant bulk. If specifying\n'
                                  'FASTQ, separate paired-end files with a comma,\n'
                                  'e.g., -b fastq1,fastq2. This option can be\n'
                                  'used multiple times.'),
                            metavar='')

        parser.add_argument('-n',
                            '--N-bulk',
                            action='store',
                            required=True,
                            type=int,
                            help='Number of individuals in the mutant bulk.',
                            metavar='')

        parser.add_argument('-o',
                            '--out',
                            action='store',
                            required=True,
                            type=str,
                            help=('Output directory. The specified directory must not\n'
                                  'already exist.'), 
                            metavar='')

        parser.add_argument('-t',
                            '--threads',
                            action='store',
                            default=2,
                            type=int,
                            help=('Number of threads. If a value less than 1 is specified,\n'
                                  'MutMap will use the maximum available threads. [2]'),
                            metavar='')

        parser.add_argument('-w',
                            '--window',
                            action='store',
                            default=2000,
                            type=int,
                            help='Window size in kilobases (kb). [2000]',
                            metavar='')

        parser.add_argument('-s',
                            '--step',
                            action='store',
                            default=100,
                            type=int,
                            help='Step size in kilobases (kb). [100]',
                            metavar='')

        parser.add_argument('-D',
                            '--max-depth',
                            action='store',
                            default=250,
                            type=int,
                            help=('Maximum depth of variants to be used. This cutoff\n'
                                  'applies to both the cultivar and the bulk. [250]'),
                            metavar='')

        parser.add_argument('-d',
                            '--min-depth',
                            action='store',
                            default=8,
                            type=int,
                            help=('Minimum depth of variants to be used. This cutoff\n'
                                  'applies to both the cultivar and the bulk. [8]'),
                            metavar='')

        parser.add_argument('-N',
                            '--N-rep',
                            action='store',
                            default=5000,
                            type=int,
                            help=('Number of replicates for simulations to generate\n'
                                  'null distribution. [5000]'),
                            metavar='')

        parser.add_argument('-T',
                            '--trim',
                            action='store_true',
                            default=False,
                            help='Trim FASTQ files using Trimmomatic.')

        parser.add_argument('-a',
                            '--adapter',
                            action='store',
                            type=str,
                            help=('FASTA file containing adapter sequences. This option\n'
                                  'is used when "-T" is specified for trimming.'),
                            metavar='')

        parser.add_argument('--trim-params',
                            action='store',
                            type=str,
                            help=('Parameters for Trimmomatic. Input parameters\n'
                                  'must be comma-separated in the following order:\n'
                                  'Phred score, ILLUMINACLIP, LEADING, TRAILING,\n'
                                  ' SLIDINGWINDOW, MINLEN. To remove Illumina adapters,\n'
                                  'specify the adapter FASTA file with "--adapter".\n'
                                  'If not specified, adapter trimming will be skipped.\n'
                                  '[33,<ADAPTER_FASTA>:2:30:10,20,20,4:15,75]'),
                            metavar='')

        parser.add_argument('-e',
                            '--snpEff',
                            action='store',
                            type=str,
                            help=('Predict causal variants using SnpEff. Check\n'
                                  'available databases in SnpEff.'),
                            metavar='')

        parser.add_argument('--line-colors',
                            action='store',
                            default='#D55D00,#009E72,#FDB003',
                            type=str,
                            help=('Colors for threshold lines in plots. Specify a\n'
                                  'comma-separated list in the order of SNP-index,\n'
                                  'p95, and p99. ["#D55D00,#009E72,#FDB003"]'),
                            metavar='')

        parser.add_argument('--dot-color',
                            action='store',
                            default='#B3B8DD',
                            type=str,
                            help=('Color of the dots in plots. ["#B3B8DD"]'),
                            metavar='')


        parser.add_argument('--mem',
                            action='store',
                            default='1G',
                            type=str,
                            help=('Maximum memory per thread when sorting BAM files;\n'
                                  'suffixes K/M/G are recognized. [1G]'),
                            metavar='')

        parser.add_argument('-q',
                            '--min-MQ',
                            action='store',
                            default=40,
                            type=int,
                            help='Minimum mapping quality for mpileup. [40]',
                            metavar='')

        parser.add_argument('-Q',
                            '--min-BQ',
                            action='store',
                            default=18,
                            type=int,
                            help='Minimum base quality for mpileup. [18]',
                            metavar='')

        parser.add_argument('-C',
                            '--adjust-MQ',
                            action='store',
                            default=50,
                            type=int,
                            help=('Adjust the mapping quality for mpileup. The default\n'
                                  'setting is optimized for BWA. [50]'),
                            metavar='')

        # set version
        parser.add_argument('-v',
                            '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        return parser

    def mutplot_options(self):
        parser = argparse.ArgumentParser(description='MutPlot version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('mutplot -v <VCF> -o <OUT_DIR> -n <INT> [-w <INT>] [-s <INT>]\n'
                        '               [-D <INT>] [-d <INT>] [-N <INT>] [-m <FLOAT>]\n'
                        '               [-S <INT>] [-e <DATABASE>] [--igv] [--indel]\n')

        # set options
        parser.add_argument('-v',
                            '--vcf',
                            action='store',
                            required=True,
                            type=str,
                            help=('VCF file which contains cultivar and mutant bulk.\n'
                                  'in this order. This VCF file must have AD field.'),
                            metavar='')

        parser.add_argument('-o',
                            '--out',
                            action='store',
                            required=True,
                            type=str,
                            help=('Output directory. The specified directory can already\n'
                                  'exist.'),
                            metavar='')

        parser.add_argument('-n',
                            '--N-bulk',
                            action='store',
                            required=True,
                            type=int,
                            help='Number of individuals in the mutant bulk.',
                            metavar='')

        parser.add_argument('-w',
                            '--window',
                            action='store',
                            default=2000,
                            type=int,
                            help='Window size in kilobases (kb). [2000]',
                            metavar='')

        parser.add_argument('-s',
                            '--step',
                            action='store',
                            default=100,
                            type=int,
                            help='Step size in kilobases (kb). [100]',
                            metavar='')

        parser.add_argument('-D',
                            '--max-depth',
                            action='store',
                            default=250,
                            type=int,
                            help=('Maximum depth of variants to be used. This cutoff\n'
                                  'applies to both the cultivar and the bulk. [250]'),
                            metavar='')

        parser.add_argument('-d',
                            '--min-depth',
                            action='store',
                            default=8,
                            type=int,
                            help=('Minimum depth of variants to be used. This cutoff\n'
                                  'applies to both the cultivar and the bulk. [8]'),
                            metavar='')

        parser.add_argument('-N',
                            '--N-rep',
                            action='store',
                            default=5000,
                            type=int,
                            help=('Number of replicates for simulations to generate\n'
                                  'null distribution. [5000]'),
                            metavar='')

        parser.add_argument('-m',
                            '--min-SNPindex',
                            action='store',
                            default=0.3,
                            type=float,
                            help='Cutoff of minimum SNP-index for clear results. [0.3]',
                            metavar='')

        parser.add_argument('-S',
                            '--strand-bias',
                            action='store',
                            default=7,
                            type=int,
                            help=('Filter out spurious homozygous genotypes in the cultivar\n'
                                  'based on strand bias. If ADF (or ADR) is higher than\n'
                                  'this cutoff when ADR (or ADF) is 0, that SNP will be\n'
                                  'filtered out. If you want to disable this filtering,\n'
                                  'set this cutoff to 0. [7]\n'),
                            metavar='')

        parser.add_argument('-e',
                            '--snpEff',
                            action='store',
                            type=str,
                            help=('Predict causal variants using SnpEff. Check\n'
                                  'available databases in SnpEff.'),
                            metavar='')

        parser.add_argument('--igv',
                            action='store_true',
                            default=False,
                            help='Output IGV format file to check results on IGV.')

        parser.add_argument('--indel',
                            action='store_true',
                            default=False,
                            help='Plot SNP-index with INDEL.')

        parser.add_argument('--line-colors',
                            action='store',
                            default='#D55D00,#009E72,#FDB003',
                            type=str,
                            help=('Colors for threshold lines in plots. Specify a\n'
                                  'comma-separated list in the order of SNP-index,\n'
                                  'p95, and p99. ["#D55D00,#009E72,#FDB003"]'),
                            metavar='')

        parser.add_argument('--dot-color',
                            action='store',
                            default='#B3B8DD',
                            type=str,
                            help=('Color of the dots in plots. ["#B3B8DD"]'),
                            metavar='')

        parser.add_argument('--fig-width',
                            action='store',
                            default=7.5,
                            type=float,
                            help='Width allocated in chromosome figure. [7.5]',
                            metavar='')

        parser.add_argument('--fig-height',
                            action='store',
                            default=4.0,
                            type=float,
                            help='Height allocated in chromosome figure. [4.0]',
                            metavar='')

        parser.add_argument('--white-space',
                            action='store',
                            default=0.6,
                            type=float,
                            help=('White space between figures. (This option\n'
                                  'only affects vertical direction.) [0.6]'),
                            metavar='')

        parser.add_argument('-f',
                            '--format',
                            action='store',
                            default='png',
                            type=str,
                            help=('Specify the format of an output image.\n'
                                  'eps/jpeg/jpg/pdf/pgf/png/rgba/svg/svgz/tif/tiff'),
                            metavar='')

        # set version
        parser.add_argument('--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))

        return parser

    def check_max_threads(self, args):
        max_cpu = multi.cpu_count()
        if max_cpu < args.threads:
            sys.stderr.write(('!!WARNING!! You can use up to {0} threads. '
                              'The program will now use {0} threads.\n').format(max_cpu))
            sys.stderr.flush()
            args.threads = max_cpu
        elif args.threads < 1:
            args.threads = max_cpu
        return args

    def check_args(self, args):
        args.out = args.out.strip()
        args.ref = args.ref.strip()

        if args.out == '':
            sys.stderr.write(('  Error: Output directory path is empty.\n\n'))
            sys.exit(1)

        if args.ref == '':
            sys.stderr.write(('  Error: Reference FASTA path is empty.\n\n'))
            sys.exit(1)

        if os.path.exists(args.out):
            sys.stderr.write(('  Error: Output path already exists.\n'
                              '  Please rename the output path or '
                              'remove the existing path.\n\n'))
            sys.exit(1)

        if not os.path.isfile(args.ref):
            sys.stderr.write(('  Error: Reference FASTA was not found.\n'
                              '  Please check "{}".\n\n').format(args.ref))
            sys.exit(1)

        if args.adapter is not None and not args.trim:
            sys.stderr.write(('!!WARNING!! You specified "--adapter" but did not specify "--trim".\n'
                              '!!WARNING!! "--trim" has been automatically added.\n'))
            sys.stderr.flush()
            args.trim = True

        if args.trim_params is not None and not args.trim:
            sys.stderr.write(('!!WARNING!! You specified "--trim-params" but did not specify "--trim".\n'
                              '!!WARNING!! "--trim" has been automatically added.\n'))
            sys.stderr.flush()
            args.trim = True

        if args.trim and args.trim_params is None:
            args.trim_params = '33,<ADAPTER_FASTA>:2:30:10,20,20,4:15,75'

        N_fastq = 0

        for i, input_name in enumerate(args.cultivar):
            input_name = input_name.strip()
            args.cultivar[i] = input_name
            n_comma = input_name.count(',')
            if n_comma == 0:
                root, ext = os.path.splitext(input_name)
                if ext.lower() != '.bam':
                    sys.stderr.write(('  Error: Please check "{}".\n'
                                      '  The file extension is not "bam".\n'
                                      '  If you intended to specify FASTQ files, please provide them as paired-end reads separated by a comma. e.g., -c fastq1,fastq2\n\n').format(input_name))
                    sys.exit(1)
                if not os.path.isfile(input_name):
                    sys.stderr.write(('  Error: Input file was not found.\n'
                                      '  Please check "{}".\n\n').format(input_name))
                    sys.exit(1)
            elif n_comma == 1:
                fastqs = [fastq.strip() for fastq in input_name.split(',')]
                if '' in fastqs:
                    sys.stderr.write(('  Error: Please check "{}".\n'
                                      '  FASTQ file path is empty around a comma.\n\n').format(input_name))
                    sys.exit(1)
                args.cultivar[i] = ','.join(fastqs)
                for fastq in fastqs:
                    root, ext = os.path.splitext(fastq)
                    if ext.lower() == '.bam':
                        sys.stderr.write(('  Error: Please check "{}".\n'
                                          '  The file extension should not be "bam".\n'
                                          '  If you intended to specify BAM files, please input them separately. e.g., -c bam1 -c bam2\n\n').format(input_name))
                        sys.exit(1)
                    if not os.path.isfile(fastq):
                        sys.stderr.write(('  Error: Input file was not found.\n'
                                          '  Please check "{}".\n\n').format(fastq))
                        sys.exit(1)
                N_fastq += 1
            else:
                sys.stderr.write(('  Error: Please check "{}".\n'
                                  '  You have specified too many files, or the file name contains a comma.\n\n').format(input_name))
                sys.exit(1)

        for i, input_name in enumerate(args.bulk):
            input_name = input_name.strip()
            args.bulk[i] = input_name
            n_comma = input_name.count(',')
            if n_comma == 0:
                root, ext = os.path.splitext(input_name)
                if ext.lower() != '.bam':
                    sys.stderr.write(('  Error: Please check "{}".\n'
                                      '  The file extension is not "bam".\n'
                                      '  If you intended to specify FASTQ files, please provide them as paired-end reads separated by a comma. e.g., -b fastq1,fastq2\n\n').format(input_name))
                    sys.exit(1)
                if not os.path.isfile(input_name):
                    sys.stderr.write(('  Error: Input file was not found.\n'
                                      '  Please check "{}".\n\n').format(input_name))
                    sys.exit(1)
            elif n_comma == 1:
                fastqs = [fastq.strip() for fastq in input_name.split(',')]
                if '' in fastqs:
                    sys.stderr.write(('  Error: Please check "{}".\n'
                                      '  FASTQ file path is empty around a comma.\n\n').format(input_name))
                    sys.exit(1)
                args.bulk[i] = ','.join(fastqs)
                for fastq in fastqs:
                    root, ext = os.path.splitext(fastq)
                    if ext.lower() == '.bam':
                        sys.stderr.write(('  Error: Please check "{}".\n'
                                          '  The file extension should not be "bam".\n'
                                          '  If you intended to specify BAM files, please input them separately. e.g., -b bam1 -b bam2\n\n').format(input_name))
                        sys.exit(1)
                    if not os.path.isfile(fastq):
                        sys.stderr.write(('  Error: Input file was not found.\n'
                                          '  Please check "{}".\n\n').format(fastq))
                        sys.exit(1)
                N_fastq += 1
            else:
                sys.stderr.write(('  Error: Please check "{}".\n'
                                  '  You have specified too many files, or the file name contains a comma.\n\n').format(input_name))
                sys.exit(1)

        if N_fastq == 0 and args.trim:
            sys.stderr.write(('  Error: "--trim" can only be used when FASTQ files are specified.\n\n'))
            sys.exit(1)

        if args.N_bulk < 1:
            sys.stderr.write(('  Error: "--N-bulk" must be >= 1.\n\n'))
            sys.exit(1)

        if args.window < 1:
            sys.stderr.write(('  Error: "--window" must be >= 1 (kb).\n\n'))
            sys.exit(1)

        if args.step < 1:
            sys.stderr.write(('  Error: "--step" must be >= 1 (kb).\n\n'))
            sys.exit(1)

        if args.N_rep < 1:
            sys.stderr.write(('  Error: "--N-rep" must be >= 1.\n\n'))
            sys.exit(1)

        if args.min_depth < 1:
            sys.stderr.write(('  Error: "--min-depth" must be >= 1.\n\n'))
            sys.exit(1)

        if args.max_depth < args.min_depth:
            sys.stderr.write(('  Error: "--max-depth" must be greater than or equal to "--min-depth".\n\n'))
            sys.exit(1)

        if args.min_MQ < 0:
            sys.stderr.write(('  Error: "--min-MQ" must be >= 0.\n\n'))
            sys.exit(1)

        if args.min_BQ < 0:
            sys.stderr.write(('  Error: "--min-BQ" must be >= 0.\n\n'))
            sys.exit(1)

        if args.adjust_MQ < 0:
            sys.stderr.write(('  Error: "--adjust-MQ" must be >= 0.\n\n'))
            sys.exit(1)

        mem_value = args.mem.strip()
        if not re.match(r'^\d+[KMGkmg]$', mem_value):
            sys.stderr.write(('  Error: "--mem" must be specified as an integer followed by K, M, or G.\n'
                              '  Example: --mem 1G\n\n'))
            sys.exit(1)
        args.mem = mem_value.upper()

        line_colors = [color.strip() for color in args.line_colors.split(',')]
        if len(line_colors) != 3 or '' in line_colors:
            sys.stderr.write(('  Error: "--line-colors" must have exactly 3 comma-separated values.\n'
                              '  Order: SNP-index,p95,p99\n\n'))
            sys.exit(1)
        args.line_colors = ','.join(line_colors)

        args.dot_color = args.dot_color.strip()
        if args.dot_color == '':
            sys.stderr.write(('  Error: "--dot-color" must not be empty.\n\n'))
            sys.exit(1)

        return N_fastq
