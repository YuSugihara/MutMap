#!/usr/bin/env python3

from mutmap.params import Params
pm = Params('mutplot')
args = pm.set_options()


import os
import sys
import pandas as pd
import subprocess as sbp
from multiprocessing import set_start_method
try:
    set_start_method('fork')
except (RuntimeError, ValueError):
    pass
if os.environ.get('JAVA_TOOL_OPTIONS') == None:
    os.environ['JAVA_TOOL_OPTIONS'] = '-Xmx4000m'
from mutmap.vcf2index import Vcf2Index
from mutmap.plot import Plot
from mutmap.utils import time_stamp, clean_cmd, check_executables


class MutPlot(object):

    def __init__(self, args):
        self.out = args.out
        self.vcf = args.vcf
        self.snpEff = args.snpEff
        self.args = args
        self.validate_args()

        if not os.path.isfile(self.vcf):
            sys.stderr.write(('{} Error: VCF file was not found.\n'
                              '{} Please check "{}".\n').format(time_stamp(),
                                                                time_stamp(),
                                                                self.vcf))
            sys.exit(1)

        if self.snpEff is not None:
            check_executables(['snpEff'], program_name='MutPlot')

    def validate_args(self):
        if self.args.N_bulk < 1:
            sys.stderr.write(('{} Error: "--N-bulk" must be >= 1.\n').format(time_stamp()))
            sys.exit(1)

        if self.args.window < 1:
            sys.stderr.write(('{} Error: "--window" must be >= 1 (kb).\n').format(time_stamp()))
            sys.exit(1)

        if self.args.step < 1:
            sys.stderr.write(('{} Error: "--step" must be >= 1 (kb).\n').format(time_stamp()))
            sys.exit(1)

        if self.args.N_rep < 1:
            sys.stderr.write(('{} Error: "--N-rep" must be >= 1.\n').format(time_stamp()))
            sys.exit(1)

        if self.args.min_depth < 1:
            sys.stderr.write(('{} Error: "--min-depth" must be >= 1.\n').format(time_stamp()))
            sys.exit(1)

        if self.args.max_depth < self.args.min_depth:
            sys.stderr.write(('{} Error: "--max-depth" must be greater than or equal to "--min-depth".\n').format(time_stamp()))
            sys.exit(1)

        if self.args.min_SNPindex < 0 or self.args.min_SNPindex > 1:
            sys.stderr.write(('{} Error: "--min-SNPindex" must be between 0 and 1.\n').format(time_stamp()))
            sys.exit(1)

        if self.args.strand_bias < 0:
            sys.stderr.write(('{} Error: "--strand-bias" must be >= 0.\n').format(time_stamp()))
            sys.exit(1)

        if self.args.fig_width <= 0 or self.args.fig_height <= 0:
            sys.stderr.write(('{} Error: "--fig-width" and "--fig-height" must be > 0.\n').format(time_stamp()))
            sys.exit(1)

        if self.args.white_space < 0:
            sys.stderr.write(('{} Error: "--white-space" must be >= 0.\n').format(time_stamp()))
            sys.exit(1)

    def check_outdir(self):
        if os.path.exists(self.out):
            if not os.path.isdir(self.out):
                sys.stderr.write(('{} Error: Output path exists but is not a directory.\n'
                                  '{} Please check "{}".\n').format(time_stamp(),
                                                                    time_stamp(),
                                                                    self.out))
                sys.exit(1)
            print(time_stamp(),
                  'Output directory already exists.',
                  flush=True)
        else:
            os.mkdir(self.out)

    def check_required_tables(self):
        required = ['snp_index.tsv', 'sliding_window.tsv']
        missing = []
        for filename in required:
            path = '{}/{}'.format(self.out, filename)
            if not os.path.isfile(path):
                missing.append(path)
        if missing:
            raise FileNotFoundError(('Required intermediate file(s) were not found:\n'
                                     '{}').format('\n'.join(missing)))

    def run_snpEff(self):
        cmd = 'snpEff ann -Xmx4g \
                          -s {0}/snpEff_summary.html \
                          {1} \
                          {2} \
                          1> {0}/mutmap.snpEff.vcf \
                          2> {0}/snpEff.log'.format(self.out,
                                                    self.snpEff,
                                                    self.vcf)
        cmd = clean_cmd(cmd)

        print(time_stamp(), 'Running SnpEff annotation.', flush=True)
        try:
            sbp.run(cmd, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        except sbp.CalledProcessError:
            sys.stderr.write(('{} Error: SnpEff annotation failed.\n'
                              '{} Please check "{}/snpEff.log".\n').format(time_stamp(),
                                                                           time_stamp(),
                                                                           self.out))
            sys.exit(1)
        print(time_stamp(), 'SnpEff annotation completed successfully.', flush=True)

        self.args.vcf = '{}/mutmap.snpEff.vcf'.format(self.out)

    def get_outlier_SNPindex(self):
        if self.snpEff is None:
            snp_index = pd.read_csv('{}/snp_index.tsv'.format(self.out),
                                    sep='\t',
                                    names=['CHROM',
                                           'POSI',
                                           'variant',
                                           'depth',
                                           'p99',
                                           'p95',
                                           'SNPindex'])
        else:
            snp_index = pd.read_csv('{}/snp_index.tsv'.format(self.out),
                                    sep='\t',
                                    names=['CHROM',
                                           'POSI',
                                           'variant',
                                           'impact',
                                           'depth',
                                           'p99',
                                           'p95',
                                           'SNPindex'])

        snp_index[abs(snp_index['p99']) <= \
                  abs(snp_index['SNPindex'])].to_csv('{}/snp_index.p99.tsv'.format(self.out),
                                                     sep='\t',
                                                     index=False)

        snp_index[abs(snp_index['p95']) <= \
                  abs(snp_index['SNPindex'])].to_csv('{}/snp_index.p95.tsv'.format(self.out),
                                                     sep='\t',
                                                     index=False)

    def get_outlier_windows(self):
        sliding_window = pd.read_csv('{}/sliding_window.tsv'.format(self.out),
                                     sep='\t',
                                     names=['CHROM',
                                            'POSI',
                                            'mean_p99',
                                            'mean_p95',
                                            'mean_SNPindex'])

        sliding_window[abs(sliding_window['mean_p99']) <= \
                       abs(sliding_window['mean_SNPindex'])].to_csv('{}/sliding_window.p99.tsv'.format(self.out),
                                                                    sep='\t',
                                                                    index=False)

        sliding_window[abs(sliding_window['mean_p95']) <= \
                       abs(sliding_window['mean_SNPindex'])].to_csv('{}/sliding_window.p95.tsv'.format(self.out),
                                                                    sep='\t',
                                                                    index=False)

    def make_igv_file(self):
        if self.snpEff is None:
            snp_index = pd.read_csv('{}/snp_index.tsv'.format(self.out),
                                    sep='\t',
                                    names=['CHROM',
                                           'POSI',
                                           'variant',
                                           'depth',
                                           'p99',
                                           'p95',
                                           'SNPindex'])
        else:
            snp_index = pd.read_csv('{}/snp_index.tsv'.format(self.out),
                                    sep='\t',
                                    names=['CHROM',
                                           'POSI',
                                           'variant',
                                           'impact',
                                           'depth',
                                           'p99',
                                           'p95',
                                           'SNPindex'])

        snp_index['Start'] = snp_index['POSI'] - 1
        snp_index['End'] = snp_index['POSI']
        snp_index['Feature'] = snp_index['CHROM'].astype('str') \
                               + ':' + snp_index['POSI'].astype('str') \
                               + ':' + snp_index['variant'].astype('str') \
                               + ':' + snp_index['p99'].astype('str') \
                               + ':' + snp_index['p95'].astype('str')

        snp_index = snp_index[['CHROM',
                               'Start',
                               'End',
                               'Feature',
                               'SNPindex']]

        snp_index.to_csv('{}/snp_index.igv'.format(self.out),
                         sep='\t',
                         index=False)

        sliding_window = pd.read_csv('{}/sliding_window.tsv'.format(self.out),
                                     sep='\t',
                                     names=['CHROM',
                                            'POSI',
                                            'mean_p99',
                                            'mean_p95',
                                            'mean_SNPindex'])

        sliding_window['Start'] = sliding_window['POSI'] - 1
        sliding_window['End'] = sliding_window['POSI']
        sliding_window['Feature'] = sliding_window['CHROM'].astype('str') \
                                    + ':' + sliding_window['POSI'].astype('str') \
                                    + ':' + sliding_window['mean_p99'].astype('str') \
                                    + ':' + sliding_window['mean_p95'].astype('str')

        sliding_window = sliding_window[['CHROM',
                                         'Start',
                                         'End',
                                         'Feature',
                                         'mean_SNPindex']]

        sliding_window.to_csv('{}/sliding_window.igv'.format(self.out),
                              sep='\t',
                              index=False)

    def run(self):
        self.check_outdir()

        if self.snpEff is not None:
            self.run_snpEff()

        v2i = Vcf2Index(self.args)
        try:
            v2i.run()
        except SystemExit:
            raise
        except Exception as err:
            sys.stderr.write(('{} Error: Failed while converting VCF to MutPlot intermediate tables.\n'
                              '{} Reason: {}\n').format(time_stamp(),
                                                        time_stamp(),
                                                        err))
            sys.exit(1)

        print(time_stamp(), 'Generating plot...', flush=True)
        pt = Plot(self.args)
        try:
            pt.run()
        except (RuntimeError, ValueError) as err:
            sys.stderr.write('{} Error: {}\n'.format(time_stamp(), err))
            sys.exit(1)
        except Exception as err:
            sys.stderr.write(('{} Error: Plot generation failed unexpectedly.\n'
                              '{} Reason: {}\n').format(time_stamp(),
                                                        time_stamp(),
                                                        err))
            sys.exit(1)

        try:
            self.check_required_tables()
            self.get_outlier_SNPindex()
            self.get_outlier_windows()
        except Exception as err:
            sys.stderr.write(('{} Error: Failed while generating summary tables (p95/p99).\n'
                              '{} Reason: {}\n').format(time_stamp(),
                                                        time_stamp(),
                                                        err))
            sys.exit(1)

        if self.args.igv:
            try:
                self.make_igv_file()
            except Exception as err:
                sys.stderr.write(('{} Error: Failed while writing IGV output files.\n'
                                  '{} Reason: {}\n').format(time_stamp(),
                                                            time_stamp(),
                                                            err))
                sys.exit(1)

def main():
    print(time_stamp(), 'Starting MutPlot.', flush=True)
    mp = MutPlot(args)
    mp.run()
    print(time_stamp(), 'MutPlot completed successfully.', flush=True)

if __name__ == '__main__':
    main()
