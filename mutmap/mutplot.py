#!/usr/bin/env python3

from mutmap.params import Params
pm = Params('mutplot')
args = pm.set_options()


import os
import pandas as pd
import subprocess as sbp
from multiprocessing import set_start_method
set_start_method('fork')
from mutmap.vcf2index import Vcf2Index
from mutmap.plot import Plot
from mutmap.utils import time_stamp, clean_cmd


class MutPlot(object):

    def __init__(self, args):
        self.out = args.out
        self.vcf = args.vcf
        self.snpEff = args.snpEff
        self.args = args

    def check_outdir(self):
        if os.path.exists(self.out):
            print(time_stamp(),
                  'Output directory already exists.'.format(self.out),
                  flush=True)
        else:
            os.mkdir(self.out)

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
        sbp.run(cmd, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
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
        v2i.run()

        print(time_stamp(), 'Generating plot...', flush=True)
        pt = Plot(self.args)
        pt.run()

        self.get_outlier_SNPindex()
        self.get_outlier_windows()

        if self.args.igv:
            self.make_igv_file()

def main():
    print(time_stamp(), 'Starting MutPlot.', flush=True)
    mp = MutPlot(args)
    mp.run()
    print(time_stamp(), 'MutPlot completed successfully.', flush=True)

if __name__ == '__main__':
    main()
