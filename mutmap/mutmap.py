#!/usr/bin/env python3

from mutmap.params import Params
pm = Params('mutmap')
args, config = pm.set_options()


import os
import sys
import glob
import subprocess as sbp
from mutmap.utils import time_stamp
from mutmap.utils import clean_cmd
from mutmap.trim import Trim
from mutmap.alignment import Alignment
from mutmap.bamfilt import BamFilt
from mutmap.mpileup import Mpileup
from mutmap.vcf2index import Vcf2Index


class MutMap(object):

    def __init__(self, args, config):
        args = pm.check_max_threads(args)
        self.N_fastq = pm.check_args(args)
        self.args = args
        self.config = config
        self.out = args.out
        self.cultivar_fastq, self.cultivar_bam = self.get_files(args.cultivar)
        self.bulk_fastq, self.bulk_bam = self.get_files(args.bulk)

        self.mkdir()
        self.link_ref()
        self.link_bam('cultivar', self.cultivar_bam)
        self.link_bam('bulk', self.bulk_bam)

    def get_files(self, input_names):
        fastq = []
        bam = []
        for input_name in input_names:
            if ',' in input_name:
                fastq.append(input_name)
            else:
                bam.append(input_name)
        return fastq, bam

    def mkdir(self):
        os.mkdir('{}'.format(self.out))
        os.mkdir('{}/log'.format(self.out))
        if self.N_fastq > 0 and self.args.trim:
            os.mkdir('{}/00_fastq'.format(self.out))
        os.mkdir('{}/10_ref'.format(self.out))
        os.mkdir('{}/20_bam'.format(self.out))

    def link_ref(self):
        path_to_ref = os.path.abspath(self.args.ref)
        ref = os.path.basename(self.args.ref)
        sym_ref = '{}/10_ref/{}'.format(self.out, ref)
        os.symlink(path_to_ref, sym_ref)
        self.args.ref = sym_ref

    def link_bam(self, label, input_files):
        for i, input_name in enumerate(input_files):
            path_to_bam = os.path.abspath(input_name)
            index = '{}.1{:0>3}'.format(label, i)
            os.symlink(path_to_bam,
                       '{}/20_bam/{}.bam'.format(self.out, index))

    def mkindex(self):
        print(time_stamp(),
              'start to index reference fasta by BWA.',
              flush=True)
        cmd = 'bwa index {} &>> {}/log/bwa.log'.format(self.args.ref,
                                                       self.out)
        sbp.run(cmd, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        print(time_stamp(),
              'indexing of reference successfully finished.',
              flush=True)

    def trimming(self):
        tm = Trim(args, config)
        for i, fastq in enumerate(self.cultivar_fastq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'cultivar.0{:0>3}'.format(i)
            tm.run(fastq1, fastq2, index)

        for i, fastq in enumerate(self.bulk_fastq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'bulk.0{:0>3}'.format(i)
            tm.run(fastq1, fastq2, index)

    def alignment(self):
        aln = Alignment(args, config)
        for i, fastq in enumerate(self.cultivar_fastq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'cultivar.0{:0>3}'.format(i)
            aln.run(fastq1, fastq2, index)

        for i, fastq in enumerate(self.bulk_fastq):
            fastq1 = fastq.split(',')[0]
            fastq2 = fastq.split(',')[1]
            index = 'bulk.0{:0>3}'.format(i)
            aln.run(fastq1, fastq2, index)

    def bamfilt(self):
        bt = BamFilt(self.args, self.config)
        bt.run()

    def mpileup(self):
        os.mkdir('{}/30_vcf'.format(self.out))
        mp = Mpileup(self.args, self.config)
        mp.run()

    def mutplot(self):
        if args.snpEff is None:
            cmd = 'mutplot -v {0}/30_vcf/mutmap.vcf.gz \
                                  -n {1} \
                                  -o {0}/40_mutmap \
                                  &>> {0}/log/mutplot.log'.format(self.out,
                                                                  self.args.N_bulk)
        else:
            cmd = 'mutplot -v {0}/30_vcf/mutmap.vcf.gz \
                                  -n {1} \
                                  -o {0}/40_mutmap \
                                  -e {2}\
                                  &>> {0}/log/mutplot.log'.format(self.out,
                                                                  self.args.N_bulk,
                                                                  self.args.snpEff)

        cmd = clean_cmd(cmd)
        sbp.run(cmd,
                stdout=sbp.DEVNULL,
                stderr=sbp.DEVNULL,
                shell=True,
                check=True)

    def run(self):
        self.mkindex()
        if self.args.trim:
            self.trimming()
        else:
            self.alignment()
        self.bamfilt()
        self.mpileup()
        self.mutplot()

def main():
    print(time_stamp(), 'start to run MutMap.', flush=True)
    MutMap(args, config).run()
    print(time_stamp(), 'MutMap successfully finished.', flush=True)

if __name__ == '__main__':
    main()
