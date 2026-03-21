#!/usr/bin/env python3

from mutmap.params import Params
pm = Params('mutmap')
args = pm.set_options()


import os
import sys
import subprocess as sbp
from multiprocessing import set_start_method
try:
    set_start_method('fork')
except (RuntimeError, ValueError):
    pass
from mutmap.refindex import RefIndex
from mutmap.trim import Trim
from mutmap.alignment import Alignment
from mutmap.mpileup import Mpileup
from mutmap.utils import time_stamp, check_executables


class MutMap(object):

    def __init__(self, args):
        args = pm.check_max_threads(args)
        self.N_fastq = pm.check_args(args)
        self.args = args
        self.out = args.out
        self.cultivar_fastq, self.cultivar_bam = self.get_files(args.cultivar)
        self.bulk_fastq, self.bulk_bam = self.get_files(args.bulk)

        self.check_dependencies()
        self.set_runtime_environment()
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

    def check_dependencies(self):
        commands = ['bwa', 'samtools', 'bcftools', 'tabix']
        if self.args.trim:
            commands.append('trimmomatic')
        if self.args.snpEff is not None:
            commands.append('snpEff')
        check_executables(commands, program_name='MutMap')

    def set_runtime_environment(self):
        # These environment variables can avoid potential threading/runtime issues.
        if os.environ.get('OMP_NUM_THREADS') is None:
            os.environ['OMP_NUM_THREADS'] = str(self.args.threads)
        if os.environ.get('USE_SIMPLE_THREADED_LEVEL3') is None:
            os.environ['USE_SIMPLE_THREADED_LEVEL3'] = str(self.args.threads)
        if os.environ.get('JAVA_TOOL_OPTIONS') is None:
            os.environ['JAVA_TOOL_OPTIONS'] = '-Xmx4000m'

    def mkdir(self):
        os.makedirs('{}'.format(self.out))
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

    def refindex(self):
        ri = RefIndex(self.args)
        ri.run()

    def trimming(self):
        tm = Trim(self.args)
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
        aln = Alignment(self.args)
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

    def mpileup(self):
        os.mkdir('{}/30_vcf'.format(self.out))
        mp = Mpileup(self.args)
        mp.run()

    def mutplot(self):
        cmd = [sys.executable,
               '-m',
               'mutmap.mutplot',
               '-v',
               '{}/30_vcf/mutmap.vcf.gz'.format(self.out),
               '-n',
               str(self.args.N_bulk),
               '-w',
               str(self.args.window),
               '-s',
               str(self.args.step),
               '-N',
               str(self.args.N_rep),
               '-D',
               str(self.args.max_depth),
               '-d',
               str(self.args.min_depth),
               '--dot-color',
               self.args.dot_color,
               '--line-colors',
               self.args.line_colors,
               '-o',
               '{}/40_mutmap'.format(self.out)]

        if self.args.snpEff is not None:
            cmd.extend(['-e', self.args.snpEff])

        p = sbp.Popen(cmd,
                      stdout=sbp.PIPE,
                      stderr=sbp.STDOUT,
                      text=True)

        for line in p.stdout:
            print(line.rstrip(), flush=True)

        return_code = p.wait()
        if return_code != 0:
            sys.stderr.write(('{} Error: MutPlot failed during the post-processing step.\n').format(time_stamp()))
            sys.exit(return_code)

    def run(self):
        self.refindex()
        if self.args.trim:
            self.trimming()
        else:
            self.alignment()
        self.mpileup()
        self.mutplot()

def main():
    print(time_stamp(), 'Starting MutMap analysis.', flush=True)
    MutMap(args).run()
    print(time_stamp(), 'MutMap analysis completed successfully.\n', flush=True)

if __name__ == '__main__':
    main()
