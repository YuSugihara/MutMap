import sys
import re
import os
import glob
import subprocess as sbp
from multiprocessing import Pool
from mutmap.vcf2index import Vcf2Index
from mutmap.utils import time_stamp, clean_cmd, call_log


class Mpileup(object):

    def __init__(self, args):
        self.out = args.out
        self.args = args

    def get_bams(self, label):
        bams = glob.glob('{}/20_bam/{}*.bam'.format(self.out, label))
        return bams

    def merge(self):
        for label in ['cultivar', 'bulk']:
            bams = self.get_bams(label)
            if len(bams) == 1:
                path_to_bam = os.path.abspath(bams[0])
                cmd1 = 'ln -s {} {}/20_bam/{}.unsorted.bam'.format(path_to_bam,
                                                                   self.out,
                                                                   label)
            else:
                cmd1 = 'samtools merge -f {0}/20_bam/{1}.unsorted.bam \
                                          {0}/20_bam/{1}*.bam \
                                          >> {0}/log/samtools.log \
                                          2>&1'.format(self.out, label)

            cmd2 = 'samtools sort -m {0} \
                                  -@ {1} \
                                  -o {2}/20_bam/{3}.bam \
                                  {2}/20_bam/{3}.unsorted.bam \
                                  >> {2}/log/samtools.log \
                                  2>&1'.format(self.args.mem,
                                               self.args.threads,
                                               self.out,
                                               label)

            cmd3 = 'samtools index {0}/20_bam/{1}.bam \
                                   >> {0}/log/samtools.log \
                                   2>&1'.format(self.out, label)

            cmd4 = 'samtools index -c \
                                   {0}/20_bam/{1}.bam \
                                   >> {0}/log/samtools.log \
                                   2>&1'.format(self.out, label)

            cmd5 = 'rm -f {}/20_bam/{}.*.bam'.format(self.out, label)

            cmd1 = clean_cmd(cmd1)
            cmd2 = clean_cmd(cmd2)
            cmd3 = clean_cmd(cmd3)
            cmd4 = clean_cmd(cmd4)
            cmd5 = clean_cmd(cmd5)

            try:
                sbp.run(cmd1, 
                        stdout=sbp.DEVNULL, 
                        stderr=sbp.DEVNULL, 
                        shell=True, 
                        check=True)
            except sbp.CalledProcessError:
                call_log(self.out, 'samtools', cmd1)
                sys.exit(1)

            try:
                sbp.run(cmd2, 
                        stdout=sbp.DEVNULL, 
                        stderr=sbp.DEVNULL, 
                        shell=True, 
                        check=True)
            except sbp.CalledProcessError:
                call_log(self.out, 'samtools', cmd2)
                sys.exit(1)

            try:
                sbp.run(cmd3, 
                        stdout=sbp.DEVNULL, 
                        stderr=sbp.DEVNULL, 
                        shell=True, 
                        check=True)
            except sbp.CalledProcessError:
                call_log(self.out, 'samtools', cmd3)
                sys.exit(1)

            try:
                sbp.run(cmd4, 
                        stdout=sbp.DEVNULL, 
                        stderr=sbp.DEVNULL, 
                        shell=True, 
                        check=True)
            except sbp.CalledProcessError:
                call_log(self.out, 'samtools', cmd4)
                sys.exit(1)

            sbp.run(cmd5, 
                    stdout=sbp.DEVNULL, 
                    stderr=sbp.DEVNULL, 
                    shell=True, 
                    check=True)

    def get_header(self):
        ref = open(self.args.ref, 'r')
        pattern = re.compile('>')
        chr_names = []
        for line in ref:
            if pattern.match(line):
                line = line.rstrip('\n')
                chr_name = re.split('[> ]',line)[1]
                chr_names.append(chr_name)
        return chr_names

    def mpileup(self, chr_name):
        cmd = 'bcftools mpileup -a AD,ADF,ADR \
                                -B \
                                -q {0}\
                                -Q {1} \
                                -C {2} \
                                -O u \
                                -r {3} \
                                -f {4} \
                                {5}/20_bam/cultivar.bam \
                                {5}/20_bam/bulk.bam | \
               bcftools call -vm \
                             -f GQ,GP \
                             -O u | \
               bcftools filter -i "INFO/MQ>={0}" \
                               -O z \
                               -o {5}/30_vcf/mutmap.{3}.vcf.gz \
                               >> {5}/log/bcftools.{3}.log \
                               2>&1'.format(self.args.min_MQ,
                                            self.args.min_BQ,
                                            self.args.adjust_MQ,
                                            chr_name,
                                            self.args.ref,
                                            self.out)

        cmd = clean_cmd(cmd)

        try:
            sbp.run(cmd, 
                    stdout=sbp.DEVNULL, 
                    stderr=sbp.DEVNULL, 
                    shell=True, 
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'bcftools', cmd)
            sys.exit(1)
    
    def concat(self):
        cmd1 = 'bcftools concat -O z \
                                -o {0}/30_vcf/mutmap.vcf.gz \
                                {0}/30_vcf/mutmap.*.vcf.gz \
                                >> {0}/log/bcftools.log \
                                2>&1'.format(self.out)
        cmd2 = 'cat {0}/log/bcftools.*.log > {0}/log/bcftools.log'.format(self.out)
        cmd3 = 'rm -f {}/30_vcf/mutmap.*.vcf.gz'.format(self.out)
        cmd4 = 'rm -f {}/log/bcftools.*.log'.format(self.out)
        cmd5 = 'rm -f {}/log/tabix.*.log'.format(self.out)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)
        cmd3 = clean_cmd(cmd3)
        cmd4 = clean_cmd(cmd4)
        cmd5 = clean_cmd(cmd5)

        try:
            sbp.run(cmd1, 
                    stdout=sbp.DEVNULL, 
                    stderr=sbp.DEVNULL, 
                    shell=True, 
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'bcftools', cmd3)
            sys.exit(1)

        sbp.run(cmd2, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd3, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd4, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd5, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

    def mkindex(self):
        cmd1 = 'tabix -f \
                      -p vcf \
                      {0}/30_vcf/mutmap.vcf.gz \
                      >> {0}/log/tabix.log \
                      2>&1'.format(self.out)
        
        cmd2 = 'tabix -f \
                      -C \
                      -p vcf \
                      {0}/30_vcf/mutmap.vcf.gz \
                      >> {0}/log/tabix.log \
                      2>&1'.format(self.out)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)

        try:
            sbp.run(cmd1, 
                    stdout=sbp.DEVNULL, 
                    stderr=sbp.DEVNULL, 
                    shell=True, 
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'tabix', cmd1)
            sys.exit(1)

        try:
            sbp.run(cmd2, 
                    stdout=sbp.DEVNULL, 
                    stderr=sbp.DEVNULL, 
                    shell=True, 
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'tabix', cmd2)
            sys.exit(1)

    def run(self):
        print(time_stamp(), 'start to merge BAMs.', flush=True)
        self.merge()
        print(time_stamp(), 'merging process successfully finished.', flush=True)

        print(time_stamp(), 'start to call variants.', flush=True)
        chr_names = self.get_header()

        p = Pool(self.args.threads)
        p.map(self.mpileup, chr_names)
        p.close()

        self.concat()
        print(time_stamp(), 'variant calling successfully finished.', flush=True)

        print(time_stamp(), 'start to index VCF.', flush=True)
        self.mkindex()
        print(time_stamp(), 'indexing VCF successfully finished.', flush=True)
