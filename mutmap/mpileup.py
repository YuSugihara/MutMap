import sys
import re
import os
import glob
import subprocess as sbp
from multiprocessing import Pool
from mutmap.utils import time_stamp, clean_cmd, call_log


class Mpileup(object):

    def __init__(self, args):
        self.out = args.out
        self.args = args

    def get_bams(self, label):
        bams = glob.glob('{}/20_bam/{}*.bam'.format(self.out, label))
        bams = [bam for bam in bams if not bam.endswith('.unsorted.bam')]
        return bams

    def sort(self):
        for label in ['cultivar', 'bulk']:
            bams = self.get_bams(label)
            unsorted_bam = '{}/20_bam/{}.unsorted.bam'.format(self.out, label)
            if len(bams) == 0:
                sys.stderr.write(('{} Error: No input BAM files were found for "{}".\n'
                                  '{} Please check the files under "{}/20_bam".\n').format(time_stamp(),
                                                                                             label,
                                                                                             time_stamp(),
                                                                                             self.out))
                sys.exit(1)
            elif len(bams) == 1:
                path_to_bam = os.path.abspath(bams[0])
                if os.path.lexists(unsorted_bam):
                    os.remove(unsorted_bam)
                try:
                    os.symlink(path_to_bam, unsorted_bam)
                except OSError as err:
                    sys.stderr.write(('{} Error: Failed to prepare an intermediate BAM symlink.\n'
                                      '{} Source: "{}"\n'
                                      '{} Target: "{}"\n'
                                      '{} Reason: {}\n'
                                      '{} Please check file permissions and available disk space.\n').format(time_stamp(),
                                                                                                                time_stamp(),
                                                                                                                path_to_bam,
                                                                                                                time_stamp(),
                                                                                                                unsorted_bam,
                                                                                                                time_stamp(),
                                                                                                                err,
                                                                                                                time_stamp()))
                    sys.exit(1)
            else:
                cmd1 = 'samtools merge -f {0}/20_bam/{1}.unsorted.bam \
                                          {0}/20_bam/{1}*.bam \
                                          >> {0}/log/samtools.log \
                                          2>&1'.format(self.out, label)
                cmd1 = clean_cmd(cmd1)

            cmd2 = 'samtools collate -O \
                                     -@ {3} \
                                     {0}/20_bam/{1}.unsorted.bam \
                                     2>> {0}/log/samtools.log | \
                    samtools fixmate -m \
                                     - \
                                     - \
                                     2>> {0}/log/samtools.log | \
                    samtools sort -m {2} \
                                  -@ {3} \
                                  2>> {0}/log/samtools.log | \
                    samtools markdup -r \
                                     - \
                                     - \
                                     2>> {0}/log/samtools.log | \
                    samtools view -b \
                                  -f 2 \
                                  -F 2048 \
                                  -o {0}/20_bam/{1}.bam \
                                  >> {0}/log/samtools.log \
                                  2>&1'.format(self.out,
                                               label,
                                               self.args.mem,
                                               self.args.threads)

            cmd3 = 'samtools index {0}/20_bam/{1}.bam \
                                   >> {0}/log/samtools.log \
                                   2>&1'.format(self.out, label)

            cmd4 = 'samtools index -c \
                                   {0}/20_bam/{1}.bam \
                                   >> {0}/log/samtools.log \
                                   2>&1'.format(self.out, label)

            cmd5 = 'rm -f {}/20_bam/{}.*.bam'.format(self.out, label)

            cmd2 = clean_cmd(cmd2)
            cmd3 = clean_cmd(cmd3)
            cmd4 = clean_cmd(cmd4)
            cmd5 = clean_cmd(cmd5)

            if len(bams) > 1:
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
        cmd_idx1 = 'tabix -f -p vcf {0}/30_vcf/mutmap.{1}.vcf.gz \
                    >> {0}/log/bcftools.{1}.log \
                    2>&1'.format(self.out, chr_name)
        cmd_idx2 = 'tabix -f -C -p vcf {0}/30_vcf/mutmap.{1}.vcf.gz \
                    >> {0}/log/bcftools.{1}.log \
                    2>&1'.format(self.out, chr_name)
        cmd_idx1 = clean_cmd(cmd_idx1)
        cmd_idx2 = clean_cmd(cmd_idx2)

        try:
            sbp.run(cmd, 
                    stdout=sbp.DEVNULL, 
                    stderr=sbp.DEVNULL, 
                    shell=True, 
                    check=True)
        except sbp.CalledProcessError as err:
            raise RuntimeError(('bcftools mpileup/call failed for contig "{}". '
                                'Please check "{}/log/bcftools.{}.log".').format(chr_name,
                                                                                  self.out,
                                                                                  chr_name)) from err

        try:
            sbp.run(cmd_idx1,
                    stdout=sbp.DEVNULL,
                    stderr=sbp.DEVNULL,
                    shell=True,
                    check=True)
        except sbp.CalledProcessError:
            try:
                sbp.run(cmd_idx2,
                        stdout=sbp.DEVNULL,
                        stderr=sbp.DEVNULL,
                        shell=True,
                        check=True)
            except sbp.CalledProcessError as err:
                raise RuntimeError(('tabix indexing failed for contig "{}". '
                                    'Please check "{}/log/bcftools.{}.log".').format(chr_name,
                                                                                      self.out,
                                                                                      chr_name)) from err
    
    def concat(self):
        vcf_files = glob.glob('{}/30_vcf/mutmap.*.vcf.gz'.format(self.out))
        if len(vcf_files) == 0:
            sys.stderr.write(('{} Error: No per-contig VCF files were found for concatenation.\n'
                              '{} Variant calling may have failed earlier. Please check "{}/log/bcftools.*.log".\n').format(time_stamp(),
                                                                                                                             time_stamp(),
                                                                                                                             self.out))
            sys.exit(1)

        cmd1 = 'bcftools concat -O z \
                                -o {0}/30_vcf/mutmap.vcf.gz \
                                {0}/30_vcf/mutmap.*.vcf.gz \
                                >> {0}/log/bcftools.log \
                                2>&1'.format(self.out)
        cmd2 = 'cat {0}/log/bcftools.*.log > {0}/log/bcftools.log'.format(self.out)
        cmd3 = 'rm -f {}/30_vcf/mutmap.*.vcf.gz'.format(self.out)
        cmd4 = 'rm -f {}/log/bcftools.*.log'.format(self.out)

        cmd1 = clean_cmd(cmd1)
        cmd2 = clean_cmd(cmd2)
        cmd3 = clean_cmd(cmd3)
        cmd4 = clean_cmd(cmd4)

        try:
            sbp.run(cmd1, 
                    stdout=sbp.DEVNULL, 
                    stderr=sbp.DEVNULL, 
                    shell=True, 
                    check=True)
        except sbp.CalledProcessError:
            call_log(self.out, 'bcftools', cmd1)
            sys.exit(1)

        sbp.run(cmd2, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd3, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)
        sbp.run(cmd4, stdout=sbp.DEVNULL, stderr=sbp.DEVNULL, shell=True, check=True)

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
        print(time_stamp(), 'Sorting BAM files.', flush=True)
        self.sort()
        print(time_stamp(), 'BAM file sorting completed successfully.', flush=True)

        print(time_stamp(), 'Calling variants.', flush=True)
        chr_names = self.get_header()
        if len(chr_names) == 0:
            sys.stderr.write(('{} Error: No contig/chromosome entries were found in the reference FASTA header.\n'
                              '{} Please check "{}".\n').format(time_stamp(),
                                                                time_stamp(),
                                                                self.args.ref))
            sys.exit(1)

        p = Pool(self.args.threads)
        try:
            p.map(self.mpileup, chr_names)
        except RuntimeError as err:
            sys.stderr.write('{} Error: {}\n'.format(time_stamp(), err))
            sys.exit(1)
        except Exception as err:
            sys.stderr.write(('{} Error: Variant calling failed in multiprocessing workers.\n'
                              '{} Reason: {}\n'
                              '{} Please check "{}/log/bcftools.*.log".\n').format(time_stamp(),
                                                                                    time_stamp(),
                                                                                    err,
                                                                                    time_stamp(),
                                                                                    self.out))
            sys.exit(1)
        finally:
            p.close()
            p.join()

        self.concat()
        print(time_stamp(), 'Variant calling completed successfully.', flush=True)

        print(time_stamp(), 'Indexing the VCF file.', flush=True)
        self.mkindex()
        print(time_stamp(), 'VCF indexing completed successfully.', flush=True)
