
import os
import re
import sys
import gzip
import numpy as np
import pandas as pd
from functools import lru_cache
from mutmap.snpfilt import SnpFilt
from mutmap.smooth import Smooth
from mutmap.utils import time_stamp


class Vcf2Index(object):

    def __init__(self, args):
        self.out = args.out
        self.vcf = args.vcf
        self.snpEff = args.snpEff
        self.N_bulk = args.N_bulk
        self.N_replicates = args.N_rep
        self.min_SNPindex = args.min_SNPindex
        self.snp_index = '{}/snp_index.tsv'.format(self.out)
        self.args = args
        if self.snpEff is not None:
            self.ANN_re = re.compile(';ANN=(.*);*')

        self.p99_index = int(0.99*self.N_replicates) - 1
        self.p95_index = int(0.95*self.N_replicates) - 1

    def get_field(self):
        root, ext = os.path.splitext(self.vcf)
        if ext == '.gz':
            vcf = gzip.open(self.vcf, 'rt')
        else:
            vcf = open(self.vcf, 'r')
        found_variant = False
        for line in vcf:
            if re.match(r'[^#]', line):
                found_variant = True
                cols = line.split('\t')
                if len(cols) < 11:
                    sys.stderr.write(('{} Error: VCF must contain at least 11 columns (up to bulk sample column).\n'
                                      '{} Please check "{}".\n').format(time_stamp(),
                                                                        time_stamp(),
                                                                        self.vcf))
                    sys.exit(1)
                fields = cols[8].split(':')
                try:
                    GT_pos = fields.index('GT')
                except ValueError:
                    sys.stderr.write(('{} Error: No GT (Genotype) field found in your VCF file!\n').format(time_stamp()))
                    sys.exit(1)
                try:
                    AD_pos = fields.index('AD')
                except ValueError:
                    sys.stderr.write(('{} Error: No AD (Allelic Depth) field found in your VCF file!\n').format(time_stamp()))
                    sys.exit(1)

                if 'ADF' in fields and 'ADR' in fields:
                    ADF_pos = fields.index('ADF')
                    ADR_pos = fields.index('ADR')
                else:
                    ADF_pos = None
                    ADR_pos = None
                    sys.stderr.write(('{} Warning: No ADF (Allelic Depth Forward) or ADR (Allelic Depth Reverse) field found in your VCF file.\n').format(time_stamp()))
                    sys.stderr.write(('{} Warning: Strand bias filtering will be skipped.\n').format(time_stamp()))
                break
        vcf.close()
        if not found_variant:
            sys.stderr.write(('{} Error: No variant records were found in your VCF file.\n'
                              '{} Please check "{}".\n').format(time_stamp(),
                                                                time_stamp(),
                                                                self.vcf))
            sys.exit(1)
        return GT_pos, AD_pos, ADF_pos, ADR_pos

    def get_variant_impact(self, annotation):
        ANN_fields = self.ANN_re.findall(annotation)
        if len(ANN_fields) == 0 or ANN_fields[0] == '':
            raise ValueError('ANN field is missing in the INFO column.')
        ANN = ANN_fields[0]
        genes = ANN.split(',')
        # The impact of the variant is determined by the most severe impact
        # of the genes that the variant is located in.
        # The impact order is HIGH > MODERATE > LOW > MODIFIER
        # The impacts of variants are oredered by severity in the ANN field.
        impacts = []
        for gene in genes:
            gene_fields = gene.split('|')
            if len(gene_fields) > 2 and gene_fields[2] != '':
                impacts.append(gene_fields[2])
        if len(impacts) == 0:
            raise ValueError('ANN field exists but does not contain valid impact annotations.')
        if 'HIGH' in impacts:
            impact = 'HIGH'
        elif 'MODERATE' in impacts:
            impact = 'MODERATE'
        elif 'LOW' in impacts:
            impact = 'LOW'
        else:
            impact = 'MODIFIER'
        return impact

    def check_variant_type(self, REF, ALT):
        if (len(REF) == 1) and (len(ALT) == 1):
            variant = 'snp'
        else:
            variant = 'indel'
        return variant

    @lru_cache(maxsize=256)
    def F2_simulation(self, depth):
        GT = [0, 1]
        replicates = []
        n = 0
        while n < self.N_replicates:
            F2_GT = np.random.choice(GT, 2*self.N_bulk, p=[0.5, 0.5])
            mut_AD = np.random.choice(F2_GT, depth, replace=True)
            SNPindex = mut_AD.sum()/depth
            if SNPindex < self.min_SNPindex:
                continue
            else:
                replicates.append(SNPindex)
                n += 1
        replicates.sort()
        p99 = replicates[self.p99_index]
        p95 = replicates[self.p95_index]
        return p99, p95

    def calc_SNPindex(self, field_pos):
        root, ext = os.path.splitext(self.vcf)
        if ext == '.gz':
            vcf = gzip.open(self.vcf, 'rt')
        else:
            vcf = open(self.vcf, 'r')
        snp_index = open(self.snp_index, 'w')

        sf = SnpFilt(self.args)
        n_variants = 0
        n_kept = 0
        for i, line in enumerate(vcf, start=1):
            if re.match(r'[^#]', line):
                n_variants += 1
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 11:
                    raise ValueError(('VCF record has fewer than 11 columns at line {}. '
                                      'Please ensure the VCF contains cultivar and bulk sample columns.').format(i))
                CHR = cols[0]
                POS = cols[1]
                REF = cols[3]
                ALT = cols[4]
                annotation = cols[7]
                GT_pos = field_pos[0]
                AD_pos = field_pos[1]
                ADF_pos = field_pos[2]
                ADR_pos = field_pos[3]

                cultivar_fields = cols[9].split(':')
                bulk_fields = cols[10].split(':')
                if len(cultivar_fields) <= max(GT_pos, AD_pos):
                    raise ValueError(('FORMAT/sample mismatch at line {} in cultivar column. '
                                      'Required fields (GT/AD) were not found in expected positions.').format(i))
                if len(bulk_fields) <= AD_pos:
                    raise ValueError(('FORMAT/sample mismatch at line {} in bulk column. '
                                      'Required AD field was not found in expected position.').format(i))

                cultivar_GT = cultivar_fields[GT_pos]
                cultivar_AD = cultivar_fields[AD_pos]
                bulk_AD = bulk_fields[AD_pos]

                if ADF_pos != None and ADR_pos != None:
                    if len(cultivar_fields) <= max(ADF_pos, ADR_pos):
                        raise ValueError(('FORMAT/sample mismatch at line {} in cultivar column. '
                                          'ADF/ADR were expected but not found.').format(i))
                    cultivar_ADF = cultivar_fields[ADF_pos]
                    cultivar_ADR = cultivar_fields[ADR_pos]
                    ADFR = (cultivar_ADF, cultivar_ADR)
                else:
                    ADFR = None

                try:
                    record = sf.filt(cultivar_GT, cultivar_AD, bulk_AD, ADFR)
                except (ValueError, IndexError, ZeroDivisionError) as err:
                    raise ValueError(('Failed to parse sample fields at VCF line {}. '
                                      'Please check GT/AD/ADF/ADR values.\n'
                                      'Reason: {}').format(i, err)) from err

                if record['type'] == 'keep':
                    variant = self.check_variant_type(REF, ALT)

                    p99, p95 = self.F2_simulation(record['bulk_depth'])

                    if self.snpEff is None:
                        snp_index.write(('{}\t{}\t{}\t{}\t'
                                         '{:.4f}\t{:.4f}\t{:.4f}\n').format(CHR,
                                                                            POS,
                                                                            variant,
                                                                            record['bulk_depth'],
                                                                            p99,
                                                                            p95,
                                                                            record['SNPindex']))
                        n_kept += 1
                    else:
                        try:
                            impact = self.get_variant_impact(annotation)
                        except ValueError as err:
                            raise ValueError(('Failed to parse ANN annotation at VCF line {}.\n'
                                              'Reason: {}\n'
                                              'If you are using "-e/--snpEff", please confirm this VCF was annotated by SnpEff ann.').format(i, err)) from err
                        snp_index.write(('{}\t{}\t{}\t{}\t{}\t'
                                         '{:.4f}\t{:.4f}\t{:.4f}\n').format(CHR,
                                                                            POS,
                                                                            variant,
                                                                            impact,
                                                                            record['bulk_depth'],
                                                                            p99,
                                                                            p95,
                                                                            record['SNPindex']))
                        n_kept += 1

        snp_index.close()
        vcf.close()
        if n_variants == 0:
            raise RuntimeError('No variant records were found in the VCF body.')
        if n_kept == 0:
            raise RuntimeError(('No variants passed filtering and SNP-index criteria.\n'
                                'Please review filter parameters: -d/--min-depth, -D/--max-depth,\n'
                                '-m/--min-SNPindex, and -S/--strand-bias.'))

    def run(self):
        print(time_stamp(), 'Calculating SNP-index.', flush=True)
        field_pos = self.get_field()
        try:
            self.calc_SNPindex(field_pos)
        except (RuntimeError, ValueError) as err:
            sys.stderr.write('{} Error: {}\n'.format(time_stamp(), err))
            sys.exit(1)
        print(time_stamp(), 'SNP-index calculation completed successfully.', flush=True)

        print(time_stamp(), 'Smoothing SNP-index.', flush=True)
        sm = Smooth(self.args)
        try:
            sm.run()
        except Exception as err:
            sys.stderr.write(('{} Error: Failed while creating the sliding-window table.\n'
                              '{} Reason: {}\n').format(time_stamp(),
                                                        time_stamp(),
                                                        err))
            sys.exit(1)
        print(time_stamp(), 'SNP-index smoothing completed successfully.', flush=True)
