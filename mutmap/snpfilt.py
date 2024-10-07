class SnpFilt(object):

    def __init__(self, args):
        self.min_SNPindex = args.min_SNPindex
        self.maxDP = args.max_depth
        self.minDP = args.min_depth
        self.strand_bias = args.strand_bias

    def filt_cultivar_gt(self, cultivar_GT, bulk_AD):
        record = {}
        # only use homo in cultivar
        if cultivar_GT in ['0/0', '0|0', '1/1', '1|1']:
            ADs = bulk_AD.split(',')
            # check biallele or multi-allele
            if len(ADs) == 2:
                # filter missing
                if not '.' in ADs:
                    # check whether REF homo or ALT homo.
                    if cultivar_GT in ['0/0', '0|0']:
                        record['bulk_ref_AD'] = int(ADs[0])
                        record['bulk_alt_AD'] = int(ADs[1])
                        # if depth of ALT is zero in bulk,
                        # it will be discarded because SNP-index will be zero.
                        if record['bulk_alt_AD'] != 0:
                            record['type'] = 'keep'
                            record['cultivar_GT'] = '0/0'
                        else:
                            record['type'] = 'discard'
                    else:
                        record['bulk_ref_AD'] = int(ADs[1])
                        record['bulk_alt_AD'] = int(ADs[0])
                        # if depth of REF is zero in bulk,
                        # it will be discarded because SNP-index will be zero.
                        if record['bulk_ref_AD'] != 0:
                            record['type'] = 'keep'
                            record['cultivar_GT'] = '1/1'
                        else:
                            record['type'] = 'discard'
                else:
                    record['type'] = 'discard'
            # check ALT homo in cultivar
            elif len(ADs) == 3:
                # filter missing
                if not '.' in ADs:
                    if cultivar_GT in ['1/1', '1|1']:
                        if int(ADs[0]) == 0:
                            record['bulk_ref_AD'] = int(ADs[1])
                            record['bulk_alt_AD'] = int(ADs[2])
                            record['type'] = 'keep'
                            record['cultivar_GT'] = '1/1'
                        else:
                            record['type'] = 'discard'
                    else:
                        record['type'] = 'discard'
                else:
                    record['type'] = 'discard'
            else:
                record['type'] = 'discard'
        else:
            record['type'] = 'discard'
        return record

    def filt_depth(self, record, cultivar_AD):
        record['cultivar_depth'] = sum([int(AD) for AD in cultivar_AD.split(',')])
        record['bulk_depth'] = record['bulk_ref_AD'] + record['bulk_alt_AD']

        if record['cultivar_depth'] < self.minDP or self.maxDP < record['cultivar_depth']:
            record['type'] = 'discard'
        elif record['bulk_depth'] < self.minDP or self.maxDP < record['bulk_depth']:
            record['type'] = 'discard'
        return record

    def filt_index(self, record):
        record['SNPindex'] = record['bulk_alt_AD']/record['bulk_depth']
        if record['SNPindex'] < self.min_SNPindex:
            record['type'] = 'discard'
        return record

    def filt_strand_bias(self, record, ADFR):
        ADF = ADFR[0]
        ADR = ADFR[1]

        if record['cultivar_GT'] == '0/0':
            cultivar_ADF = int(ADF.split(',')[0])
            cultivar_ADR = int(ADR.split(',')[0])
        else:
            cultivar_ADF = int(ADF.split(',')[1])
            cultivar_ADR = int(ADR.split(',')[1])

        if cultivar_ADF == 0 and cultivar_ADR > self.strand_bias:
            record['type'] = 'discard'
        elif cultivar_ADR == 0 and cultivar_ADF > self.strand_bias:
            record['type'] = 'discard'
        return record


    def filt(self, cultivar_GT, cultivar_AD, bulk_AD, ADFR):
        record = self.filt_cultivar_gt(cultivar_GT, bulk_AD)
        if record['type'] == 'keep':
            record = self.filt_depth(record, cultivar_AD)
            if record['type'] == 'keep':
                record = self.filt_index(record)
                if record['type'] == 'keep':
                    if ADFR != None:
                        record = self.filt_strand_bias(record, ADFR)

        return record
