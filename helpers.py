import sys
import re
import data
from constants import *
from sortedcontainers import SortedList, SortedSet
from intervaltree import Interval, IntervalTree

class Filter (object):
    
    def passesFilter(self):
        return False

class MaxSizeFilter(Filter):
    def __init__(self, INT_maxSize):
        self._maxSize = INT_maxSize

    def passesFilter(self, IOC_variant):
        return IOC_variant.getSize() <= self._maxSize

class MinSizeFilter(Filter):
    def __init__(self, INT_minSize):
        self._minSize = INT_minSize

    def passesFilter(self, IOC_variant):
        return IOC_variant.getSize() >= self._minSize          

class File(object):
    def __init__(self, STR_file):
        self._file = STR_file
        self._fstream = open(STR_file, 'r')

class PedigreeFile(File):
    affected    = '2'
    nonaffected    = '1'

    def __init__(self, STR_file):
        File.__init__(self, STR_file)

    def parse(self):
        self._fstream.seek(0)
        for STR_line in self._fstream.readlines():
            [STR_fam, STR_samp, STR_father, STR_mother, STR_sex, STR_status] = STR_line.rstrip().split('\t')
            IOC_family = data.Family.Get(STR_fam) if data.Family.Get(STR_fam) else data.Family.GetSet(data.Family(STR_fam))
            IOC_sample = data.Sample(IOC_family, STR_samp, STR_father, STR_mother, STR_sex, STR_status)
            IOC_family.addSample(IOC_sample)

class PopSVFile(File):
    def __init__(self, STR_file, STR_sourceBatch=''):
        """
            this class is used to process raw variant call data from a single caller, in a single callset (with multiple samples)
            at this stage of processing, all variants only have a *single* sample associated
            it is only at later processing stages that multiple sample-variants are merged into single variants
        """
        File.__init__(self, STR_file)
        self._sourceBatch = STR_sourceBatch

    def parse(self, SET_Filters):
        
        self._fstream.seek(0)
        # skip header line
        DIC_chr_interval_tree = dict()

        # prepare to handle skipped samples
        SET_skipped = set()
        
        # parse values from input file (skip header)
        for STR_line in self._fstream.readlines()[0:]:
            [   STR_sample, 
                STR_chr,  STR_start, STR_end,
                STR_nb_bin_cons, STR_z, STR_fc,
                STR_mean_cov, STR_pv, STR_qv,
                STR_cn2_dev, STR_cn, STR_prop_single_bin, STR_status
            ] = STR_line.rstrip().split('\t')
            # determine variant type based on copy number (skip cn==2 i.e. wildtype)          
            if int(STR_cn) < 2:
                STR_varType=Deletion
            elif int(STR_cn) > 2:
                STR_varType=Duplication
            else: continue
            
            # instantiate Variant with parsed data
            IOC_sample = data.Sample.Get(STR_sample)
            # if sample not in data.Sample (i.e. not in pedFile loaded in segregation.py), skip this variant
            if IOC_sample is None:
                if STR_sample not in SET_skipped:
                    sys.stderr.write('\tWARNING: skipping variants from sample {0} (sample not found in ped file)\n'.format(STR_sample))
                    SET_skipped.add(STR_sample)
                continue
            IOC_var = data.RawVariant(STR_chr, int(STR_start), STR_chr, int(STR_end), STR_varType, IOC_sample, popsv, {},
                # {
                #     'chr':STR_chr,
                #     'start':STR_start,
                #     'end':STR_end,
                #     'nb_bin_cons':STR_nb_bin_cons,
                #     'z':STR_z,
                #     'fc':STR_fc,
                #     'mean_cov':STR_mean_cov,
                #     'pv':STR_pv,
                #     'qv':STR_qv,
                #     'cn2_dev':STR_cn2_dev,
                #     'cn':STR_cn,
                #     'prop_single_bin':STR_prop_single_bin
                # },
                self._sourceBatch )
                 
            # if variant does not satisfy the filter rules, flush it                  
            passes_filters = True
            for my_filter in SET_Filters:
                if not my_filter.passesFilter(IOC_var):
                    passes_filters = False
                    break
            if not passes_filters: 
                continue

            # if variants already present, die hard with a vengeance (NOTE: used to simply skip the second instance of the variant and keep going
            if IOC_var in IOC_sample._contained_variants:
                raise RuntimeError('variant {0} found more than once in input dataset "{1}". Clean up your dataset before running SV segregation analysis'.format(IOC_var, self._sourceBatch))
            
            #variant is valid, adding it to sample
            IOC_sample.addVariantToSample(IOC_var)
            #adding keys to the dict
            if IOC_var._chrom1 not in DIC_chr_interval_tree: DIC_chr_interval_tree[IOC_var._chrom1] = dict()
            if IOC_sample not in DIC_chr_interval_tree[IOC_var._chrom1]: DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample] = dict()
            if STR_varType not in DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample]: DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample][STR_varType] = IntervalTree()
            #adding variant to tree
            DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample][STR_varType][IOC_var._pos1:IOC_var._pos2] = IOC_var                

        return DIC_chr_interval_tree


class LumpyFile(File):
    
    def __init__(self, STR_file, STR_sourceBatch=''):
        File.__init__(self, STR_file)
        self._sourceBatch = STR_sourceBatch
          
    def parse(self, SET_Filters):

        self._fstream.seek(0)
        # skip header line
        DIC_chr_interval_tree = dict()
        
        # prepare to handle skipped samples
        SET_skipped = set()
        
        for STR_line in self._fstream.readlines():
            [
                STR_sample,                             #1
                STR_chromosome_1,                       #2
                STR_interval_1_start,                   #3
                STR_interval_1_end,                     #4
                STR_chromosome_2,                       #5
                STR_interval_2_start,                   #6
                STR_interval_2_end,                     #7
                STR_id,                                 #8
                STR_evidence_set_score,                 #9
                STR_strand_1,                           #10
                STR_strand_2,                           #11
                STR_lumpy_type,                         #12
                STR_evidence_sets_supporting_reads,     #13
                STR_strand_configuration_counts,        #14
                STR_maximum_probability_positions,      #15
                STR_top95_percent_confidence_intervals  #16
            ] = STR_line.rstrip().split('\t')
         
            LIS_breakpoints = STR_maximum_probability_positions[4:].split(';')
            [INT_break1Pos, INT_break2Pos] = [int(STR_brk.split(':')[1]) for STR_brk in LIS_breakpoints]
            STR_varType=LUMPY_TYPES[STR_lumpy_type]
             
            
            #Modify calls with start = stop
            if STR_interval_1_start == STR_interval_1_end:
                STR_interval_1_end = int(STR_interval_1_end) + 1
            if STR_interval_2_start == STR_interval_2_end:
                STR_interval_2_end = int(STR_interval_2_end) + 1    
            
            # instantiate Variant with parsed data
            IOC_sample = data.Sample.Get(STR_sample)
            # if sample not in data.Sample (i.e. not in pedFile loaded in segregation.py), skip this variant
            if IOC_sample is None:
                if STR_sample not in SET_skipped:
                    sys.stderr.write('\tWARNINGwarning: skipping variants from sample {0} (sample not found in ped file)\n'.format(STR_sample))
                    SET_skipped.add(STR_sample)
                continue
            
            IOC_var = data.RawVariant(
                STR_chromosome_1, 
                INT_break1Pos, 
                STR_chromosome_2, 
                INT_break2Pos, 
                STR_varType, 
                IOC_sample, 
                lumpy, 
                {},
                # {
                #     'chromosome_1':STR_chromosome_1,
                #     'interval_1_start':STR_interval_1_start,
                #     'interval_1_end':STR_interval_1_end,
                #     'chromosome_2':STR_chromosome_2,
                #     'interval_2_start':STR_interval_2_start,
                #     'interval_2_end':STR_interval_2_end,
                #     'id':STR_id,
                #     'evidence_set_score':STR_evidence_set_score,
                #     'strand_1':STR_strand_1,
                #     'strand_2':STR_strand_2,
                #     'type':STR_lumpy_type,
                #     'evidence_sets_supporting_reads':STR_evidence_sets_supporting_reads,
                #     'strand_configuration_counts':STR_strand_configuration_counts,
                #     'maximum_probability_positions':STR_maximum_probability_positions,
                #     'top95_percent_confidence_intervals':STR_top95_percent_confidence_intervals
                # }, 
                self._sourceBatch,               
                int (STR_interval_1_start), 
                int (STR_interval_1_end), 
                int (STR_interval_2_start), 
                int (STR_interval_2_end),
            )              
               
            # if variant does not satisfy the filter rules, flush it          
            passes_filters = True
            if IOC_var._type == Translocation:
                passes_filters = True
            else:
                for my_filter in SET_Filters:
                    if not my_filter.passesFilter(IOC_var):
                        passes_filters = False
                        break
            if not passes_filters: 
                continue

            # if variants already present, die hard with a vengeance (NOTE: used to simply skip the second instance of the variant and keep going
            if IOC_var in IOC_sample._contained_variants:
                raise RuntimeError('variant {0} found more than once in input dataset "{1}". Clean up your dataset before running SV segregation analysis'.format(IOC_var, self._sourceBatch))

            #variant is valid, adding it to sample
            IOC_sample.addVariantToSample(IOC_var) 
            #adding keys to the dict            
            if IOC_var._chrom1 not in DIC_chr_interval_tree:                    
                DIC_chr_interval_tree[IOC_var._chrom1] = dict()            
            if STR_varType == Translocation and IOC_var._chrom2 not in DIC_chr_interval_tree:
                DIC_chr_interval_tree[IOC_var._chrom2] = dict()
            if IOC_sample not in DIC_chr_interval_tree[IOC_var._chrom1]:
                DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample] = dict()
            if STR_varType == Translocation and IOC_sample not in DIC_chr_interval_tree[IOC_var._chrom2]:
                DIC_chr_interval_tree[IOC_var._chrom2][IOC_sample] = dict()           
            if STR_varType not in DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample]:
                DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample][STR_varType] = IntervalTree()
            if STR_varType == Translocation and STR_varType not in DIC_chr_interval_tree[IOC_var._chrom2][IOC_sample]:
                DIC_chr_interval_tree[IOC_var._chrom2][IOC_sample][STR_varType] = IntervalTree()            
            
            if STR_varType == Translocation:
                DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample][STR_varType][IOC_var._INT_pos1_brkp1:IOC_var._INT_pos1_brkp2] = IOC_var
                DIC_chr_interval_tree[IOC_var._chrom2][IOC_sample][STR_varType][IOC_var._INT_pos2_brkp1:IOC_var._INT_pos2_brkp2] = IOC_var
            else:
                DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample][STR_varType][IOC_var._pos1:IOC_var._pos2] = IOC_var

        return DIC_chr_interval_tree

class ConiferFile(File):

    def __init__(self, STR_file, STR_sourceBatch=''):
        """
            Conifer files handled by this class must have a header, and 5 tab-separated columns in this specific format:
            sampleID	chromosome	start	stop	state
            S28108	chr1	147954794	148017703	dup
            S28108	chr1	196748266	196759385	del
        """
        File.__init__(self, STR_file)
        self._sourceBatch = STR_sourceBatch

    def parse(self, SET_Filters):

        self._fstream.seek(0)
        # skip header line
        self._fstream.readline()

        DIC_fixed_raw_variants = dict()

        # prepare to handle skipped samples
        SET_skipped = set()

        # parse values from input file (skip header)       
        for STR_line in self._fstream.readlines():
            LIS_line = STR_line.rstrip().split('\t')
            # sanity check: conifer format MUST have 5 tab-separated values
            if len(LIS_line) != 5: raise RuntimeError('incorrect number of columns in file {0} for line {1} (expected 5 (sampleID, chromosome, start, stop, state), got {2}'.format(self._file, STR_line.rstrip(), len(LIS_line)))   
            [STR_sample, STR_chr, STR_start, STR_end, STR_conifer_type] = LIS_line
            #removing the 'chr' in the chromosome id prefix
            STR_chr=STR_chr.replace('chr','')
            # little check to make sure the type fits with expected values
            if STR_conifer_type not in CONIFER_TYPES.keys(): raise RuntimeError('unrecognizable SV variant type (5th column of input file {0}) for line {1}'.format(self._file, STR_line.rstrip()))
            STR_varType = CONIFER_TYPES[STR_conifer_type]

            IOC_sample = data.Sample.Get(STR_sample)
            # if sample not in data.Sample (i.e. not in pedFile loaded in segregation.py), skip this variant
            if IOC_sample is None:
                if STR_sample not in SET_skipped:
                    sys.stderr.write('\tWARNING: skipping variants from sample {0} (sample not found in ped file)\n'.format(STR_sample))
                    SET_skipped.add(STR_sample)
                continue
            # instantiate Variant with parsed data
            IOC_var = data.RawVariant(STR_chr, int(STR_start), STR_chr, int(STR_end), STR_varType, IOC_sample, conifer, {}, self._sourceBatch )

            # if variant does not satisfy the filter rules, flush it                  
            passes_filters = True
            for my_filter in SET_Filters:
                if not my_filter.passesFilter(IOC_var):
                    passes_filters = False
                    break
            if not passes_filters: continue

            # if variants already present, die hard with a vengeance (NOTE: used to simply skip the second instance of the variant and keep going
            if IOC_var in IOC_sample._contained_variants:
                raise RuntimeError('variant {0} found more than once in input dataset "{1}". Clean up your dataset before running SV segregation analysis'.format(IOC_var, self._sourceBatch))
            #variant is valid, adding it to sample
            if not DIC_fixed_raw_variants.has_key(IOC_sample):
                DIC_fixed_raw_variants[IOC_sample] = dict()
            if not DIC_fixed_raw_variants[IOC_sample].has_key(STR_varType):
                DIC_fixed_raw_variants[IOC_sample][STR_varType] = SortedList(key=data.BaseVariant.sort_by_start)

            DIC_fixed_raw_variants[IOC_sample][STR_varType].add(IOC_var)

        return DIC_fixed_raw_variants


class GenericBedFile(File):
    
    def __init__(self, STR_file, STR_sourceBatch=''):
        """
            generic BED files handled by this class must have no header, and 5 tab-separated columns in this specific format:
            chrom, start, end, SV_type, sample_id
            NOTE: this class is not appropriate for translocation variants, which must be provided in BEDPE format; 
            a separate class will handle raw variant data files that include translocations
        """
        File.__init__(self, STR_file)
        self._sourceBatch = STR_sourceBatch

    def parse(self, SET_Filters):
        
        # DIC_fixed_raw_variants = dict()
        DIC_chr_interval_tree = dict()

        # prepare to handle skipped samples
        SET_skipped = set()
        
        # parse values from input file (skip header)       
        for STR_line in self._fstream.readlines():
            LIS_line = STR_line.rstrip().split('\t')
            # sanity check: genericBed format MUST have 5 tab-separated values
            if len(LIS_line) != 5: raise RuntimeError('incorrect number of columns in file {0} for line {1} (expected 5 (chrom, start, end, SV_type, sample), got {2}'.format(self._file, STR_line.rstrip(), len(LIS_line)))
            [STR_chr, STR_start, STR_end, STR_generic_type, STR_sample] = LIS_line
            # will convert any non-standard SV type strings into usable format
            STR_key = re.findall('(del|dup|inv|tra)', STR_generic_type.lower())[0].upper()
            # little check to make sure the type fits with expected values
            if STR_key not in GENERIC_TYPES.keys(): raise RuntimeError('unrecognizable SV variant type (4th column of input file {0}) for line {1}'.format(self._file, STR_line.rstrip())) 
            STR_varType = GENERIC_TYPES[STR_key]
            
            IOC_sample = data.Sample.Get(STR_sample)
            # if sample not in data.Sample (i.e. not in pedFile loaded in segregation.py), skip this variant
            if IOC_sample is None:
                if STR_sample not in SET_skipped:
                    sys.stderr.write('\tWARNING: skipping variants from sample {0} (sample not found in ped file)\n'.format(STR_sample))
                    SET_skipped.add(STR_sample)
                continue
            # instantiate Variant with parsed data
            IOC_var = data.RawVariant(STR_chr, int(STR_start), STR_chr, int(STR_end), STR_varType, IOC_sample, genericBed, {}, self._sourceBatch )
            
            # if variant does not satisfy the filter rules, flush it                  
            passes_filters = True
            for my_filter in SET_Filters:
                if not my_filter.passesFilter(IOC_var):
                    passes_filters = False
                    break
            if not passes_filters: 
                continue       

            # if variants already present, die hard with a vengeance (NOTE: used to simply skip the second instance of the variant and keep going
            if IOC_var in IOC_sample._contained_variants:
                raise RuntimeError('variant {0} found more than once in input dataset "{1}". Clean up your dataset before running SV segregation analysis'.format(IOC_var, self._sourceBatch))
            
            #variant is valid, adding it to sample
            IOC_sample.addVariantToSample(IOC_var)
            #adding keys to the dict            
            if IOC_var._chrom1 not in DIC_chr_interval_tree: DIC_chr_interval_tree[IOC_var._chrom1] = dict()
            if IOC_sample not in DIC_chr_interval_tree[IOC_var._chrom1]: DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample] = dict()
            if STR_varType not in DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample]: DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample][STR_varType] = IntervalTree()
            #adding variant to tree
            DIC_chr_interval_tree[IOC_var._chrom1][IOC_sample][STR_varType][IOC_var._pos1:IOC_var._pos2] = IOC_var
            # print (DIC_chr_interval_tree)

        return DIC_chr_interval_tree

class CnmopsFile(File):

    def __init__(self, STR_file, STR_sourceBatch=''):
        """
            Cnmops files handled by this class must have a header, and 9 tab-separated columns in this specific format:
            seqnames  start  end  width  strand  sampleName  median  mean  CN
            chrX  54334342  54335731 1390  *  XI001  -0.852460064002923  -0.852460064002923  CN1
        """
        File.__init__(self, STR_file)
        self._sourceBatch = STR_sourceBatch

    def parse(self, SET_Filters):

        self._fstream.seek(0)
        # skip header line 
        self._fstream.readline()

        DIC_fixed_raw_variants = dict()

        # prepare to handle skipped samples
        SET_skipped = set()

        # parse values from input file (skip header)       
        for STR_line in self._fstream.readlines():
            LIS_line = STR_line.rstrip().split(' ')
            # sanity check: conifer format MUST have 5 tab-separated values
            if len(LIS_line) != 9: raise RuntimeError('incorrect number of columns in file {0} for line {1} (expected 9 (seqnames start end width strand sampleName median mean CN ), got {2}'.format(self._file, STR_line.rstrip(), len(LIS_line)))   
            [STR_chr, STR_start, STR_end, STR_width, STR_strand, STR_sample, STR_median, STR_mean, STR_cnmops_type] = LIS_line
            #removing the 'chr' in the chromosome id prefix
            STR_chr=STR_chr.replace('chr','')
            # little check to make sure the type fits with expected values
            if STR_cnmops_type not in CNMOPS_TYPES.keys(): raise RuntimeError('unrecognizable SV variant type (9th column of input file {0}) for line {1}'.format(self._file, STR_line.rstrip()))
            STR_varType = CNMOPS_TYPES[STR_cnmops_type]

            IOC_sample = data.Sample.Get(STR_sample)
            # if sample not in data.Sample (i.e. not in pedFile loaded in segregation.py), skip this variant
            if IOC_sample is None:
                if STR_sample not in SET_skipped:
                    sys.stderr.write('\tWARNING: skipping variants from sample {0} (sample not found in ped file)\n'.format(STR_sample))
                    SET_skipped.add(STR_sample)
                continue
            # instantiate Variant with parsed data
            IOC_var = data.RawVariant(STR_chr, int(STR_start), STR_chr, int(STR_end), STR_varType, IOC_sample, cnmops, {}, self._sourceBatch )

            # if variant does not satisfy the filter rules, flush it                  
            passes_filters = True
            for my_filter in SET_Filters:
                if not my_filter.passesFilter(IOC_var):
                    passes_filters = False
                    break
            if not passes_filters: continue

            # if variants already present, die hard with a vengeance (NOTE: used to simply skip the second instance of the variant and keep going
            if IOC_var in IOC_sample._contained_variants:
                raise RuntimeError('variant {0} found more than once in input dataset "{1}". Clean up your dataset before running SV segregation analysis'.format(IOC_var, self._sourceBatch))
            #variant is valid, adding it to sample
            if not DIC_fixed_raw_variants.has_key(IOC_sample):
                DIC_fixed_raw_variants[IOC_sample] = dict()
            if not DIC_fixed_raw_variants[IOC_sample].has_key(STR_varType):
                DIC_fixed_raw_variants[IOC_sample][STR_varType] = SortedList(key=data.BaseVariant.sort_by_start)

            DIC_fixed_raw_variants[IOC_sample][STR_varType].add(IOC_var)

        return DIC_fixed_raw_variants


    # def getVariants(self):
    #     return self._contained_variants


class XhmmFile(File):

    def __init__(self, STR_file, STR_sourceBatch=''):
        """
        Xhmm files handled by this class must have a header, and 15 tab-separated columns in this specific format:
        SAMPLE  CNV  INTERVAL  KB  CHR  MID_BP  TARGETS  NUM_TARG  Q_EXACT  Q_SOME  Q_NON_DIPLOID  Q_START  Q_STOP  MEAN_RD  MEAN_ORIG_RD
        XI001  DEL  chrX:20028957-20044146  15.19  chrX 20036551  730..738  9  4  74  77  5  14  -2.50  14.12

        """
        File.__init__(self, STR_file)
        self._sourceBatch = STR_sourceBatch

    def parse(self, SET_Filters):

        self._fstream.seek(0)
        # skip header line  
        self._fstream.readline()

        DIC_fixed_raw_variants = dict()

        # prepare to handle skipped samples
        SET_skipped = set()

        # parse values from input file (skip header)       
        for STR_line in self._fstream.readlines():
            LIS_line = STR_line.rstrip().split('\t')
            # sanity check: conifer format MUST have 5 tab-separated values
            if len(LIS_line) != 15: raise RuntimeError('incorrect number of columns in file {0} for line {1} (expected 15 (SAMPLE  CNV  INTERVAL  KB  CHR  MID_BP  TARGETS  NUM_TARG  Q_EXACT  Q_SOME  Q_NON_DIPLOID  Q_START  Q_STOP  MEAN_RD  MEAN_ORIG_RD ), got {2}'.format(self._file, STR_line.rstrip(), len(LIS_line)))   
            [STR_sample, STR_xhmm_type, STR_interval, STR_kb, STR_chr, STR_mid_bp, STR_targets, STR_num_targets, STR_qexact, STR_qsome, STR_non_diploid, STR_qstart, STR_qstop, STR_mean_rd, STR_mean_origrd] = LIS_line
            #removing the 'chr' in the chromosome id prefix
            STR_chr=STR_chr.replace('chr','')
            #getting the start and end from the interval
            STR_interval = STR_interval.split(':')[1]
            STR_interval = STR_interval.split('-')
            STR_start = STR_interval[0]
            STR_end = STR_interval[1]
            # little check to make sure the type fits with expected values
            if STR_xhmm_type not in XHMM_TYPES.keys(): raise RuntimeError('unrecognizable SV variant type (2nd column of input file {0}) for line {1}'.format(self._file, STR_line.rstrip()))
            STR_varType = XHMM_TYPES[STR_xhmm_type]

            IOC_sample = data.Sample.Get(STR_sample)
            # if sample not in data.Sample (i.e. not in pedFile loaded in segregation.py), skip this variant
            if IOC_sample is None:
                if STR_sample not in SET_skipped:
                    sys.stderr.write('\tWARNING: skipping variants from sample {0} (sample not found in ped file)\n'.format(STR_sample))
                    SET_skipped.add(STR_sample)
                continue
            # instantiate Variant with parsed data
            IOC_var = data.RawVariant(STR_chr, int(STR_start), STR_chr, int(STR_end), STR_varType, IOC_sample, xhmm, {}, self._sourceBatch )

            # if variant does not satisfy the filter rules, flush it                  
            passes_filters = True
            for my_filter in SET_Filters:
                if not my_filter.passesFilter(IOC_var):
                    passes_filters = False
                    break
            if not passes_filters: continue

            # if variants already present, die hard with a vengeance (NOTE: used to simply skip the second instance of the variant and keep going
            if IOC_var in IOC_sample._contained_variants:
                raise RuntimeError('variant {0} found more than once in input dataset "{1}". Clean up your dataset before running SV segregation analysis'.format(IOC_var, self._sourceBatch))
            #variant is valid, adding it to sample
            if not DIC_fixed_raw_variants.has_key(IOC_sample):
                DIC_fixed_raw_variants[IOC_sample] = dict()
            if not DIC_fixed_raw_variants[IOC_sample].has_key(STR_varType):
                DIC_fixed_raw_variants[IOC_sample][STR_varType] = SortedList(key=data.BaseVariant.sort_by_start)

            DIC_fixed_raw_variants[IOC_sample][STR_varType].add(IOC_var)

        return DIC_fixed_raw_variants


    # def getVariants(self):
    #     return self._contained_variants



class ExomeDepthFile(File):

    def __init__(self, STR_file, STR_sourceBatch=''):
        """
            ExomeDepth files handled by this class must have a header, and 13 tab-separated columns in this specific format:
            start.p  end.p  type  nexons  start  end  chromosome  id  BF  reads.expected  reads.observed  reads.ratio  sample
            819827  "deletion"  920028956  20044146  "chrX"  "chrchrX:20028956-20044146"  10.2  439  275  0.626  "XI001"

        """
        File.__init__(self, STR_file)
        self._sourceBatch = STR_sourceBatch

    def parse(self, SET_Filters):

        self._fstream.seek(0)
        # skip header line
        self._fstream.readline()

        DIC_fixed_raw_variants = dict()

        # prepare to handle skipped samples
        SET_skipped = set()

        # parse values from input file (skip header)       
        for STR_line in self._fstream.readlines():
            LIS_line = STR_line.rstrip().split('\t')
            # sanity check: conifer format MUST have 5 tab-separated values
            if len(LIS_line) != 13: raise RuntimeError('incorrect number of columns in file {0} for line {1} (expected 13 (start.p  end.p  type  nexons  start  end  chromosome  id  BF  reads.expected  reads.observed  reads.ratio  sample), got {2}'.format(self._file, STR_line.rstrip(), len(LIS_line)))   
            [STR_startp, STR_endp, STR_exomedepth_type, STR_nexons, STR_start, STR_end, STR_chr, STR_id, STR_BF, STR_reads_expected, STR_reads_obs, STR_reads_ratio, STR_sample ] = LIS_line
            #removing the "" from type, sample and chr 
            STR_exomedepth_type = STR_exomedepth_type[1:-1]
            STR_sample = STR_sample[1:-1]
            STR_chr = STR_chr[1:-1]
           
            #removing the 'chr' in the chromosome id prefix
            STR_chr=STR_chr.replace('chr','')
            # little check to make sure the type fits with expected values
            if STR_exomedepth_type not in EXOMEDEPTH_TYPES.keys(): raise RuntimeError('unrecognizable SV variant type (3rd column of input file {0}) for line {1}'.format(self._file, STR_line.rstrip()))
            STR_varType = EXOMEDEPTH_TYPES[STR_exomedepth_type]
            
            IOC_sample = data.Sample.Get(STR_sample)
            # if sample not in data.Sample (i.e. not in pedFile loaded in segregation.py), skip this variant
            if IOC_sample is None:
                if STR_sample not in SET_skipped:
                    sys.stderr.write('\tWARNING: skipping variants from sample {0} (sample not found in ped file)\n'.format(STR_sample))
                    SET_skipped.add(STR_sample)
                continue
            # instantiate Variant with parsed data
            IOC_var = data.RawVariant(STR_chr, int(STR_start), STR_chr, int(STR_end), STR_varType, IOC_sample, exomedepth, {}, self._sourceBatch )

            # if variant does not satisfy the filter rules, flush it                  
            passes_filters = True
            for my_filter in SET_Filters:
                if not my_filter.passesFilter(IOC_var):
                    passes_filters = False
                    break
            if not passes_filters: continue

            # if variants already present, die hard with a vengeance (NOTE: used to simply skip the second instance of the variant and keep going
            if IOC_var in IOC_sample._contained_variants:
                raise RuntimeError('variant {0} found more than once in input dataset "{1}". Clean up your dataset before running SV segregation analysis'.format(IOC_var, self._sourceBatch))
            #variant is valid, adding it to sample
            if not DIC_fixed_raw_variants.has_key(IOC_sample):
                DIC_fixed_raw_variants[IOC_sample] = dict()
            if not DIC_fixed_raw_variants[IOC_sample].has_key(STR_varType):
                DIC_fixed_raw_variants[IOC_sample][STR_varType] = SortedList(key=data.BaseVariant.sort_by_start)

            DIC_fixed_raw_variants[IOC_sample][STR_varType].add(IOC_var)

        return DIC_fixed_raw_variants


    # def getVariants(self):
    #     return self._contained_variants




#================================================

class CodexFile(File):

    def __init__(self, STR_file, STR_sourceBatch=''):
        """
            Codex files handled by this class must have a header, and 13 tab-separated columns in this specific format:
sample_name  chr  cnv  st_bp  ed_bp  length_kb  st_exon  ed_exon  raw_cov  norm_cov  copy_no  lratio  mBIC
XI001 chrX  del  54321013 54360197  39.185  2048 2054  186  347  1  44.584  23.019

        """
        File.__init__(self, STR_file)
        self._sourceBatch = STR_sourceBatch

    def parse(self, SET_Filters):

        self._fstream.seek(0)
        # skip header line
        self._fstream.readline()

        DIC_fixed_raw_variants = dict()

        # prepare to handle skipped samples
        SET_skipped = set()

        # parse values from input file (skip header)       
        for STR_line in self._fstream.readlines():
            LIS_line = STR_line.rstrip().split('\t')
            # sanity check: conifer format MUST have 5 tab-separated values
            if len(LIS_line) != 13: raise RuntimeError('incorrect number of columns in file {0} for line {1} (expected 13 (sample_name  chr  cnv  st_bp  ed_bp  length_kb  st_exon  ed_exon  raw_cov  norm_cov  copy_no  lratio  mBIC), got {2}'.format(self._file, STR_line.rstrip(), len(LIS_line)))   
            [STR_sample, STR_chr, STR_codex_type, STR_start, STR_end, STR_length, STR_start_exon, STR_end_exon, STR_raw_cov, STR_norm_cov, STR_copy_no, STR_lratio, STR_mBIC] = LIS_line
            #removing the 'chr' in the chromosome id prefix
            STR_chr=STR_chr.replace('chr','')
            # little check to make sure the type fits with expected values
            if STR_codex_type not in CODEX_TYPES.keys(): raise RuntimeError('unrecognizable SV variant type (3rd column of input file {0}) for line {1}'.format(self._file, STR_line.rstrip()))
            STR_varType = CODEX_TYPES[STR_codex_type]

            IOC_sample = data.Sample.Get(STR_sample)
            # if sample not in data.Sample (i.e. not in pedFile loaded in segregation.py), skip this variant
            if IOC_sample is None:
                if STR_sample not in SET_skipped:
                    sys.stderr.write('\tWARNING: skipping variants from sample {0} (sample not found in ped file)\n'.format(STR_sample))
                    SET_skipped.add(STR_sample)
                continue
            # instantiate Variant with parsed data
            IOC_var = data.RawVariant(STR_chr, int(STR_start), STR_chr, int(STR_end), STR_varType, IOC_sample, codex, {}, self._sourceBatch )

            # if variant does not satisfy the filter rules, flush it                  
            passes_filters = True
            for my_filter in SET_Filters:
                if not my_filter.passesFilter(IOC_var):
                    passes_filters = False
                    break
            if not passes_filters: continue

            # if variants already present, die hard with a vengeance (NOTE: used to simply skip the second instance of the variant and keep going
            if IOC_var in IOC_sample._contained_variants:
                raise RuntimeError('variant {0} found more than once in input dataset "{1}". Clean up your dataset before running SV segregation analysis'.format(IOC_var, self._sourceBatch))
            #variant is valid, adding it to sample
            if not DIC_fixed_raw_variants.has_key(IOC_sample):
                DIC_fixed_raw_variants[IOC_sample] = dict()
            if not DIC_fixed_raw_variants[IOC_sample].has_key(STR_varType):
                DIC_fixed_raw_variants[IOC_sample][STR_varType] = SortedList(key=data.BaseVariant.sort_by_start)

            DIC_fixed_raw_variants[IOC_sample][STR_varType].add(IOC_var)

        return DIC_fixed_raw_variants


    # def getVariants(self):
    #     return self._contained_variants



# class GenericBedFile(File):
    
#     def __init__(self, STR_file, STR_sourceBatch=''):
#         """
#             generic BED files handled by this class must have no header, and 5 tab-separated columns in this specific format:
#             chrom, start, end, SV_type, sample_id
#             NOTE: this class is not appropri
# ate for translocation variants, which must be provided in BEDPE format; 
#             a separate class will handle raw variant data files that include translocations
#         """
#         File.__init__(self, STR_file)
#         self._sourceBatch = STR_sourceBatch

#     def parse(self, SET_Filters):
        
#         DIC_fixed_raw_variants = dict()
        
#         # prepare to handle skipped samples
#         SET_skipped = set()
        
#         # parse values from input file (skip header)       
#         for STR_line in self._fstream.readlines():
#             LIS_line = STR_line.rstrip().split('\t')
#             # sanity check: genericBed format MUST have 5 tab-separated values
#             if len(LIS_line) != 5: raise RuntimeError('incorrect number of columns in file {0} for line {1} (expected 5 (chrom, start, end, SV_type, sample), got {2}'.format(self._file, STR_line.rstrip(), len(LIS_line)))
#             [STR_chr, STR_start, STR_end, STR_generic_type, STR_sample] = LIS_line
#             # will convert any non-standard SV type strings into usable format
#             STR_key = re.findall('(del|dup|inv|tra)', STR_generic_type)[0].upper()
#             # little check to make sure the type fits with expected values
#             if STR_key not in GENERIC_TYPES.keys(): raise RuntimeError('unrecognizable SV variant type (4th column of input file {0}) for line {1}'.format(self._file, STR_line.rstrip())) 
#             STR_varType = GENERIC_TYPES[STR_key]
            
#             IOC_sample = data.Sample.Get(STR_sample)
#             # if sample not in data.Sample (i.e. not in pedFile loaded in segregation.py), skip this variant
#             if IOC_sample is None:
#                 if STR_sample not in SET_skipped:
#                     sys.stderr.write('\tWARNING: skipping variants from sample {0} (sample not found in ped file)\n'.format(STR_sample))
#                     SET_skipped.add(STR_sample)
#                 continue
#             # instantiate Variant with parsed data
#             IOC_var = data.RawVariant(STR_chr, int(STR_start), STR_chr, int(STR_end), STR_varType, IOC_sample, genericBed, {}, self._sourceBatch )
            
#             # if variant does not satisfy the filter rules, flush it                  
#             passes_filters = True
#             for my_filter in SET_Filters:
#                 if not my_filter.passesFilter(IOC_var):
#                     passes_filters = False
#                     break
#             if not passes_filters: continue

#             # if variants already present, die hard with a vengeance (NOTE: used to simply skip the second instance of the variant and keep going
#             if IOC_var in IOC_sample._contained_variants:
#                 raise RuntimeError('variant {0} found more than once in input dataset "{1}". Clean up your dataset before running SV segregation analysis'.format(IOC_var, self._sourceBatch))
# #                 del IOC_var
# #                 continue
            
#             #variant is valid, adding it to sample
#             if not DIC_fixed_raw_variants.has_key(IOC_sample):
#                 DIC_fixed_raw_variants[IOC_sample] = dict()
#             if not DIC_fixed_raw_variants[IOC_sample].has_key(STR_varType):
#                 DIC_fixed_raw_variants[IOC_sample][STR_varType] = SortedList(key=data.BaseVariant.sort_by_start)

#             DIC_fixed_raw_variants[IOC_sample][STR_varType].add(IOC_var)
            
#         return DIC_fixed_raw_variants

    
#     def getVariants(self):
#         return self._contained_variants


#===============================================
class GenericBedpeFile(File):
    
    def __init__(self, STR_file, STR_sourceBatch):
        """
            generic BEDPE files handled by this class must have no header, and 8 tab-separated columns in this specific format:
            chrom_1, start_1, end_1, chrom_2, start_2, end_2, SV_type, sample_id
            NOTE: for translocations, the first set of 3 fields define one breakpoint, and the second set of 3 fields defines the other
            if the calling algorithm provides breakpoints as a single point (chr1, pos1, chr2, pos2), then this should be translated into bedpe format as
            chr1, pos1-1, pos1, chr2, pos2-1, pos2 
            this class can still process non-translocation variants; these must however be fit into the bedpe format;
            the conversion from standard bed (chrom, start, end) to bedpe should be: chrom, start, start+1, chrom, end-1, end
        """
        File.__init__(self, STR_file)
        self._sourceBatch = STR_sourceBatch

    def parse(self, SET_Filters):
        
        DIC_fixed_raw_variants = dict()
    
        # prepare to handle skipped samples
        SET_skipped = set()
        
        # parse values from input file (skip header)
        for STR_line in self._fstream.readlines():
            LIS_line = STR_line.rstrip().split('\t')
            # sanity check: genericBed format MUST have 5 tab-separated values
            if len(LIS_line) != 8: raise RuntimeError('incorrect number of columns in file {0} for line {1} (expected 8 (chrom1, start1, end1, chrom2, start2, end2, SV_type, sample), got {2}'.format(self._file, STR_line.rstrip(), len(LIS_line)))
            [STR_chrom1, STR_start1, STR_end1, STR_chrom2, STR_start2, STR_end2, STR_generic_type, STR_sample] = LIS_line
            # will convert any non-standard SV type strings into usable format
            STR_key = re.findall('(del|dup|inv|tra)', STR_generic_type.lower())[0].upper()
            # little check to make sure the type fits with expected values
            if STR_key not in GENERIC_TYPES.keys(): raise RuntimeError("""unrecognizable SV variant type (4th column of input file {file}) at line {line}""".format(file=self._file, line=STR_line.rstrip())) 
            STR_varType = GENERIC_TYPES[STR_key]

            IOC_sample = data.Sample.Get(STR_sample)
            # if sample not in data.Sample (i.e. not in pedFile loaded in segregation.py), skip this variant
            if IOC_sample is None:
                if STR_sample not in SET_skipped:
                    sys.stderr.write('\tWARNING: skipping variants from sample {0} (sample not found in ped file)\n'.format(STR_sample))
                    SET_skipped.add(STR_sample)
                continue
            # get single breakpoints from midpoint of bedpe entries - NOTE: use of int() will always round down any fractional midpoints
            INT_breakpoint1 = (int(STR_start1)+int(STR_end1))/2
            INT_breakpoint2 = (int(STR_start2)+int(STR_end2))/2
            # instantiate Variant with parsed data
            IOC_var = data.RawVariant(
                STR_chrom1, INT_breakpoint1, 
                STR_chrom2, INT_breakpoint2, 
                STR_varType, IOC_sample, 
                genericBedpe, {}, self._sourceBatch, 
                int(STR_start1), int(STR_end1), int(STR_start2), int(STR_end2)
            )
            
            # if variant does not satisfy the filter rules, flush it                  
            passes_filters = True
            for my_filter in SET_Filters:
                if not my_filter.passesFilter(IOC_var):
                    passes_filters = False
                    break
            if not passes_filters: continue

            # if variants already present, die hard with a vengeance (NOTE: used to simply skip the second instance of the variant and keep going
            if IOC_var in IOC_sample._contained_variants:
                raise RuntimeError('variant {0} found more than once in input dataset "{1}". Clean up your dataset before running SV segregation analysis'.format(IOC_var, self._sourceBatch))
            
            #variant is valid, adding it to sample
            IOC_sample.addVariantToSample(IOC_var)
            if not DIC_fixed_raw_variants.has_key(IOC_sample):
                DIC_fixed_raw_variants[IOC_sample] = dict()
            if not DIC_fixed_raw_variants[IOC_sample].has_key(STR_varType):
                DIC_fixed_raw_variants[IOC_sample][STR_varType] = SortedList(key=data.BaseVariant.sort_by_start)

            DIC_fixed_raw_variants[IOC_sample][STR_varType].add(IOC_var)
            
        return DIC_fixed_raw_variants

 
class AnnotationFile(File):

    def __init__(self, STR_file, STR_sourceType):
        File.__init__(self, STR_file)
        self._sourceType = STR_sourceType
    
#    def parseBed(self):
#        """ a function to parse generic 4-column annotation BED file """
#        LIS_annotations = []
#        self._fstream.seek(0)
#        for STR_line in self._fstream.readlines()[0:]:
#            [   STR_chr,  STR_start, STR_end, 
#                STR_annotation
#            ] = STR_line.rstrip().split('\t')
#            LIS_annotations.append(data.Annotation(self._sourceType, STR_chr, int(STR_start), STR_chr, int(STR_end), STR_annotation, self._minRO))
#        return LIS_annotations

    def parseBed(self):
        """ a function to parse generic 4-column annotation BED file """
        DIC_annotations = dict()
        self._fstream.seek(0)
        for STR_line in self._fstream.readlines()[0:]:
            [   STR_chr,  STR_start, STR_end, 
                STR_annotation
            ] = STR_line.rstrip().split('\t')
            if not STR_chr in DIC_annotations: DIC_annotations[STR_chr] = IntervalTree()
            new_interval = Interval(int(STR_start), int(STR_end), data.Annotation(self._sourceType, STR_chr, int(STR_start), STR_chr, int(STR_end), STR_annotation))
            DIC_annotations[STR_chr].add(new_interval)
        return DIC_annotations

    
    # def parseBedPE(self):
    #     """ a function to parse generic 7-column annotation BEDPE file """
    #     LIS_annotations = []
    #     self._fstream.seek(0)
    #     for STR_line in self._fstream.readlines()[0:]:
    #         [   STR_chrom1,  STR_start1, STR_end1,
    #             STR_chrom2,  STR_start2, STR_end2, 
    #             STR_annotation
    #         ] = STR_line.rstrip().split('\t')
    #         LIS_annotations.append(data.Annotation(self._sourceType, STR_chrom1, int(STR_start1), STR_chrom2, int(STR_end2), STR_annotation))
    #     return LIS_annotations

def flatten(l, ltypes=(list, tuple)):
    """ generic function to flatten a list of lists into a single list """
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def findOverlappingLeftMost(variant, leftSortedList):
    """ helper function which is used (in conjunction with findOverlappingRightMost) 
        to grow a list of overlapping variants starting from a single reference variant """
    currentIndex=leftSortedList.bisect_left(variant)
    currentVariant=leftSortedList[currentIndex]
    testIndex=currentIndex-1
    testVariant=leftSortedList[testIndex]
    while (testVariant._pos2 >= currentVariant._pos1) & (testIndex >= 0):
        currentIndex -= 1
        currentVariant=leftSortedList[currentIndex]
        testIndex -= 1
        testVariant=leftSortedList[testIndex]
    return leftSortedList[currentIndex]
    
def findOverlappingRightMost(variant, rightSortedList):
    """ helper function which is used (in conjunction with findOverlappingLeftMost) 
        to grow a list of overlapping variants starting from a single reference variant """
    currentIndex=rightSortedList.bisect_left(variant)
    currentVariant=rightSortedList[currentIndex]
    testIndex=currentIndex-1
    testVariant=rightSortedList[testIndex] 
    while (testVariant._pos1 <= currentVariant._pos2) & (testIndex >= 0):
        currentIndex -= 1
        currentVariant=rightSortedList[currentIndex]
        testIndex -= 1
        testVariant=rightSortedList[testIndex]
    return rightSortedList[currentIndex]

def findSingleLeftMost(variant, leftSortedList):
    """ helper function which is used (in conjunction with findSingleRightMost)
        to find all variants that overlap a reference variant """
    testIndex=leftSortedList.bisect_left(variant)
    testVariant=leftSortedList[testIndex]
    while (testVariant._pos2 >= variant._pos1) & (testIndex > 0):
        testIndex -= 1
        testVariant=leftSortedList[testIndex]
    return leftSortedList[testIndex]

def findSingleRightMost(variant, rightSortedList):
    """ helper function which is used (in conjunction with findSingleLeftMost)
        to find all variants that overlap a reference variant """
    testIndex=rightSortedList.bisect_left(variant)
    testVariant=rightSortedList[testIndex] 
    while (testVariant._pos1 <= variant._pos2) & (testIndex > 0):
        testIndex -= 1
        testVariant=rightSortedList[testIndex]
    return rightSortedList[testIndex]

def getVariantsFromBorders(leftMost, rightMost, sortedList):
    return sortedList[sortedList.index(leftMost):sortedList.index(rightMost)+1]



