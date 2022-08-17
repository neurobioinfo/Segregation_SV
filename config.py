import sys
import re
import os
import configparser
import argparse
from helpers import flatten
from constants import *
from datetime import *

class Config(object):
    
    # fundamental methods
    def __init__(self, LIS_inputArgs):
        # order of priority: command-line args > custom .cfg file > default .ini file
        # step 1: parse input args with argparse: 
        # populate name/file attributes first, preserve any override arguments (for later) to local TMP_ variables
        # "TMP_" prefix to local attributes to make them easier to spot in the code
        # TMP variables will be properly set as attributes only after parsing ini + cfg files
        self._inputArgs = LIS_inputArgs
        (
            self.pedFile, self.dataFiles, self.dataSources, self.dataNames, self.runName, self.cfgFile, self.outFile, self.customAnnotationFiles, self.customAnnotationNames,self.customAnnotationMinROs, self.excludeDefaultAnnotations, TMP_INT_rawDataFilters_minVariantSize, TMP_INT_rawDataFilters_maxVariantSize, TMP_STR_otherRulesAndParameters_emitStatus, self.verbose
        ) = self.parseInputArgs()
        # sanity check: make sure all input files exist!
        for infile in self.getInputFiles():
            if not os.path.isfile(infile): raise RuntimeError("input file {0} does not exist\n".format(infile))
        # step 2: parse default .ini file with ConfigParser to set all run parameters to default values
        CParser = configparser.SafeConfigParser()
        CParser.optionxform = str   # preserve case of variables to be parsed from cfg_file 
        CParser.read(FILE_defaultIni)
        try:
            self.FLO_minRO_cleanupRaw_Deletion = CParser.getfloat('DEFAULT', 'FLO_minRO_cleanupRaw_Deletion')
            self.FLO_minRO_cleanupRaw_Duplication = CParser.getfloat('DEFAULT', 'FLO_minRO_cleanupRaw_Duplication')
            self.FLO_minRO_cleanupRaw_Inversion = CParser.getfloat('DEFAULT', 'FLO_minRO_cleanupRaw_Inversion')
            self.FLO_minRO_cleanupRaw_Translocation = CParser.getfloat('DEFAULT', 'FLO_minRO_cleanupRaw_Translocation')
            self.FLO_minRO_mergeCallsets_Deletion = CParser.getfloat('DEFAULT', 'FLO_minRO_mergeCallsets_Deletion')
            self.FLO_minRO_mergeCallsets_Duplication = CParser.getfloat('DEFAULT', 'FLO_minRO_mergeCallsets_Duplication')
            self.FLO_minRO_mergeCallsets_Inversion = CParser.getfloat('DEFAULT', 'FLO_minRO_mergeCallsets_Inversion')
            self.FLO_minRO_mergeCallsets_Translocation = CParser.getfloat('DEFAULT', 'FLO_minRO_mergeCallsets_Translocation')
            self.FLO_minRO_segregation_Deletion = CParser.getfloat('DEFAULT', 'FLO_minRO_segregation_Deletion')
            self.FLO_minRO_segregation_Duplication = CParser.getfloat('DEFAULT', 'FLO_minRO_segregation_Duplication')
            self.FLO_minRO_segregation_Inversion = CParser.getfloat('DEFAULT', 'FLO_minRO_segregation_Inversion')
            self.FLO_minRO_segregation_Translocation = CParser.getfloat('DEFAULT', 'FLO_minRO_segregation_Translocation')
            # self.FLO_minRO_annotation_DGV_inclusive = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_DGV_inclusive')
            # self.FLO_minRO_annotation_DGV_stringent = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_DGV_stringent')
            # self.FLO_minRO_annotation_KG_SV = CParser.getfloat('DEFAULT','FLO_minRO_annotation_KG_SV')
            # self.FLO_minRO_annotation_Segmental_Dups = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_Segmental_Dups')
            # self.FLO_minRO_annotation_micro_exons = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_micro_exons')
            # self.FLO_minRO_annotation_repeatMasker = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_repeatMasker')
            # self.FLO_minRO_annotation_RefSeq_genes = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_RefSeq_genes')
            # self.FLO_minRO_annotation_RefSeq_genes_exons = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_RefSeq_genes_exons')
            self.INT_rawDataFilters_minVariantSize = CParser.getint('DEFAULT', 'INT_rawDataFilters_minVariantSize')
            self.INT_rawDataFilters_maxVariantSize = CParser.getint('DEFAULT', 'INT_rawDataFilters_maxVariantSize')
            self.STR_otherRulesAndParameters_emitStatus = CParser.get('DEFAULT', 'STR_otherRulesAndParameters_emitStatus')
        except Exception as e:
            sys.stderr.write('parsing problem:\n')
            print (type(e), e)
            sys.exit(1)
        # step 3: parse .cfg file with ConfigParser to modify any run parameters
        if self.cfgFile:
            CParser = configparser.SafeConfigParser()
            CParser.optionxform = str   # preserve case of variables to be parsed from cfg_file 
            CParser.read(self.cfgFile)
            try:
                if CParser.has_option('DEFAULT', 'FLO_minRO_cleanupRaw_Deletion'): self.FLO_minRO_cleanupRaw_Deletion = CParser.getfloat('DEFAULT', 'FLO_minRO_cleanupRaw_Deletion')
                if CParser.has_option('DEFAULT', 'FLO_minRO_cleanupRaw_Duplication'): self.FLO_minRO_cleanupRaw_Duplication = CParser.getfloat('DEFAULT', 'FLO_minRO_cleanupRaw_Duplication')
                if CParser.has_option('DEFAULT', 'FLO_minRO_cleanupRaw_Inversion'): self.FLO_minRO_cleanupRaw_Inversion = CParser.getfloat('DEFAULT', 'FLO_minRO_cleanupRaw_Inversion')
                if CParser.has_option('DEFAULT', 'FLO_minRO_cleanupRaw_Translocation'): self.FLO_minRO_cleanupRaw_Translocation = CParser.getfloat('DEFAULT', 'FLO_minRO_cleanupRaw_Translocation')
                if CParser.has_option('DEFAULT', 'FLO_minRO_mergeCallsets_Deletion'): self.FLO_minRO_mergeCallsets_Deletion = CParser.getfloat('DEFAULT', 'FLO_minRO_mergeCallsets_Deletion')
                if CParser.has_option('DEFAULT', 'FLO_minRO_mergeCallsets_Duplication'): self.FLO_minRO_mergeCallsets_Duplication = CParser.getfloat('DEFAULT', 'FLO_minRO_mergeCallsets_Duplication')
                if CParser.has_option('DEFAULT', 'FLO_minRO_mergeCallsets_Inversion'): self.FLO_minRO_mergeCallsets_Inversion = CParser.getfloat('DEFAULT', 'FLO_minRO_mergeCallsets_Inversion')
                if CParser.has_option('DEFAULT', 'FLO_minRO_mergeCallsets_Translocation'): self.FLO_minRO_mergeCallsets_Translocation = CParser.getfloat('DEFAULT', 'FLO_minRO_mergeCallsets_Translocation')
                if CParser.has_option('DEFAULT', 'FLO_minRO_segregation_Deletion'): self.FLO_minRO_segregation_Deletion = CParser.getfloat('DEFAULT', 'FLO_minRO_segregation_Deletion')
                if CParser.has_option('DEFAULT', 'FLO_minRO_segregation_Duplication'): self.FLO_minRO_segregation_Duplication = CParser.getfloat('DEFAULT', 'FLO_minRO_segregation_Duplication')
                if CParser.has_option('DEFAULT', 'FLO_minRO_segregation_Inversion'): self.FLO_minRO_segregation_Inversion = CParser.getfloat('DEFAULT', 'FLO_minRO_segregation_Inversion')
                if CParser.has_option('DEFAULT', 'FLO_minRO_segregation_Translocation'): self.FLO_minRO_segregation_Translocation = CParser.getfloat('DEFAULT', 'FLO_minRO_segregation_Translocation')
                # if CParser.has_option('DEFAULT', 'FLO_minRO_annotation_DGV_inclusive'): self.FLO_minRO_annotation_DGV_inclusive = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_DGV_inclusive')
                # if CParser.has_option('DEFAULT', 'FLO_minRO_annotation_DGV_stringent'): self.FLO_minRO_annotation_DGV_stringent = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_DGV_stringent')
                # if CParser.has_option('DEFAULT','FLO_minRO_annotation_KG_SV'): self.FLO_minRO_annotation_KG_SV = CParser.getfloat('DEFAULT','FLO_minRO_annotation_KG_SV')
                # if CParser.has_option('DEFAULT', 'FLO_minRO_annotation_Segmental_Dups'): self.FLO_minRO_annotation_Segmental_Dups = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_Segmental_Dups')
                # if CParser.has_option('DEFAULT', 'FLO_minRO_annotation_micro_exons'): self.FLO_minRO_annotation_micro_exons = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_micro_exons')
                # if CParser.has_option('DEFAULT', 'FLO_minRO_annotation_repeatMasker'): self.FLO_minRO_annotation_repeatMasker = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_repeatMasker')
                # if CParser.has_option('DEFAULT', 'FLO_minRO_annotation_RefSeq_genes'): self.FLO_minRO_annotation_RefSeq_genes = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_RefSeq_genes')
                # if CParser.has_option('DEFAULT', 'FLO_minRO_annotation_RefSeq_genes_exons'): self.FLO_minRO_annotation_RefSeq_genes_exons = CParser.getfloat('DEFAULT', 'FLO_minRO_annotation_RefSeq_genes_exons')
                if CParser.has_option('DEFAULT', 'INT_rawDataFilters_minVariantSize'): self.INT_rawDataFilters_minVariantSize = CParser.getint('DEFAULT', 'INT_rawDataFilters_minVariantSize')
                if CParser.has_option('DEFAULT', 'INT_rawDataFilters_maxVariantSize'): self.INT_rawDataFilters_maxVariantSize = CParser.getint('DEFAULT', 'INT_rawDataFilters_maxVariantSize')
                if CParser.has_option('DEFAULT', 'STR_otherRulesAndParameters_emitStatus'): self.STR_otherRulesAndParameters_emitStatus = CParser.get('DEFAULT', 'STR_otherRulesAndParameters_emitStatus')
            except Exception as e:
                sys.stderr.write('parsing problem:\n')
                print (type(e), e)
                sys.exit(1)
        # step 4: update any override options that were provided with script flags/arguments (stored locally as TMP_ variables)
        if TMP_INT_rawDataFilters_minVariantSize: self.INT_rawDataFilters_minVariantSize = TMP_INT_rawDataFilters_minVariantSize
        if TMP_INT_rawDataFilters_maxVariantSize: self.INT_rawDataFilters_maxVariantSize = TMP_INT_rawDataFilters_maxVariantSize 
        if TMP_STR_otherRulesAndParameters_emitStatus: self.STR_otherRulesAndParameters_emitStatus = TMP_STR_otherRulesAndParameters_emitStatus
        self.outFH = self.outFile if self.outFile == sys.stdout else open(self.outFile, 'w')
        
    def __repr__(self):
        return '\n'.join([':'.join(map(str, [attr, self.__dict__[attr]])) for attr in sorted(self.__dict__.keys())])
    
    def __str__(self):
        return self.__repr__()
    
    def parseInputArgs(self):
        # script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
        parser = argparse.ArgumentParser(description='run segregation analysis on structural variants dataset')
        parser.add_argument('-i', '--input_data', dest='inputData', action='append', help='/path/to/input_data_file:input_type:[batch_name] (can be specified multiple times; accepted values are "lumpy", "popsv", "conifer", "codex", "exomedepth", "cnmops", "xhmm", "generic_bed" and "generic_bedpe")', required=True)
        parser.add_argument('-p', '--ped_file',   dest='pedFile',  help='/path/to/ped_file', required=True)
        parser.add_argument('-r', '--run_name', dest='runName', help='run name prefix for all output files', default='SV_segregation')
        parser.add_argument('-a','--custom_annotation_file',dest='customAnnoData', action='append', help='/path/to/custom_annotation_file:anno_type:minimum_reciprocal_overlap_fraction ;  (can be specified multiple times; file must be .bed or .bedpe format)', default=[])
        parser.add_argument('-o', '--out_file',   dest='outFile',  help='/path/to/output/file; where outputs will be written', default=sys.stdout)
        parser.add_argument('-cfg', '--config_file', dest='cfgFile', help='path/to/cfg_file (used to modify run parameters)', default=None)
        parser.add_argument('-eda', '--exclude_default_annotations', dest='excludeDefaultAnnotations', action='store_true', help='exclude default annotations', default=False)
        parser.add_argument('-min', '--min_size_filter', dest='minSizeFilter', help='minimal size of raw variant to keep (smaller will be excluded from any processing)', default=None)
        parser.add_argument('-max', '--max_size_filter', dest='maxSizeFilter', help='maximal size of raw variant to keep (larger will be excluded from any processing', default=None)  
        parser.add_argument('-e', '--emit_vars_from_samples_with_status', dest='emitStatus', help="emit variants from samples with clinical status of [" + affected + "," + nonaffected + "," + allsamples + "]", default=None )
        parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='verbose output', default=False)
        args = parser.parse_args(self._inputArgs)
        # parse inputData to identify file(s), type(s) and batchName(s), and raise errors if not specified properly


#self.pedFile, self.dataFiles, self.dataSources, self.dataNames, self.runName, self.cfgFile, self.outFile, self.customAnnotationFiles, self.customAnnotationNames,self.customAnnotationMinROs, self.excludeDefaultAnnotations, 

        dataFiles, dataSources, dataNames = ([],[],[])
        for STR_input in args.inputData:
            LIS_input = STR_input.split(':')
            dataFile = LIS_input[0]
            # die if data type not specified
            if len(LIS_input) < 2: raise RuntimeError("input error: data file {0} has no associated input type\n".format(dataFile))
            # die if data type not among allowed options
            dataSource = LIS_input[1].lower()
            if dataSource not in DATA_PARSERS: raise RuntimeError("input error: specified data type ({1}) for data file {0} not among accepted options ({2})\n".format(dataFile, dataSource, DATA_PARSERS))
            # if batch name not given, default to data type
            dataName = LIS_input[2] if len(LIS_input)==3 else dataSource
            dataFiles.append(dataFile)
            dataSources.append(dataSource)
            dataNames.append(dataName)
        # make sure no batch names are duplicated
        if len(dataNames) != len(set(dataNames)): raise RuntimeError('cannot have >1 batch with the same name. Make sure each batch has a unique name')
        customAnnotationFiles, customAnnotationNames, customAnnotationMinROs = ([],[], [])
        for STR_input in args.customAnnoData:
            LIS_input = STR_input.split(':')
            # if any customAnnoType or customAnnoMinRO are missing or improperly formatted, die; otherwise populate variables and move on
            try:
                customAnnotationFile, customAnnotationName, customAnnotationMinRO = LIS_input
            except: raise RuntimeError('custom annotations MUST be supplied with an anno_type and a minimum_reciprocal_overlap_fraction value (see help for details). You have provided {0} of those values for customAnnoFile "{1}": {2}'.format(len(LIS_input[1:]), LIS_input[0], LIS_input[1:]))
            if not customAnnotationName: raise RuntimeError('please specify anno_type for custom annotation {0} (see help for details)'.format(customAnnotationFile))
            if not customAnnotationMinRO: raise RuntimeError('please specify minimum reciprocal overlap for custom annotation {0} (see help for details'.format(customAnnotationFile))
            customAnnotationFiles.append(customAnnotationFile)
            customAnnotationNames.append(customAnnotationName)
            try: 
                customAnnotationMinROs.append(float(customAnnotationMinRO))
            except: 
                raise RuntimeError('minimum reciprocal overlap for custom annotation {0} must be an integer or floating value. (you gave a "{1}")'.format(customAnnotationFile, customAnnotationMinRO))
        runName = args.runName
        pedFile, cfgFile, outFile, excludeDefaultAnnotations = (args.pedFile, args.cfgFile, args.outFile, args.excludeDefaultAnnotations)
        minSizeFilter = int(args.minSizeFilter) if args.minSizeFilter else None 
        maxSizeFilter = int(args.maxSizeFilter) if args.maxSizeFilter else None
        emitStatus, verbose = (args.emitStatus, args.verbose)
        return [
            pedFile, dataFiles, dataSources, dataNames,
            runName, cfgFile, outFile,
            customAnnotationFiles, customAnnotationNames, customAnnotationMinROs, excludeDefaultAnnotations,
            minSizeFilter, maxSizeFilter,
            emitStatus, verbose
        ]
    
    def printConfig(self):
        with open(self.runName + '.' + str(datetime.now().isoformat()) + '.setup', 'w') as f:
            pat = re.compile("^(STR|INT|FLO)_")
            DIC_config = {attr:self.__dict__[attr] for attr in self.__dict__ if re.match(pat, attr)}
            DIC_switch = {attr:self.__dict__[attr] for attr in self.__dict__ if attr not in DIC_config}
            STR_config = '\n'.join([':'.join(map(str, [k, DIC_config[k]])) for k in sorted(DIC_config.keys())])
            STR_switch = '\n'.join([':'.join(map(str, [k, DIC_switch[k]])) for k in sorted(DIC_switch.keys())])
            f.write('[settings from .ini and/or .cfg]')
            f.write(STR_config + '\n')
            f.write('[settings from segregation.py args]')
            f.write(STR_switch)
            # print >> f, '[settings from .ini and/or .cfg]'
            # print >> f, STR_config + '\n'
            # print >> f, '[settings from segregation.py args]'
            # print >> f, STR_switch
                
    def getInputFiles(self):
        return [File for File in flatten([self.pedFile, self.dataFiles, self.cfgFile, self.customAnnotationFiles]) if File is not None]
    
    def getOutputFH(self):
        return self.outFH
    
    def getOutputDir(self):
        return '.' if self.outFile == sys.stdout or not self.outFile else os.path.dirname(os.path.realpath(self.outFile))
