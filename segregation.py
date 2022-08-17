#!/usr/bin/env python
import sys
import os
import argparse
if (sys.version_info > (3, 0)):
    import pickle as pickle
else:
    import cPickle as pickle

import helpers
from data import *
from constants import *
from intervaltree import Interval, IntervalTree
import config
from sortedcontainers import SortedSet, SortedList
import logging
import random
import threading
import time

if __name__ == "__main__":

    # STEP 1: process script inputs using config.Config class and print to .setup file
    sys.stderr.write("1. Parsing input arguments\n")
    Config = config.Config(LIS_inputArgs=sys.argv[1:])
    Config.printConfig()
    # disable in-depth stacktrace for error reporting if run is not verbose
    sys.tracebacklimit = None if Config.verbose else 0
    # NOTE: creating Config above will update config.py's global variables with parsed values; 
    # open the outfile (if it was given as argument) and return a file handle, which will be used for all subsequent printing
    outFH = Config.getOutputFH()
    
    SET_Filter = set([helpers.MaxSizeFilter(Config.INT_rawDataFilters_maxVariantSize), helpers.MinSizeFilter(Config.INT_rawDataFilters_minVariantSize)])
    
    # STEP 2: read input files and perform primary data parsing
    # read pedigree file
    ped = helpers.PedigreeFile(Config.pedFile)
    ped.parse()
    # read + process raw variants
    DIC_consensusVariants = dict()
    for i in range(0,len(Config.dataFiles)):
        sys.stderr.write('2. Processing raw variants from input file ' + Config.dataFiles[i] + '\n')
        dataFile, dataSource, dataName = (Config.dataFiles[i], Config.dataSources[i], Config.dataNames[i])
        dataFileParser = None
        if   dataSource == popsv: dataFileParser = helpers.PopSVFile(dataFile, dataName)
        elif dataSource == lumpy: dataFileParser = helpers.LumpyFile(dataFile, dataName)
        elif dataSource == conifer: dataFileParser = helpers.ConiferFile(dataFile, dataName)
        elif dataSource == genericBed: dataFileParser = helpers.GenericBedFile(dataFile, dataName)
        elif dataSource == genericBedpe: dataFileParser = helpers.GenericBedpeFile(dataFile, dataName)
        elif dataSource == codex: dataFileParser = helpers.CodexFile(dataFile, dataName)
        elif dataSource == exomedepth: dataFileParser = helpers.ExomeDepthFile(dataFile, dataName)
        elif dataSource == xhmm: dataFileParser = helpers.XhmmFile(dataFile, dataName)
        elif dataSource == cnmops: dataFileParser = helpers.CnmopsFile(dataFile, dataName)
        DIC_variantsInDataSource = dataFileParser.parse(SET_Filter)
        
        # STEP 3 : merge calls per type per sample from each raw callset (in case callset gives overlapping or fragmented calls)
        sys.stderr.write('3. Creating consensus for input file ' + Config.dataFiles[i] + '\n')
        for chrom in DIC_variantsInDataSource:
            #initialise the dict
            if chrom not in DIC_consensusVariants: DIC_consensusVariants[chrom] = dict()
            for IOC_sample in DIC_variantsInDataSource[chrom]:
                #initialise the dict
                if IOC_sample not in DIC_consensusVariants[chrom]: DIC_consensusVariants[chrom][IOC_sample] = dict()                             
                for STR_varType in DIC_variantsInDataSource[chrom][IOC_sample]:  
                    #initialise the dict
                    if STR_varType not in DIC_consensusVariants[chrom][IOC_sample]: 
                        DIC_consensusVariants[chrom][IOC_sample][STR_varType] = IntervalTree()

                    if STR_varType == Deletion : FLO_minRO = Config.FLO_minRO_cleanupRaw_Deletion
                    elif STR_varType == Duplication : FLO_minRO = Config.FLO_minRO_cleanupRaw_Duplication
                    elif STR_varType == Inversion : FLO_minRO = Config.FLO_minRO_cleanupRaw_Inversion
                    elif STR_varType == Translocation : 
                        DIC_consensusVariants[chrom][IOC_sample][STR_varType].update(DIC_variantsInDataSource[chrom][IOC_sample][STR_varType])        
                        continue
                    
                    while ( DIC_variantsInDataSource[chrom][IOC_sample][STR_varType] ) :
                        IOC_smallest = SortedList(DIC_variantsInDataSource[chrom][IOC_sample][STR_varType], key=Interval.length)[0]
                        DIC_variantsInDataSource[chrom][IOC_sample][STR_varType].remove(IOC_smallest)                    
                        LIS_consensusVariants = [IOC_smallest.data]
                        BOOL_mergedIntervals = False
                        for IOC_compareInterval in SortedList(DIC_variantsInDataSource[chrom][IOC_sample][STR_varType][IOC_smallest], key=Interval.length):
                            if (IOC_smallest.data.getPairwiseReciprocalOverlap(IOC_compareInterval.data) > FLO_minRO): 
                                #merging two intervals                            
                                DIC_variantsInDataSource[chrom][IOC_sample][STR_varType].remove(IOC_compareInterval)
                                LIS_consensusVariants.append(IOC_compareInterval.data)
                                BOOL_mergedIntervals = True
                        if BOOL_mergedIntervals: #merging variants
                            IOC_callerConsensus = ConsensusVariant (LIS_consensusVariants)
                            DIC_variantsInDataSource[chrom][IOC_sample][STR_varType][IOC_callerConsensus._pos1:IOC_callerConsensus._pos2] = IOC_callerConsensus 
                        else: #nothing to merge                                
                            DIC_consensusVariants[chrom][IOC_sample][STR_varType].add(IOC_smallest)

    # STEP 4: merge filtered calls from multiple callsets
    sys.stderr.write("4. Merging calls from multiple callsets\n")    
    DIC_mergedVariants = dict()
    for chrom in DIC_consensusVariants:
        if chrom not in DIC_mergedVariants: DIC_mergedVariants[chrom] = dict()
        for IOC_sample in DIC_consensusVariants[chrom]:
            for STR_varType in DIC_consensusVariants[chrom][IOC_sample]:
                if STR_varType not in DIC_mergedVariants[chrom]: DIC_mergedVariants[chrom][STR_varType] = IntervalTree()                     
                if STR_varType == Deletion : FLO_minRO = Config.FLO_minRO_mergeCallsets_Deletion
                elif STR_varType == Duplication : FLO_minRO = Config.FLO_minRO_mergeCallsets_Duplication
                elif STR_varType == Inversion : FLO_minRO = Config.FLO_minRO_mergeCallsets_Inversion
                elif STR_varType == Translocation : 
                    DIC_mergedVariants[chrom][STR_varType] = DIC_consensusVariants[chrom][IOC_sample][STR_varType]
                    continue
                
                while ( DIC_consensusVariants[chrom][IOC_sample][STR_varType] ) :
                    IOC_smallest = SortedList(DIC_consensusVariants[chrom][IOC_sample][STR_varType], key=Interval.length)[0]
                    DIC_consensusVariants[chrom][IOC_sample][STR_varType].remove(IOC_smallest)
                    BOOL_mergedIntervals = False
                    LIS_mergedVariants = [IOC_smallest.data]
                    for IOC_compareInterval in SortedList(DIC_consensusVariants[chrom][IOC_sample][STR_varType][IOC_smallest], key=Interval.length):
                        if (IOC_smallest.data.getPairwiseReciprocalOverlap(IOC_compareInterval.data) > FLO_minRO): 
                            #merging two intervals                              
                            LIS_mergedVariants.append(IOC_compareInterval.data)
                            DIC_consensusVariants[chrom][IOC_sample][STR_varType].remove(IOC_compareInterval)
                            BOOL_mergedIntervals = True                   
                    if  BOOL_mergedIntervals: #merging
                        IOC_callerConsensus = MergedVariant (LIS_mergedVariants)
                        DIC_consensusVariants[chrom][IOC_sample][STR_varType][IOC_callerConsensus._pos1:IOC_callerConsensus._pos2] = IOC_callerConsensus      
                    else: #nothing to merge                         
                        DIC_mergedVariants[chrom][STR_varType].add(IOC_smallest)

    # STEP 5: segregation in pedigrees
    sys.stderr.write("5. Segregating merged calls\n")

    LIS_segregated_variants = SortedList(key=BaseVariant.sort_by_start)          
    for chrom in DIC_mergedVariants:
        for STR_varType in DIC_mergedVariants[chrom]:
            if STR_varType == Deletion : FLO_minRO = Config.FLO_minRO_segregation_Deletion
            elif STR_varType == Duplication : FLO_minRO = Config.FLO_minRO_segregation_Duplication
            elif STR_varType == Inversion : FLO_minRO = Config.FLO_minRO_segregation_Inversion
            elif STR_varType == Translocation : 
                FLO_minRO = Config.FLO_minRO_segregation_Translocation
                for INTERVAL_var in DIC_mergedVariants[chrom][STR_varType]:
                    LIS_segregated_variants.add(SegregatedVariant([INTERVAL_var.data]))
                continue        
            
            while ( DIC_mergedVariants[chrom][STR_varType] ) : 
                IOC_smallest = SortedList(DIC_mergedVariants[chrom][STR_varType], key=Interval.length)[0]
                DIC_mergedVariants[chrom][STR_varType].remove(IOC_smallest)
                BOOL_mergedIntervals = False
                LIS_localSegregatedVariants = [IOC_smallest.data]
                for IOC_compareInterval in SortedList(DIC_mergedVariants[chrom][STR_varType][IOC_smallest], key=Interval.length):
                    if (IOC_smallest.data.getPairwiseReciprocalOverlap(IOC_compareInterval.data) > FLO_minRO): 
                        DIC_mergedVariants[chrom][STR_varType].remove(IOC_compareInterval)
                        BOOL_mergedIntervals = True           
                        LIS_localSegregatedVariants.append(IOC_compareInterval.data)
                if BOOL_mergedIntervals: #merging
                    IOC_callerConsensus = SegregatedVariant (LIS_localSegregatedVariants)
                    DIC_mergedVariants[chrom][STR_varType][IOC_callerConsensus._pos1:IOC_callerConsensus._pos2] = IOC_callerConsensus  
                else: #nothing to merge            
                    if not isinstance(IOC_smallest.data, SegregatedVariant):
                        LIS_segregated_variants.add(SegregatedVariant([IOC_smallest.data]))
                    else: 
                        LIS_segregated_variants.add(IOC_smallest.data)

    # STEP 6: annotation with external databases
    sys.stderr.write('6. Annotate variants\n')
    DIC_annotationMinROs = {}
    DIC_annotations = {}
    LIS_annotationNames = []
    
    if not Config.excludeDefaultAnnotations:
        TREE_defaultAnnotations, LIS_defaultAnnotationNames, DIC_annotationMinROs = pickle.load(open(FILE_resources_pickle, 'rb'))
        LIS_annotationNames.extend(LIS_defaultAnnotationNames)
        for chrom in TREE_defaultAnnotations:            
            if not chrom in DIC_annotations: DIC_annotations[chrom] = TREE_defaultAnnotations[chrom]
            else: DIC_annotations[chrom] |= TREE_defaultAnnotations[chrom]

    if len(Config.customAnnotationFiles) > 0:
        for i in range(0,len(Config.customAnnotationFiles)):
            STR_annoFile = Config.customAnnotationFiles[i]
            STR_annotationName = Config.customAnnotationNames[i]
            FLO_annotationMinRO = Config.customAnnotationMinROs[i]
            # sanity check: if custom annotation is already one of the default types, raise error asking user to rename his custom anno
            if STR_annotationName in LIS_annotationNames: 
                raise RuntimeError('custom annotation "{0}" has the same anno_type as a default annotation. Please use a different anno_type for custom annotation file {1}'.format(STR_annotationName, STR_annoFile))
            DIC_annotationMinROs.update( {STR_annotationName : float(FLO_annotationMinRO)} )
            LIS_annotationNames.append(STR_annotationName)
            IOC_annoFile = helpers.AnnotationFile(STR_annoFile, STR_annotationName)
            STR_extension  = STR_annoFile.split('/')[-1].split('.')[-1]
            if STR_extension == 'bed': 
                annotationFile_tree = IOC_annoFile.parseBed()
                for chrom in annotationFile_tree:            
                    if not chrom in DIC_annotations: DIC_annotations[chrom] = annotationFile_tree[chrom]
                    else: DIC_annotations[chrom] |= annotationFile_tree[chrom]
            # elif STR_ext == 'bedpe':
            #     raise RuntimeError("bedpe not implemented yet for custom annotations")
            #     LIS_customAnnotations.extend(IOC_annoFile.parseBedPE())
            else: raise RuntimeError("input error: annotation file {0} must be either .bed or .bedpe".format(STR_annoFile))


    # print(DIC_annotationMinROs)

    # STEP 7: print outputs
    sys.stderr.write('7. Printing variants\n')
    # one IOC_family at a time...
        
    INT_counter=0          
    #construct header line
    LIS_header = [
        'family', 'type',
        'chrom_1','pos_1','chrom_2','pos_2',
        'count_Fam_Aff', 'count_Fam_NAff', 'count_nonFam_Aff', 'count_nFam_NAff',
        'ids_Fam_Aff','ids_Fam_NAff','ids_nonFam_Aff','ids_nonFam_NAff',
        'source_batches'
    ]
    LIS_header.extend(LIS_annotationNames)
    #print header line
    outFH.write('\t'.join(LIS_header)+"\n")

    #Print each variant
    for IOC_currentVariant in LIS_segregated_variants:

        #Get annotation for this variant
        LIS_annotationsToPrint = list()
        if DIC_annotations and IOC_currentVariant._chrom1 in DIC_annotations:
            #Get all annotations included in the current variant using the tree search operator[]
            LIS_overlappingAnnotationsInVariant = DIC_annotations[IOC_currentVariant._chrom1][IOC_currentVariant._pos1:IOC_currentVariant._pos2]
            for currentAnnotationName in LIS_annotationNames:
                SET_annotationsWithOverlap = set()    
                #Iterate annotations by name (type)
                for annotationFromInterval in (interval.data for interval in LIS_overlappingAnnotationsInVariant if interval.data._type == currentAnnotationName):
                    #Check to see if reciprical overlap rule is satisfied
                    if (IOC_currentVariant.getPairwiseReciprocalOverlap(annotationFromInterval) > float(DIC_annotationMinROs[currentAnnotationName])):
                        SET_annotationsWithOverlap.add(annotationFromInterval._annotation)
                LIS_annotationsToPrint.append (",".join(sorted(list(SET_annotationsWithOverlap))))

        #printing variants in affected
        if Config.STR_otherRulesAndParameters_emitStatus == affected:
            for IOC_family in IOC_currentVariant._families_with_affected:
                STR_segInFam = IOC_currentVariant.segregationInFamily(IOC_family)
                STR_sourceBatches = ','.join(IOC_currentVariant._sourceBatches)
                LIS_output = [IOC_family, IOC_currentVariant, STR_segInFam, STR_sourceBatches]
                LIS_output.extend(LIS_annotationsToPrint)
                outFH.write('\t'.join(map(str, LIS_output)) + '\n')                

        #printing variants in non-affected
        elif Config.STR_otherRulesAndParameters_emitStatus == nonaffected:
            for IOC_family in IOC_currentVariant._families_with_nonaffected:
                STR_segInFam = IOC_currentVariant.segregationInFamily(IOC_family)
                STR_sourceBatches = ','.join(IOC_currentVariant._sourceBatches)
                LIS_output = [IOC_family, IOC_currentVariant, STR_segInFam, STR_sourceBatches]
                LIS_output.extend(LIS_annotationsToPrint)
                outFH.write('\t'.join(map(str, LIS_output)) + '\n')
            
        #print all variants
        elif Config.STR_otherRulesAndParameters_emitStatus == allsamples:
            for IOC_family in IOC_currentVariant._families_with_nonaffected | IOC_currentVariant._families_with_affected:
                STR_segInFam = IOC_currentVariant.segregationInFamily(IOC_family)
                STR_sourceBatches = ','.join(IOC_currentVariant._sourceBatches)
                LIS_output = [IOC_family, IOC_currentVariant, STR_segInFam, STR_sourceBatches]
                LIS_output.extend(LIS_annotationsToPrint)
                outFH.write('\t'.join(map(str, LIS_output))+ '\n')
    
    # STEP 8: print run statistics (optional)
    if Config.verbose:
        sys.stderr.write('8. printing run statistics (optional)')
        DIC_null = {Deletion:[],Duplication:[],Inversion:[],Translocation:[]}
        sys.stderr.write('raw variants by sample\n')
        for IOC_sample in Sample.getInstances():
            DIC_raw = {}
            if IOC_sample in DIC_variantsInDataSource:
                for STR_type in DIC_null.keys():
                    DIC_tmp = {STR_type:DIC_variantsInDataSource[IOC_sample][STR_type]} if STR_type in DIC_variantsInDataSource[IOC_sample] else {STR_type:[]}
                    DIC_raw.update(DIC_tmp)
            else:
                DIC_raw = DIC_null 
            sys.stderr.write('\t'.join([str(IOC_sample), '\t'.join([k + ':' + str(len(DIC_raw[k])) for k in sorted(DIC_raw.keys())])]) + '\n')  
