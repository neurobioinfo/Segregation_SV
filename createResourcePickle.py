#!/usr/bin/env python

import sys
import os
import argparse
import config
import data
import helpers
import sortedcontainers
if (sys.version_info > (3, 0)):
    import pickle as pickle
else:
    import cPickle as pickle
from helpers import *
from config import *
from sortedcontainers.sortedset import SortedSet

if __name__ == "__main__":

    # STEP 1: process script inputs
    parser = argparse.ArgumentParser(description='create picke object storing all desired annotations for future segregation runs')
    parser.add_argument('-a','--custom_annotation_file',dest='customAnnoData', action='append', help='/path/to/annotation_file:anno_type:min_reciprocal_overlap ([.bed or .bedpe]:[see config.py for accepted list of types])')
    parser.add_argument('-o', '--out_file',   dest='outFile',  help='/path/to/output/annotation.picklefile; where outputs will be written', default=sys.stdout)
    args = parser.parse_args()
    LIS_annotationFiles = []
    LIS_annotationNames = []
    DIC_annotationMinROs = {}
    for STR_input in args.customAnnoData:
        LIS_input = STR_input.split(':')
        STR_annotationFile = LIS_input[0]
        # if input not specified, parse from name and hope for the best
        STR_annotationName = LIS_input[1] if len(LIS_input) >= 2 else '.'.join(STR_annotationFile.split('/')[-1].split('.')[0:-1])
        FLO_annotationMinRO = LIS_input[2] if len(LIS_input) == 3 else 0.0
        LIS_annotationFiles.append(STR_annotationFile)
        LIS_annotationNames.append(STR_annotationName)
        DIC_annotationMinROs.update( {STR_annotationName : float(FLO_annotationMinRO) } ) #DIC_annotationMinROs.update({STR_annotationName : float(FLO_annotationMinRO)})
    STR_outFile = args.outFile

    # STEP 2: parse external database annotations into Annotation objects
    sys.stderr.write('1. Parsing annotation databases\n')
    DIC_annotations = {}
    for i in range(0,len(LIS_annotationFiles)):
        STR_annotationFile, STR_annotationName = (LIS_annotationFiles[i], LIS_annotationNames[i])
        print('Loading ' + str(STR_annotationFile))        
        STR_ext  = STR_annotationFile.split('/')[-1].split('.')[-1]
        file = helpers.AnnotationFile(STR_annotationFile, STR_annotationName)
        if STR_ext == 'bed':
            TREE_current = file.parseBed()
            for chrom in TREE_current:            
                if not chrom in DIC_annotations: DIC_annotations[chrom] = TREE_current[chrom]
                else: DIC_annotations[chrom] |= TREE_current[chrom]            
        # elif STR_ext == 'bedpe':
        #     TREE_current = file.parseBedPE()
        else: raise RuntimeError("input error: annotation file {0} must be either .bed or .bedpe".format(STR_annotationFile))
        # LIS_annotations.extend(LIS_current)
    # STEP 3: save annotation objects to python pickle for easier re-loading
    sys.stderr.write('2. Saving annotations to pickle file ' + str(STR_outFile) + '\n')
    pickle.dump( (DIC_annotations, LIS_annotationNames, DIC_annotationMinROs), open( STR_outFile, 'wb' ) )
