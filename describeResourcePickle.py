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
    parser = argparse.ArgumentParser(description='descibe the content of picke object storing all desired annotations')
    parser.add_argument('-f','--file',dest='annotationFile', help='/path/to/annotation_file')
    args = parser.parse_args()
    
    print ("Loading \"" + args.annotationFile + "\"")

    TREE_defaultAnnotations, LIS_defaultAnnotationNames, DIC_defaultAnnotationMinRO = pickle.load(open(args.annotationFile, 'rb'))
    print ("Annotation names: " + str(LIS_defaultAnnotationNames))
    print ("Annotation min overlaps: " + str(DIC_defaultAnnotationMinRO))

