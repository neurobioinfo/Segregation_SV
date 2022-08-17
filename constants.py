import sys
import os

# program constants that are NOT to be modified by segregation.ConfigParser
# variant types
Deletion = 'Deletion'
Duplication = 'Duplication'
Inversion = 'Inversion'
Translocation = 'Translocation'

# variant processing steps
cleanupRaw = 'cleanupRaw'
mergeCallsets = 'mergeCallsets'
segregation = 'segregation'
annotation = 'annotation'

# raw SV caller names
lumpy = 'lumpy'
popsv = 'popsv'
conifer = 'conifer'
genericBed = 'generic_bed'
genericBedpe = 'generic_bedpe'
codex = 'codex'
exomedepth = 'exomedepth'
cnmops = 'cnmops'
xhmm = 'xhmm'

# variant merge types - not yet implemented!
Union = 'Union'
Intersection = 'Intersection'
BestScoring = 'BestScoring'
Average = 'Average'
CleanestBreakpoint = 'CleanestBreakpoint'

# annotation types
DGV_inclusive = 'DGV_inclusive'
DGV_stringent = 'DGV_stringent'
KG_SV = 'KG_SV'
micro_exons = 'micro_exons'
repeatMasker = 'repeatMasker'
Segmental_Dups = 'Segmental_Dups'
RefSeq_genes = 'RefSeq_genes'
RefSeq_genes_exons = 'RefSeq_genes_exons'

# miscellaneous parameters
minOverlapRatio = 'minOverlapRatio'
affected = 'affected'
nonaffected = 'nonaffected'
allsamples = 'all'

# conversion dict for lumpy variant types
LUMPY_TYPES = { 'TYPE:DELETION' : Deletion, 'TYPE:DUPLICATION' : Duplication, 'TYPE:INVERSION' : Inversion, 'TYPE:INTERCHROM' : Translocation }
CONIFER_TYPES = { 'del' : Deletion, 'dup' : Duplication }
GENERIC_TYPES = { 'DEL' : Deletion, 'DUP' : Duplication, 'INV' : Inversion, 'TRA' : Translocation }
CODEX_TYPES = { 'del' : Deletion, 'dup' : Duplication }
EXOMEDEPTH_TYPES = { 'deletion' : Deletion, 'duplication' : Duplication }
XHMM_TYPES = { 'DEL' : Deletion, 'DUP' : Duplication }
#CN2 is missing. its meaning is normal 
CNMOPS_TYPES = { 'CN0' : Deletion, 'CN1' : Deletion, 'CN3' : Duplication, 'CN4' : Duplication, 'CN5' : Duplication, 'CN6' : Duplication, 'CN7' : Duplication, 'CN8' : Duplication, 'CN16' : Duplication, 'CN32' : Duplication, 'CN64' : Duplication } 

# allowed SV callers
DATA_PARSERS = set([lumpy, popsv, genericBed, genericBedpe, conifer, exomedepth, codex, xhmm, cnmops])

# default (hidden) files 
FILE_defaultIni = os.path.dirname(os.path.realpath(sys.argv[0])) + '/resources/default.cfg'
FILE_resources_pickle = os.path.dirname(os.path.realpath(sys.argv[0])) + '/resources/default_annotations.pickle'
