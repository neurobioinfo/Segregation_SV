import sys
import itertools
import helpers
from operator import neg
from constants import *
from sortedcontainers import *
from intervaltree import Interval, IntervalTree

def cmp(a, b):
    return (a > b) - (a < b)

class BaseVariant(object):
        
    def __init__(self, STR_chrom1, INT_pos1, STR_chrom2, INT_pos2, STR_varType, INT_pos1_brkp1=None, INT_pos1_brkp2=None, INT_pos2_brkp1=None, INT_pos2_brkp2=None):
        self._chrom1                 = STR_chrom1
        self._pos1                   = INT_pos1
        self._chrom2                 = STR_chrom2
        self._pos2                   = INT_pos2
        self._type                   = STR_varType
        self._INT_pos1_brkp1         = INT_pos1_brkp1
        self._INT_pos1_brkp2         = INT_pos1_brkp2
        self._INT_pos2_brkp1         = INT_pos2_brkp1
        self._INT_pos2_brkp2         = INT_pos2_brkp2
        self._hash                   = hash((STR_chrom1, INT_pos1, STR_chrom2, INT_pos2))
                
    def __repr__(self):
        """ a more detailed representation of the variant, mainly used for debugging """
        return str(self.__class__) + ":" + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2) + ":" + str(self._type)
    
    def __str__(self):
        """ a more general representation of the variant, mainly used for printing final outputs """
        return str(self.__class__) + ":" + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2) + ":" + str(self._type)
    
    def __hash__(self):
        return self._hash
     
    def __cmp__(self, IOC_other):
        if self.__class__ != IOC_other.__class__ : return cmp(self.__class__, IOC_other.__class__)
        if self._type != IOC_other._type: return cmp(self._type, IOC_other._type)
        if self._chrom1 != IOC_other._chrom1: return cmp(self._chrom1, IOC_other._chrom1)
        if self._pos1 != IOC_other._pos1: return cmp(self._pos1, IOC_other._pos1)
        if self._chrom2 != IOC_other._chrom2: return cmp(self._chrom2, IOC_other._chrom2)
        if self._pos2 != IOC_other._pos2: return cmp(self._pos2, IOC_other._pos2)
        return 0
    
    def __lt__(self, IOC_other):
        return self.__cmp__(IOC_other) < 0
    
    # getters    
    def getSize(self):
        if self._chrom1 != self._chrom2: raise RuntimeError("can't get size from translocation")
        return self._pos2 - self._pos1 +1
    
    def getType(self):
        return self._type
     
    def sort_by_start(self):
        return (self._chrom1,self._chrom2,self._pos1,self._pos2,self._type)    
     
    def sort_by_stop(self):
        return (self._chrom1,self._chrom2,-self._pos2,-self._pos1,self._type)

    def sort_by_size(self):
        if self._chrom1 != self._chrom2: raise RuntimeError("Size unavailable for variant with 2 different chroms")
        return (self._chrom1,self._chrom2,(self._pos2 - self._pos1),self._pos1,self._pos2,self._type)
    
    # doers
    def getPairwiseReciprocalOverlap(self, IOC_other):
        # different chromosomes? FAIL
        if self._chrom1 != IOC_other._chrom1 or self._chrom2 != IOC_other._chrom2: return 0.0
        FLO_RO=None        
        if self._type == Translocation:
            INT_overlap_break1 = min(self._INT_pos1_brkp2, IOC_other._INT_pos1_brkp2) - max(self._INT_pos1_brkp1, IOC_other._INT_pos1_brkp1)
            INT_overlap_break2 = min(self._INT_pos2_brkp2, IOC_other._INT_pos2_brkp2) - max(self._INT_pos2_brkp1, IOC_other._INT_pos2_brkp1)            
            FLO_break1_self     = float(INT_overlap_break1) / (self._INT_pos1_brkp2 - self._INT_pos1_brkp1)
            FLO_break1_other    = float(INT_overlap_break1) / (IOC_other._INT_pos1_brkp2 - IOC_other._INT_pos1_brkp1)
            FLO_break2_self     = float(INT_overlap_break2) / (self._INT_pos2_brkp2 - self._INT_pos2_brkp1)
            FLO_break2_other    = float(INT_overlap_break2) / (IOC_other._INT_pos2_brkp2 - IOC_other._INT_pos2_brkp1)    
            FLO_RO = min([FLO_break1_self,FLO_break1_other,FLO_break2_self,FLO_break2_other])
        else:
            # non-overlapping intervals? FAIL
            if self._pos2 + 1 < IOC_other._pos1 or self._pos1 > IOC_other._pos2 + 1: return 0.0
            (INT_selfSize, INT_otherSize) = (self.getSize(), IOC_other.getSize())
            INT_overlap = min(self._pos2 +1, IOC_other._pos2 +1) - max(self._pos1, IOC_other._pos1) 
            FLO_RO = min(float(INT_overlap)/INT_selfSize, float(INT_overlap)/INT_otherSize)
        return FLO_RO
                   
class RawVariant(BaseVariant):
    def __init__(self, STR_chrom1, INT_pos1, STR_chrom2, INT_pos2, STR_varType, IOC_sample, STR_caller, DIC_callerDetails, STR_sourceBatch, INT_pos1_brkp1=None, INT_pos1_brkp2=None, INT_pos2_brkp1=None, INT_pos2_brkp2=None):
        BaseVariant.__init__(self, STR_chrom1, INT_pos1, STR_chrom2, INT_pos2, STR_varType, INT_pos1_brkp1, INT_pos1_brkp2, INT_pos2_brkp1, INT_pos2_brkp2)
        self._sample                 = IOC_sample
        self._caller                 = STR_caller
        self._caller_details         = DIC_callerDetails
        self._sourceBatch            = STR_sourceBatch
        self._hash                   = hash((IOC_sample, STR_caller, STR_varType, STR_sourceBatch, STR_chrom1, INT_pos1, STR_chrom2, INT_pos2))

    def __hash__(self):
        return self._hash

    def __repr__(self):
        """ a more detailed representation of the variant, mainly used for debugging """
        return "RawVariant::" + str(self._sample) + '::' + self._type + '::' + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2) + '::' + str(self._sourceBatch) 
    
    def __str__(self):
        """ a more general representation of the variant, mainly used for printing final outputs """
        return "RawVariant::[" + str(self._sample) + '::' + self._type + ']::' + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2) + ':' + str(self._sourceBatch) 
     
    def __cmp__(self, IOC_other):
        if self.__class__ != IOC_other.__class__ : return cmp(self.__class__, IOC_other.__class__)
        if self._type != IOC_other._type: return cmp(self._type, IOC_other._type)
        if self._sourceBatch !=  IOC_other._sourceBatch: return cmp(self._sourceBatch, IOC_other._sourceBatch)
        if self._chrom1 != IOC_other._chrom1: return cmp(self._chrom1, IOC_other._chrom1)
        if self._pos1 != IOC_other._pos1: return cmp(self._pos1, IOC_other._pos1)
        if self._chrom2 != IOC_other._chrom2: return cmp(self._chrom2, IOC_other._chrom2)
        if self._pos2 != IOC_other._pos2: return cmp(self._pos2, IOC_other._pos2)
        if self._sample != IOC_other._sample: return cmp(self._sample, IOC_other._sample)
        return 0    

    def getType(self):
        return self._type

    def getSample(self):
        return self._sample
        
    def getCaller(self):
        return self._caller
    
    def getSource(self):
        return self._sourceBatch

class TempCombinedVariant (BaseVariant):         
    
    def __init__(self, LIS_input_variants):
        
        SET_chrom1     = SortedSet()
        SET_chrom2     = SortedSet()
        SET_pos1       = SortedSet()
        SET_pos2       = SortedSet()
        SET_pos1_brkp1 = SortedSet()
        SET_pos1_brkp2 = SortedSet()        
        SET_pos2_brkp1 = SortedSet()
        SET_pos2_brkp2 = SortedSet()
        
        self._type = None
        self._contained_variants = SortedSet(key=BaseVariant.sort_by_start)      
        self._raw_variants = SortedSet(key=BaseVariant.sort_by_start)
        
        for IOC_var in LIS_input_variants:
            SET_chrom1.add(IOC_var._chrom1)
            SET_chrom2.add(IOC_var._chrom2)
            SET_pos1.add(IOC_var._pos1)
            SET_pos2.add(IOC_var._pos2)
            if IOC_var._INT_pos1_brkp1 is not None: SET_pos1_brkp1.add(IOC_var._INT_pos1_brkp1)
            else: SET_pos1_brkp1.add(IOC_var._pos1)
            if IOC_var._INT_pos1_brkp2 is not None: SET_pos1_brkp2.add(IOC_var._INT_pos1_brkp2)
            else: SET_pos1_brkp2.add(IOC_var._pos1)            
            if IOC_var._INT_pos2_brkp1 is not None: SET_pos2_brkp1.add(IOC_var._INT_pos2_brkp1)
            else: SET_pos2_brkp1.add(IOC_var._pos2) 
            if IOC_var._INT_pos2_brkp2 is not None: SET_pos2_brkp2.add(IOC_var._INT_pos2_brkp2)
            else: SET_pos2_brkp2.add(IOC_var._pos2)
            
            if   isinstance(IOC_var, RawVariant) : 
                self._raw_variants.add(IOC_var)
                self._contained_variants.add(IOC_var)
            elif isinstance(IOC_var, ConsensusVariant) : 
                self._raw_variants.update(IOC_var._raw_variants)
                self._contained_variants.add(IOC_var)
            elif isinstance(IOC_var, MergedVariant) : 
                self._raw_variants.update(IOC_var._raw_variants)
                self._contained_variants.add(IOC_var)
            elif isinstance(IOC_var, TempCombinedVariant) : 
                self._raw_variants.update(IOC_var._raw_variants)
                self._contained_variants.update(IOC_var._contained_variants)
            else: raise RuntimeError("new TempCombinedVariant with unknown class " + IOC_var.__class__.__name__)
            
            if self._type is not None:
                if (self._type!= IOC_var._type):
                    raise RuntimeError("Multiple types in TempCombinedVariant")
            else: self._type = IOC_var._type              
      
        # sanity check: only one set of chrom1/chrom2 allowed!
        if len(SET_chrom1)>1 or len(SET_chrom2)>1: raise RuntimeError('only one set of chrom1/chrom2 allowed')
        STR_chrom1 = SET_chrom1.pop()
        STR_chrom2 = SET_chrom2.pop()
                   
        INT_pos1 = min(SET_pos1)
        INT_pos2 = max(SET_pos2)
        INT_pos1_brkp1 = min(SET_pos1_brkp1)
        INT_pos1_brkp2 = max(SET_pos1_brkp2)
        INT_pos2_brkp1 = min(SET_pos2_brkp1)
        INT_pos2_brkp2 = max(SET_pos2_brkp2)

        self._hash = hash(( self._type, tuple(self._contained_variants), STR_chrom1, INT_pos1, STR_chrom2, INT_pos2))
        BaseVariant.__init__(self, STR_chrom1, INT_pos1, STR_chrom2, INT_pos2, self._type, INT_pos1_brkp1, INT_pos1_brkp2, INT_pos2_brkp1, INT_pos2_brkp2)

    def __str__(self):
        """ a more general representation of the variant, mainly used for printing final outputs """
        return "TempCombinedVariant::[" + self._type + ']::' + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2) + "::" + str(self._contained_variants)

class ConsensusVariant(BaseVariant):
    def __init__(self, LIS_input_raw_variants):
        # parse coordinates of all variants into sets for subsequent steps
        SET_chrom1     = SortedSet()
        SET_chrom2     = SortedSet()
        SET_pos1       = SortedSet()
        SET_pos2       = SortedSet()
        SET_pos1_brkp1 = SortedSet()
        SET_pos1_brkp2 = SortedSet()        
        SET_pos2_brkp1 = SortedSet()
        SET_pos2_brkp2 = SortedSet()

        LIS_raw_from_input = list()
        for in_var in LIS_input_raw_variants:
            if  isinstance(in_var, ConsensusVariant):
                LIS_raw_from_input.extend(list(in_var._raw_variants))
            else: LIS_raw_from_input.append(in_var)

        self._sourceBatch = None
        self._sample = None
        self._type = None
        self._contained_variants = SortedSet(iterable=LIS_raw_from_input)
        self._raw_variants = self._contained_variants
        
        for IOC_var in LIS_raw_from_input:
            SET_chrom1.add(IOC_var._chrom1)
            SET_chrom2.add(IOC_var._chrom2)
            SET_pos1.add(IOC_var._pos1)
            SET_pos2.add(IOC_var._pos2)
            if IOC_var._INT_pos1_brkp1 is not None: SET_pos1_brkp1.add(IOC_var._INT_pos1_brkp1)
            else: SET_pos1_brkp1.add(IOC_var._pos1)
            if IOC_var._INT_pos1_brkp2 is not None: SET_pos1_brkp2.add(IOC_var._INT_pos1_brkp2)
            else: SET_pos1_brkp2.add(IOC_var._pos1)            
            if IOC_var._INT_pos2_brkp1 is not None: SET_pos2_brkp1.add(IOC_var._INT_pos2_brkp1)
            else: SET_pos2_brkp1.add(IOC_var._pos2) 
            if IOC_var._INT_pos2_brkp2 is not None: SET_pos2_brkp2.add(IOC_var._INT_pos2_brkp2)
            else: SET_pos2_brkp2.add(IOC_var._pos2)
            
            if not isinstance(IOC_var, RawVariant):
                raise RuntimeError("ConsensusVariant with non RawVariant")
            
            if self._sourceBatch is not None:
                if self._sourceBatch!= IOC_var._sourceBatch:
                    raise RuntimeError("Multiple batches in ConsensusVariant")
            else: self._sourceBatch = IOC_var._sourceBatch       
            
            if self._type is not None:
                if self._type!= IOC_var._type:
                    raise RuntimeError("Multiple types in ConsensusVariant")
            else: self._type = IOC_var._type  
            
            if self._sample is not None:
                if self._sample!= IOC_var._sample:
                    raise RuntimeError("Multiple samples in ConsensusVariant")
            else: self._sample = IOC_var._sample                            

        # sanity check: only one set of chrom1/chrom2 allowed!
        if len(SET_chrom1)>1 or len(SET_chrom2)>1: raise RuntimeError('only one set of chrom1/chrom2 allowed')
        STR_chrom1 = SET_chrom1.pop()
        STR_chrom2 = SET_chrom2.pop()
 
        INT_pos1 = min(SET_pos1)
        INT_pos2 = max(SET_pos2)
        INT_pos1_brkp1 = min(SET_pos1_brkp1)
        INT_pos1_brkp2 = max(SET_pos1_brkp2)
        INT_pos2_brkp1 = min(SET_pos2_brkp1)
        INT_pos2_brkp2 = max(SET_pos2_brkp2)
                              
        self._hash = hash((self._sample, self._type, self._sourceBatch, tuple(self._contained_variants), STR_chrom1, INT_pos1, STR_chrom2, INT_pos2))
        BaseVariant.__init__(self, STR_chrom1, INT_pos1, STR_chrom2, INT_pos2, self._type, INT_pos1_brkp1, INT_pos1_brkp2, INT_pos2_brkp1, INT_pos2_brkp2)
        
    def __repr__(self):
        """ a more detailed representation of the variant, mainly used for debugging """
        return "ConsensusVariant::[" + str(self._sample) + '::' + self._type + '::' + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2) + "::" + str(self._sourceBatch) 
    
    def __str__(self):
        """ a more general representation of the variant, mainly used for printing final outputs """
        return "ConsensusVariant::[" + self._type + ']::' + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2) + ':' + str(self._sourceBatch) 
    
    def __hash__(self):
        return self._hash
     
    def __cmp__(self, IOC_other):
        if self.__class__ != IOC_other.__class__ : return cmp(self.__class__, IOC_other.__class__)
        if self._type != IOC_other._type: return cmp(self._type, IOC_other._type)
        if self._sourceBatch != IOC_other._sourceBatch: return cmp(self._sourceBatch, IOC_other._sourceBatch)
        if self._chrom1 != IOC_other._chrom1: return cmp(self._chrom1, IOC_other._chrom1)
        if self._pos1 != IOC_other._pos1: return cmp(self._pos1, IOC_other._pos1)
        if self._chrom2 != IOC_other._chrom2: return cmp(self._chrom2, IOC_other._chrom2)
        if self._pos2 != IOC_other._pos2: return cmp(self._pos2, IOC_other._pos2)
        if self._sample != IOC_other._sample: return cmp(self._sample, IOC_other._sample)
        if self._contained_variants != IOC_other._contained_variants: return cmp(tuple(self._contained_variants), tuple(IOC_other._contained_variants))
        return 0
             
class MergedVariant(BaseVariant):
    def __init__(self, LIS_input_variants):
        """ create a new variant from a list of variants """
        # parse coordinates of all variants into sets for subsequent steps
        SET_chrom1     = SortedSet()
        SET_chrom2     = SortedSet()
        SET_pos1       = SortedSet()
        SET_pos2       = SortedSet()
        SET_pos1_brkp1 = SortedSet()
        SET_pos1_brkp2 = SortedSet()        
        SET_pos2_brkp1 = SortedSet()
        SET_pos2_brkp2 = SortedSet()

        self._type = None
        self._sourceBatches = SortedSet()
        self._samples = SortedSet()
        self._raw_variants = SortedSet(key=BaseVariant.sort_by_start)
        self._contained_variants = SortedSet(iterable=LIS_input_variants, key=BaseVariant.sort_by_start)
        
        for IOC_var in LIS_input_variants:
            SET_chrom1.add(IOC_var._chrom1)
            SET_chrom2.add(IOC_var._chrom2)
            SET_pos1.add(IOC_var._pos1)
            SET_pos2.add(IOC_var._pos2)

            if IOC_var._INT_pos1_brkp1 is not None: SET_pos1_brkp1.add(IOC_var._INT_pos1_brkp1)
            else: SET_pos1_brkp1.add(IOC_var._pos1)
            if IOC_var._INT_pos1_brkp2 is not None: SET_pos1_brkp2.add(IOC_var._INT_pos1_brkp2)
            else: SET_pos1_brkp2.add(IOC_var._pos1)
            if IOC_var._INT_pos2_brkp1 is not None: SET_pos2_brkp1.add(IOC_var._INT_pos2_brkp1)
            else: SET_pos2_brkp1.add(IOC_var._pos2) 
            if IOC_var._INT_pos2_brkp2 is not None: SET_pos2_brkp2.add(IOC_var._INT_pos2_brkp2)
            else: SET_pos2_brkp2.add(IOC_var._pos2)      
                       
            if isinstance(IOC_var, RawVariant):
                self._raw_variants.add(IOC_var)
                self._samples.add(IOC_var._sample)
                self._sourceBatches.add(IOC_var._sourceBatch)
            elif isinstance(IOC_var, ConsensusVariant):
                self._raw_variants.update(IOC_var._raw_variants)
                self._samples.add(IOC_var._sample)
                self._sourceBatches.add(IOC_var._sourceBatch)                                  
            elif isinstance(IOC_var, MergedVariant) :
                self._raw_variants.update(IOC_var._raw_variants)
                self._samples.update(IOC_var._samples)
                self._sourceBatches.update(IOC_var._sourceBatches)          
            else: raise RuntimeError("Non variant in LIS_variants")
            
            if self._type is not None:
                if (self._type!= IOC_var._type):
                    raise RuntimeError("Multiple types in MergedVariant")
            else: self._type = IOC_var._type  

        # sanity check: only one set of chrom1/chrom2 allowed!
        if len(SET_chrom1)>1 or len(SET_chrom2)>1: raise RuntimeError('only one set of chrom1/chrom2 allowed')
        STR_chrom1 = SET_chrom1.pop()
        STR_chrom2 = SET_chrom2.pop()
 
        INT_pos1 = min(SET_pos1)
        INT_pos2 = max(SET_pos2)
        INT_pos1_brkp1 = min(SET_pos1_brkp1)
        INT_pos1_brkp2 = max(SET_pos1_brkp2)
        INT_pos2_brkp1 = min(SET_pos2_brkp1)
        INT_pos2_brkp2 = max(SET_pos2_brkp2)          
              
        self._hash = hash((tuple(self._samples), self._type, tuple(self._sourceBatches), tuple(self._contained_variants), tuple(self._raw_variants), STR_chrom1, INT_pos1, STR_chrom2, INT_pos2))
        BaseVariant.__init__(self, STR_chrom1, INT_pos1, STR_chrom2, INT_pos2, self._type, INT_pos1_brkp1, INT_pos1_brkp2, INT_pos2_brkp1, INT_pos2_brkp2)

    def __repr__(self):
        """ a more detailed representation of the variant, mainly used for debugging """
        return "MergedVariant::" + str(self._samples) + '::' + self._type + '::' + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2)
    
    def __str__(self):
        """ a more general representation of the variant, mainly used for printing final outputs """
        return "MergedVariant::[" + self._type + ']::' + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2)
    
    def __hash__(self):
        return self._hash
     
    def __cmp__(self, IOC_other):
        if self.__class__ != IOC_other.__class__ : return cmp(self.__class__, IOC_other.__class__)
        if self._type != IOC_other._type: return cmp(self._type, IOC_other._type)
        if self._sourceBatches != IOC_other._sourceBatches: return cmp(self._sourceBatches, IOC_other._sourceBatches)
        if self._chrom1 != IOC_other._chrom1: return cmp(self._chrom1, IOC_other._chrom1)
        if self._pos1 != IOC_other._pos1: return cmp(self._pos1, IOC_other._pos1)
        if self._chrom2 != IOC_other._chrom2: return cmp(self._chrom2, IOC_other._chrom2)
        if self._pos2 != IOC_other._pos2: return cmp(self._pos2, IOC_other._pos2)
        if self._samples != IOC_other._samples: return cmp(self._samples, IOC_other._samples)
        if self._contained_variants != IOC_other._contained_variants: return cmp(tuple(self._contained_variants), tuple(IOC_other._contained_variants))
        # if self._raw_variants != IOC_other._contained_raw_variants: return cmp(tuple(self._contained_raw_variants), tuple(IOC_other._contained_raw_variants))
        return 0 

class SegregatedVariant(BaseVariant):
    def __init__(self, LIS_input_variants):
        """ create a new variant from a list of variants """
        # parse coordinates of all variants into sets for subsequent steps
        SET_chrom1     = SortedSet()
        SET_chrom2     = SortedSet()
        SET_pos1       = SortedSet()
        SET_pos2       = SortedSet()
        SET_type       = SortedSet()
        SET_pos1_brkp1 = SortedSet()
        SET_pos1_brkp2 = SortedSet()        
        SET_pos2_brkp1 = SortedSet()
        SET_pos2_brkp2 = SortedSet()

        self._sourceBatches = SortedSet()
        self._samples = SortedSet()
        self._contained_variants = SortedSet(LIS_input_variants, key=BaseVariant.sort_by_start)
        self._raw_variants = SortedSet(key=BaseVariant.sort_by_start)
        self._families_with_affected = SortedSet()
        self._families_with_nonaffected = SortedSet()
                
        for IOC_var in LIS_input_variants:
            SET_chrom1.add(IOC_var._chrom1)
            SET_chrom2.add(IOC_var._chrom2)
            SET_pos1.add(IOC_var._pos1)
            SET_pos2.add(IOC_var._pos2)
            SET_type.add(IOC_var._type)

            if IOC_var._INT_pos1_brkp1 is not None: SET_pos1_brkp1.add(IOC_var._INT_pos1_brkp1)
            else: SET_pos1_brkp1.add(IOC_var._pos1)
            if IOC_var._INT_pos1_brkp2 is not None: SET_pos1_brkp2.add(IOC_var._INT_pos1_brkp2)
            else: SET_pos1_brkp2.add(IOC_var._pos1)
            if IOC_var._INT_pos2_brkp1 is not None: SET_pos2_brkp1.add(IOC_var._INT_pos2_brkp1)
            else: SET_pos2_brkp1.add(IOC_var._pos2) 
            if IOC_var._INT_pos2_brkp2 is not None: SET_pos2_brkp2.add(IOC_var._INT_pos2_brkp2)
            else: SET_pos2_brkp2.add(IOC_var._pos2)    
                       
            if isinstance(IOC_var, RawVariant):
                self._raw_variants.add(IOC_var)
                self._contained_variants.add(IOC_var)
                self._samples.add(IOC_var._sample)
                self._sourceBatches.add(IOC_var._sourceBatch)
            elif isinstance(IOC_var, ConsensusVariant):
                self._raw_variants.update(IOC_var._raw_variants)
                self._samples.add(IOC_var._sample)
                self._sourceBatches.add(IOC_var._sourceBatch)                                  
            elif isinstance(IOC_var, MergedVariant) :
                self._contained_variants.update(IOC_var._raw_variants)
                self._raw_variants.update(IOC_var._raw_variants)
                self._samples.update(IOC_var._samples)
                self._sourceBatches.update(IOC_var._sourceBatches)   
            elif isinstance(IOC_var, SegregatedVariant) :
                # print "seg"
                self._contained_variants.update(IOC_var._contained_variants)
                self._raw_variants.update(IOC_var._raw_variants)
                self._samples.update(IOC_var._samples)
                self._sourceBatches.update(IOC_var._sourceBatches)
            else: raise RuntimeError("Non variant in LIS_variants")

        # sanity check: only one set of chrom1/chrom2 allowed!
        if len(SET_chrom1)>1 or len(SET_chrom2)>1: raise RuntimeError('only one set of chrom1/chrom2 allowed')
        STR_chrom1 = SET_chrom1.pop()
        STR_chrom2 = SET_chrom2.pop()
                
        if len(SET_type)>1: raise RuntimeError('only one type per SegregatedVariant')
        STR_varType = SET_type.pop()    
                 
        INT_pos1 = min(SET_pos1)
        INT_pos2 = max(SET_pos2)
        INT_pos1_brkp1 = min(SET_pos1_brkp1)
        INT_pos1_brkp2 = max(SET_pos1_brkp2)
        INT_pos2_brkp1 = min(SET_pos2_brkp1)
        INT_pos2_brkp2 = max(SET_pos2_brkp2)          
              
        for sample in self._samples:
            if sample._isAffected:
                self._families_with_affected.add(sample._family)
            else:
                self._families_with_nonaffected.add(sample._family)
              
        self._hash = hash((tuple(self._samples), STR_varType, tuple(self._sourceBatches), tuple(self._contained_variants), tuple(self._contained_variants), STR_chrom1, INT_pos1, STR_chrom2, INT_pos2))
        BaseVariant.__init__(self, STR_chrom1, INT_pos1, STR_chrom2, INT_pos2, STR_varType, INT_pos1_brkp1, INT_pos1_brkp2, INT_pos2_brkp1, INT_pos2_brkp2)

    def __repr__(self):
        """ a more detailed representation of the variant, mainly used for debugging """
        return "SegregatedVariant::" + str(self._samples) + '::' + self._type + '::' + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2)
    
    def __str__(self):
        """ a more general representation of the variant, mainly used for printing final outputs """
        #''.join(list1)
        return self._type + '\t' + str(self._chrom1) + '\t' + str(self._pos1) + '\t' + str(self._chrom2) + '\t' + str(self._pos2)
    
    def __hash__(self):
        return self._hash
     
    def __cmp__(self, IOC_other):
        if self.__class__ != IOC_other.__class__ : return cmp(self.__class__, IOC_other.__class__)
        if self._type != IOC_other._type: return cmp(self._type, IOC_other._type)
        if self._sourceBatches != IOC_other._sourceBatches: return cmp(self._sourceBatches, IOC_other._sourceBatches)
        if self._chrom1 != IOC_other._chrom1: return cmp(self._chrom1, IOC_other._chrom1)
        if self._pos1 != IOC_other._pos1: return cmp(self._pos1, IOC_other._pos1)
        if self._chrom2 != IOC_other._chrom2: return cmp(self._chrom2, IOC_other._chrom2)
        if self._pos2 != IOC_other._pos2: return cmp(self._pos2, IOC_other._pos2)
        if self._samples != IOC_other._samples: return cmp(self._samples, IOC_other._samples)
        if self._contained_variants != IOC_other._contained_variants: return cmp(tuple(self._contained_variants), tuple(IOC_other._contained_variants))
        if self._contained_variants != IOC_other._contained_variants: return cmp(tuple(self._contained_variants), tuple(IOC_other._contained_variants))
        return 0
    
    def segregationInFamily(self, IOC_fam, INT_numSamplesToPrintLimit=10):
        """ returns a print-ready string detailing counts and ids of samples in/out of family with the variant """       
        INT_famAff=0
        INT_famNAff=0
        INT_nFamAff=0
        INT_nFamNAff=0
        SET_famAff = SortedSet()
        SET_famNAff = SortedSet()
        SET_nFamAff = SortedSet()
        SET_nFamNAff = SortedSet()
        
        # iterate all samples associated with variant, and increment/update DIC_ref values according to aff/naff and in/out-fam
        for IOC_sample in self._samples:            
            if IOC_sample.inFamily(IOC_fam):
                if IOC_sample.isAffected():
                    INT_famAff+=1
                    SET_famAff.add(str(IOC_sample))
                else:
                    INT_famNAff+=1
                    SET_famNAff.add(str(IOC_sample))
            else:
                if IOC_sample.isAffected():
                    INT_nFamAff+=1
                    SET_nFamAff.add(str(IOC_sample))
                else:
                    INT_nFamNAff+=1
                    SET_nFamNAff.add(str(IOC_sample))

        if len(SET_famAff) > INT_numSamplesToPrintLimit: 
            SET_famAff=set(['>' + str(INT_numSamplesToPrintLimit) + ' samples'])
        if len(SET_famNAff) > INT_numSamplesToPrintLimit: 
            SET_famNAff=set(['>' + str(INT_numSamplesToPrintLimit) + ' samples'])
        if len(SET_nFamAff) > INT_numSamplesToPrintLimit: 
            SET_nFamAff=set(['>' + str(INT_numSamplesToPrintLimit) + ' samples'])
        if len(SET_nFamNAff) > INT_numSamplesToPrintLimit: 
            SET_nFamNAff=set(['>' + str(INT_numSamplesToPrintLimit) + ' samples'])                        
            
        STR_counts = '\t'.join(map(str, [INT_famAff, INT_famNAff, INT_nFamAff, INT_nFamNAff]))            
        STR_ids = '\t'.join(map(lambda x:','.join(x), [SET_famAff, SET_famNAff, SET_nFamAff, SET_nFamNAff]))   
            
        return '\t'.join([STR_counts, STR_ids])

class ConsensusVariantGenerator(object):

    """ a class which takes a list of Variants, performs a grouping operation, and returns a list of Variants """
    def __init__(self, LIS_input):
        self._contained_variants = LIS_input     
    
    def variant_consensus(self, FLO_minRO=0.0):
        """ takes a list of variants and groups into a new list of variants, based on Union of coordinates """
   
        remaining_variants = SortedList(self._contained_variants, key=BaseVariant.sort_by_size)
        size = remaining_variants[0].getSize()       
        
        while (len(remaining_variants) > 1):
            matrix = SortedDict(neg)      
            for IOC_var1 in remaining_variants:
                for IOC_var2 in remaining_variants:                    
                    if IOC_var1 == IOC_var2 : break
                    if IOC_var1.getSize() != size and IOC_var2.getSize() != size : break
                    matrix[IOC_var1.getPairwiseReciprocalOverlap(IOC_var2)] = (IOC_var1,IOC_var2)     
            
            (overlap,(IOC_var1,IOC_var2)) = matrix.peekitem(0)
                  
            if overlap < FLO_minRO:                
                remaining_variants.discard(IOC_var1)
                size = remaining_variants[0].getSize()
            else:                                  
                remaining_variants.discard(IOC_var1)
                remaining_variants.discard(IOC_var2)                
                IOC_new = TempCombinedVariant ([IOC_var1, IOC_var2]) 
                remaining_variants.add(IOC_new)
                size = remaining_variants[0].getSize()                
     
        if isinstance(remaining_variants[0], TempCombinedVariant):
            IOC_new = ConsensusVariant (list(remaining_variants[0]._raw_variants))
        else :
            IOC_new = ConsensusVariant ([self._contained_variants[0]])
        return IOC_new

    
    def variant_combine(self, FLO_minRO=0.0):
        """ takes a list of variants and groups into a new list of variants, based on Union of coordinates """
        remaining_variants = SortedList(self._contained_variants, key=BaseVariant.sort_by_size)
        size = remaining_variants[0].getSize()       
        
        while (len(remaining_variants) > 1):    
            matrix = SortedDict(neg)      
            for IOC_var1 in remaining_variants:
                for IOC_var2 in remaining_variants:
                    if IOC_var1 == IOC_var2 : break                    
                    if IOC_var2.getSize() != size and IOC_var1.getSize() != size : break                    
                    matrix[IOC_var1.getPairwiseReciprocalOverlap(IOC_var2)] = (IOC_var1,IOC_var2)     
                    
            (overlap,(IOC_var1,IOC_var2)) = matrix.peekitem(0)
            
            if overlap < FLO_minRO:                
                remaining_variants.discard(IOC_var1)
                size = remaining_variants[0].getSize()
                continue
            else :    
                #we have found a match              
                remaining_variants.discard(IOC_var1)
                remaining_variants.discard(IOC_var2)
                
                IOC_new = TempCombinedVariant ([IOC_var1, IOC_var2]) 
                remaining_variants.add(IOC_new)
                size = remaining_variants[0].getSize()                
        
        if isinstance(remaining_variants[0], TempCombinedVariant):
            IOC_returnVar = MergedVariant (list(remaining_variants[0]._contained_variants))
        else :
            IOC_returnVar = MergedVariant ([self._contained_variants[0]])
        
        return IOC_returnVar

    
    def variant_segregate(self, FLO_minRO=0.0):
        """ takes a list of variants and groups into a new list of variants, based on Union of coordinates """
        # initialize starting lists (inputs, current progress, final outputs, overlapFunction)
        LIS_varIN = self._contained_variants
        LIS_current = SortedList([LIS_varIN.pop()], key=BaseVariant.sort_by_start)
        LIS_varOUT = SortedList(key=BaseVariant.sort_by_start)
               
        # iterate position-sorted list of variants, building new variants on the fly based on Union of variants
        while LIS_varIN:
            BOOL_hits = False
            # remove current variant from list; will eventually shrink list down to zero
            IOC_varCurrent = LIS_varIN.pop()

            for IOC_varPrevious in LIS_current:
                # only one hit is needed to see if IOC_varCurrent belongs with LIS_current
                if IOC_varCurrent.getPairwiseReciprocalOverlap(IOC_varPrevious) > FLO_minRO:
                    BOOL_hits = True
                    break
            # if any hits, add IOC_varCurrent to LIS_current and goto next variant
            if BOOL_hits:
                LIS_current.add(IOC_varCurrent)
                continue
            # if no hits, add last variant (or Union of last overlapping variants) to outputs, and prepare for next iteration
            else:
                IOC_new = SegregatedVariant (list(set(LIS_current)))# if len(LIS_current)>1 else LIS_current[0]
                LIS_varOUT.add(IOC_new)
                LIS_current = SortedList([IOC_varCurrent],key=BaseVariant.sort_by_start)
        # populate remaining variant(s) from inputs (since while loop exits before populating last batch
        IOC_new = SegregatedVariant (list(set(LIS_current)))# if len(LIS_current)>1 else LIS_current[0]
        LIS_varOUT.add(IOC_new)
        return LIS_varOUT        
    
class Annotation(BaseVariant):
    # class variants to easily retrieve all processed annotation types and associated minRO values
    # Types = SortedSet()
    # minROs = {}
    
    # fundamental methods
    def __init__(self, STR_annoType, STR_chrom1, INT_pos1, STR_chrom2, INT_pos2, STR_anno):
        BaseVariant.__init__(self, STR_chrom1, INT_pos1, STR_chrom2, INT_pos2, STR_annoType)
        self._annotation = STR_anno
        # Annotation.Types.add(STR_annoType)
        # Annotation.minROs.update({STR_annoType:FLO_minRO})
    
    def __repr__(self):
        """ a more detailed representation of the variant, mainly used for debugging """
        return str(self.__class__) + ":" + str(self._chrom1) + ':' + str(self._pos1) + '__' + str(self._chrom2) + ':' + str(self._pos2) + ":" + str(self._type) + ":" + self._annotation
    
    def __cmp__(self, IOC_other):
        if self.__class__ != IOC_other.__class__ : return cmp(self.__class__, IOC_other.__class__)
        if self._chrom1 != IOC_other._chrom1: return cmp(self._chrom1, IOC_other._chrom1)
        if self._pos1 != IOC_other._pos1: return cmp(self._pos1, IOC_other._pos1)
        if self._chrom2 != IOC_other._chrom2: return cmp(self._chrom2, IOC_other._chrom2)
        if self._pos2 != IOC_other._pos2: return cmp(self._pos2, IOC_other._pos2)
        return 0
    
    # getters
    def getAnnotation(self):
        return self._annotation
    
    # def getMinRO(self):
    #     return Annotation.minROs[self._type]
    
    # @classmethod
    # def getTypes(cls):
    #     return list(cls.Types)

# class VariantAnnotationGenerator(object):
#     """ a class which takes a Variant object, a list of annotations,  
#         and returns a string with all annotations for that Variant """
#     def __init__(self, IOC_variant, LIS_annotationsLeftSorted, LIS_annotationsRightSorted):
#         self._left = LIS_annotationsLeftSorted
#         self._right = LIS_annotationsRightSorted
#         self._variant = IOC_variant
    
#     def annotate_variant(self):
#         LIS_annoTypes = Annotation.getTypes()
#         DIC_annotations = { STR_type:SortedSet() for STR_type in LIS_annoTypes}
#         leftMostAnnotation = helpers.findSingleLeftMost(self._variant, self._left)
#         rightMostAnnotation = helpers.findSingleRightMost(self._variant, self._right)
#         allPossibleAnnotations = helpers.getVariantsFromBorders(leftMostAnnotation, rightMostAnnotation, self._left)
#         for IOC_anno in allPossibleAnnotations:
#             STR_type = IOC_anno.getType()
#             STR_errmsg = """
#                 custom annotation of type "{0}" has no associated minimum RO value. 
#                 Be sure to include this information in a properly-formatted .cfg file (provided with "--config_file [argument]") 
#             """.format(STR_type)
#             # if IOC_anno overlaps IOC_variant with sufficient RO, add to outputs
#             # sanity check: does current annotation have minRO value?
#             try: 
#                 FLO_minRO = IOC_anno.getMinRO()
#             except: raise RuntimeError(STR_errmsg, Annotation.minROs)
#             FLO_RO = self._variant.getPairwiseReciprocalOverlap(IOC_anno)
#             if FLO_RO > FLO_minRO:
#                 DIC_annotations[STR_type].add(IOC_anno.getAnnotation())
#         STR_output = '\t'.join([';'.join(DIC_annotations[STR_type]) for STR_type in LIS_annoTypes])
#         return STR_output

class Family(object):
    Instances = SortedDict()
    Names = SortedDict()
    # fundamental methods
    def __init__(self, name):
        self._name        = name
        self._samples     = SortedList() # to be filled with a list of instances of sample objects
        Family.Names[self._name]    = self
        self._hash = hash(name)

    def __hash__(self):
        return self._hash
#         return hash(self._name)
    
    def __repr__(self):
        return self._name

    def __str__(self):
        return self._name

    def __cmp__(self, IOC_other):
        # print ("tron2")
        # print (self._name, IOC_other._name)
        return cmp(self._name, IOC_other._name)
        # return cmp (self._hash, IOC_other._hash)
    
    def __lt__(self, IOC_other):
        # print ("tron")
        return self.__cmp__(IOC_other) < 0

    # getters
    def getSamples(self):
        return self._samples
    
    def getVariants(self, STR_fromSamplesWithStatus):
        """ returns the list of all variants from samples in the family, by default for affecteds only """
        if STR_fromSamplesWithStatus == affected:
            LIS_samples = filter(lambda x:x.isAffected(), self.getSamples())
        elif STR_fromSamplesWithStatus == nonaffected:
            LIS_samples = filter(lambda x:not x.isAffected(), self.getSamples())
        elif STR_fromSamplesWithStatus == allsamples:
            LIS_samples = self.getSamples()
        else: raise RuntimeError("error: Family.getVariants() only accepts the following arguments: [" + affected + "," + nonaffected + "," + allsamples + "]")
        return [var for sample in LIS_samples for var in sample.getVariants()]

    # setters
    def addSample(self, sample):
        self._samples.add(sample)
    
    # doers
    def getVariantsFromListInFamilyWithStatus(self, LIS_vars, STR_fromSamplesWithStatus):
        if STR_fromSamplesWithStatus == affected: SET_samples = set(filter(lambda x:x.isAffected(), self.getSamples()))
        elif STR_fromSamplesWithStatus == nonaffected: SET_samples = set(filter(lambda x:not x.isAffected(), self.getSamples()))
        elif STR_fromSamplesWithStatus == allsamples: SET_samples = set(self.getSamples())
        else: raise RuntimeError("error: Family.getVariantsFromListInFamilyWithStatus(LIS_vars) requires one of the following arguments: [" + affected + "," + nonaffected + "," + allsamples + "]")
        
        LIS_output = SortedList()        
        for IOC_var in LIS_vars:
            samples_in_variant = SortedSet()
            if isinstance(IOC_var, RawVariant): 
                if len(SET_samples.intersection([IOC_var.getSample()]))>0: 
                    LIS_output.add(IOC_var)
            else:
                for raw_variant in IOC_var._contained_variants:
                    samples_in_variant.add(raw_variant.getSample())
                if len(SET_samples.intersection(samples_in_variant))>0: 
                    LIS_output.add(IOC_var)
                    
        return LIS_output

    @classmethod
    def GetSet(cls, IOC_target):
        """return the Instance IOC_target if it exists; create it before returning if it doesn't"""
        if IOC_target.__hash__() not in cls.Instances: cls.Instances[IOC_target] = IOC_target
        return cls.Instances.get(IOC_target, None)

    @classmethod
    def Get(cls, STR_target):
        """return the Instance whose name is STR_target"""
        return cls.Names.get(STR_target, None)
    
    @classmethod
    def getInstances(cls):
        return cls.Instances.values()

class Sample(object):  
    Instances = {}
    # fundamental methods
    def __init__(self, IOC_family, STR_name, STR_father_id, STR_mother_id, INT_sex, STR_clinical_status):
        self._family = IOC_family
        self._name = STR_name
        self._father_id = STR_father_id
        self._mother_id = STR_mother_id
        self._sex = INT_sex
        self._isAffected = True if STR_clinical_status == helpers.PedigreeFile.affected else False
        self._contained_variants = SortedSet()        # will be populated with a list of variant objects
        Sample.Instances[self._name] = self        # try this out to see if it speeds up the overall performance of the script (to be used in Sample.Get function)
        self._hash = hash((IOC_family, STR_name, STR_father_id, STR_mother_id, INT_sex, self._isAffected))

    def __repr__(self):
        return self._name

    def __hash__(self):
        return self._hash

    def __cmp__(self, IOC_other):
        if self.__class__ != IOC_other.__class__ : return cmp(self.__class__, IOC_other.__class__)
        if self._family != IOC_other._family: return cmp(self._family, IOC_other._family)
        if self._name != IOC_other._name: return cmp(self._name, IOC_other._name)
        if self._father_id != IOC_other._father_id: return cmp(self._father_id, IOC_other._father_id)
        if self._mother_id != IOC_other._mother_id: return cmp(self._mother_id, IOC_other._mother_id)
        if self._sex != IOC_other._sex: return cmp(self._sex, IOC_other._sex)
        if self._isAffected != IOC_other._isAffected: return cmp(self._isAffected, IOC_other._isAffected)
        if self._contained_variants != IOC_other._contained_variants: return cmp(tuple(self._contained_variants), tuple(IOC_other._contained_variants))            
        return 0

    def __lt__(self, IOC_other):
        return self.__cmp__(IOC_other) < 0

    # getters
    def getName(self):
        return self._name
    
    def getFamily(self):
        return Family.Get(self._family)

    def isAffected(self):
        return self._isAffected

    def getRawVariants(self):
        return [var for var in self._contained_variants if var.passesFilters()]
    
    # doers
    def addVariantToSample(self, IOC_variant):
        self._contained_variants.add(IOC_variant)
    
    def inFamily(self, IOC_fam):
        return self._family == IOC_fam

    @classmethod
    def Get(cls, STR_target):
        """return the Instance whose name is STR_target"""
        return cls.Instances.get(STR_target, None)
    
    @classmethod
    def getInstances(cls):
        return (value for value in cls.Instances.values())
