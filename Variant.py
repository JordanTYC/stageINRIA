#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    my class Variant
"""

class Variant:
    '''
        Splits a string (vcf format) to create a variant.
        
        line -- string
    '''
    def __init__(self,line):
        line = line.split()
        self.info = self.createDict(line)
        self.chrom = line[0]
        self.pos = self.get_pos(line)
        self.id = line[2]
        self.ref = line[3]
        self.alt = line[4]
        self.qual = line[5]
        self.filter = line[6]
        
    def createDict(self,line):
        '''
            Creates a dictionnary via all the info.
        '''
        info = line[7].split(";")
        d = {}
        for field in info:
            e = field.split('=')
            if len(e) == 1: # exception for keywords without value
                d[e[0]] = ""
            else:
                d[e[0]] = e[1] 
        return d
        
    def get_pos(self,line):
        '''
            Returns the variant's start position.
        '''
        if self.get_svtype() == "INS" and "LEFT_SVINSSEQ" in self.info:
            return int(line[1]) - len(self.info["LEFT_SVINSSEQ"])
        return int(line[1])

    def get_end(self):
        '''
            Returns the variant's end position.
        '''
        svtype = self.get_svtype()
        if svtype == "INS":
            if "RIGHT_SVINSSEQ" not in self.info:
                return self.pos + int(self.info["SVLEN"])
            else:
                return self.pos + len(self.info["RIGHT_SVINSSEQ"])
        else:
            return int(self.info["END"])

    def get_svtype(self):
        '''
            Returns the variant's type.
        '''
        return self.info["SVTYPE"]
    
    def get_svlen(self):
        '''
            Returns variant's length.
        '''
        if "SVLEN" in self.info:
            return int(self.info["SVLEN"])
