#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    Comparison of the number of isolated barcodes between real variants and false one.
"""

########## RAJOUTER EXCEPTION POUR LES GET


import argparse, pysam, xlsxwriter


N_GAP = 100


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

def trueSV(file):
    '''
        Returns a list of variants within a file.
        Registers first chromosome, start and end position.
        
        file -- file
    '''
    truth = []
    with open(file,"r") as filin:
        for line in filin:
            sv = line.split()
            truth.append((sv[0],int(sv[1]),int(sv[3])))
    return truth

    
def isValid(variant,L,m):
    '''
        Returns True if a Variant is valid, else False.
        
        variant -- Variant object
        L -- list of real variants
        m -- int
    '''
    for v in L:
        if variant.chrom == v[0] and abs(v[1] - variant.pos) <= m and abs(v[2] - variant.get_end()) <= m:
            return True
    return False


def isValid_bnd(variant,L,m):
    '''
        Returs True if a BND variant is valid, else False.
        
        variant -- a list as [chrom,start,end]
        L -- list of real variants
        m -- int
    '''
    for v in L:
        if variant[0] == v[0] and abs(variant[1] - v[1]) <= m and abs(variant[2] - v[2]) <= m:
            return True
    return False


def get_chrom_bnd(v):
    '''
        Returns the chromosome name of a BND variant.
        
        v -- Variant object
    '''
    c = v.alt.split(':')
    if '[' in c[0]:
        c = c[0].split('[')
        return c[1]
    else:
        c = c[0].split(']')
        return c[1]


def get_pos_bnd(v):
    '''
        Returns the position in ALT attribute of a BND variant.
        
        v -- Variant object
    '''
    c = v.alt.split(':')
    try:
        c = c[1].split(']')
        return int(c[0])
    except ValueError:
        c = c[0].split('[')
        return int(c[0])
    
    
def get_all_Bx(file,chrom,start,end):
    '''
        Returns all the barcodes and their position from a region (set).
        
        file -- a samfile
        chrom -- chromosome name
        start -- region's start position
        end -- region's end position
    '''
    all_bx = set()
    if start > end:
        start1 = end
        end1 = end
    else:
        start1 = start
        end1 = end
    for read in file.fetch(chrom,start1,end1):
        if read.has_tag('BX'):
            bx = read.get_tag('BX')
            bx = bx[:-2]
            pos = read.reference_start
            all_bx.add((bx,pos))
    return all_bx


def get_beg_bx(s):
    '''
        Returns the position of a linked-read.   
        
        s -- string (chrom:pos:length)
    '''
    s = s.split(':')
    return int(s[1])


def get_len_bx(s):
    '''
        Returns the length of a linked-read.

        s -- string (chrom:pos:length)
    '''
    s = s.split(':')
    return int(s[2])


def get_chrom_bx(s):
    '''
        Returns the chromosome of a linked-read.
        
        s -- string (chrom:pos:length)
    '''
    s = s.split(':')
    return s[0]


def isIsolated(pos,P,gap=N_GAP):
    '''
        Returns True if a barcode is isolated, else False.
        
        pos -- barcode's position
        P -- list resulting from partition() + clean_P()
    '''
    for [a,b] in P:
        if a - gap <= pos <= b + gap:
            return False
    return True


def partition(D,bx,c,gap=N_GAP):
    '''
        Returns the clusters for a barcode bx and their number of barcodes (list).
        
        D -- dict containing all barcodes
        bx -- barcode (string)
        c -- chromosome (string)
    '''
    s = D[bx].split(",")
    i = 0
    chrom = get_chrom_bx(s[0])
    while chrom != c:
        i += 1
        chrom = get_chrom_bx(s[i])
    beg = get_beg_bx(s[i])
    n = get_len_bx(s[i])
    P = [[beg,beg+n,1]]
    for j in range(i+1,len(s)):
        chrom = get_chrom_bx(s[j])
        if chrom == c:
          beg = get_beg_bx(s[j])
          if beg - P[-1][1] <= gap:
              P[-1][1] = beg + get_len_bx(s[j])
              P[-1][2] += 1
          else:
              P.append([beg,beg + get_len_bx(s[j]),1])
    return P
    
 
def clean_P(P,n=2):
    '''
        Removes all the clusters that do not have at least n barcodes.

        P -- list from partition()
    '''
    F = []
    for [a,b,c] in P:
        if c >= n:
            F.append([a,b])
    return F


def store_bx(bci):
    '''
        Reads a file and stores the barcodes in a dict.

        bci -- file, each line is as bx;chrom:pos:length
    '''
    D = {}
    with open(bci,"r") as filin:
        for line in filin:
            line = line.rstrip().split(";")
            D[line[0]] = line[1]
    return D


def nb_isolated(L,bci,D,c):
    '''
        Returns the number of isolated barcodes.
    
        L -- set of barcodes
        D -- dict resulting from store_bx()
    '''
    cpt = 0
    for (bx,pos) in L:
        P = partition(D,bx,c)
        P = clean_P(P)
        if isIsolated(pos,P):
            cpt += 1
    return cpt


def sortSV(vcf,bam,bci,truth,margin):
    '''
        Creates results.xlsx containing the number of barcodes for real
        variants and false one.
        
        vcf -- vcf file with variants
        bam -- bam file with reads mapping in the genome reference
        bci -- bci file got by LRez
        truth -- file with real variants
        margin -- boolean
    '''
    cpt = 1
    L = []
    row = row_bis = 0
    m = 100 if margin else 0
    realSV = trueSV(truth)
    samfile = pysam.AlignmentFile(bam,"rb")
    workbook = xlsxwriter.Workbook('results.xlsx')
    worksheet = workbook.add_worksheet()
    D = store_bx(bci)
    # Used to store current chromosomes for BND, and to output BND when changing chromosome
    curChr1 = ""
    curChr2 = ""
    with open(vcf,"r") as filin:
        # skips file's head :
        line = filin.readline()
        while line.startswith('#'):
            line = filin.readline()
        # for each variant :
        while line != '':
            v = Variant(line)
            print("variant",cpt)
            cpt += 1
            # We keep filling L if both chromosomes correspond to current one
            # If not, this means we're not processing the same variant anymore, so we treat the BND we've read so far
            if v.get_svtype() == "BND" and ((curChr1 == "" and curChr2 == "") or (curChr1 == v.chrom and curChr2 == get_chrom_bnd(v))):
                if L ==[]:
                    L.append([v.chrom,v.pos,-1])
                    L.append([get_chrom_bnd(v),get_pos_bnd(v),-1])
                    curChr1 = v.chrom
                    curChr2 = get_chrom_bnd(v)
                else:
                    if v.pos > L[0][2]:
                        L[0][2] = v.pos
                    if get_pos_bnd(v) > L[1][2]:
                        L[1][2] = get_pos_bnd(v)
            else:
                # we treat the BND variants :
                if L != [] and L[0][2] != -1 and L[1][2] != -1:
                    all_Bx = get_all_Bx(samfile,L[0][0],L[0][1],L[0][2])
                    all_Bx_pair = get_all_Bx(samfile,L[1][0],L[1][1],L[1][2])
                    # first BND variant is valid :
                    if isValid_bnd(L[0],realSV,m):
                        worksheet.write(row,0,L[0][0]+":"+str(L[0][1])+"-"+str(L[0][2]))
                        worksheet.write(row,1,nb_isolated(all_Bx,bci,D,L[0][0]))
                        row += 1
                    # first BND variant is not valid :
                    else:
                        worksheet.write(row_bis,3,L[0][0]+":"+str(L[0][1])+"-"+str(L[0][2]))
                        worksheet.write(row_bis,4,nb_isolated(all_Bx,bci,D,L[0][0]))
                        row_bis += 1
                    # second BND variant is valid :
                    if isValid_bnd(L[1],realSV,m):
                        worksheet.write(row,0,L[1][0]+":"+str(L[1][1])+"-"+str(L[1][2]))
                        worksheet.write(row,1,nb_isolated(all_Bx_pair,bci,D,L[1][0]))
                        row += 1
                    # second BND variant is not valid :
                    else:
                        worksheet.write(row_bis,3,L[1][0]+":"+str(L[1][1])+"-"+str(L[1][2]))
                        worksheet.write(row_bis,4,nb_isolated(all_Bx_pair,bci,D,L[1][0]))
                        row_bis += 1
                    L = []
                    # Update current chromosomes and L if we read a BND, otherwise set them / leave them empty
                    if v.get_svtype() == "BND":
                        L.append([v.chrom,v.pos,-1])
                        L.append([get_chrom_bnd(v),get_pos_bnd(v),-1])
                        curChr1 = v.chrom
                        curChr2 = get_chrom_bnd(v)
                    else:
                        curChr1 = ""
                        curChr2 = ""
                # treatment of a non BND variant :
		# Only do if we didn't read a BND
                if v.get_svtype() != "BND":
                    end = v.get_end()
                    all_Bx = get_all_Bx(samfile,v.chrom,v.pos,end)
                    # variant is valid :
                    if isValid(v,realSV,m):
                        worksheet.write(row,0,v.chrom+":"+str(v.pos)+"-"+str(end))
                        worksheet.write(row,1,nb_isolated(all_Bx,bci,D,v.chrom))
                        row += 1
                    # variant is not valid :
                    else:
                        worksheet.write(row_bis,3,v.chrom+":"+str(v.pos)+"-"+str(end))
                        worksheet.write(row_bis,4,nb_isolated(all_Bx,bci,D,v.chrom))
                        row_bis += 1
            line = filin.readline()
    workbook.close()
    samfile.close()
    

####################################################


parser = argparse.ArgumentParser(description='Sort SV')
parser.add_argument('-vcf', type=str, required=True, help='vcf file')
parser.add_argument('-bam', type=str, required=True, help='bam file')
parser.add_argument('-bci', type=str, required=True, help='bci file got by LRez')
parser.add_argument('-t', type=str, required=True, help='Truth file')
parser.add_argument('-m', action='store_true', help="Allows a margin to increase variants's length")
args = parser.parse_args()

if __name__ == '__main__':
    if args.m:
        sortSV(args.vcf,args.bam,args.bci,args.t,True)
    else:
        sortSV(args.vcf,args.bam,args.bci,args.t,False)
    
