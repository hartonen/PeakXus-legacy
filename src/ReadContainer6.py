# Copyright (c) Tuomo Hartonen, 2015-2016
#
# THIS PROGRAM COMES WITH ABSOLUTELY NO WARRANTY!
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as 
# published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

import csv
from array import array
from scipy import stats
import sys
import pysam
import math
from time import clock
from operator import itemgetter

import numpy as np

class ReadContainer:

#CONSTRUCTOR
    def __init__(self):
        # intitialize an empty object
        self.reads = {} #a dictionary for reads, key=chromosome number.
        #Each chromosome has a dictionary with two keys, + and -, for positive and negative strands respectively.

        self.allreads = {} #a dictionary for saving all reads when umis are utilized
        
        self.umis = {}# a dictionary for umi-identifiers, key=chromosome number.
        #Each chromosome has a dictionary with two keys, + and -. Indices of the arrays for + and -
        #strands correspond to indices in self.reads for corresponding chromosome and strand.
        
        self.umi_histo = {}
        self.covered = {} #positions covered by reads
        self.read_number = 0 #total number of reads
        self.fragment_size = 0 #distance beween peak summits (read length)
        self.rejected = 0 #total number of rejected reads
    #end __init__

    def addRead(self, chrom, strand, fiveprime, save_umi=False, umi=None):
        #add read to dictionary
        #chrom = chromosome name
        #strand = +/-
        #fiveprime = coordinate of the 5'-end of the read
        #save_umi = if True, also the UMI is saved
        #umi = the corresponding UMI-sequence
        
        fiveprime = int(fiveprime)
        self.rejected += 1

        if umi not in self.umi_histo: self.umi_histo[umi] = 1
        else: self.umi_histo[umi] +=1

        if chrom not in self.reads:
            
            self.reads[chrom] = {}
            if save_umi!=None: self.allreads[chrom] = {}
            self.covered[chrom] = {}

            # store reads as an array of unsigned integers (4 bytes)
            self.reads[chrom]['+'] = array('i',[])
            self.reads[chrom]['-'] = array('i',[])
            if save_umi!=None: self.allreads[chrom]['+'] = array('i',[])
            if save_umi!=None: self.allreads[chrom]['-'] = array('i',[])
            self.covered[chrom]['+'] = set()
            self.covered[chrom]['-'] = set()
            
            self.covered[chrom][strand].add(int(fiveprime))
            self.reads[chrom][strand].append(fiveprime)
            if save_umi!=None: self.allreads[chrom][strand].append(fiveprime)
            
            #saving umi's if needed
            if save_umi!=None: 
                
                self.umis[chrom] = {}
                self.umis[chrom]['+'] = {}
                self.umis[chrom]['-'] = {}
                self.umis[chrom][strand][fiveprime] = set()
                self.umis[chrom][strand][fiveprime].add(umi)
                
        else:

            if save_umi!=None: self.allreads[chrom][strand].append(fiveprime)
            #first read of the strand is always saved
            if len(self.reads[chrom][strand])<1:
                self.reads[chrom][strand].append(fiveprime)
                self.covered[chrom][strand].add(int(fiveprime))
                self.read_number += 1
                if save_umi!=None:
                    self.umis[chrom][strand][fiveprime] = set()
                    self.umis[chrom][strand][fiveprime].add(umi)
                    
            else:
                if save_umi==None:
                    self.reads[chrom][strand].append(fiveprime)
                    self.covered[chrom][strand].add(fiveprime)
                    self.read_number += 1
                else:
                    #if there are no reads mapped to the same position, read is saved automatically
                    
                    if fiveprime not in self.covered[chrom][strand]:
                        self.reads[chrom][strand].append(fiveprime)
                        self.covered[chrom][strand].add(fiveprime)
                        self.umis[chrom][strand][fiveprime] = set()
                        self.umis[chrom][strand][fiveprime].add(umi)
                        self.read_number += 1
                    else:
                        #here we check if the same umi is already used for this position
                        if umi not in self.umis[chrom][strand][fiveprime]:
                            self.reads[chrom][strand].append(fiveprime)
                            self.covered[chrom][strand].add(fiveprime)
                            self.umis[chrom][strand][fiveprime].add(umi)
                            self.read_number += 1
        
    #end addRead

    def readSam(self,path,save_umi,outdir):
        #read a sam file and add contents to self
        #path = path to the input sam/bam file

        #fragment size is estimated from q_max first reads
        q_max = 1000
        frag_sizes = []

        #Reading in the reference umi-sequences from file
        if save_umi!=None:
            ref_umis = set()
            with open(save_umi,'rb') as csvfile:#f = open(save_umi,'r')
                r = csv.reader(csvfile,delimiter='\t')
                for line in r: ref_umis.add(line[0])
            
        ext = path.split('.')[-1]
        if ext=='sam': samfile = pysam.Samfile(path,'r')
        else: samfile = pysam.Samfile(path,'rb')

        iter = samfile.fetch() #iterator over the reads in sam-file

        for read in iter:
            if read.is_unmapped: continue
            chrom = samfile.getrname(read.rname)

            #start = read.pos
            #end = start+read.rlen
            strand = '+'
            if read.is_reverse: strand = '-'

            if save_umi!=None:
                index = [i for i,v in enumerate(read.tags) if v[0].count('BC')>0]
                umi = read.tags[index[0]][1]
                if umi not in ref_umis: continue
            else: umi = None

            if len(frag_sizes)<q_max: frag_sizes.append(read.pos+read.rlen-read.pos)

            #determine five prime end
            if strand == '+': fiveprime = read.aend-read.rlen+1
            else: fiveprime = read.aend#read.aend+read.rlen

            if fiveprime<0: continue
            self.addRead(chrom,strand,fiveprime,save_umi,umi)

        self.fragment_size = int(sum(frag_sizes)/q_max)
        
    #end readSam


    def readBed(self, path):
        #read a bed file and add contents to self
        #path = path to the input bed file

        #fragment size is estimated from q_max first reads
        q_max = 1000
        frag_sizes = []
        for row in csv.reader(open(path,'r'), delimiter='\t'):
            try:
                chrom = row[0]
                start = int(row[1])
                end = int(row[2])
                strand = row[5]

                if len(frag_sizes)<q_max: frag_sizes.append(end-start)
                    
                #determine five prime end
                if strand == '+': fiveprime = start
                else: fiveprime = end

                #add read to the container
                self.addRead(chrom, strand, fiveprime)
            except Exception as detail:
                sys.stderr.write("Input file is not in BED format!\n")
                sys.exit(1)

        #estimating fragment size  
        self.fragment_size = int(sum(frag_sizes)/q_max)    

    #end readBed

    def getCounts(self,chrom,start,end,count_type):
        #chrom = chromosome name
        #start = starting coordinate
        #end = ending coordinate
        #count_type = + if counts are |neg| + |pos|, - if |pos|-|neg|

        if start==end:
            c_pos = len(np.where(self.reads[chrom]['+']==start)[0])
            c_neg = len(np.where(self.reads[chrom]['-']==start)[0])                  
            
            if count_type=='-': return c_pos-c_neg
            return c_pos+c_neg

        c_pos = np.zeros(end-start)
        c_neg = np.zeros(end-start)
        index = 0
        for i in range(start,end):
            c_pos[index] = len(np.where(self.reads[chrom]['+']==i)[0])
            c_neg[index] = len(np.where(self.reads[chrom]['-']==i)[0])
            index += 1

        if count_type == '-': return c_pos-c_neg
        return  c_pos+c_neg
    #end getCounts

    def getAllCounts():
        #chrom = chromosome name
        #start = starting coordinate
        #end = ending coordinate
        #count_type = + if counts are |neg| + |pos|, - if |pos|-|neg|

        if start==end:
            c_pos = len(np.where(self.allreads[chrom]['+']==start)[0])
            c_neg = len(np.where(self.allreads[chrom]['-']==start)[0])                  
            
            if count_type=='-': return c_pos-c_neg
            return c_pos+c_neg

        c_pos = np.zeros(end-start)
        c_neg = np.zeros(end-start)
        index = 0
        for i in range(start,end):
            c_pos[index] = len(np.where(self.allreads[chrom]['+']==i)[0])
            c_neg[index] = len(np.where(self.allreads[chrom]['-']==i)[0])
            index += 1

        if count_type == '-': return c_pos-c_neg
        return  c_pos+c_neg
    #end getAllCounts
        
    def sortReads(self):
        #sort all reads in ascending order by the coordinate of the 5'-end while preserving the array data structure
        #also converts the arrays to numpy arrays for speed
        for chrom in self.reads.keys():
            #as sorted returns conversion back to array is required
            self.reads[chrom]['+'] = np.array(sorted(self.reads[chrom]['+']))
            self.reads[chrom]['-'] = np.array(sorted(self.reads[chrom]['-']))
    #end sortReads

    def getChromSize(self, chrom):
        #chromosome size to consider for scanning of both strands
        if (len(self.reads[chrom]['+'])>0) and len(self.reads[chrom]['-'])>0:
            return max([self.reads[chrom]['-'][-1],self.reads[chrom]['+'][-1]])+1, min([self.reads[chrom]['-'][0],self.reads[chrom]['+'][0]])
        else: return 0,0
    #end getChromSize

    def getGenomeSize(self):
        #genome size to consider for scanning of both strands
        genome_size = 0
        for chrom in self.reads.keys(): genome_size += self.getChromSize(chrom)
        return genome_size
    #end getGenomeSize 

    def getReads(self, chrom, strand):
        #returns all reads for chromosome chrom
        #READS ARE RETURNED AS NUMPY ARRAY FOR FASTER ACCESSING
        if strand==None:
            p = set(self.reads[chrom]['+'])
            n = set(self.reads[chrom]['-'])
            u = p.union(n)
            return sorted(list(u))
        if chrom in self.reads: return np.array(self.reads[chrom][strand])
        else:
            sys.stderr.write("Chromosome "+chrom+" not found from data!\n")
            sys.exit(1)
    #end getReads

    def getReadHisto(self,chrom,hist_type,s,nproc):
        #this is for getting histogram of reads when UMIs are not used
        #chrom = chromosome name
        #s = sliding window size used

        N,n = self.getChromSize(chrom)
        if hist_type=='++' or hist_type=='-' or hist_type=='+':
            h_pos = np.bincount(self.reads[chrom]['+'],minlength=N)
            if s>1: h_pos = np.convolve(h_pos,np.ones(s),mode='same')
            if hist_type=='++':return h_pos

        h_neg = np.bincount(self.reads[chrom]['-'],minlength=N)
        if s>1: h_neg = np.convolve(h_neg,np.ones(s),mode='same')

        if hist_type=='-': return h_pos-h_neg
        if hist_type=='--': return h_neg
        return h_pos+h_neg
    #end getReadHisto

    def getAllReadHisto(self,chrom,hist_type,s,nproc):
        #this is for getting histogram of all reads when umis are used
        #meaning these histograms are not filtered for duplicates
        #chrom = chromosome name
        #s = sliding window size used
        
        if len(self.allreads.keys())<1: return self.getReadHisto(chrom,hist_type,s,nproc)
        
        N,n = self.getChromSize(chrom)
        if hist_type=='++' or hist_type=='-' or hist_type=='+':
            h_pos = np.bincount(self.allreads[chrom]['+'],minlength=N)
            if s>1: h_pos = np.convolve(h_pos,np.ones(s),mode='same')
            if hist_type=='++': return h_pos

        h_neg = np.bincount(self.allreads[chrom]['-'],minlength=N)
        if s>1: h_neg = np.convolve(h_neg,np.ones(s),mode='same')

        if hist_type=='-': return h_pos-h_neg
        if hist_type=='--': return h_neg
        return h_pos+h_neg
    #end getAllReadHisto

    def getReadCount(self, chrom, strand):
        #find out how many reads are mapped to a particular chromosomome and strand
        return len(self.reads[chrom][strand])
    #end getReadCount

    def getChromNames(self):
        #retrieve a sorted list of all chromosome names
        return sorted(self.reads.keys())
    #end getChromNames

    def getFragmentSize(self): return self.fragment_size

    def getUmis(self,chrom,strand):
        if len(self.umis.keys())==0: return None
        else: return np.array(self.umis[chrom][strand])

    def zeroOrderLeastSquares(self,data):
        #Fits a zeroth order polynomial to data, calculates the MLE
        #for the parameter theta_0, and returns the P-value

        N = float(len(data)) #number of data points
        m = 1.0 #number of parameters

        #MLE for the constant parameter
        theta_0 = float(sum(data))/N
        #test statistics:
        chi2_min = sum([(y-theta_0)*(y-theta_0) for y in data])
        
        P = 1-stats.chi2.cdf(chi2_min,N-m)
        if math.isnan(P):
            return 1
        return P
