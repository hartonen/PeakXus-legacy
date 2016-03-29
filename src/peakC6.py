#!/usr/bin/env python
# Copyright (c) Tuomo Hartonen, 2015-2016
#
# THIS PROGRAM COMES WITH ABSOLUTELY NO WARRANTY!
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License version 2 as 
# published by the Free Software Foundation.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see http://www.gnu.org/licenses/.

import argparse
import sys
import csv
import numpy as np

from time import time
from os import listdir
from glob import glob

from ReadContainer6 import ReadContainer
from analyzeTransitions6_symmetric import tPoints

def peakC():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    #MANDATORY PARAMETERS
    #infiles = path to a text file that contains full paths of all input files in separate rows
    parser.add_argument("infiles", help="full paths to input files. File type is determined by the file name extension. Compatible file types are .bam and .sam.", type=str,nargs='+')
    #outdir = directory for output
    parser.add_argument("outdir", help="full path to output directory", type=str,nargs=1)

    #OPTIONAL PARAMETERS
    #if these are not given, defaults are used

    parser.add_argument("-m","--mindist",help="minimum distance between two peaks (default=30)",type=int,default=30)
    parser.add_argument("-w", "--winsize", help="width of the window around 5'-end (default=5)", type=int, default=5)
    parser.add_argument("-l_l", "--l_limit", help="l_limit+w is the largest allowed width for a peak (default=60)", type=int, default=60)
    parser.add_argument("-u_l", "--u_limit", help="u_limit is the smallest allowed peak width (default=5)", type=int,default=5)
    parser.add_argument("-b", "--b_threshold", help="Threshold p-value for background model. Default is 0.05", type=float, default=0.05)
    parser.add_argument("-r", "--saveregions", help="If 1, read densities for all peak regions are output into a file, if 0, not (default=0).",type=int, choices=[0,1], default=0)
    parser.add_argument("-u", "--readumis", type=str, help="Full path to the file containing the UMI-sequences (1 per row, default=None).",default=None)
    parser.add_argument("-v", "--verbosity", help="verbosity level, 1=verbose on, 0=verbose off", type=int, choices=[0,1], default=1)
    parser.add_argument("-t", "--timing", help="Timing on(1)/off(0)", type=int, choices=[0,1], default=1)
    parser.add_argument("-y","--yates",help="1, if chi2-test with Yates' correction is used, 0 meaning G-test is used. 2=KS-test (default=0)",type=int,choices=[0,1,2],default=0)
    parser.add_argument("-d","--allpairs",help="1, if all positive peaks are considered for each candidate site, 0 if only the highest (default=0)",type=int,choices=[0,1],default=0)
    parser.add_argument("-s","--sliding_window",help="Width of a sliding window used in averaging density of 5' ends of reads. This should be an odd number (default=1).",type=int,default=1)
    parser.add_argument("-n","--nproc",help="Number of parallel processes used in calculating significance of peaks (default=1).",type=int,default=1)
    parser.add_argument("-o","--onlysignal",help="Only report the top 1000 regions with highest read count for each chromosome if o=1 (default=0)",type=int,choices=[0,1],default=0)
    parser.add_argument("-p","--pseudocount",help="Pseudocount added to each position when testing significance of candidate binding sites (default=5.0).",type=float,default=5.0)

    args = parser.parse_args()
    logfile = args.outdir[0]+"/log_short.txt"

    #################################
    #Testing input, setting defaults#
    #################################


    w = args.winsize
    l_limit = args.l_limit
    u_limit = args.u_limit
    b = args.b_threshold
    mindist = args.mindist
    s = args.sliding_window

    #counter for total number of test conducted
    test_counter = 0
    save_umi = args.readumis

    if args.verbosity==1:
        f = open(logfile,'a')
        f.write("\n")
        f.write("\n")
        f.write("****************\n")
        f.write("*Parameters are*\n")
        f.write("****************\n")
        f.write("* minimum distance between peaks:   mindist  ="+str(mindist)+"\n")
        f.write("* window width:                     w        = "+str(w)+"\n")
        f.write("* lower limit of 5' end of + reads: l_limit  = "+str(l_limit)+"\n")
        f.write("* upper limit of 5' end of - reads: u_limit  = "+str(u_limit)+"\n")
        f.write("* background threshold:             b        = "+str(b)+"\n")
        f.write("* width of sliding window:          s        = "+str(s)+"\n")
        f.write("* number of CPU's used:             nproc    = "+str(args.nproc)+"\n")
        f.write("\n")
        f.write("\n")
        f.write("Input files: "+str(args.infiles)+"\n")
        f.write("Output directory: "+args.outdir[0]+"\n")
        f.write("\n")
        f.write("\n")
        if args.verbosity==1: f.write("* Verbosity is on.\n")
        else: f.write("* Verbosity is off.\n")
        if args.timing==1: f.write("* Timing is on.\n")
        else: f.write("* Timing is off.\n")
        if save_umi: f.write("* Using UMIs to detect artefacts.\n")
        else: f.write("* Not using UMIs.\n")
        f.write("\n")
        f.write("\n")
        f.close()

    #testing that the input files open

    if args.verbosity==1:
		f = open(logfile,'a')	
		f.write("Testing that input files exist...")
		f.close()

    try:
        infilenames = args.infiles
        #testing that all the individual files open
        for name in infilenames:
            f = open(name,'r')
            f.close()
        #testing that the output directory is writable
        f = open(args.outdir[0]+"all_transition_points.igv",'w')
        f.close()
        if args.verbosity==1:
			f = open(logfile,'a')
			f.write("done!\n")
			f.close()

    except Exception as detail:
        if args.verbosity==1:
            sys.stderr.write("Input file doesn't exist or output directory is not writable!")
            f = open(logfile,'a')
            f.write(str(detail))
            f.write("\n")
            f.close()
        sys.exit(1)

    
    for infile in infilenames:

        if args.verbosity==1:
			f = open(logfile,'a')
			f.write("Processing "+infile+"\n")
			f.close()

        #################
        #Reading in data#
        #################

        #determining if input file is .bed or .sam
        ext = infile.split('.')
        ext = ext[-1]

        if args.timing==1:
            times = []
            start = time()

        if args.verbosity==1:
			f = open(logfile,'a')		
			f.write("Reading in data...")
			f.close()

        reads = ReadContainer()
        if ext.count('bed')>0:
            reads.readBed(infile)
            reads.sortReads()
        elif ext.count('sam')>0: 
            reads.readSam(infile,save_umi,args.outdir[0])
            reads.sortReads()
        elif ext.count('bam')>0:
            reads.readSam(infile,save_umi,args.outdir[0])
            reads.sortReads()
        else:
            f = open(logfile,'a')
            f.write("Input file not in .sam or .bed format! Check the extension of the file.\n")
            f.close()
            continue

        if args.verbosity==1:
			f = open(logfile,'a')
			f.write("done!\n")
			f.close()

        if args.timing==1:
            end = time()
            times.append(end-start)
            if args.verbosity==1:
				f = open(logfile,'a')
				f.write("Reading data took "+str(end-start)+" s\n")
				f.close()

        #Running the rest of the code just for one chromosome at a time saves memory

        if True:
            #running for one chromosome at a time
            for chrom in reads.getChromNames():
                print chrom+" ",
                if args.verbosity==1:
                    f = open(logfile,'a')
                    f.write("**************************************\n")
                    f.write("CHROMOSOME: "+chrom+"\n")
                    f.close()

                #######################
                #Creating peak regions#
                #######################

                if args.verbosity==1:
					f = open(logfile,'a')
					f.write("Creating peak regions...")
					f.close()
                if args.timing==1: start = time()

                peak_regions = tPoints(w,l_limit,u_limit,mindist,s)
                if args.onlysignal==1: peak_regions.calcTransitions_counts(chrom,reads,args.allpairs,args.nproc)
                else: peak_regions.calcTransitions(chrom,reads,args.allpairs,args.nproc)

                if args.verbosity==1:
					f = open(logfile,'a')
					f.write("done!\n")
					f.close()
                if args.timing==1:
                    end = time()
                    times.append(end-start)
                    if args.verbosity==1:
						f = open(logfile,'a')
						f.write("Creating peak regions took "+str(end-start)+" s\n")
						f.close()

                #######################
                #Estimating background#
                #######################
 
                if b>0:
        
                    if args.verbosity==1:
                                                f = open(logfile,'a')
                                                f.write("Calculating background...")
                                                f.close()
                    if args.timing==1: start = time()

                    if args.onlysignal==0:
                        peak_regions.testBackground(chrom,b,args.yates,args.nproc,args.pseudocount)
                        test_counter += peak_regions.getNumTests()
                    if args.verbosity==1:
                                                f = open(logfile,'a')
                                                f.write("done!\n")
                                                f.close()
                    if args.timing==1:
                        end = time()
                        times.append(end-start)
                        if args.verbosity==1:
                                                        f = open(logfile,'a')
                                                        f.write("Estimating background took "+str(end-start)+" s\n")
                                                        f.close()
                tRegions = peak_regions.getTRegions(chrom)
            
                #####################
                #Saving called peaks#
                #####################

                if args.saveregions==1: regionfile = args.outdir[0]+"region_w="+str(w)+"_l_limit="+str(l_limit)+"_k="+str(k)+"_u_limit="+str(u_limit)+"_"+chrom+".csv"
                else: regionfile = None
                
                f = open(args.outdir[0]+chrom+"_transition_points.igv",'w')
                if args.onlysignal==1:
                    f.write("chromosome\tstart\tend\tid\tsignal\n")
                    for t in range(0,len(tRegions)): f.write(chrom+"\t"+str(tRegions[t,0]-1)+"\t"+str(tRegions[t,1])+"\t"+str(tRegions[t,0]-1)+"-"+str(tRegions[t,1])+"\t"+str(tRegions[t,2])+"\n")
                else:
                    f.write("chromosome\tstart\tend\tid\tsignal\tp-value\tscore\n")
                    for t in range(0,len(tRegions)): f.write(chrom+"\t"+str(tRegions[t,0]-1)+"\t"+str(tRegions[t,1])+"\t"+str(tRegions[t,0]-1)+"-"+str(tRegions[t,1])+"\t"+str(sum(tRegions[t,3:-1]))+"\t"+str(tRegions[t,-2])+"\t"+str(tRegions[t,-1]*sum(tRegions[t,3:-1]))+"\n")
                f.close()

                if regionfile != None:
                    if args.verbosity==1:
                        f = open(logfile,'a')				
                        f.write("done!\n")
                        f.write("Saving regions...")
                        f.close()

                    Regions = np.zeros((tRegions.shape[0],2*(l_limit+2*w)+1))
                    for r in range(0,len(tRegions)): Regions[r,:] = tRegions[r,2:-1]
                    np.savetxt(regionfile,Regions,delimiter=',')

                    if args.verbosity==1:
                        f = open(logfile,'a')
                        f.write("done!\n")
                        f.close()

            #analysis done for one chromosome
            
        else:
            #this means analysis is done for all chromosomes simultaneously
			f = open(logfile,'a')            
			f.write("Not implemented\n")
			f.close()

    #total time
    if args.verbosity==1 and args.timing==1:
		f = open(logfile,'a')
		f.write("Total duration "+str(sum(times))+" s")
		f.close()
    #saving the total number of tests conducted
    f = open(args.outdir[0]+"numtests.txt",'w')
    f.write(str(test_counter))
    f.close()

#end

peakC()
