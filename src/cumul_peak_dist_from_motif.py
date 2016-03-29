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
import string
import csv

import numpy as np

from operator import itemgetter
from matplotlib import pyplot as plt

def cumul_peak_dist_from_motif():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    
    #MANDATORY PARAMETERS
    parser.add_argument("motifs",help="Full paths to input motif-files.",type=str,nargs=1)
    parser.add_argument("outname",help="Full path to output file (image)",type=str,nargs=1)
    parser.add_argument("peaks",help="Full paths to input peak-files",type=str,nargs='+')

    #OPTIONAL PARAMETERS
    parser.add_argument("-p","--colors",help="Histogram colors, default=r",type=str,default='r',nargs='+')
    parser.add_argument("-l","--legend",type=str,nargs='+')
    parser.add_argument("-m","--maxdist",help="maximum distance between peak and motif",type=int,default=50,nargs=1)
    parser.add_argument("-N",help="Number of top peaks considered (default is all peaks, assumes input file is sorted in desired order.)",type=int,default=None,nargs=1)
    parser.add_argument("--center",help="If 1, only center to center distances calculated, otherwise all (default=0)",type=int,default=0,nargs=1)
    parser.add_argument("--log",help="If 1, the x-axis scale is log. If 0 (default), it is not.",type=int,choices=[0,1],default=0)

    args = parser.parse_args()

    center = {}
    left = {}
    right = {}

    for i in range(0,len(args.peaks)):
        #print "reading in motifs from ",
        print args.motifs[0]

        center[args.legend[i]] = [] #distances between motif center and peak center
        left[args.legend[i]] = [] #distances between motif left and peak left
        right[args.legend[i]] = [] #distances between motif right and peak right

        motifs = {}
        motiflen = 0
        #key = chromosome
        #value = [motifs center 1, motif center 2,...]
        with open(args.motifs[0],'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r:
                chrom = row[0]
                start = int(row[3])
                end = int(row[4])
                motiflen = end-start
                loc = (end-start)/2+start
                if chrom not in motifs: motifs[chrom] = set([loc])
                else: motifs[chrom].add(loc)

        #now reading in overlapping peaks
        #print "reading in peaks from ",
        print args.peaks[i]

        with open(args.peaks[i],'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            rowcount = 0
            for row in r:
                rowcount += 1
                if args.N!=None:
                    if rowcount>args.N[0]: break
                if row[0]=='chromosome': continue
                chrom = row[0]
                start = int(float(row[1]))
                end = int(float(row[2]))
                loc = (end-start)/2+start
                
                for j in range(0,args.maxdist[0]):
                    if (loc-j in motifs[chrom]) or (loc+j in motifs[chrom]):
                        #j is the distance between motif center and peak center
                        center[args.legend[i]].append(j)
                        isleft = True
                        isright = True
                        k = 0
                        while isleft and isright:
                            if k>1000:
                                if not isleft: left[args.legend[i]].append(k)
                                if not isright: right[args.legend[i]].append(k)
                                break
                            if ((start-k-motiflen/2) in motifs[chrom]) or (start+k-motiflen/2) in motifs[chrom]:
                                #k is the distance between left edge of motif and left edge of peak
                                left[args.legend[i]].append(k)
                                isleft = False
                            if ((end-k+motiflen/2) in motifs[chrom]) or (end+k+motiflen/2) in motifs[chrom]:
                                #k is the distance between right edge of motif and right edge of peak
                                right[args.legend[i]].append(k)
                                isright = False
                            k += 1
                        break
        #all peaks are read in
    #plotting the histograms

    #print "Plotting...",

    if args.center[0]==1:

        i = 0
        for leg in args.legend:
            #calculating the histogram of middle distances
            N = float(len(center[leg]))
            x = np.sort(np.array(center[leg]))
            y = np.array(range(int(N)))/N
            plt.plot(x,y,args.colors[i],label=leg,linewidth=2)
            i += 1
            #saving data to file
            with open(args.outname[0][:-4]+"_"+leg+".txt",'wb') as csvfile:
                w = csv.writer(csvfile,delimiter="\t")
                for j in range(0,len(x)): w.writerow([x[j],y[j]])

        plt.xlabel("Distance between centers")
        plt.ylabel("Cumulative fraction")
        if args.log==1: plt.xscale('log')
        plt.legend(loc='lower right')

    else:
        f,ax = plt.subplots(1,3,sharey=True)

        i = 0
        for leg in args.legend:
            #calculating the histogram of left edge distances
            N = float(len(left[leg]))
            x = np.sort(np.array(left[leg]))
            y = np.array(range(int(N)))/N
            ax[0].plot(x,y,args.colors[i],label=leg,linewidth=2)

            #calculating the histogram of middle distances
            N = float(len(center[leg]))
            x = np.sort(np.array(center[leg]))
            y = np.array(range(int(N)))/N
            ax[1].plot(x,y,args.colors[i],label=leg+"\n("+str(int(N))+" peaks)",linewidth=2)         

            #calculating the histogram of right edge distances
            N = float(len(right[leg]))
            x = np.sort(np.array(right[leg]))
            y = np.array(range(int(N)))/N
            ax[2].plot(x,y,args.colors[i],label=leg,linewidth=2)

            i += 1
        ax[0].set_xlabel("Distance between left edges")
        ax[0].set_ylabel("Cumulative fraction")
        if args.log==1: ax[0].set_xscale('log')

        ax[1].set_xlabel("Distance between centers")
        if args.log==1: ax[1].set_xscale('log')

        ax[2].set_xlabel("Distance between right edges")
        if args.log==1: ax[2].set_xscale('log')

        plt.subplots_adjust(bottom=0.08,top=0.87)
        ax[1].legend(loc='upper center', bbox_to_anchor=(0.5,1.17),ncol=3, fancybox=True, shadow=True)
    plt.savefig(args.outname[0])
    plt.clf()

    #print "done!" 
#end

cumul_peak_dist_from_motif()
