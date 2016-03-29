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
import csv
import string

import numpy as np
from operator import itemgetter
from glob import glob
import matplotlib
from matplotlib import pyplot as plt
from random import randint

def consMotifHitsPerPeak():
    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    
    #MANDATORY PARAMETERS
    parser.add_argument("indirs",help="full path to the file containing the full paths to input directories. Each directory contains input from one run. First column is the directory, second column is the legend text.",type=str)
    parser.add_argument("outname",help="Full path for the output file.",type=str)
    
    #OPTIONAL PARAMETERS
    #if these are not given, defaults are used
    parser.add_argument("-N","--numpeaks",help="How many of the top peaks are analyzed, default=1000.",type=int,default=1000)
    parser.add_argument("-m","--maxdist",help="Maximum distance of peak summit from motif border, default=20.",type=int,default=20)
    parser.add_argument("-s","--sorting",help="Sorting order,1=signal,2=p-value, 3=random, 4=score (default=1),0=already sorted, no sorting is done.",type=int,default=1,choices=[0,1,2,3,4])
    parser.add_argument("-o","--motif_once",help="If 1, each motif is only counted matching to one peak, if 0, there is no limit (default=0).",type=int,default=0)
    parser.add_argument('--title',type=str,default=None)

    #parameters for plotting the random slope
    parser.add_argument("-S",help="Motif size (default=24bp).",type=float,default=24.0)
    parser.add_argument("-t",help="Matching tolerance (default=40bp).",type=float,default=40.0)
    parser.add_argument("-M",help="M top motifs considered (default=300000).",type=int,default=300000)

    args = parser.parse_args()

    matplotlib.rcParams.update({'font.size': 18})

    N = args.numpeaks
    m = args.maxdist
    slope = args.S/(2867188341.0/args.M-args.t-args.S)

    indirs = []

    #print "Reading in input directories...",
    #with open(args.indirs,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter=',')
        for line in r: indirs.append(line)
    #print "done!"
    #print "Reading in peaks...",

    peaks = {}
    motif_sizes = {}
    for key in indirs:
        peaks[key[0]] = []
        motif_sizes[key[0]] = 0.0
        if args.sorting==0: inname = glob(key[0]+"*all_transition_points_sorted.igv")[0]
        else: inname = glob(key[0]+"*all_transition_points.igv")[0]
        print inname+"...",
        with open(inname,'rb') as csvfile:
            r = csv.reader(csvfile,delimiter='\t')
            for row in r:
                if row[0].count('chromosome')>0: continue
                if len(row)==5: signal = 0.0
                start = int(float(row[1]))
                end = int(float(row[2]))
                motif_sizes[key[0]] += end-start
                if end-start>1: summit = start+(end-start)/2
                else: summit = start
                
                if len(row)==6:
                    signal = float(row[-2])
                    p_val_score = float(row[-1])

                if len(row)==7:
                    signal = float(row[-3])
                    p_val_score = float(row[-1])

                if signal<0.001: peaks[key[0]].append([row[0],summit,float(row[-1])]) #MACE gives only p-value
                else: peaks[key[0]].append([row[0],summit,p_val_score,signal])
        #print "# of peaks: "+str(len(peaks[key[0]]))

    #print "done!"
    #print "Sorting the peaks...",

    for key in peaks.keys():
        if args.sorting==3:
            #random peaks
            used = set()
            new_peaks = []
            for i in range(0,N):
                r = randint(0,len(peaks[key])-1)
                while r in used: r = randint(0,len(peaks[key])-1)
                new_peaks.append(peaks[key][r])
            
            peaks[key] = new_peaks
            continue

        if len(peaks[key][0])==3: peaks[key] = sorted(peaks[key],key=itemgetter(2))[:N]
        elif args.sorting!=0:
            if args.sorting==1: peaks[key] = sorted(peaks[key],key=lambda x: (-x[3],x[2]))[:N]
            elif args.sorting==2: peaks[key] = sorted(peaks[key],key=lambda x: (x[2],-x[3]))[:N]
            else:
                #sorting by score
                if key.find("peakzilla")>-1: peaks[key] = sorted(peaks[key],key=itemgetter(2),reverse=True)[:N]
                else: peaks[key] = sorted(peaks[key],key=lambda x: (-x[2],-x[3]))[:N]
        else: peaks[key] = peaks[key][:N]
    #print "done!"
    #print "Plotting...",

    x_axis = range(1,N,10)
    y_axes = []
    y_rand = []
    sizes = []

    for ind in indirs:
        y_axes.append([])
        y_rand.append([])
        S = motif_sizes[ind[0]]
        S = 40.0
        #conserved motif-file is the fourth column of the input file
        motif_coords = {}
        with open(ind[3],'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r:
                chrom = row[0]
                loc = [int(float(row[3])),int(float(row[4]))]
                if chrom not in motif_coords: motif_coords[chrom] = set(range(loc[0]-m,loc[1]+m))
                else:
                    for s in range(loc[0]-m,loc[1]+m): motif_coords[chrom].add(s)
            
        #motifs read in

        for x in x_axis:
            if args.motif_once!=0:
                deleted = {}
                for c in motif_coords: deleted[c] = set([])

            if x>len(peaks[ind[0]]):
                y_axes[-1].append(0.0)
                y_rand[-1].append(slope*x)
                continue
            #counting the hits
            count = 0.0
            peak_set = peaks[ind[0]][:x]
            for peak in peak_set:
                chrom = peak[0]
                loc = peak[1]

                if chrom not in motif_coords: continue
                if args.motif_once!=0:
                    if loc in deleted[chrom]: continue
                if loc in motif_coords[chrom]:
                    if args.motif_once!=0:
                        #deleting the coordinates that belong to corresponding motif
                        i = 1
                        while loc+i in motif_coords[chrom]:
                            deleted[chrom].add(loc+i)
                            i += 1
                        i = 1
                        while loc-i in motif_coords[chrom]:
                            deleted[chrom].add(loc-i)
                            i += 1
                        deleted[chrom].add(loc)
                    #hit found
                    count += 1.0
            y_axes[-1].append(count)
            y_rand[-1].append(slope*x)
            
    styles = ['-r','-b','-g','-m','-k','-c','-y','--r','--b','--g','--m','--k','--c','--y','-.r','-.b','-.g','-.m','-.k','-.c','-.y',':r',':b',':g',':m',':k',':c',':y','-*r','-*b','-*g','-*m','-*k','-*c','-*y','-Dr','-Db','-Dg','-Dm','-Dk','-Dc','-Dy','-or','-ob','-og','-om','-ok','-oc','-oy']
    styles_rand = [':k',':k',':k',':m',':k',':c',':y']
    for i in range(0,len(y_axes)):
        #saving the results as a text file
        f = open(args.outname[:-4]+"_"+indirs[i][1]+".txt",'w')
        for j in range(0,len(x_axis)): f.write(str(x_axis[j])+"\t"+str(y_axes[i][j])+"\n")
        f.close()
        plt.plot(x_axis,y_axes[i],styles[i],linewidth=2,label=indirs[i][1])
        plt.plot(x_axis,y_rand[i],styles_rand[0],linewidth=2)
    plt.xlim((0,max(x_axis)))
    plt.ylim((0,max([max(y) for y in y_axes])))
    plt.xlabel('N strongest peaks')
    plt.ylabel('Number of motifs found')
    plt.legend(loc=4)

    plt.suptitle(args.title,y=0.93)
    plt.tight_layout()
    plt.savefig(args.outname)
    plt.clf()
    #print "done!"
#end

consMotifHitsPerPeak()
