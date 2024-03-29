#!/usr/bin/env python
import argparse
import csv

import numpy as np

from operator import itemgetter

def FDR():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    
    #MANDATORY PARAMETERS
    parser.add_argument("peakfile",help="full path to the file containing called peaks in igv-format.",type=str)
    parser.add_argument("outfile",help="Full path to output file (.igv)",type=str)

    #OPTIONAL PARAMETERS
    #if these are not given, defaults are used
    parser.add_argument("-t","--input_type",help="1=peakC (default), 2=MACE",type=int,choices=[1,2],default=1)
    parser.add_argument("-m","--num_tests",help="Total number of tests conducted. If not given, the input file is assumed to contain all the test results.",type=int,default=None)
    args = parser.parse_args()

    #print "Reading in peaks...",
    data = [] # = [[row1],[row2],...]
    #p-value is always the last entry of the row

    with open(args.peakfile,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter="\t")
        for row in r:
            if row[0]=="chromosome": continue
            for i in range(1,len(row)):
                if i==3: continue
                row[i] = float(row[i])
            data.append(row)

    #print "done!"

    #sorting the list according to p-value (secondary by score)
    if args.input_type==2: data = sorted(data,key=lambda x: x[-1])
    else: data = sorted(data,key=lambda x: (x[5],-x[6]))

    #print data[:10]


    #print "Calculating the FDR's...",

    if args.num_tests==None: m = len(data)
    else: m = args.num_tests

    c = 1#sum([1/i for i in range(1,m)])
    for k in range(1,len(data)+1):
        P_k = data[k-1][-2]
        alpha = P_k*m*c/k
        data[k-1].append(alpha)


    #saving the results as an igv-file sorted by chromosome and location
    data = sorted(data,key=lambda x: (x[0],x[1]))

    with open(args.outfile,'wb') as csvfile:
        w = csv.writer(csvfile,delimiter="\t")
        
        if args.input_type==1: w.writerow(["chromosome","start","end","id","signal","p-value","score","FDR"])
        else: w.writerow(["chromosome","start","end","name","p-value","FDR"])
        for row in data: w.writerow(row)
#FDR ends

FDR()
