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
import pysam

import matplotlib
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
from operator import itemgetter
from os import system

from ReadContainer6 import ReadContainer
import clusterRegions6 as cr
#Functions for plotting fancy results from peak calling

def fancyResults():

    ########################
    #command line arguments#
    ########################

    parser = argparse.ArgumentParser()
    
    #MANDATORY PARAMETERS
    parser.add_argument("peakfile",help="full path to an igv-file containing called peaks.",type=str)
    parser.add_argument("readfile",help="full path to the bam-file used as input for peak calling.",type=str)
    parser.add_argument("motiffile",help="full path to a tsv-file containing motif locations.",type=str)
    parser.add_argument("outdir",help="full path to output directory",type=str)
    parser.add_argument("peakseqs",help="full path to a fasta-file containing the sequences under top N peaks.",type=str)
    
    #OPTIONAL PARAMETERS
    #if these are not given, defaults are used
    parser.add_argument("-N","--numpeaks",help="Number of top peaks considered",type=int,default=1000)
    parser.add_argument("-c","--comparison",help="0=calculating distance between middle of the peak and middle of the motif, 1=calculating distance between right end of the motif and right end of the peak.",type=int,choices=[0,1],default=0)
    parser.add_argument("-m","--maxdist",help="maximum distance between peak and motif (default=50).",type=int,default=50)
    parser.add_argument("-w","--peakwidth",help="width of region around the peak in heatmaps (default=200bp's)",type=int,default=200)
    parser.add_argument("-e","--expname",help="experiment name for plot titles.",type=str,default="ChIP-exo")
    parser.add_argument("-s","--sorting",help="Sorting order of peaks, 0=score(default), 1=statistics, 2=signal, 3=p-value, -1=already sorted.",type=int,choices=[-1,0,1,2,3],default=0)
    parser.add_argument("-n","--nproc",help="Number of parallel processes used in clustering (default=1).",type=int,default=1)
    parser.add_argument("-k",help="Number of clusters (default=4).",type=int,default=4)
    parser.add_argument("-p","--npass",help="number of times clustering is run (default=1).",type=int,default=1)
    parser.add_argument("-t","--test",help="Statistical test used, 0=G-test, 2=KS-test (default=0).",type=int,choices=[0,2],default=0)
    parser.add_argument("-H","--skipheatmap",help="1, if plotting heatmaps is omitted,0 (=default) otherwise.",type=int,choices=[0,1],default=0)
    parser.add_argument("-S","--smoothing",help="Window size used in calculating read distributions (default=1).",type=int,default=1)
    parser.add_argument("-ms","--motif_start",help="Start position of the motif, if middle of the motif is at 0",type=int,default=None)
    parser.add_argument("-ml","--motif_len",help="Motif length.",type=int,default=None)
    parser.add_argument("--truereads",help="1, if only reads pointing towards the peak summit are plotted, 0 if all reads are plotted (default=0).",type=int,choices=[0,1],default=0)
    parser.add_argument("--pseudo",help="Pseudocount value for g-test.",type=float,default=10.0)
    parser.add_argument("--ctest",help="g if g-test for clustering, euc if euclidian distance",type=str,choices=['g','euc'],default='g')

    #score-columns in peakfile:
    #row[-1]=score
    #row[-2]=statistics
    #row[-3]=signal

    args = parser.parse_args()
    N = args.numpeaks
    c = args.comparison
    m = args.maxdist
    w = args.peakwidth
    h = args.skipheatmap
    S = args.smoothing
    sortindex = -1-args.sorting
    if args.sorting==3: sortindex = -1

    matplotlib.rcParams.update({'font.size': 26})
    matplotlib.rcParams.update({'font.weight': 'bold'})

    #################
    #READING IN DATA#
    #################

    #reading in peaks
    #print "Reading in peaks...",
    peaks = []
    chroms = set()
    with open(args.peakfile,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        for row in r:
            if row[0].count('chromosome')>0: continue
            peaks.append(row)
            if row[0] not in chroms: chroms.add(row[0])
    print "done!"
    
    if args.sorting==-1:
        topPeaks = peaks[:N]
    if args.sorting==3:
        peaks = sorted(peaks,key=itemgetter(sortindex),reverse=False)
        topPeaks = peaks[:N]
    else:
        peaks = sorted(peaks,key=itemgetter(sortindex),reverse=True)
        topPeaks = peaks[:N]

    #reading in motifs
    #print "Reading in motifs...",
    motifs = {}
    motif_size = 0
    #key = chrom
    #value = [[+ strand motifs],[- strand motifs]]
    with open(args.motiffile,'rb') as csvfile:
        r = csv.reader(csvfile,delimiter='\t')
        for row in r:
            chrom = row[0]
            if chrom not in chroms: continue
            start = int(float(row[3]))
            end = int(float(row[4]))
            strand = row[5]
            if c==0: loc = (end-start)/2+start
            elif c==1: loc = end
            elif c==2: loc = start
            
            motif_size = end-start
            if chrom not in motifs:
                if strand=='+': motifs[chrom] = [set([loc]),set()]
                else: motifs[chrom] = [set(),set([loc])]
            else:
                if strand=='+': motifs[chrom][0].add(loc)
                else: motifs[chrom][1].add(loc)
    #print "done!"

    #reading in sequences
    #print "Reading in reads...",
    #determining if input file is .bed or .sam
    samfile = pysam.Samfile(args.readfile,'rb')
    #print "done!"

    ########################
    #HEATMAP OF TOP N PEAKS#
    ########################
    
    #print "Plotting heatmap of top "+str(N)+" peaks...",

    #creating matrix of the read densities of top peaks
    E = np.zeros((2*N,w+1))
    ind = 0

    for p in topPeaks:
        chrom = p[0]
        start = int(float(p[1]))
        end = int(float(p[2]))
        summit = (end-start)/2+start
        for read in samfile.fetch(chrom,summit-w/2-100,summit+2/2+1+100):
            if read.is_unmapped: continue
            strand = '+'
            if read.is_reverse: strand = '-'
            
            #determine five prime end
            if strand == '+': fiveprime = read.aend-read.rlen+1
            else: fiveprime = read.aend

            if fiveprime in range(summit-w/2,summit+w/2+1):
                if strand=='+':
                    if args.truereads==1:
                        #only reads starting from the left side of the summit are considered
                        if fiveprime<summit: E[ind,fiveprime-(summit-w/2)] += 1.0
                    else: E[ind,fiveprime-(summit-w/2)] += 1.0
                else:
                    if args.truereads==1:
                        #only reads starting from the right side of the summit are considered
                        if fiveprime>summit: E[ind+1,fiveprime-(summit-w/2)] -= 1.0
                    else: E[ind+1,fiveprime-(summit-w/2)] -= 1.0
        
        ind += 2
    #plotting

    #Smoothing using np.convolve
    if S>1:
        for i in range(0,N): E[i,:] = np.convolve(E[i,:],np.ones(S),'same')/float(S)

    #print "max count="+str(np.max(E))
    ##########################
    #CLUSTERING THE TOP PEAKS#
    ##########################

    x = np.array([i for i in range(-w/2,w/2+1)])

    y = np.array([i for i in range(0,int(np.shape(E)[0]))])

    if args.skipheatmap==2:
        E_clustered,centroids,clusters = cr.clusterRegions(E[range(0,len(E),2),:],np.abs(E[range(1,len(E),2),:]),args.k,args.npass,args.pseudo,args.nproc,args.ctest)
        #plotting cluster centroids

        ind = 0
        for c in centroids:
            ind += 1
            plot1d(x,[c[1],c[2]],args.outdir+"Centroid"+str(ind)+"_top_"+str(N)+".png",xlabel="Distance from peak summit",ylabel="Average ChIP-exo reads",colors=['r','b'],title="Centroid="+str(ind)+", "+str(len(np.where(clusters==c[0])[0]))+" peaks")
        #plotting heatmap of each cluster
        ind = 0
        for c in set(clusters):
            ind += 1
            E_aux = E_clustered[np.where(E_clustered[:,-1]==c)]
            plotHeatMap(E_aux[:,:-1],x,np.array([i for i in range(0,int(np.shape(E_aux[:,:-1])[0]))]),args.outdir+"Heatmap_top_"+str(N)+"cluster"+str(ind)+".png",xlabel="Distance from peak summit (bp's)",ylabel="Cluster"+str(ind),xticks=[-w/2,w/2])

        plotHeatMap(E,x,y,args.outdir+"Heatmap_top_"+str(N)+".png",xlabel="Distance from peak summit (bp's)",ylabel="Top "+str(N)+" peaks",xticks=[-w/2,w/2])

    #print "done!"

    #############################################
    #AVERAGE READ DENSITY PROFILE OF TOP N PEAKS#
    #############################################

    #print "Plotting read density profile of top "+str(N)+" peaks...",
    plot1d(x,[np.mean(np.abs(E[range(0,len(E),2)]),axis=0),np.mean(np.abs(E[range(1,len(E),2)]),axis=0)],args.outdir+"Avg_top_"+str(N)+".png",colors=['r','b'],xlabel="Distance from peak summit",ylabel="Average 5'-end count",xticks=[-w/2,0,w/2])   
    #print "done!"

    return 0
    ###########################################
    #HEATMAP OF TOP N PEAKS WITH A MOTIF MATCH#
    ###########################################

    #print "Plotting heatmap of top"+str(N)+" peaks overlapping with motif...",
    #selecting top N peaks overlapping with a motif
    topPeaks = []
    dists_from_summit = []
    strands = []
    hitcount = 0
    allowed_dists = [i for i in range(0,m+1)]
    matching_motifs = []

    for peak in peaks:
        if hitcount>=N: break
        chrom = peak[0]
        if chrom not in motifs: continue
        peak_start = int(float(peak[1]))
        peak_end = int(float(peak[2]))
        if args.comparison==0: peak_loc = (peak_end-peak_start)/2+peak_start
        elif args.comparison==1: peak_loc = peak_end
        elif args.comparison==2: peak_loc = peak_start
        
        for d in allowed_dists:
            if (peak_loc+d in motifs[chrom][0]) or (peak_loc-d in motifs[chrom][0]) or (peak_loc+d in motifs[chrom][1]) or (peak_loc-d in motifs[chrom][1]):
                #we have a hit
                hitcount += 1
                topPeaks.append(peak)
                #removing the used motif
                if peak_loc+d in motifs[chrom][0]:
                    dists_from_summit.append(d)
                    motifs[chrom][0].remove(peak_loc+d)
                    strands.append('+')
                    matching_motifs.append(peak_loc+d)
                elif peak_loc+d in motifs[chrom][1]:
                    dists_from_summit.append(d)
                    motifs[chrom][1].remove(peak_loc+d)
                    strands.append('-')
                    matching_motifs.append(peak_loc+d)
                elif peak_loc-d in motifs[chrom][0]:
                    dists_from_summit.append(-d)
                    motifs[chrom][0].remove(peak_loc-d)
                    strands.append('+')
                    matching_motifs.append(peak_loc-d)
                else:
                    dists_from_summit.append(-d)
                    motifs[chrom][1].remove(peak_loc-d)
                    strands.append('-')
                    matching_motifs.append(peak_loc-d)
                break

    #print str(hitcount)+" peaks overlapping a motif"
    #creating matrix of the read densities of top peaks
    E = np.zeros((2*N,w+1))
    ind = 0

    for p in range(0,len(topPeaks)):
        chrom = topPeaks[p][0]
        start = int(float(topPeaks[p][1]))
        end = int(float(topPeaks[p][2]))
        summit = (end-start)/2+start
        loc = summit+dists_from_summit[p]
        for read in samfile.fetch(chrom,loc-w/2-100,loc+2/2+1+100):
            if read.is_unmapped: continue
            strand = '+'
            if read.is_reverse: strand = '-'
            
            #determine five prime end
            if strand == '+': fiveprime = read.aend-read.rlen+1
            else: fiveprime = read.aend

            if fiveprime in range(loc-w/2,loc+w/2+1):
                if strand=='+':
                    if args.truereads==1:
                        #true reads are on the left side of the motif center
                        if fiveprime<loc: E[ind,fiveprime-(loc-w/2)] += 1.0
                    else: E[ind,fiveprime-(loc-w/2)] += 1.0
                else:
                    if args.truereads==1:
                        #true reads are on the right side of the motif center
                        if fiveprime>loc: E[ind+1,fiveprime-(loc-w/2)] -= 1.0
                    else: E[ind+1,fiveprime-(loc-w/2)] -= 1.0
        ind += 2
    #plotting

    #Smoothing using np.convolve
    if S>1:
        for i in range(0,N): E[i,:] = np.convolve(E[i,:],np.ones(S),'same')/float(S)

    #calculating the average peak width
    widths = []
    for i in range(0,2*N,2):
        d = 0.0
        left = E[i,0:w/2][::-1]
        right = E[i+1,w/2+1:]
        if len(np.where(left>0)[0])>0: d += np.where(left>0)[0][0]
        else: continue
        if len(np.where(right<0)[0])>0: d += np.where(right<0)[0][0]
        else: continue
        widths.append(d)

    x = np.array([i for i in range(-w/2,w/2+1)])
    y = np.array([i for i in range(0,np.shape(E)[0])])
    plotHeatMap(E,x,y,args.outdir+"Heatmap_motif_match_top_"+str(N)+".png",xlabel="Distance from motif center (bp's)",ylabel="Top "+str(N)+" peaks with motif match",xticks=[-w/2,w/2],yticks=[])
    #print "done!"

    #############################################
    #CLUSTERING THE TOP PEAKS WITH A MOTIF MATCH#
    #############################################
    if args.skipheatmap==2:
        E_clustered,centroids,clusters = cr.clusterRegions(E[range(0,len(E),2),:],np.abs(E[range(1,len(E),2),:]),args.k,args.npass,args.ctest,args.nproc)
        #plotting cluster centroids
        x = np.array([i for i in range(-w/2,w/2+1)])
        y = np.array([i for i in range(0,int(np.shape(E)[0]))])
        ind = 0
        for c in centroids:
            ind += 1
            plot1d(x,[c[1],c[2]],args.outdir+"Centroid"+str(ind)+"_top_"+str(N)+"motif_match.png",xlabel="Distance from peak summit",ylabel="Average ChIP-Nexus reads",colors=['r','b'],title="Centroid="+str(ind)+", "+str(len(np.where(clusters==c[0])[0]))+" peaks")
        #plotting heatmap of each cluster
        ind = 0
        for c in set(clusters):
            ind += 1
            E_aux = E_clustered[np.where(E_clustered[:,-1]==c)]
            plotHeatMap(E_aux[:,:-1],x,np.array([i for i in range(0,int(np.shape(E_aux[:,:-1])[0]))]),args.outdir+"Heatmap_top_"+str(N)+"cluster"+str(ind)+"motif_match.png",xlabel="Distance from peak summit (bp's)",ylabel="Cluster"+str(ind),xticks=[-w/2,w/2])

    #######################################################
    #HISTOGRAM OF DISTANCES BETWEEN MOTIF AND PEAK SUMMITS#
    #######################################################

    bins = [i for i in range(-m,m+1)]
    plotHistogram(np.array(dists_from_summit),bins,args.outdir+"Dist_from_motif_center_top_"+str(N)+".png",color='g',xlabel="Distance from motif center (bp's)",ylabel="Hit count",legend=args.expname)
    

    ################################################################
    #AVERAGE READ DENSITY PROFILE OF TOP N PEAKS WITH A MOTIF MATCH#
    ################################################################

    #print "Plotting read density profile of top "+str(N)+" peaks with a motif match...",
    plot1d(x,[np.mean(np.abs(E[range(0,len(E),2)]),axis=0),np.mean(np.abs(E[range(1,len(E),2)]),axis=0)],args.outdir+"Avg_motif_match_top_"+str(N)+".png",colors=['r','b'],xlabel="Distance from motif center",ylabel="Average 5'-end count")   
    #print "done!"

    ##########################################################################
    #AVERAGE READ DENSITY PROFILE OF TOP N PEAKS WITH AN ORIENTED MOTIF MATCH#
    ##########################################################################

    #for motifs on + strand: left=start, right=end
    #for motifs on - strand: left=end, right=start
    #print "Plotting heatmap of top peaks with an oriented motif match...",

    #creating matrix of the read densities of top peaks
    #list matching_motifs contains the middle positions of motifs and
    #variable motif_size the width of motif

    E = np.zeros((2*N,w+1))
    ind = 0
    dists_from_start = []
    dists_from_end = []
    for p in range(0,len(topPeaks)):
        chrom = topPeaks[p][0]
        start = int(float(topPeaks[p][1]))
        end = int(float(topPeaks[p][2]))
        summit = (end-start)/2+start
        motif_loc = matching_motifs[p]

        if strands[p]=='+':
            loc = summit+dists_from_summit[p]
            #start of motif is paired with start of peak
            dists_from_start.append(abs(start-(motif_loc-motif_size/2)))
            dists_from_end.append(abs(end-(motif_loc+motif_size/2)))
        else:
            loc = summit-dists_from_summit[p]
            #start of motif is paired with end of peak
            dists_from_start.append(abs(end-(motif_loc+motif_size/2)))
            dists_from_end.append(abs(start-(motif_loc-motif_size/2)))

        for read in samfile.fetch(chrom,motif_loc-w/2-100,motif_loc+2/2+1+100):
            if read.is_unmapped: continue
            strand = '+'
            if read.is_reverse: strand = '-'
            
            #determine five prime end
            if strand == '+': fiveprime = read.aend-read.rlen+1
            else: fiveprime = read.aend

            if fiveprime in range(motif_loc-w/2,motif_loc+w/2+1):
                if strand=='+':
                    if args.truereads==1:
                        #true reads are on the left side of the motif center
                        if fiveprime<motif_loc: E[ind,fiveprime-(motif_loc-w/2)] += 1.0
                    else: E[ind,fiveprime-(motif_loc-w/2)] += 1.0
                else:
                    if args.truereads==1:
                        #true reads are on the right side of the motif center
                        if fiveprime>motif_loc: E[ind+1,fiveprime-(motif_loc-w/2)] -= 1.0
                    else: E[ind+1,fiveprime-(motif_loc-w/2)] -= 1.0

        #now orienting the read density with respect to motif direction
        if strands[p]=='-':
            aux = -1*E[ind+1,:][::-1]
            E[ind+1,:] = -1*E[ind,:][::-1]
            E[ind,:] = aux

        ind += 2
            
    #plotting

    #Smoothing using np.convolve
    if S>1:
        for i in range(0,N): E[i,:] = np.convolve(E[i,:],np.ones(S),'same')/float(S)

    #calculating average peak size and its standard deviation
    std_peak_size = np.std(np.array(widths))
    avg_peak_size = np.mean(np.array(widths))
    med_peak_size = np.median(np.array(widths))

    #zooming E-matrix to the 30 bp's around summit
    x = np.array([i for i in range(-100,101)])
    y = np.array([i for i in range(0,np.shape(E)[0])])
    
    newx = np.array([i for i in range(-15,16)])
    newy = np.array([i for i in range(0,np.shape(E)[0])])

    x = np.array([i for i in range(-w/2,w/2+1)])
    y = np.array([i for i in range(0,np.shape(E)[0])])
    
    plotHeatMap(E,x,y,args.outdir+"Heatmap_motif_match_oriented_top"+str(N)+".png",xlabel="Distance from motif center (bp's)",ylabel="Top "+str(N)+" peaks with motif match",xticks=[-w/2,w/2+1],yticks=[])
    #print "done!"  

    #######################################################
    #CLUSTERING THE TOP PEAKS WITH AN ORIENTED MOTIF MATCH#
    #######################################################
    if args.skipheatmap==0:
        E_clustered,centroids,clusters = cr.clusterRegions(E[range(0,len(E),2),:],np.abs(E[range(1,len(E),2),:]),args.k,args.npass,args.pseudo,args.nproc,args.ctest)
        #plotting cluster centroids
        x = np.array([i for i in range(-w/2,w/2+1)])
        y = np.array([i for i in range(0,int(np.shape(E)[0]))])
        ind = 0
        for c in centroids:
            ind += 1
            plot1d(x,[c[1],c[2]],args.outdir+"Centroid"+str(ind)+"_top_"+str(N)+"_oriented_motif_match.png",xlabel="Distance from peak summit",ylabel="Average ChIP-Nexus reads",colors=['r','b'],title="Centroid="+str(ind)+", "+str(len(np.where(clusters==c[0])[0]))+" peaks")
        #plotting heatmap of each cluster
        ind = 0
        for c in set(clusters):
            ind += 1
            E_aux = E_clustered[np.where(E_clustered[:,-1]==c)]
            plotHeatMap(E_aux[:,:-1],x,np.array([i for i in range(0,int(np.shape(E_aux[:,:-1])[0]))]),args.outdir+"Heatmap_top_"+str(N)+"cluster"+str(ind)+"_oriented_motif_match.png",xlabel="Distance from peak summit (bp's)",ylabel="Cluster"+str(ind),xticks=[-w/2,w/2])  

    ##########################################################################
    #AVERAGE READ DENSITY PROFILE OF TOP N PEAKS WITH AN ORIENTED MOTIF MATCH#
    ##########################################################################

    auxtitle = "Mean peak size="+str(avg_peak_size)+", std="+str(std_peak_size)+", median="+str(med_peak_size)
    #print "Plotting read density profile of top "+str(N)+" peaks with oriented motif match...",
    plot1d(x,[np.mean(np.abs(E[range(0,len(E),2)]),axis=0),np.mean(np.abs(E[range(1,len(E),2)]),axis=0)],args.outdir+"Avg_motif_match_top_oriented_"+str(N)+".png",colors=['r','b'],xlabel="Distance from motif center",ylabel="Average 5'-end count",m_left=args.motif_start,m_len=args.motif_len,opacity=0.3,errorbar=[np.std(E[range(0,len(E),2)],axis=0),np.std(E[range(0,len(E),1)],axis=0)])   
    #print "done!"    

    ########################################################
    #HISTOGRAMS OF DISTANCES FROM MOTIF START AND MOTIF END#
    ########################################################

    plotHistogram(np.array(dists_from_start),bins,args.outdir+"Dist_from_motif_start_top_"+str(N)+".png",color='r',xlabel="Distance from motif start (bp's)",ylabel="Hit count",legend=args.expname)

    plotHistogram(np.array(dists_from_end),bins,args.outdir+"Dist_from_motif_end_top_"+str(N)+".png",color='b',xlabel="Distance from motif end (bp's)",ylabel="Hit count",legend=args.expname)

    return

    ################################################
    #RUNNING MEME TO LOOK FOR MOTIFS FROM TOP PEAKS#
    ################################################

    #parameters for meme
    #meme -revcomp -dna -minw 5 -maxw 25 -evt 0.05 -oc meme_out_1000_p-value/ first_1000_p-value.fasta
    minw = 3
    maxw = 30
    evt = 0.05
    oc = args.outdir+"/meme_out_"+str(N)+"/"
    p = args.nproc
    mod = "anr"
    nmotifs = 5

    #print "meme -revcomp -dna -minw "+str(minw)+" -maxw "+str(maxw)+" -evt "+str(evt)+" -p "+str(p)+" -mod "+mod+" -nmotifs "+str(nmotifs)+" -oc "+oc+" "+args.peakseqs
    system("meme -revcomp -dna -minw "+str(minw)+" -maxw "+str(maxw)+" -evt "+str(evt)+" -mod "+mod+" -nmotifs "+str(nmotifs)+" -oc "+oc+" "+args.peakseqs)
    
    
#end

def plotHeatMap(E,x,y,outname,color='seismic',title=None,xlabel=None,ylabel=None,xticks=None,xticklabels=None,yticks=None,yticklabels=None,figsize=None,resolution=900):
    #plots a heat map using pcolormesh from matplotlib

    fig = plt.figure()
    gs = gridspec.GridSpec(100,100,top=0.96,bottom=0.04,left=0.04,right=0.96)
    ax = fig.add_subplot(gs[:,0:96])
    axC = fig.add_subplot(gs[:,97:99])
    quadmesh = ax.pcolormesh(x,y,E,cmap=color,vmin=-1*np.abs(E).max(),vmax=np.abs(E).max())
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(quadmesh,ax=ax,cax=axC)
    if title!=None: ax.set_title(title)
    if ylabel!=None: ax.set_ylabel(ylabel)
    if xlabel!=None: ax.set_xlabel(xlabel)
    if xticks!=None: ax.set_xticks(xticks)
    if yticks!=None: ax.set_yticks(yticks)
    if figsize!=None: fig.set_size_inches((8,6))

    plt.savefig(outname,dpi=resolution)
    plt.clf()
    plt.close()

def plot1d(x,y,outname,legend=None,colors=['r','b','g','c','m','k'],title=None,xlabel=None,xticks=None,xticklabels=None,ylabel=None,yticks=None,yticklabels=None,figsize=None,resolution=None,opacity=1,m_left=None,m_len=None,errorbar=None):
    #plot normal 1d-plot using plot from matplotlib

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    for i in range(0,len(y)):
        if legend==None: ax.plot(x,y[i],color=colors[i],linewidth=4)
        else: ax.plot(x,y[i],color=colors[i],linewidth=4,label=legend[i])
    ax.axis([min(x), max(x), 0, max([max(y[i]) for i in range(0,len(y))])])
    maxy = max([max(y[i]) for i in range(0,len(y))])
    #plt.xlim([min(x),max(x)])
    #error bars
    if errorbar!=None:
        for i in range(0,len(errorbar)): ax.fill_between(x,y[i]-errorbar[i]/4.0,y[i]+errorbar[i]/4.0,alpha=0.2,edgecolor=colors[i],facecolor=colors[i])

    if opacity<1:
        x = []
        i = 0
        while i<(m_len-1):
            x.append(m_left+i)
            i += 1
        y = [maxy for i in x]
        ax.bar(x,y,width=1,color='g',alpha=opacity,edgecolor=None,linewidth=0)

    if title!=None: ax.set_title(title)
    if ylabel!=None: ax.set_ylabel(ylabel,fontsize=26,weight='bold')
    if xlabel!=None: ax.set_xlabel(xlabel,fontsize=26,weight='bold')
    if xticks!=None: ax.set_xticks(xticks)
    if yticks!=None: ax.set_yticks(yticks)
    if figsize!=None: fig.set_size_inches((8,6))
    if legend!=None: ax.legend(loc='upper right')

    plt.locator_params(axis='y',nbins=3,tight=True,trim=False)
    plt.tight_layout()
    plt.savefig(outname,dpi=resolution)
    plt.clf()
    plt.close()

def plotHistogram(data,bins,outname,color='r',legend=None,xlabel=None,ylabel=None,title=None,figsize=None,resolution=None):

    hist,edges = np.histogram(data,bins=bins,density=False)
    plt.bar(edges[:-1],hist,color=color,width=edges[1]-edges[0],label=legend+"\nmean="+str(np.mean(np.array(data)))+"\nmedian="+str(np.median(np.array(data))))
    if xlabel!=None: plt.xlabel(xlabel)
    if ylabel!=None: plt.ylabel(ylabel)
    plt.margins(0.05)
    plt.subplots_adjust(bottom=0.08,top=0.87)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5,1.17),ncol=1, fancybox=True, shadow=True)
    
    plt.savefig(outname)
    plt.clf()
    plt.close()

fancyResults()
