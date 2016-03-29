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
from os import system
from sys import exit
import csv
import textwrap

import multiprocessing as mp

def PeakXus():

    ########################
    #COMMAND LINE ARGUMENTS#
    ########################

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,description=textwrap.dedent('''\

    ***  *****   *   *  * *   * *   *  ****
    *  * *      * *  * *   * *  *   * * 
    ***  *****  ***  **     *   *   *  ****
    *    *     *   * * *   * *  *   *      *
    *    ***** *   * *  * *   *  ***  ****

    Copyright (c) Tuomo Hartonen, 2015-2016

    THIS PROGRAM COMES WITH ABSOLUTELY NO WARRANTY!
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License version 2 as 
    published by the Free Software Foundation.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see http://www.gnu.org/licenses/.
    
    ----------
    |OVERVIEW|
    ----------

    This script performs the following functionality if a fastq-file is given as an input:
    
      1) Trimming adapters from ends of reads using Cutadapt-tool.
      2) Aligning the reads in the input fastq-file to a specified reference genome using bwa aln-algorithm.
      3) Peak calling using PeakXus.
      4-6) Plotting different kind of figures for graphical inspection of the peak calling results.
    
    Steps 1-2 are skipped if the input is already in bam-format. Steps 4-6 are optional.
    Below is a detailed description of all the possible input parameters for different stages of the pipeline.

    ----------------
    |USAGE EXAMPLES|
    ----------------

    Majority of the optional input parameter values seldom need tweaking. Below are listed some example
    calls of the pipeline for conducting different tasks. Input file locations are assumed to be:

    ref_genome_path/genome.idx = reference genome index
    ref_genome_path/wg.fasta   = reference genome in fasta-format
    exp_path/exp.fastq         = reads from ChIP-Nexus experiment in fastq-format
    exp_path/exp.bam           = aligned reads from ChIP-Nexus experiment in bam-format
    exp_path/UMIs.txt          = UMI-labels used in the experiment
    exp_path/chroms.txt        = chromosome names in the ChIP-Nexus experiment
    motif_path/matrix_hits.gff = hits of the corresponding TF-pwm to the used reference genome in gff-format

      i) PEAK CALLING ONLY FROM A BAM-FILE
      --------------------------------------

        without UMIs:
        PeakXus.py exp_path/exp.bam out/ exp_path/chroms.txt

        with UMIs:
        PeakXus.py exp_path/exp.bam out/ exp_path/chroms.txt --UMIs exp_path/UMIs.txt

      ii) PEAK CALLING ONLY FROM A FASTQ-FILE
      ---------------------------------------

        aligning with 8 cores
        PeakXus.py exp_path/exp.fastq out/ exp_path/chroms.txt --genome2 ref_genome_path/genome.idx --t2 8

        UMIs of length 4 and 5 used
        PeakXus.py exp_path/exp.fastq out/ exp_path/chroms.txt --genome2 ref_genome_path/genome.idx --UMIs exp_path/UMIs.txt --B2 4 5

      iii) PEAK CALLING FROM A BAM-FILE, PLOTTING GRAPHICAL RESULTS
      -------------------------------------------------------------

        assuming the binding motif width is 8 bp's
        PeakXus.py exp_path/exp.bam out/ exp_path/chroms.txt --matrixhits motif_path/matrix_hits.gff --ml4 8 --ms4 -4

    ----------------------------------
    |DETAILED LIST OF INPUT ARGUMENTS|
    ----------------------------------
    '''))
    
    #INPUT FILES, MANDATORY
    parser.add_argument("fastq",help="Full path to the input fastq-file. If input is a bam-file, alignment-step is skipped and the pipeline starts straight from peak calling. File type is determined based on the file name extension. If input is in bam-format, also multiple files are accepted.",type=str,nargs='+')
    parser.add_argument("outdir",help="Full path to output directory.",type=str,nargs=1)
    parser.add_argument("chromnames",help="Full path to a file containing chromosome names each on its on line.",type=str,nargs=1)
    
    #INPUT FILES, OPTIONAL
    parser.add_argument("--UMIs",help="Full path to file(s) containing all UMIs. If UMIs of different length were used, they should be provided in different input files. Files should have two tab-separated columns, first is the UMI name and second is the sequence. If left blank, UMIs are not used in the analysis.",type=str,nargs='+',default=None)
    parser.add_argument("--matrixhits",help="Full path to .gff-file containing the hits to binding matrix. If left unused, part of the results are not calculated.",type=str,nargs=1,default=None)
    
    #1 CUTADAPT
    group1 = parser.add_argument_group('1. CUTADAPT','These arguments are passed to Cutadapt-tool for trimming adapters ends of reads in the input fastq-file.')
    group1.add_argument("--adapters1",help="Adapters trimmed from the ends of reads using cutadatp.",type=str,nargs='+',default=None)

    #2 ALIGNER
    group2 = parser.add_argument_group('2. BWA','These arguments are passed to bwa aln-algorithm for aligning the reads in fastq-file.')
    group2.add_argument("--genome2",help="Full path to the reference genome index.",type=str,default=None)
    group2.add_argument("--B2",help="Length(s) of UMI-label(s). Should be in the same order as the UMI-files.",type=int,default=None,nargs='*')
    group2.add_argument("--q2",help="Quality threshold for base calls (default=20).",type=int,default=20)
    group2.add_argument("--t2",help="Number of parallel processes used by the aligner (default=1)",type=int,default=1)

    #3 PEAK CALLING
    group3 = parser.add_argument_group('3. PEAKXUS','These arguments are passed to PeakXus-peak caller.')
    #group3.add_argument("--m3",help="minimum distance between two peaks (default=30)",type=int,default=30)
    group3.add_argument("--w3", help="width of the window around 5'-end (default=5)", type=int, default=5)
    group3.add_argument("--l_l3", help="l_limit+w is the largest allowed width for a peak (default=60)", type=int, default=60)
    group3.add_argument("--u_l3", help="u_limit is the smallest allowed peak width (default=5)", type=int,default=5)
    group3.add_argument("--b3",help="Threshold p-value for background model. (default=0.05)", type=float, default=0.05)
    group3.add_argument("--s3",help="Width of a sliding window used in averaging density of 5' ends of reads. This should be an odd number (default=1).",type=int,default=1)
    group3.add_argument("--n3",help="Number of parallel processes used in calculating significance of peaks (default=1).",type=int,default=1)
    group3.add_argument("--p3",help="Pseudocount added to each position when testing significance of candidate binding sites (default=10.0).",type=float,default=10.0)

    #4 READ DENSITIES AROUND PEAKS
    group4 = parser.add_argument_group('4. READ DENSITY FIGURES','These arguments are passed to a script that plots figures detailing the read density profiles around called peaks. These results are not calculated if --matrixhits argument is not given.')
    group4.add_argument("--wg4",help="Full path to reference genome file in fasta-format",type=str,default=None)
    group4.add_argument("--N4",help="Number of top peaks considered (default=1000)",type=int,default=1000)
    group4.add_argument("--m4",help="maximum distance between peak and motif (default=10).",type=int,default=10)
    group4.add_argument("--w4",help="width of region around the peak in heatmaps (default=100bp's)",type=int,default=100)
    group4.add_argument("--e4",help="experiment name for plot titles.",type=str,default="ChIP-exo")
    group4.add_argument("--s4",help="Sorting order of peaks, 0=peak score(default), 1=test statistics, 2=read count.",type=int,choices=[0,1,2],default=0)
    group4.add_argument("--n4",help="Number of parallel processes used in clustering (default=1).",type=int,default=1)
    group4.add_argument("--k4",help="Number of clusters (default=4).",type=int,default=4)
    group4.add_argument("--4p",help="number of times clustering is run (default=1).",type=int,default=1)
    group4.add_argument("--S4",help="Window size used in calculating read distributions (default=1).",type=int,default=1)
    group4.add_argument("--ms4",help="Start position of the motif, if middle of the motif is at 0 (default=-8, CTCF)",type=int,default=-8)
    group4.add_argument("--ml4",help="Motif length (default=17, CTCF).",type=int,default=17)
    
    #5 BINDING MOTIF HITS PER PEAK
    group5 = parser.add_argument_group('5. FIGURES OF MOTIFS OVERLAPPING A PEAK','These arguments are passed to a script that calculates the number of peaks that overlap with a TF-specific high-affinity binding sequence. These results are not calculated if --matrixhits argument is not given.')
    group5.add_argument("--N5",help="Number of top peaks are analyzed. (default=1000).",type=int,default=1000)
    group5.add_argument("--m5",help="Maximum distance from peak summit to motif border. (default=20).",type=int,default=20)
    group5.add_argument("--s5",help="Sorting order,1=signal, 3=random, 4=score (default=4),0=already sorted, no sorting is done.",type=int,default=4,choices=[0,1,2,3,4])
    group5.add_argument("--o5",help="If 1, each motif is only counted matching to one peak, if 0, there is no limit (default=1).",type=int,default=1)
    group5.add_argument('--title5',type=str,default="PeakXus")

    #6 PEAK DISTANCE FROM MOTIF
    group6 = parser.add_argument_group('6. FIGURES OF PEAK LOCALIZATION ACCURACY','These arguments are passed to a script that calculates the distances of called peaks from TF-specific high-affinity binding sequences. These results are not calculated if --matrixhits argument is not given.')
    group6.add_argument("--l6","--legend",type=str,default="PeakXus")
    group6.add_argument("--m6",help="maximum distance between peak and motif (default=100)",type=int,default=100)
    group6.add_argument("--N6",help="Number of top peaks considered (default is all peaks, assumes input file is sorted in desired order.)",type=int,default=None)
    group6.add_argument("--center6",help="If 1, only center to center distances calculated, otherwise all (default=0)",type=int,default=0)

    parser.add_argument("-v","--verbosity",help="1 (default) if progress is printed to screen, 0 otherwise.",type=int,default=1)
    a = parser.parse_args()

    if a.verbosity==1: print "Testing input parameters...",
    if not testinput(a): exit(0)
    if a.verbosity==1: print "succesful!"

    if a.fastq[-1][-6:]==".fastq":
        #1) TRIMMING ADAPTERS FROM THE ENDS OF READS

        aln_inname = a.fastq[0]
        if a.adapters1!=None:
            if a.verbosity==1: print "Trimming adapters...",
            ca_call = "~/.local/bin/cutadapt -m 22 -O 4 -e 0.2 "
            for adapter in a.adapters1: ca_call += "-a "+adapter+" "
            ca_call += a.fastq[0]
            ca_call += " > "
            ca_call += a.fastq[0][:-6]+"_trimmed.fastq"
            system(ca_call)
            aln_inname = a.fastq[0][:-6]+"_trimmed.fastq"
            if a.verbosity==1: print "succesful!"

        #2) ALIGNING
        if a.verbosity==1: print "Aligning...",
        
        full_bam = a.outdir[0]+aln_inname.split('/')[-1][:-6]+"_sorted_filtered.bam"

        if a.B2==None:
            aln_call = "bwa_align.py "+aln_inname+" "+a.genome2+" "+a.outdir[0]+" -t "+str(a.t2)+" -q "+str(a.q2)
            system(aln_call)

        elif len(a.B2)==1:
            #print aln_inname
            #print a.genome2
            #print a.outdir[0]
            #print a.t2
            #print a.q2
            #print a.B2
            aln_call = "bwa_align.py "+aln_inname+" "+a.genome2+" "+a.outdir[0]+" -t "+str(a.t2)+" -q "+str(a.q2)+" -B "+str(a.B2[0])
            system(aln_call)

        else:
            for B in a.B2:
                aln_call = "bwa_align.py "+aln_inname+" "+a.genome2+" "+a.outdir[0]+"BC"+str(B)+"_ -t "+str(a.t2)+" -q "+str(a.q2)+" -B "+str(B)
                system(aln_call)
            system("samtools merge "+full_bam+" "+a.outdir[0]+"BC*_"+aln_inname.split('/')[-1][:-6]+"_sorted_filtered.bam")
            system("samtools index "+full_bam)
        if a.verbosity==1: print "succesful!"

        #separating chromosomes to distinct files
        chroms = []
        with open(a.chromnames[0],'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r: chroms.append(row[0])

        #processing the files in parallel if possible
        if a.t2>1:
            pool = mp.Pool(processes=a.t2)
            res = [pool.apply_async(samtools,args=("samtools view -b "+a.outdir[0]+aln_inname.split('/')[-1][:-6]+"_sorted_filtered.bam"+" "+i+" > "+a.outdir[0]+i+".bam && samtools index "+a.outdir[0]+i+".bam",a.verbosity)) for i in chroms]
            for r in res: r.get()
            pool.close()
            pool.terminate()
            pool.join()
        else:
            for i in chroms:
                if a.verbosity==1: print "samtools view -b "+a.outdir[0]+aln_inname.split('/')[-1][:-6]+"_sorted_filtered.bam"+" "+i+" > "+a.outdir[0]+i+".bam && samtools index "+a.outdir[0]+i+".bam"
                system("samtools view -b "+a.outdir[0]+aln_inname.split('/')[-1][:-6]+"_sorted_filtered.bam"+" "+i+" > "+a.outdir[0]+i+".bam && samtools index "+a.outdir[0]+i+".bam")
            
    #3 PEAK CALLING

    #creating a UMI-input file for peak calling and later applications (only one column)
    if a.UMIs!=None:
        with open(a.UMIs[0],'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            with open(a.outdir[0]+"UMI.bc",'wb') as outfile:
                w = csv.writer(outfile,delimiter="\t")
                for row in r: w.writerow([row[1]])

    if a.verbosity==1: print "Calling peaks..."
    numtests = 0
    if a.fastq[0][-3:]=="bam":
        for f in a.fastq:
            pc_call = "peakC6.py "+f+" "+a.outdir[0]
            pc_call += " -w "+str(a.w3)+" -l_l "+str(a.l_l3)+" -u_l "+str(a.u_l3)+" -b "+str(a.b3)+" -s "+str(a.s3)+" -n "+str(a.n3)+" -p "+str(a.p3)
            if a.UMIs!=None: pc_call += " -u "+a.outdir[0]+"UMI.bc"

            if a.verbosity==1: print pc_call
            system(pc_call)
            with open(a.outdir[0]+"numtests.txt",'r') as csvfile:
                r = csv.reader(csvfile,delimiter="\t")
                for row in r:
                    numtests += int(float(row[0]))
                    break
    else:
        pc_call ="peakC6.py "+a.outdir[0]+"c*.bam "+a.outdir[0]

        pc_call += " -w "+str(a.w3)+" -l_l "+str(a.l_l3)+" -u_l "+str(a.u_l3)+" -b "+str(a.b3)+" -s "+str(a.s3)+" -n "+str(a.n3)+" -p "+str(a.p3)
        if a.UMIs!=None: pc_call += " -u "+a.outdir[0]+"UMI.bc"

        if a.verbosity==1: print pc_call
        system(pc_call)
        with open(a.outdir[0]+"numtests.txt",'r') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r:
                numtests += int(float(row[0]))
                break

    #Merging the results from different chromosomes
    system("mergeIGV.py "+a.outdir[0]+"*.igv "+a.outdir[0]+"all_transition_points_nofdr.igv")

    #Calculating the false discovery rates with Benjamini-Hochberg
    #with open(a.outdir[0]+"numtests.txt",'r') as csvfile:
    #    r = csv.reader(csvfile,delimiter="\t")
    #    for row in r: numtests = int(float(row[0]))

    system("FDR.py "+a.outdir[0]+"all_transition_points_nofdr.igv "+a.outdir[0]+"all_transition_points.igv -m "+str(numtests))

    if a.verbosity==1: print "succesful!"

    #4 READ DENSITIES AROUND PEAKS

    if a.verbosity==1: print "Plotting results...",

    #fetching the sequences under top peaks
    if a.matrixhits!=None and a.wg4!=None:
        system("seqs_under_peaks.py "+a.wg4+" "+a.outdir[0]+"all_transition_points.igv "+a.outdir[0]+"top"+a.N4+"_peaks.fasta -N "+a.N4+" -w 50 -s 4")
        fr_call = "fancyResults_lowMem.py "+a.outdir[0]+"all_transition_points.igv "+full_bam+" "+a.matrixhits[0]+" "+a.outdir[0]+" "+a.wg4+" -N "+a.N4+" -m "+a.m4+" -w "+a.w4+" -e "+a.e4+" -s "+a.s4+" -n "+a.n4+" -k "+a.k4+" -p "+a.p4+" -S "+a.S4+" -m_s "+a.ms4+" -m_l "+a.ml4
        system(fr_call)

    #5 BINDING MOTIF HITS PER PEAK
    if a.matrixhits!=None:
        #creating input file for motif hit analysis
        with open(a.outdir[0]+"motif.input",'wb') as csvfile:
            r = csv.writer(csvfile,delimiter=',')
            r.writerow([a.outdir[0],"PeakXus","aux",a.matrixhits[0]])

        mh_call = "consMotifHitsPerPeak.py "+a.outdir[0]+"motif.input "+a.outdir[0]+"motifHitsPerPeak_N="+str(a.N5)+"_m="+str(a.m5)+"_o="+str(a.o5)+".png -N "+str(a.N5)+" -m "+str(a.m5)+" -o "+str(a.o5)+" -s "+str(a.s5)+" --title "+a.title5
        system(mh_call)

    #6 DISTANCE FROM MOTIF
    if a.matrixhits!=None:
        #sorting the peaks based on peak score
        system("sort -k 7,7 -n -r "+a.outdir[0]+"all_transition_points.igv > "+a.outdir[0]+"all_transition_points_sorted.igv")
        
        if a.N6==None: cd_call = "cumul_peak_dist_from_motif.py "+a.matrixhits[0]+" "+a.outdir[0]+"dist_from_motif_m="+str(a.m6)+"_N=all.png "+a.outdir[0]+"all_transition_points_sorted.igv -p g -l "+a.l6+" -m "+str(a.m6)+" --center "+str(a.center6)
        else: cd_call = "cumul_peak_dist_from_motif.py "+a.matrixhits[0]+" "+a.outdir[0]+"dist_from_motif_m="+str(a.m6)+"_N="+str(a.N6)+".png "+a.outdir[0]+"all_transition_points_sorted.igv -p g -l "+a.l6+" -m "+str(a.m6)+" -N "+str(a.N6)+" --center "+str(a.center6)
        system(cd_call)
    if a.verbosity==1: print "succesful!"

#end

def samtools(call_string,verbosity):
    
    if verbosity==1: print call_string
    system(call_string)
    return 1

def testinput(a):
    #this function test that the most crucial input parameters are sensible
    
    #1 Reference genome
    #print a.fastq
    #print a.fastq[0][-5:]
    if a.genome2==None and a.fastq[0][-5:]=='fastq':
        print "Reference genome index not given! Terminating..."
        return False
         
    if a.wg4==None:
        print "** Reference genome file in fasta-format not given, some of the graphical output is not produced."
        #z = input("1=yes, otherwise terminating:")
        #if z!=1: return False
    else:
        try:
            f = open(a.wg4,'rb')
            f.close()
        except Exception:
            print e
            print "Reference genome file cannot be opened! Terminating..."
            return False

    #2 fastq/bam input
    try:
        for fil in a.fastq:
            f = open(fil,'rb')
            f.close()

    except Exception as e:
        print e
        print "Input file(/s) cannot be opened! Terminating..."
        return False

    #3 UMI-file
    if a.UMIs!=None:
        try:
            aux = []
            with open(a.UMIs[0],'rb') as csvfile:
                r = csv.reader(csvfile,delimiter="\t")
                for row in r: aux = [row[0],row[1]]
        except Exception as e:
            print e
            print "Error with the UMI-file! Terminating..."
            return False
        
    #4 chromosome names-file
    try:
        aux = ""
        with open(a.chromnames[0],'rb') as csvfile:
            r = csv.reader(csvfile,delimiter="\t")
            for row in r: aux = row[0]

    except Exception as e:
        print e
        print "Problem with chromnames-file! Terminating..."
        return False

    #5 matrix hits-file
    if a.matrixhits!=None:
        try:
            f = open(a.matrixhits[0],'rb')
            f.close()
        except Exception as e:
            print e
            print "matrixhits-file could not be opened! Terminating..."
            return False
    else:
        print "** List of binding matrix hits to genome not provided, graphical results are not printed."
        #z = input("1=yes, otherwise terminating:")
        #if z!=1: return False

    return True
        
            

PeakXus()
