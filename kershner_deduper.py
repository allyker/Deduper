#!/usr/bin/env python3
#This program takes a SAM file that is sorted (at least chromosomes must be grouped together)
#And outputs a SAM file with PCR duplicated removed
#Doesn't do paired end files right now

import stringdist
import re
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-u", "--umi", type=str, help='file with list of UMIs')
parser.add_argument("-e", "--error", action="store_true", help='error corrects UMIs by 1 character, StringDist package required')
parser.add_argument("-f","--file", type=str, required=True, help='SORTED sam file')
parser.add_argument("-p", "--paired", action="store_true", help='put this flag if your file is paired end')
args = parser.parse_args()

def FLAG_check(FLAG):
    """looks at FLAG to see if forward or reverse/checks for paired/unmapped"""
    if((4 & FLAG) == 4):
        return "UNMAPPED"
    if((1 & FLAG) == 1):
        print("NO PAIRED END FUNCTIONALITY AT THIS TIME")
        quit()
    if((16 & FLAG) != 16):
        return "F"
    else:
        return "R"
#FLAG_check(16) should return R

def UMI_check(UMI):
    '''sees if umi is in reference list'''
    if UMI in UMI_dict:
        return(UMI)
    elif args.error ==True:
        return(UMI_correct(UMI))
    else:
        return("BAD")

def UMI_correct(UMI):
    """Corrects UMI by up to one letter"""
    for item in UMI_dict:
        if (stringdist.levenshtein(item, UMI)) <= 1:
            return(item)
        else:
            return("BAD")

def UMI_check2(UMI):
    '''N check for randomers'''
    if "N" in UMI:
        return("BAD")
    else:
        return(UMI)

def CIGAR_check(CIGAR, POS, dire):
    '''looks at CIGAR string to determine actual 5' POS'''
    import re
    adjf=re.findall(r'^\d+S', CIGAR)
    adjr=re.findall(r'\d+', CIGAR)
    if adjf==[]:
        if dire=="R":
            rev=0
            for item in adjr:
                rev+=int(item)
            POS+=rev
    else:
        FSadj=int(adjf[0].rstrip("S"))
        if dire=="F":
            POS-=FSadj
        if dire=="R":
            rev=0
            for item in adjr:
                rev+=int(item)
            rev-=FSadj
            POS+=rev
    return(POS)
#EX: CIGAR(3S10M, 100, forward) would return 97

#Check if paired option used, don't have paired functionality
if args.paired == True:
    print("NO PAIRED END FUNCTIONALITY AT THIS TIME")
    quit()

#Make dictionary of UMIs if given file with UMIs
if args.umi is not None:
    UMI_dict={}
    with open(args.umi, "rt") as fh1:
        for line in fh1:
            line=line.strip()
            UMI_dict[line]="v"
#Set variables that will be used
Store_dict={}
CHR=0
infile=args.file.rstrip(".sam")
#Open input and output sam files
with open(args.file, "rt") as fh, open(infile+"_deduped.sam", "wt") as outfile:
    while True:
#Read in file one line at a time, split into list
        line=fh.readline().strip().split()
#Empty any remaiing lines in storage dict at end of file
        if not line:
            for key in Store_dict:
                print(*Store_dict[key], file=outfile)
            break
#Print header lines
        if "@" in line[0]:
            print(*line, file=outfile)
#For all lines in each chromosome, check for duplicates, store unique in storage dict
        else:
            umilist=line[0].split(":")
#Print out lines and empty storage dict after each chromosome
            if CHR != int(line[2]):
                for key in Store_dict:
                    print(*Store_dict[key], file=outfile)
                CHR=int(line[2])
                Store_dict={}
#Compare POS, FLAG, and UMI if given UMI file
            if args.umi is not None:
                UMI=UMI_check(umilist[7])
                if UMI != "BAD":
                    FLAG=FLAG_check(int(line[1]))
                    if FLAG != "UNMAPPED":
                        POS=CIGAR_check(line[5], line[3], FLAG)
                        Store_dict.setdefault((POS,FLAG,UMI), line)
#Compare POS, FLAG, and UMI if not given UMI file
            if args.umi is None:
                UMI=UMI_check2(umilist[7])
                if UMI != "BAD":
                    FLAG=FLAG_check(int(line[1]))
                    if FLAG != "UNMAPPED":
                        POS=CIGAR_check(line[5], line[3], FLAG)
                        Store_dict.setdefault((POS,FLAG,UMI), line)
