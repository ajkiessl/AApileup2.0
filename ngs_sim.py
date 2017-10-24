"""
This app takes a fastq file from sangers sequencing and partitions it to replicate NGS fragments.
The output of this app is meant to be used with our alignment/pileup pipeline (alignvariants.sh) to test if it works with both sanger and NGS sequences.
The first section of this app filters the sanger sequence by length and phred score.
The only required arguments are path/to/fastqfile.fastq and the output file id.
Optional arguments are the length to filter at, phred score minimum, number of allowed poor quality reads, the amount of partitions in the simulated NGS file, and the type of output file.
"""

from Bio import SeqIO
import argparse
import random

#Argument parser
parser = argparse.ArgumentParser(description="NGS Simulator")
parser.add_argument('-fastqinput', required=True)
parser.add_argument('-lengthfilter', type=int, default=168)
parser.add_argument('-phredmin', type=int, default=10)
parser.add_argument('-allowance', type=int, default=3)
parser.add_argument('-partitions', type=int, default=3)
parser.add_argument('-outputtype',type=str, default='fastq-sanger')
parser.add_argument('-outputid', required=True)
args = parser.parse_args()

fastqfile = open(args.fastqinput,'rU')

#Defining functions

#Number rounding function
def myround(x, base=3):
    return int(base * round(float(x)/base))

#Sequence partitioning function
def seq_part(seq,parts):
    portions=(len(seq))/parts
    randsec = myround(random.randint(int(round(portions*0.8)),int((round(portions/0.8)))))
    primarysec=[]
    parted_seq=[]
    parted_seq.append(seq[0:(randsec+(random.randint(0,2)))])
    primarysec.append(seq[0:randsec])
    for number in range(2,(parts+1)):
        portions2=((len(seq))/parts)*number
        length_cut=0
        randsecvar = myround(random.randint(int(round(portions2*0.8)),int(round(portions2/0.8))))
        for item in primarysec:
            length_cut+=len(item)
        if number == parts:
            parted_seq.append(seq[length_cut+(random.randint(-2,0)):len(seq)])
        else:    
            parted_seq.append(seq[length_cut+(random.randint(-2,0)):(randsecvar+(random.randint(0,2)))])
            primarysec.append(seq[length_cut:randsecvar])
    return parted_seq

#Filter reads

#Filter sanger seqs by length
filter1_fastq = [] 
for record in SeqIO.parse(fastqfile, "fastq"):
    if int(args.lengthfilter) == len(record.seq):
        filter1_fastq.append(record)

#Filter again by phred quality
filter2_fastq=[]
for record in filter1_fastq:    
    bad_quality_count=0
    for phred in record.letter_annotations["phred_quality"]:
        if phred <= int(args.phredmin):
            bad_quality_count+=1
    if bad_quality_count < int(args.allowance):
        filter2_fastq.append(record)

#NGS Simulator

#Cutting up the sequences
pseudoNGS=[]
for record in filter2_fastq:
    for item in seq_part(record, args.partitions):
        pseudoNGS.append(item)
    
#Changing sequence IDs and appending to final list
pseudoNGSfinal=[]
for number in range(1,(args.partitions+1)):
    for record in pseudoNGS[number-1::args.partitions]:
        record.id+= '(%s)' %(number)
        record.name+= '(%s)' %(number)
        pseudoNGSfinal.append(record)

#write to new file
file = open((args.outputid)+".fastq",'w')
    
for record in pseudoNGSfinal:
    file.write(str(record.format(args.outputtype)))

file.close()
fastqfile.close()
