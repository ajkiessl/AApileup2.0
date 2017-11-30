"""
Module for the parts of the AApileup pipeline.
"""

import pysam
import argparse
from Bio import SeqIO
import MySQLdb
from Bio.Seq import Seq
from pysam import AlignedSegment
import sys

class ab1_to_fastq:
    def __init__(self):
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('inputs',nargs='+')
        self.parser.add_argument('-output',required=True)
        self.parser.add_argument('-rc',action='store_true')
        self.parser.add_argument('-clip5prime')
        self.parser.add_argument('-clip3prime')
        self.args = self.parser.parse_args()

    def convert(self):
        with open(self.args.output,'wb') as out:
            for inp in self.args.inputs:
                record = SeqIO.parse(open(inp,'rb'),'abi').next()
                if self.args.rc:
                    record = record.reverse_complement(id=True,name=True,description=True)
                if self.args.clip5prime:
                    start = record.seq.find(self.args.clip5prime)
                    if start > 0:
                        record = record[start+len(self.args.clip5prime):]
                if self.args.clip3prime:
                    stop = record.seq.find(self.args.clip3prime)
                    if stop > 0:
                        record = record[:stop]
                SeqIO.write(record,out,'fastq')


class fastq_filter:
    def __init__(self):
        self.parser = argparse.Arumentparser()
        self.parser.add_argument('-fastqinput', required=True)
        self.parser.add_argument('-lengthfilter', type=int, default=0)
        self.parser.add_argument('-phredmin', type=int, default=10)
        self.parser.add_argument('-allowance', type=int, default=3)
        self.parser.add_argument('-outputtype',type=str, default='fastq-sanger')
        self.parser.add_argument('-outputid', required=True)
        self.args = self.parser.parse_args()

    def fq_filter(self)
        fastqfile = open(self.args.fastqinput,'rU')
        filter1_fastq = []
        for record in SeqIO.parse(fastqfile, "fastq"):
            if int(self.args.lengthfilter) == len(record.seq):
                filter1_fastq.append(record)
            if int(self.args.lengthfilter) == 0:
                filter1_fastq.append(record)
        filter2_fastq=[]
        for record in filter1_fastq:
            bad_quality_count=0
            for phred in record.letter_annotations["phred_quality"]:
                if phred <= int(self.args.phredmin):
                    bad_quality_count+=1
            if bad_quality_count < int(self.args.allowance):
                filter2_fastq.append(record)
        file = open((self.args.outputid),'w')
        for record in filter2_fastq:
            file.write(str(record.format(self.args.outputtype)))
        file.close()
        fastqfile.close()


class bowtie2_align:
    def __init__(self):
        self.parser = argparse.Arumentparser()
        self.parser.add_argument('-fa', required = True, help='path to fasta file')
        self.parser.add_argument('-fq', required = True, help='path to fastq file')
        self.args = self.parser.parse_args()

    def align(self):
        indexBaseName = str(self.args.fa)[:-6]
        baseName = str(self.args.fq)[:-6]
        subprocess.call(['bowtie2-build', self.args.fa, indexBaseName])
        subprocess.call(['bowtie2', '-x', indexBaseName, '-U', self.args.fq, '-S', baseName + '.sam'])


class sam_to_bam:
    def __init__(self):
        self.parser = argparse.Argumentparser()
        self.parser.add_argument('-inputsam', required=True)
        self.args = self.parser.parse_args()

    def convert(self):
        file1 = open(os.path.splitext(self.args.inputsam)[0] + '.bam','w')
        rows = pysam.view("-Sb", self.args.inputsam  )
        for r in rows:
            file1.write(r)
        file1.close()

class bam_sort:
    def __init__(self):
        self.parser = argparse.Arumentparser()
        self.parser.add_argument('-inputbam', required=True)
        self.args = self.parser.parse_args()
    
    def sort(self):
        rows = pysam.sort(self.args.inputbam, os.path.splitext(self.args.inputbam)[0] + '.sorted')


class bam_index:
    def __init__(self):
        self.parser = argparse.Arumentparser()
        self.parser.add_argument('-inputbam', required=True)
        self.args = self.parser.parser_args()

    def index(self):
        index = pysam.index(self.args.inputbam)


class aa_pileup:
    def __init__(self):
        self.parser = argparse.Arumentparser()
        self.parser.add_argument('-BAMinput',required=True)
        self.parser.add_argument('-fastainput',required=True)
        self.parser.add_argument('-sslpath')
        self.parser.add_argument('-server', required=True)
        self.parser.add_arguemnt('-user', reqired=True)
        self.parser.add_argument('-password', required=True)
        self.parser.add_argument('-database', required=True)
        self.args = self.parser.parser_args()

    def pileup(self):
        fastaFile = pysam.FastaFile(self.args.fastainput)
        bamFile = pysam.AlignmentFile(self.args.BAMinput, "rb")
        ssl_settings = {'ca':self.args.sslpath}
        con = MySQLdb.connect(self.args.server, self.args.user, self.args.password, self.args.database, ssl=ssl_settings)
        with con:
            cur = con.cursor()
            cur.execute("USE " + self.args.database)
        def batch_gen(data, batch_size):
            for i in range(0, len(data), batch_size):
                    yield data[i:i+batch_size]
        references = sorted(set(bamFile.getrname(read.tid) for read in samfile.fetch()))
        referencesLeng = sorted(set(len(fastaFile.fetch(reference=str(item)))for item in references))
        for ref, leng in zip(references, referencesLeng):
            print ref, leng
            cur.execute('INSERT INTO templates(protein, length) VALUES(%s, %s)' ,(ref, leng))
        for reference in references:
            returned_position_lines=[]
            length=0
            refcodonpos=0
            counter=0
            for codon in batch_gen(fastaFile.fetch(reference=str(reference)),3):
                length+=3
                markerlist=[]
                referenceid = str(reference)+ ' '
                refnucpos1=0 +(3*refcodonpos)
                refnucpos2=1 +(3*refcodonpos)
                refnucpos3=2 +(3*refcodonpos)
                if 1 <= (refcodonpos+1) <= 9:
                    refcodonposid = str(refcodonpos+1)+ " "
                else:
                    refcodonposid = str(refcodonpos+1)
                refAAid = str(Seq(codon).translate()[0])
                marker_list=[]
                for read in samfile.fetch():
                    read_codon=[]
                    for seq, pos in zip(read.seq,AlignedSegment.get_reference_positions(read)):
                        if pos == refnucpos1:
                            read_codon.append(seq)
                        if pos == refnucpos2:
                            read_codon.append(seq)
                        if pos == refnucpos3:
                            read_codon.append(seq)
                    if any(read_codon) is True:
                        if len(read_codon) == 3:
                            counter+=1
                            if ''.join(read_codon) == codon:
                                marker_list.append('.')
                            else:
                                marker_list.append(str(Seq("".join(read_codon)).translate()[0]))
                print (referenceid, refcodonposid, refAAid, counter, ''.join(str(item)for item in marker_list))
                returned_position_lines.append(''.join(str(item)for item in marker_list))
                cur.execute("INSERT INTO sites(template_id, position, wild_type_AA) VALUES((SELECT id from templates WHERE protein=%s), %s, %s)" ,(reference, refcodonposid, refAAid))
                counter=0
                refcodonpos+=1
            print returned_position_lines
            AAs = ('A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','*')
            for AA in AAs:
                position=0
                for line in returned_position_lines:
                    position+=1
                    count=0
                    for readAA in line:
                        if readAA==AA:
                            count+=1
                    if (count >= 1):
                        print count, AA, position
                        cur.execute("INSERT INTO substitutions(site_id, substitution, count) VALUES((SELECT id from sites WHERE position=%s AND template_id=(SELECT id from templates WHERE protein=%s)), %s, %s)" ,(position, reference, AA, count))
        con.commit()
        fastaFile.close()
        bamFile.close()

