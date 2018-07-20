import os 
import pandas as pd
import glob
from Bio import AlignIO 
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Given a msa in fasta file get a modified fasta file format acording to LDHat specifications')

parser.add_argument('-f','--msafile',required=True)

args = parser.parse_args()
msafile = args.msafile


output_modified_fasta = msafile.split('.fa')[0]+'.modFasta'

msa = AlignIO.read(msafile,'fasta')
L = msa.get_alignment_length()
nseq = len(msa[:,0])

with open(output_modified_fasta, "w") as handle:
    handle.write('%s %s %s \n'%(nseq,L,2))
    count = SeqIO.write(msa, handle, "fasta")
handle.close()





