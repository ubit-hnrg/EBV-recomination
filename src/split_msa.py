import argparse   
import pandas as pd
import os
from Bio import AlignIO 
from Bio import SeqIO
import numpy as np

def cut_function(msa,st,end):
    if(end > st):
        cutted = msa[:,int(st-1):int(end-1)] # biopython start from 0 !!
        output = msapath+'_'+str(st)+'_'+str(end)+'.fa'
        AlignIO.write(cutted,output,'fasta')

    else:
        last = msa.get_alignment_length()
        cut_first = msa[:,int(st-1):(last-1)] # biopython start from 0 !!
        cut_second =msa[:,0:int(end-1)]
        cutted = cut_first + cut_second
        output = msapath+'_'+str(st)+'_'+str(end)+'.fa'
        AlignIO.write(cutted,output,'fasta')


def pseudocircular(a,overlap):
    pseudocircular = np.pad(a, pad_width=(0, overlap-1), mode='wrap')
    return(pseudocircular)

def auto_split_points(msa,length,overlap = None):
    if overlap == None:
        overlap = int(round(length/2))

    ncols = msa.get_alignment_length()
    
    circular = np.arange(ncols)+1
    circular = pseudocircular(circular,overlap=length)
    
    ind_start = np.arange(start = 0,stop = len(circular)-1,step = overlap)
    ind_start = [i for i in ind_start if i<ncols]

    startpoints = circular[ind_start]

    ind_end = [j + length -1 for j in ind_start]
    endpoints = circular[ind_end]
    return(startpoints,endpoints)






parser = argparse.ArgumentParser(description='cut a msa file in fasta format')

if(True):
    parser.add_argument('-m','--msapath',required=True,help='msa file in fasta format')
    parser.add_argument('-s','--start',help='end coordinate',type=int,default = None)
    parser.add_argument('-e','--end',help='end coordinate',type=int,default = None)
    parser.add_argument('-l','--length',help='length for automatic cutting',type=int,default=0)
    
    args = parser.parse_args()
    msapath = args.msapath
    refseq = args.start
    outputfile = args.end
    length = args.length

msa = AlignIO.read(msapath,'fasta')
base = os.path.basename(msapath)

st = args.start
end = args.end

#if not os.path.isdir(outputpath):
#    os.system('mkdir %s'%outputpath)
if (st!=None)&(end!=None):
    cut_function(msa,st=st,end=end)
else:
    startpoints, endpoints = auto_split_points(msa,length=length)
    for i in range(len(startpoints)):
        cut_function(msa,st=startpoints[i],end=endpoints[i])
            
    

    