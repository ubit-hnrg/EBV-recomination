import os 
import pandas as pd
import glob
from Bio import SeqIO

datapath = '/data/EBV/byACCIDs/'
files = glob.glob(datapath + "*.gb") 


dfl =[]
for k,f in enumerate(files):
#    if (k & (k-1) == 0):
#        print k
    record = SeqIO.read(f, "gb")
    auxdf = pd.Series({'AccId':record.annotations['accessions'][0],'seqlen':len(record.seq)})
    dfl.append(auxdf)

seqLengths = pd.concat(dfl,axis = 1).transpose()
seqLengths.set_index('AccId',inplace = True)

seqLengths.to_csv(datapath+'seq_lengths.txt',sep=',',index = True)



