import os 
import pandas as pd
import glob
from Bio import SeqIO

path = '/data/EBV/byACCIDs/'+'seq_lengths.txt'
seqLengths = pd.read_csv(path,index_col=['AccId'])


full_seqs = seqLengths[seqLengths.seqlen>150000] 
full_seqs.to_csv('/data/EBV/seqGroups/AllseqsUpto150k.txt',sep ='\t',index = True)

refseq = 'M80517'
short_reads = seqLengths[(seqLengths.seqlen<=150000)|(seqLengths.index.isin([refseq]))] 
short_reads.to_csv('/data/EBV/seqGroups/AllseqsUnderto150k_withRef.txt',sep ='\t',index = True)


