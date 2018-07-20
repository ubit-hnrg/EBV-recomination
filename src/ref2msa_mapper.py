import argparse   
import pandas as pd
import os
from Bio import AlignIO 
from Bio import SeqIO
import numpy as np


parser = argparse.ArgumentParser(description='remapping of msa coordinates given a refeference sequence (also containded in the msa)')
if(True):
    parser.add_argument('-m','--msapath',required=True,help='msa file in fasta format')
    parser.add_argument('-r','--refseq',required=True,help='reference sequence file in fasta format')
    parser.add_argument('-o','--outputfile',required=True,help='output file *.csv')
    parser.add_argument('-R','--refgenome',dest = 'refgenome',default ='NC_007605.1',help='Character, indicating genome version of reference to use. Default: "NC_007605"')

    args = parser.parse_args()
    msapath = args.msapath
    refseq = args.refseq
    outputfile = args.outputfile
    refgenome = args.refgenome
else:
    msapath = '/data/EBV/msas/N.America/N.America_msa_gap100.fa'
    refseq = '/data/EBV/byACCIDs/NC_007605.fasta'
    refgenome = 'NC_007605.1'
    outputfile = '/data/EBV/msas/N.America/mapp_coordinates_N.America_msa_gap100.csv'

def ref2msa_mapper(ref,refInMsa):
    mismach =[]
    cumm_gaps = 0
    mapper = {}
    gapscounter = 0
    for i in range(len(ref)):
        if(ref[i] != refInMsa[i+cumm_gaps]):
            gapscounter = gapscounter+1
            mismach.append(i)
            targ = ref[i]
            while refInMsa[i+cumm_gaps]=='-':
                cumm_gaps=cumm_gaps+1

        mapper.update({i:i+cumm_gaps})
    return({'mapp':mapper,'gapscounter':gapscounter,'cumm_maps':cumm_gaps})

#######################################
## extract reference sequence from mmsa
msa = AlignIO.read(msapath,'fasta')
ids_msa = [msa[i,0:5].id for i in range(len(msa))]
ii = np.array(ids_msa) == refgenome
j = np.where(ii)[0][0]
refmsa = msa[j,:]
refmsa = refmsa.upper()
print 'selected sequence of msa %s'%refmsa.name

#######################################
## extract reference sequence from fasta file

ref=SeqIO.read(refseq,'fasta')
try:
    refgenome = ref.id
    print 'used reference genome %s '%refgenome
except:
    print 'your reference sequence do not have any Id. Please try passing with --refgenome option'
    exit()
ref=ref.upper()


#################################################
########   mapp coordinates ########
#################################################
########   mapp coordinates ########
mapa = ref2msa_mapper(ref,refmsa)
mapdict = mapa['mapp']
mapp = pd.Series(mapdict).to_frame().reset_index()
mapp.columns =['refCoord','msaCoord']
mapp = mapp+1
print 'number of gaps: %s'%mapa['gapscounter']
meanGapSize = np.round(mapa['cumm_maps']/float(mapa['gapscounter']),3)
print 'mean size of gaps: %s'%meanGapSize



#### write output ####
mapp.to_csv(outputfile,sep = '|',index=False)


