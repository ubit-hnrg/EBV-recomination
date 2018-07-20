import argparse   
import pandas as pd
import os
from Bio import AlignIO 
from Bio import SeqIO
import numpy as np

## runing examploe
## python ~/repo/TBI/EBV-db/src/trimmed_2_ref.py -m ./sin_trimear/all150_sin_n.fas -r /data/EBV/byACCIDs/NC_007605.fasta -M mapp_trimed_columns.txt -o ./whole_mapping.csv 



parser = argparse.ArgumentParser(description='[reference - msa - msa-trimmed ] mapping \n reference must be contained in msa files')
if(True):
    parser.add_argument('-m','--msa',required=True,help='msa file in fasta format')
#    parser.add_argument('-t','--trimmed_msa',required=True,help='msa file in fasta format')
    parser.add_argument('-r','--refseq',required=True,help='reference sequence file in fasta format')
    parser.add_argument('-R','--refgenome',dest = 'refgenome',default ='NC_007605.1',help='Character, indicating genome version of reference to use. Default: "NC_007605"')

    # Add trimm to msa file 
    parser.add_argument('-M','--trimm_2_msa',required=True,help='output mapp of trimal allgorithm mapping original with trimmed msa files')
    parser.add_argument('-o','--outputfile',required=True,help='output file *.csv')

    args = parser.parse_args()
    msa = args.msa
    refseq = args.refseq    
    refgenome = args.refgenome  # character indicating version of reference to use (


    trimmed_2_msa = args.trimm_2_msa
    outputfile = args.outputfile
else:
    msa = '/data/EBV/msas/splited/sin_trimear/all150_sin_n.fas'
    refseq = '/data/EBV/byACCIDs/NC_007605.fasta'
    refgenome = 'NC_007605.1'

    trimmed_2_msa = '/data/EBV/msas/splited/mapp_trimed_columns.txt'
    outputfile = '/data/EBV/msas/splited/sin_trimear/NC_007605_coordinates_in_alignment.csv'


#    msa_2_ref = '/data/EBV/msas/splited/sin_trimear/sin_trimear_mapp_coordinates_to_NC007605.csv'


### run mapping original msa 2 reference 
msa_2_ref = msa + '_to_ref_%s.csv'%refgenome
os.system(' python /home/ariel/repo/TBI/EBV-db/src/ref2msa_mapper.py -m %s -r %s -o %s'%(msa,refseq,msa_2_ref))
mapp_to_reference = pd.read_csv(msa_2_ref,sep = '|')



# load and transform  trim2msa file
trimmed_out = pd.read_csv(trimmed_2_msa,sep=',',header=None)
trim_mapper = trimmed_out.transpose()
trim_mapper.columns = ['old_msa_position']
trim_mapper[trim_mapper=='#ColumnsMap\t0'] =0
trim_mapper.reset_index(inplace = True)
trim_mapper.columns = ['trimmedCoord','position_in_original_msa']



trim_to_ref = pd.merge(mapp_to_reference,trim_mapper,how = 'inner',left_on='msaCoord',
                       right_on='position_in_original_msa').drop('position_in_original_msa',axis = 1)

print 'mapping file looks like this\n'
print trim_to_ref.tail()
#print outputfile
trim_to_ref.to_csv(outputfile,sep ='|',index =False)



