import argparse   
import pandas as pd
import os
from Bio import AlignIO 
from Bio import SeqIO
import numpy as np



# gen coordinate mapping
gentablepath = '/home/ariel/Downloads/gene_result.txt'
coordinate_mapping_file = '/data/EBV/msas/Africa/mapp_coordinates_Africa_msa_gap100.csv' ## output of ref2msa_mapper


msafile = '/data/EBV/msas/Africa/A_msa_gap100.fa'
msa = AlignIO.read(msafile,'fasta')
base = os.path.basename(msafile)


outputpath = '/data/EBV/by_gen/Africa/'
if not os.path.isdir(outputpath):
    os.system('mkdir %s'%outputpath)
    

min_len = 50 
max_len = 10000
coordinate_table = pd.read_csv(gentablepath,sep='\t',
                                   usecols=['GeneID','Symbol','Aliases','description',
                                            'start_position_on_the_genomic_accession',
                                            'end_position_on_the_genomic_accession'])

## mapp
sm = pd.read_csv(coordinate_mapping_file,sep='|')
sm.set_index(['refCoord'],inplace=True)
series_mapp = sm.msaCoord
series_mapp.tail()



coordinate_table['start_msa'] = series_mapp[coordinate_table.start_position_on_the_genomic_accession].values
coordinate_table['end_msa'] = series_mapp[coordinate_table.end_position_on_the_genomic_accession].values

coordinate_table['original_len']=coordinate_table['end_position_on_the_genomic_accession']-coordinate_table['start_position_on_the_genomic_accession']
coordinate_table['msa_len']=coordinate_table['end_msa']-coordinate_table['start_msa']

gentable = coordinate_table[(coordinate_table.original_len>=min_len)&(coordinate_table.original_len<=max_len)]
gentable.set_index(['GeneID'],inplace = True)

    

for g in gentable.index:
    st = gentable.ix[g,:]['start_msa']
    end = gentable.ix[g,:]['end_msa']
    cutted = msa[:,int(st-1):int(end-1)] # biopython start from 0 !!
    output = outputpath+str(g)+'_'+base
    AlignIO.write(cutted,output,'fasta')
    

    