import argparse   
import pandas as pd
import glob
import numpy as np
import os
from Bio import AlignIO 
from Bio import SeqIO



# input : path and patter of files to be analized
#example

#python ~/repo/TBI/EBV-db/src/unify_rdp4_outputs.py -p ~/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/ -n 1Miter -m /data/EBV/msas/ids_171_ebv/ids_171_ebv_msa_gap1.51_trimmed0.05.fas


parser = argparse.ArgumentParser(description='unify rdp4 results, remapping to absolute coordinates')

if(True):
    parser.add_argument('-p','--path',required=True,help='path')
    parser.add_argument('-n','--niters',help='number of iterations used in rdp4 runs',default = None)
    parser.add_argument('-m','--original_msa_file',required=True,help='path to trimmed msa file (or the file was used in RDP program)')
    
    args = parser.parse_args()
    path = args.path
    niters = args.niters
    msa_path = args.original_msa_file

output_file = path + 'whole_'+niters+'.tsv'


##### WARNING #####
##############################
# required file format :
#'start_end_Niters.txt'

# example 
#'60001_120000_1Miter.txt'
def rdp_absolute_coordinate_table(rdp4_patern_file):
    rdp4_files = glob.glob(rdp4_patern_file)
    basenames = [os.path.basename(f) for f in rdp4_files]
    coords = []
    for f in basenames:
        serie = pd.Series(f.split('.txt')[-2].split('_')[-3:-1])
        coords.append(serie)
    coords = pd.concat(coords,1).transpose()
    coords.columns = ['start','end']
    coords.index = basenames
    return(coords)

## extract starting coordinate from each file
coords = rdp_absolute_coordinate_table(path+'*'+niters+'.txt')

# extract total aligment length of original msa file (
msa = AlignIO.read(msa_path,'fasta')
l = msa.get_alignment_length()

whole_data = []
cols = ['Position','Mean Rho/bp.','-95% CI','+95% CI']

for f in coords.index:
    data = pd.read_csv(path + f)
    data['Position in alignment'] = data['Position in alignment'] + int(coords.loc[f].start) -1
    correct_positions = [c - l if c > l else c for c in data['Position in alignment'].values]
    data['Position in alignment'] = correct_positions
    data.columns = cols
    whole_data.append(data)
whole_data = pd.concat(whole_data)


averaged = whole_data.copy().groupby(['Position']).apply(np.mean)
averaged.sort_values(by=['Position'],inplace = True)
averaged.rename(columns = {'Position':'trimmed_position'},inplace = True)
#averaged['trimmed_position'] = [int(x) for x in averaged['trimmed_position']]
averaged.to_csv(output_file, sep = '|',index = False)








