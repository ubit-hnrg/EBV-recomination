from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
from operator import itemgetter
import re
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

mapfile = '/data/EBV/msas/ids_171_ebv/ids_171_ebv_msa_gap1.51.fa_whole_mapping.file'
outpath = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/10Miter/results/'
genbank_ref_file = '/data/EBV/byACCIDs/genBankRecord_NC_007605.gb'
whole_recom_rate = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/10Miter/whole_10Miters.tsv'
gap_windows = 10
gene_coords_file = '/data/EBV/by_gene/gene_result.txt'

def get_items(d,keylist):
    return itemgetter(*keylist)(d)


def search_fasta(pattern, file_path):
    st = []
    end = []
    pat =[]
    reverse = str(Seq(pattern).reverse_complement())
    
    for patt in [pattern,reverse]:
        for record in SeqIO.parse(open(file_path, "rU"), "fasta"):
            chrom = record.id
            for match in re.finditer(patt, str(record.seq)):
                start_pos = match.start() + 1
                end_pos = match.end() + 1
                st.append(start_pos)
                end.append(end_pos)
                pat.append(patt)
    df = pd.DataFrame({'st':st,'end':end,'pattern':pat,'call':pattern})
    df = df[['st','end','pattern','call']]
    df.sort_values(by=['st'],inplace = True)
    
    return(df)


def compute_local_obs(rrate,r1,r2):
    return(rrate[r1:r2].apply(np.mean))


# requiere map_dict  (reference vs trimed aligment mapping)
# require rrate: dataframe containing stimated recombination rate in whole genome. Its index must be a trimed position. 
# this function require a dataframe indicating regions with two coordinates columns indicating position in reference genome!
# output: tuple (same dataframe plus mapped columns, rho series)
def mapp_and_compute_rho(df,start_col = 'st',end_col = 'end'):
    df = df.dropna(subset=[start_col,end_col])
    start_aux = [x if x in mapp_dict.keys() else (x +1 ) for x in df[start_col].values]
    end_aux = [x if x in mapp_dict.keys() else (x -1 ) for x in df[end_col].values]


    df['trimmed_st'] = get_items(mapp_dict,start_aux)
    df['trimmed_end'] = get_items(mapp_dict,end_aux)

    computed_rates = df.apply(lambda x: pd.concat([x,compute_local_obs(rrate,r1 = x.trimmed_st,r2 = x.trimmed_end)]),axis =1)
    return df,computed_rates

###########################
####### mapp file #########
###########################
mapp = pd.read_csv(mapfile,sep = '|')
rrate = pd.read_csv(whole_recom_rate,sep = '|')
rrate.set_index(['trimmed_position'],inplace = True)
mapp_dict = mapp[['refCoord','trimmedCoord']]

mapp_dict =  {}
refcoords = mapp.loc[:,'refCoord']
trimcoords = mapp.loc[:,'trimmedCoord']
for i in range(mapp.shape[0]):
    mapp_dict.update({refcoords[i]:trimcoords[i]})

seq_record = SeqIO.read(genbank_ref_file,'genbank')
seq = seq_record.seq

G = seq.count('G')
C = seq.count('C')
A = seq.count('A')
T = seq.count('T')
l = len(seq_record)-kmer
proba_map = {'G':G/float(l),'C':C/float(l),'A':A/float(l),'T':T/float(l)}


full = pd.DataFrame()
full['st'] =np.arange(len(seq_record)-kmer-gap_windows)+gap_windows
full['end'] = full.st +5
full['call'] = full.apply(lambda x: seq[x.st:x.end],1)
vc = full.call.value_counts()
vc = vc.to_frame().reset_index()
vc.columns = ['motif','counts']
vc['len']=vc.motif.apply(lambda x: len(x))
vc_no_gap = vc[vc.len ==kmer]
vc_no_gap['expected'] = vc_no_gap.motif.apply(lambda x: int(l*np.prod(get_items(d=proba_map,keylist=x))))
vc_no_gap['ratio'] = vc_no_gap.counts/vc_no_gap.expected
vc_no_gap.sort_values(by=['ratio'],ascending = False,inplace = True)
vc_no_gap['mot'] = vc_no_gap.motif.apply(lambda x: str(x))
motif_statistics = vc_no_gap



#### fiter those without mapping 
ii = (full.st.isin(mapp_dict.keys())|(full.st+1).isin(mapp_dict.keys())) 
ii2 = ((full.st-gap_windows).isin(mapp_dict.keys())|(full.st-gap_windows+1).isin(mapp_dict.keys()))

jj = (full.end.isin(mapp_dict.keys())|(full.end-1).isin(mapp_dict.keys()))
jj2 = ((full.end+gap_windows).isin(mapp_dict.keys())|(full.end+gap_windows-1).isin(mapp_dict.keys()))
full_ok = full[ii&ii2&jj&jj2]
#################################

## this take ~45 minutes 
import time 
t1 = time.clock()
fullgaps = compute_gap(full_ok,window = gap_windows)
t2 = time.clock()
print t2-t1

fullgaps2 = fullgaps.copy()
mmin = fullgaps[['Mean Rho/bp._left','Mean Rho/bp._right']].apply(lambda x: np.nanmin(x)).min()
fullgaps2['Mean Rho/bp._left'].fillna(mmin,inplace = True)
fullgaps2['Mean Rho/bp._right'].fillna(mmin,inplace = True)

fullgaps2['diff'] =fullgaps2['Mean Rho/bp._right']-randgaps2['Mean Rho/bp._left']

fullgaps2['gap'] =fullgaps2.apply(lambda x:
                          (x['Mean Rho/bp._right']-x['Mean Rho/bp._left'])/np.min([x['Mean Rho/bp._right'],x['Mean Rho/bp._left']]),1)

fullgaps2.to_csv(outpath+'fullgaps_by_5mers.csv',index = False)

