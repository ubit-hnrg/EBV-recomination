import numpy as np
import pandas as pd

mapfile = '/data/EBV/msas/ids_171_ebv/ids_171_ebv_msa_gap1.51.fa_whole_mapping.file'
whole_recom_rate = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/10Miter/whole_10Miters.tsv'






####### map for all downstream analysis####
###########################################
mapp = pd.read_csv(mapfile,sep = '|')
mapp_dict = mapp[['refCoord','trimmedCoord']]

rrate = pd.read_csv(whole_recom_rate,sep = '|')
rrate.set_index(['trimmed_position'],inplace = True)

mapp_dict =  {}
refcoords = mapp.loc[:,'refCoord']
trimcoords = mapp.loc[:,'trimmedCoord']
for i in range(mapp.shape[0]):
    mapp_dict.update({refcoords[i]:trimcoords[i]})


def get_rate_by_refefrence_genome(rrate,mapp):
    rrates_reference = pd.merge(mapp[['refCoord','trimmedCoord']],rrate.reset_index(),how='inner',left_on='trimmedCoord',right_on='trimmed_position')
    rrates_reference.head()
    rrates_reference.set_index(['refCoord'],inplace = True)
    rrates_reference.drop(['trimmedCoord','trimmed_position'],1,inplace = True)
    return(rrates_reference)
    
def get_coverage(rrates_reference):
    s = rrates_reference.reset_index().refCoord
    count, division = np.histogram(s,300)
    countbin = [1 if c != 0 else 0 for c in count]
    diff = division[1::]-division[:-1]
    return(np.sum(diff*countbin)/max(division))
    

rate_ref = get_rate_by_refefrence_genome(rrate,mapp)
print get_coverage(rate_ref)