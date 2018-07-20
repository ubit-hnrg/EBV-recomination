from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
from operator import itemgetter
import re
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

brown_repeats = ['TGGTGG','CCTCCCCT','TGGAG','AGGAG','CCCAG','GGGCT']
inverted_repeats_file = '/data/EBV/repeat_regions/inverted_positions_with_complemented.txt'
file_path = '/data/EBV/byACCIDs/NC_007605.fasta'
mapfile = '/data/EBV/msas/ids_171_ebv/ids_171_ebv_msa_gap1.51.fa_whole_mapping.file'
outpath = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/10Miter/results/'
genbank_ref_file = '/data/EBV/byACCIDs/genBankRecord_NC_007605.gb'
whole_recom_rate = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/10Miter/whole_10Miters.tsv'


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



# this function allow us to extract complete list of features of each desidered region
#repeat_region = complete_features(tipo ='repeat_region' ,n=len(tipo))
#seq_record  must be a Sequence load from a genbank file

def complete_features(seq_record, TIPO='regulatory',n=5 ):
    ttyp = []
    loc = []
    get = []
    pos = []
    for feature in seq_record.features:
        ttyp.append(feature.type)

        if TIPO[0:n] in feature.type:
            get.append({k: feature.qualifiers[k][0] for k in feature.qualifiers.keys()})

            partes = []
            for p in feature.location.parts:
                partes.append(pd.Series([p.start.position,p.end.position]))
            positions = pd.concat(partes,1).transpose()    
            positions.columns = ['start','end']
            pos.append(positions)        #print feature

    get = pd.DataFrame(get)        
    pos = pd.concat(pos)

#    return(get,pos)

    pos.index = get.index
    get_features = pd.concat([get,pos],1)
        
    return(get_features)


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



##### NCBI classification ####
print 'analyzing NCBI registred repeats'

seq_record = SeqIO.read(genbank_ref_file,'genbank')
tipo = 'repeat_region'
repeat_regions = complete_features(seq_record,TIPO=tipo,n=len(tipo))
repeat_regions.rename(columns = {'start':'st'},inplace=True)
_, ncbi_repeat_rates = mapp_and_compute_rho(repeat_regions,start_col = 'st')
ncbi_repeat_rates['call'] = ncbi_repeat_rates['rpt_type'].fillna('') + ncbi_repeat_rates['rpt_family'].fillna('')


############### motifs rates #############
print 'analyzing Brown motifs'
positions = {}
for br_rep in brown_repeats:
    df = search_fasta(br_rep,file_path)
    positions.update({br_rep:df})
    
repeat_postions = pd.concat(positions)
_, motifs_rates = mapp_and_compute_rho(repeat_postions)



#### inverted repeats ###
print 'analyzing Inverted repeats'

inverted = pd.read_csv(inverted_repeats_file,header = None,sep = '\t')
inv = pd.concat([inverted.iloc[:,1].str.split('>NC_007605.1_').apply(lambda x: (x[1].split('_')[0])),
                inverted.iloc[:,1].str.split('>NC_007605.1_').apply(lambda x: (x[1].split('_')[1]))],1)
inv.columns = ['st','end']
inv['st'] = pd.to_numeric(inv['st']) 
inv['end'] = pd.to_numeric(inv['end']) 

inv_df, inverted_rates =  mapp_and_compute_rho(df=inv.copy(),start_col = 'st',end_col = 'end')
inverted_rates['call']='inverted'


######## tandem repeats analysis ########## 
print 'analyzing Tandem repeats'

tandem = pd.read_csv('/data/EBV/repeat_regions/tandem_table.csv',sep =',')
tandem.columns = ['st','end']
_ , tandem_rates = mapp_and_compute_rho(tandem)
tandem_rates['call']='tandem'



print 'merging and writing to disk'
whole_repeats = pd.concat([inverted_rates,tandem_rates,ncbi_repeat_rates])
whole_repeats.to_csv(outpath+'repeats_analysis.csv',index = False)

initiators = motifs_rates.copy()
initiators.to_csv(outpath+ '/initiators_analysis.csv',index = False)



### boxplot over repeats
thr = np.nanmean(whole_repeats['Mean Rho/bp.'])
#print thr
aux = whole_repeats[~whole_repeats['Mean Rho/bp.'].isnull()].copy()
aux['call'] = aux.call.replace('','NC.')
#aux = whole_repeats.copy()


order =  aux.groupby(['call'])['Mean Rho/bp.'].apply(lambda x:np.median(x)).sort_values(ascending=False).index
ax = sns.boxplot(x ='call',y = 'Mean Rho/bp.', data = aux,order = order)
ax.set_yscale('log')
ax.set_xlabel('repeat region')
ax.set_ylabel('4Ne.r/bp')
plt.savefig(outpath+'boxplot_repetitive_regions'+'.png',dpi = 500)
plt.clf()

#### boxplot over initiators
aux = initiators[initiators['Mean Rho/bp.']>thr].copy()
order =  aux.groupby(['call'])['Mean Rho/bp.'].apply(lambda x:np.median(x)).sort_values(ascending=False).index
ax = sns.boxplot(x ='call',y = 'Mean Rho/bp.', data = aux,order = order)
ax.set_yscale('log')
ax.set_xlabel('motifs')
ax.set_ylabel('4Ne.r/bp')
plt.savefig(outpath+'boxplot_motifs'+'.png',dpi = 500)




