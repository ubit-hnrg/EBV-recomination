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


def pos_to_trim(p):
    if type(p) in [int,float] :
        tr = [x if x in mapp_dict.keys() else (x -1 ) for x in [p]]
    else:
        tr = [x if x in mapp_dict.keys() else (x -1 ) for x in p]
    return(get_items(mapp_dict,tr))


#df provide start- end positions
def compute_gap(df,start_col='st',end_col='end',window = 10):
    mmax = mapp.refCoord.max()
    col ='Mean Rho/bp.'
    start_aux = [x if x in mapp_dict.keys() else (x +1 ) for x in df[start_col].values]
    end_aux = [x if x in mapp_dict.keys() else (x -1 ) for x  in df[end_col].values]
    trimmed_start = get_items(mapp_dict,start_aux)
    trimmed_end = get_items(mapp_dict,end_aux)

    start_gap = [x-window if x-window in mapp_dict.keys() else (x -window +1 ) for x in df[start_col].values]
    end_gap = [x +window if x+window in mapp_dict.keys() else (x +window -1 ) for x in df[end_col].values]
    end_gap = [x if (x <= mmax) else mmax for x in end_gap]
    
    trimmed_start_gap = get_items(mapp_dict,start_gap)
    trimmed_end_gap = get_items(mapp_dict,end_gap)
    
    
    df['trimmed_st'] = trimmed_start
    df['trimmed_end'] = trimmed_end

    df['trimmed_st_gap'] = trimmed_start_gap
    df['trimmed_end_gap'] = trimmed_end_gap

    rate_prev = df.apply(lambda x: pd.concat([x,compute_local_obs(rrate,r1 = x.trimmed_st_gap,r2 = x.trimmed_st)]),axis =1)
    rate_post = df.apply(lambda x: pd.concat([x,compute_local_obs(rrate,r1 = x.trimmed_end,r2 = x.trimmed_end_gap)]),axis =1)
    
    
    return(pd.merge(rate_prev[[start_col,end_col,'call',col]],rate_post[col].to_frame(),right_index = True,
                    left_index = True,suffixes=['_left','_right']))

def look_for_gene(position,gene_coords):
    res = gene_coords[(position> gene_coords.start_position_on_the_genomic_accession) & (position<gene_coords.end_position_on_the_genomic_accession)]
    m = res[['Symbol','gen_len']].drop_duplicates()
    m.sort_values(by=['gen_len'],ascending = True)
    if m.shape[0]>0:
        return(m.Symbol[0])
    else:
        return('intr')
    


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

##########################

rrate = pd.read_csv(whole_recom_rate,sep = '|')
rrate.set_index(['trimmed_position'],inplace = True)

mapp_dict =  {}
refcoords = mapp.loc[:,'refCoord']
trimcoords = mapp.loc[:,'trimmedCoord']


###########################################
######### locate brown motifs ###

print 'analyzing Brown motifs'
positions = {}
for br_rep in brown_repeats:
    df = search_fasta(br_rep,file_path)
    positions.update({br_rep:df})
    
repeat_postions = pd.concat(positions)


####### compute gaps (jumps) arround brown motifs

gaps = compute_gap(repeat_postions,window=gap_windows)
gaps2 = gaps.copy()
mmin = gaps[['Mean Rho/bp._left','Mean Rho/bp._right']].apply(lambda x: np.nanmin(x)).min()
gaps2['Mean Rho/bp._left'].fillna(mmin,inplace = True)
gaps2['Mean Rho/bp._right'].fillna(mmin,inplace = True)

gaps2['diff'] =gaps2['Mean Rho/bp._right']-gaps2['Mean Rho/bp._left']

gaps2['gap'] =gaps2.apply(lambda x:
                          (x['Mean Rho/bp._right']-x['Mean Rho/bp._left'])/np.min([x['Mean Rho/bp._right'],x['Mean Rho/bp._left']]),1)


###################################################
######   get a null hipotesys for those gaps  #####
###################################################

### prepare random 5 mers
kmer = 5
rand = np.random.choice(np.arange(mapp.refCoord.max()-(2*gap_windows+2))+gap_windows,1000,replace=False)
rand.sort()
np.random.seed(123)

rand_df = pd.DataFrame({'st':rand,'end':rand+kmer})
rand_df['call'] = 'random'


### compute random gaps
randgaps = compute_gap(rand_df,window = gap_windows)
randgaps2 = randgaps.copy()
mmin = randgaps[['Mean Rho/bp._left','Mean Rho/bp._right']].apply(lambda x: np.nanmin(x)).min()
randgaps2['Mean Rho/bp._left'].fillna(mmin,inplace = True)
randgaps2['Mean Rho/bp._right'].fillna(mmin,inplace = True)

randgaps2['diff'] =randgaps2['Mean Rho/bp._right']-randgaps2['Mean Rho/bp._left']

randgaps2['gap'] =randgaps2.apply(lambda x:
                          (x['Mean Rho/bp._right']-x['Mean Rho/bp._left'])/np.min([x['Mean Rho/bp._right'],x['Mean Rho/bp._left']]),1)

## trhesholds
right_tresh = randgaps2[['gap']].quantile(0.975)['gap']
left_tresh = randgaps2[['gap']].quantile(0.025)['gap']

right_both_side_tresh = randgaps2[['gap']].quantile(0.975)['gap']
left_both_side_tresh = randgaps2[['gap']].quantile(0.025)['gap']


###########################################################
###### analyze separately plus and minus orientation  #####
###########################################################
gene_coords = pd.read_csv(gene_coords_file,sep = '\t')
gene_coords['gen_len']= gene_coords['end_position_on_the_genomic_accession'] - gene_coords['start_position_on_the_genomic_accession']

gene_ref_coords = gene_coords[['GeneID','Symbol','Aliases','start_position_on_the_genomic_accession','end_position_on_the_genomic_accession','gen_len']]
gene_ref_coords['rang'] = gene_coords.apply(lambda x: (x.start_position_on_the_genomic_accession,
                                                   x.end_position_on_the_genomic_accession),1)
gene_ref_coords.set_index(['rang'],inplace = True)

gaps3 = gaps2.reset_index()
#gene =  gaps3.st.apply(lambda x:look_for_gene(position=x,gene_coords=gene_ref_coords[gene_ref_coords.Symbol!='LMP2']))
gene =  gaps3.st.apply(lambda x:look_for_gene(position=x,gene_coords=gene_ref_coords))

gaps3['GENE'] = gene
rep_with_gene = gaps3.copy()

#rep_with_gene.GENE.replace({'intr':'LMP2'},inplace = True)

rep_with_gene_orient = pd.merge(rep_with_gene,gene_coords[['Symbol','orientation']],left_on=['GENE'],right_on='Symbol',how = 'left').drop(['Symbol'],1)
rep_with_gene_plus = rep_with_gene_orient[rep_with_gene_orient.orientation == 'plus']
rep_with_gene_minus = rep_with_gene_orient[rep_with_gene_orient.orientation == 'minus']
rep_with_gene_intron = rep_with_gene_orient[rep_with_gene_orient.orientation.isnull()]


initiators_PLUS = rep_with_gene_plus[rep_with_gene_plus.gap>right_tresh][['st','end','call','Mean Rho/bp._left','Mean Rho/bp._right','GENE','orientation']]#.value_counts().sort_index()#/rep_with_gene_plus['call'].value_counts()

initiators_MINUS = rep_with_gene_minus[rep_with_gene_minus.gap<left_tresh][['st','end','call','Mean Rho/bp._left','Mean Rho/bp._right','GENE','orientation']]#.value_counts()#/rep_with_gene_plus['call'].value_counts()

initiators_MINUS.drop_duplicates(inplace = True)
init_minus_bubble = initiators_MINUS.groupby(['GENE','call'])['Mean Rho/bp._left'].mean().to_frame()
init_minus_bubble.to_csv(outpath+'initiators_minus_wind10.csv')


initiators_PLUS.drop_duplicates(inplace = True)
initiators_PLUS = initiators_PLUS[initiators_PLUS['Mean Rho/bp._right']<1]
init_plus_bubble = initiators_PLUS.groupby(['GENE','call'])['Mean Rho/bp._right'].mean().to_frame()
init_plus_bubble.to_csv(outpath+'initiators_plus_wind10.csv')


## logo plots:
df = pd.concat([initiators_MINUS.call.value_counts(), initiators_PLUS.call.value_counts()],1).sum(1).apply(lambda x: np.log2(x))
df.sort_values(inplace=True,ascending = False)
import matplotlib.pyplot as plt
import gact

def draw_logo(score_list):
    fig, ax = plt.subplots(figsize=(5,3))

    all_scores = score_list

    x = 1
    maxi = 0
    for scores in all_scores:
        y = 0
        for base, score in scores:
            gact.letterAt(base, x,y, score, ax)
            y += score
        x += 1
        maxi = max(maxi, y)

    plt.xticks(range(1,x))
    plt.xlim((0, x)) 
    plt.ylim((0, maxi)) 
    plt.tight_layout()      

to_plot = []
for i in range(df.shape[0]):
    r = df.index[i]
    score = ((df.values[i]+1)/float(len(r)))
    motif = []
    for s in r:
        motif.append((s,score))
    to_plot.append(motif)

draw_logo(to_plot)

plt.ylabel('# initiation events [Bits]')
plt.xlabel('Motifs')
plt.tight_layout()
plt.savefig(outpath+'motif_logo'+'.png',dpi = 400)


