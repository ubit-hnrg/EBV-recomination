# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import argparse   
import pandas as pd
import os
from Bio import AlignIO 
from Bio import SeqIO
import numpy as np
import scipy.integrate as integrate
from operator import itemgetter
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle

parser = argparse.ArgumentParser(description='analyze results of recombination rate - how this is distribuited along the genome')
if(True):
    parser.add_argument('-m','--mapfile',required=True,help='mapping file between reference genome coordinates, trimmed alignment and original alignment')
    parser.add_argument('-r','--recomb_rate',required=True,help='whole recombination rate file')
    parser.add_argument('-o','--outpath',required=True,help='reference sequence file in fasta format')
    parser.add_argument('-g','--gene_coords_file',default = '/data/EBV/by_gene/gene_result.txt',help='this file provivde gene coordinate mapp in reference genome')
    parser.add_argument('-G','--genbank_ref_file',default = '/data/EBV/byACCIDs/genBankRecord_NC_007605.gb' ,help='this file provide coordinates for region analysis')
    parser.add_argument('-s','--suffix',default = '' ,help='suffix for output files')

    args = parser.parse_args()
    mapfile = args.mapfile
    whole_recom_rate = args.recomb_rate
    gene_coords_file = args.gene_coords_file
    outpath = args.outpath
    genbank_ref_file = args.genbank_ref_file    
    suffix = args.suffix    

else:
    mapfile = '/data/EBV/msas/ids_171_ebv/ids_171_ebv_msa_gap1.51.fa_whole_mapping.file'
    whole_recom_rate = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/whole_1Miter.tsv'
    gene_coords_file = '/data/EBV/by_gene/gene_result.txt'
    outpath = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/'
    genbank_ref_file = '/data/EBV/byACCIDs/genBankRecord_NC_007605.gb'
    suffix = ''
    
def get_items(d,keylist):
    return itemgetter(*keylist)(d)


def mean_value_by_integration(serie):
    mean_value = integrate.trapz(serie.values, serie.index.astype(np.int64))/(serie.index.max()-serie.index.min())
    return(mean_value)

#r1,r2 son valores del indice
def compute_local_obs(rrate,r1,r2):
    return(rrate[r1:r2].apply(np.mean)) 

def gc_content(seq):
    gs = seq.count('G')
    cs = seq.count('C')
    return (gs+cs)/float(len(seq))


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


##################################
####### gene coordinates #########
##################################

gene_coords = pd.read_csv(gene_coords_file,sep = '\t')
gene_ref_coords = gene_coords[['GeneID','Symbol','Aliases','start_position_on_the_genomic_accession','end_position_on_the_genomic_accession']]

start_col = 'start_position_on_the_genomic_accession'
end_col = 'end_position_on_the_genomic_accession'

st = gene_ref_coords[start_col]#[0]
end = gene_ref_coords[end_col]

### muevo las coordenadas sin mapping a otras desplazadas una unidad dentro del gen para no perder el mapeo
# tomo la coordenada de la region o la siguiente (o anterior) para mapearlas al alineamiento
ref_table = gene_ref_coords.dropna(subset=[start_col,end_col])
start_aux = [x if x in mapp_dict.keys() else (x +1 ) for x in ref_table[start_col].values]
end_aux = [x if x in mapp_dict.keys() else (x -1 ) for x in ref_table[end_col].values]

ref_table['trimmed_st'] = get_items(mapp_dict,start_aux)
ref_table['trimmed_end'] = get_items(mapp_dict,end_aux)



#########################################
##############   Aalysis   ##############
#########################################


rho_mean, rho_low, rho_up = rrate.apply(mean_value_by_integration)

result_by_gene = ref_table.apply(lambda x: pd.concat([x[0:3],rrate[x.trimmed_st:x.trimmed_end].apply(mean_value_by_integration),x[3::]]),axis =1)
result_by_gene['length_aln'] = ref_table.trimmed_end-ref_table.trimmed_st
result_by_gene['length'] = ref_table.end_position_on_the_genomic_accession-ref_table.start_position_on_the_genomic_accession
result_by_gene.sort_values(by=['Mean Rho/bp.'],inplace=True,ascending = False)


print 'sin señal %s genes'%result_by_gene['Mean Rho/bp.'].isnull().sum()
print 'CON señal %s genes'% (~result_by_gene['Mean Rho/bp.'].isnull()).sum()

out1 = outpath+'result_by_gene' + suffix+ '.csv'
result_by_gene.to_csv(out1,index = False)


#########################################
####### region characterization #########
#########################################
seq_record = SeqIO.read(genbank_ref_file,'genbank')

types = []
cds = []
region =[]
for feature in seq_record.features:
    types.append(feature.type)
    if 'gene' in feature.qualifiers.keys():
        gene = feature.qualifiers['gene']
    else:
        gene = 'None'
    
    partes = []
    for p in feature.location.parts:
        partes.append(pd.Series([p.start.position,p.end.position]))
    positions = pd.concat(partes,1).transpose()    
    positions.columns = ['start','end']
    positions['type'] = feature.type
    positions['gene'] = np.repeat(gene,positions.shape[0])
    region.append(positions)
regions = pd.concat(region)

regions = regions[regions.type!='source']

# compute length of each region
regions['length'] = regions.end - regions.start
regions['length'] = pd.to_numeric(regions.length)

regions = regions[regions.length < 40000] # con esto tiramos la anomalia del primer gen. LMP2

regions['rango'] = regions.apply(lambda x: (x.start,x.end),1)
regions.set_index('rango',inplace = True)


# tomo la coordenada de la region o la siguiente (o anterior) para mapearlas al alineamiento
start_prima = [x if x in mapp_dict.keys()  else (x +1 ) for x in regions.start.values]
end_prima = [x if x in mapp_dict.keys()  else (x -1 ) for x in regions.end.values]

# mapp to trimmed alignment
regions['trimmed_st'] = get_items(mapp_dict,start_prima)
regions['trimmed_end'] = get_items(mapp_dict,end_prima)
regions['length_aln'] = regions['trimmed_end'] - regions['trimmed_st']


result_by_region = regions.apply(lambda x: pd.concat([x[2:4],rrate[x.trimmed_st:x.trimmed_end].apply(mean_value_by_integration),x[0:2],x[-4::]]),axis =1)
result_by_region.sort_values(by=['Mean Rho/bp.'],inplace=True,ascending = False)


out2 = outpath+'result_by_region'+suffix +'.csv'
result_by_region.to_csv(out2,index = False)



###############  WHOLE RECOMBINATION PLOT ###############
#########################################################
rrates_reference = pd.merge(mapp[['refCoord','trimmedCoord']],rrate.reset_index(),how='inner',left_on='trimmedCoord',right_on='trimmed_position')
rrates_reference.head()
rrates_reference.set_index(['refCoord'],inplace = True)
rrates_reference.drop(['trimmedCoord','trimmed_position'],1,inplace = True)
#rates_reference['Mean Rho/bp.'].plot()

max_scale = rrates_reference.reset_index().refCoord.max()
f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.hist(rrates_reference.reset_index().refCoord,bins = 300)
ax1.set_xlim([0,max_scale])
ax1.set_title('segregating sites coverage')
ax1.set_ylabel('counts')
ax1.grid() 

ax2.plot(rrates_reference['Mean Rho/bp.'])
ax2.set_xlabel('NC_0007605 reference possition')
ax2.set_title('population rho')
ax2.set_ylabel('4Ne.r/bp')
ax2.grid()
f.savefig(outpath+'whole_rho'+ suffix +'.png', dpi=500)
f.clf()


################  BOXPLOT ############################
#########################################################
size_thr = 10
ok = result_by_region.type[~result_by_region['Mean Rho/bp.'].isnull()].value_counts()
ok = ok[ok>=size_thr].index
filter_data = result_by_region[result_by_region.type.isin(ok)]
mm = filter_data.groupby(['type'])['Mean Rho/bp.'].median()
order =  mm.sort_values(ascending=False).index

ax = sns.boxplot(x ='type',y = 'Mean Rho/bp.', data = filter_data,order = order)
ax.set_yscale('log')
ax.set_ylabel('4Ne.r/bp')
ax.set_xlabel('genome region')

### add mean range of whole genome ##

x =ax.get_xticks() 
y,y_low,y_up = rrate.apply(mean_value_by_integration)
currentAxis = plt.gca()
currentAxis.add_patch(Rectangle((np.min(x)-1,y_low), np.max(x)-np.min(x)+2,y_up-y_low, facecolor="grey",alpha = 0.3) )   



for item in ax.get_xticklabels():
    item.set_rotation(45)
plt.gcf().subplots_adjust(bottom=0.25)
plt.savefig(outpath+'boxplot_by_region'+suffix+'.png',dpi = 500)



#################  Gen behaviour ##########################



