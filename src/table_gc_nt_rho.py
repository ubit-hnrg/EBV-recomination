#### important note: 
#this script take as input the result of Rscript nuc_div.R. 
import argparse   
import pandas as pd
from outliers import smirnov_grubbs as grubbs


parser = argparse.ArgumentParser(description='analyze results of recombination rate - how this is distribuited along the genome')
if(True):
    parser.add_argument('-i','--input_file',required=True,help='reocombination rate by table with nt and gc contennt')
    parser.add_argument('-o','--outpath',required=True,help='output path')

    args = parser.parse_args()
    f = args.input_file
    outpath = args.outpath

else:
    f = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/10Miter/results/result_by_region_10Miters.csv_with_ndiv.csv'
    outpath = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/results/'



pval_thr = 0.00001



nd_regions = pd.read_csv(f,sep=',')
ok = nd_regions.type[~nd_regions['Mean.Rho.bp.'].isnull()].value_counts()
ok = ok[ok>=10].index
filter_data = nd_regions[nd_regions.type.isin(ok)].copy()

# fill na 
mmin = filter_data.NucDiv.min()
filter_data.NucDiv.fillna(mmin,inplace = True)




################# outliers by region ##################
#found outliers by region
gb = filter_data.groupby(['type'])['Mean.Rho.bp.']
outliers_by_region=gb.apply(lambda x: grubbs.max_test_outliers(x.dropna(), alpha=pval_thr))


outliersByReg = []
for ty in outliers_by_region.index:
    ii = filter_data[(filter_data.type==ty)&(filter_data['Mean.Rho.bp.'].isin(outliers_by_region[ty]))].index.values
    if len(ii)>0:
        outliersByReg.append(ii[0])

outlier_data = filter_data.loc[outliersByReg,] 
        
filter_data = filter_data.drop(outliersByReg)        

mm = filter_data.groupby(['type'])[['NucDiv','GC','Mean.Rho.bp.', u'X.95..CI', u'X.95..CI.1','length']].mean()
mstd = filter_data.groupby(['type'])[['NucDiv','GC','length']].std()
mstd.columns = ['std-nt','std-GC','std-legth']
mm = pd.merge(mm,mstd,left_index=True,right_index=True,how = 'inner')

mm.sort_values(by='Mean.Rho.bp.',ascending = False,inplace = True)
mm = mm.applymap(lambda x: round(x,4))

outlier_data.to_csv(outpath+'dropped_outliers.csv')
mm.to_csv(outpath+'table_GC_NTdiv_rho.csv')