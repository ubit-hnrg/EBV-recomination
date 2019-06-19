import pandas as pd
import numpy as np
from matplotlib import pylab as plt


mapfile = '/data/EBV/msas/ids_171_ebv/ids_171_ebv_msa_gap1.51.fa_whole_mapping.file'
whole_recom_rate = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/10Miter/whole_10Miters.tsv'
output_path = '/home/ariel/Projects/Gutierrez/EBV-recomb/recomb/rdp4_results/ids_171_ebv/10Miter/plots/'

mapp = pd.read_csv(mapfile,sep = '|')
rrate = pd.read_csv(whole_recom_rate,sep = '|')
rrate.set_index(['trimmed_position'],inplace = True)
mapp_dict = mapp[['refCoord','trimmedCoord']]


rrates_reference = pd.merge(mapp[['refCoord','trimmedCoord']],rrate.reset_index(),how='inner',left_on='trimmedCoord',right_on='trimmed_position')
rrates_reference.head()
rrates_reference.set_index(['refCoord'],inplace = True)
rrates_reference.drop(['trimmedCoord','trimmed_position'],1,inplace = True)
rrates_reference.to_csv(output_path+'whole_rho_data_10Miter.csv', sep = '|')
#rates_reference['Mean Rho/bp.'].plot()




###############  WHOLE RECOMBINATION PLOT ###############
#########################################################
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
f.savefig(output_path+'whole_rho.png', bbox_inches='tight')


