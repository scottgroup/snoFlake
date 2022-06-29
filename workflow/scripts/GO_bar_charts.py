#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv('../../results/SNORD22/gProfiler_web/PRPF8_uniq_GO.csv',header=17)
print(df.columns)


plt.figure(figsize=(9,4))
# GO colors
#colours = {"GO:BP":"indianred","GO:MF":"limegreen","GO:CC":"steelblue"}
colours = {"GO:BP":"indianred","GO:CC":"steelblue"}
#df_mf = df[df['source']=='GO:MF'].sort_values(by=['negative_log10_of_adjusted_p_value'])
df_cc = df[df['source']=='GO:CC'].sort_values(by=['negative_log10_of_adjusted_p_value'])
df_bp = df[df['source']=='GO:BP'].sort_values(by=['negative_log10_of_adjusted_p_value'])
#df = pd.concat([df_cc,df_mf,df_bp],ignore_index=True)
df = pd.concat([df_cc,df_bp],ignore_index=True)
bars = pd.DataFrame({'source':df.source.values.tolist(),'p_val':df.negative_log10_of_adjusted_p_value.values.tolist()},index=df.term_name.values.tolist())

bars['p_val'].plot(kind="barh",color=bars['source'].replace(colours),width=0.7)
plt.title('GO Analysis of SNORD22 & PRPF8 Targets')
plt.xlabel('$-log_{10}$(p-value)')

labels = list(colours.keys())
handles = [plt.Rectangle((0,0),1,1, color=colours[label]) for label in labels]
plt.legend(handles,labels,loc='lower right')

plt.tight_layout()
plt.savefig('../../results/SNORD22/gProfiler_web/PRPF8_uniq_GO.png')