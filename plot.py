import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import numpy as np

import seaborn as sns
file1="./LAMP_results/20200428 OB O1 NA 60C for 1hr -  Quantification Amplification Results_SYBR.csv"
file2="./LAMP_results/20200504 NA vs N4 WS vs Bst2.0 65C -  Quantification Amplification Results_SYBR.csv"
file3="./LAMP_results/20200505 0.5-1-2X N4 OR 100-300-1e3 Copy in Saliva -  Quantification Amplification Results_SYBR.csv"
file4="./LAMP_results/20200506 RNA w N4 10-30-100-300-1e3 cp -  Quantification Amplification Results_SYBR.csv"

df1 = pd.read_csv(file1,usecols=lambda x:not x.startswith('Unnamed'),index_col='Cycle')
df2 = pd.read_csv(file2,usecols=lambda x:not x.startswith('Unnamed'),index_col='Cycle')
df3 = pd.read_csv(file3,usecols=lambda x:not x.startswith('Unnamed'),index_col='Cycle')
df4 = pd.read_csv(file4,usecols=lambda x:not x.startswith('Unnamed'),index_col='Cycle')

df1.head()

df3.loc[:,df4.max()>1000].plot()
df1.loc[:,df1.max()>1000].plot()

df2.loc[:,df2.max()>3000].plot()

df3['C5'].plot()
df2m = df2.loc[:,df2.max()>1000]
df2m.loc[:,df2m.max()<3000].plot()
df2.loc[:,(df2.min()<-100) & (df2.max()<1000)].plot()
df3.loc[:,df3.max()>1000].plot()
# generate SARS COV vs SARS COV 2 on N4
N4S = pd.concat([df2['C6'],df2['C5'],df3['C5']],axis=1)
N4Snorm=N4S.apply(lambda x: x/N4S.max().max()*100,axis=1)
N4Snorm.plot()
colors = ['blue','cornflowerblue','deepskyblue','cyan','limegreen','black']

fig,ax = plt.subplots(1,1,figsize=(6,4))
ax.set_xlim([0,60])
ax.set_xticks([0,20,40,60])
ax.set_xlabel('Time (mins)')
#
# for col,color in zip(N4Snorm.columns,colors):
ax.plot(N4Snorm.index[0:55],N4Snorm.loc[0:55,:],color='blue',alpha=0.7)
ax.set_title('N4 primer with Bat SARS-CoV 2015 (MG772933.1)')
ax.set_ylabel('RFU')
ax.legend(['1e4'],loc='lower right')
plt.tight_layout()
plt.savefig('Bat SARS-CoV.svg')
N4Snorm.to_csv('Bat SARS-CoV.csv ')



# generate figure 2
N4 = pd.concat([df2['B4'],df4.loc[:,['C3',]],df4['C7'],df4['C4'],df1['A1'],df1['A2']],axis=1)
N4norm=N4.apply(lambda x: x/N4.max().max()*100,axis=1)
NA = pd.concat([df2['A8'],df2['A4'],df2['D4'],df2['F3'],df2['C1'],df2['C3']],axis=1)
NAnorm=NA.apply(lambda x: x/NA.max().max()*100,axis=1)
N5 = pd.concat([df2['E4'],df2['F8'],df2['F5'],df2['F6'],df2['F7'],df2['C8']],axis=1)
N5norm = N5.apply(lambda x: x/N5.max().max()*100,axis=1)
N5norm.plot()
colors = ['blue','cornflowerblue','deepskyblue','cyan','limegreen','black']

fig,axes = plt.subplots(1,3,figsize=(12,3))
ax1,ax2,ax3=axes
for ax in axes:
    ax.set_xlim([0,60])
    ax.set_xticks([0,20,40,60])
    ax.set_xlabel('Time (mins)')

for col,color in zip(N4.columns,colors):
    ax1.plot(N4.index[0:55],N4norm[col][0:55],color=color,alpha=0.7)
ax1.set_title('N4 primer')
ax1.set_ylabel('RFU')
fig.legend(['1e3','500','300','100','30','NTC'],loc='center right')


for col,color in zip(NA.columns,colors):
    ax2.plot(NA.index,NAnorm[col],color=color,alpha=0.7)
ax2.set_title('NA primer')


for col,color in zip(N5.columns,colors):
    ax3.plot(N5.index,N5norm[col],color=color,alpha=0.7)
ax3.set_title('N5 primer')
plt.tight_layout()
plt.savefig('Fig2.svg')

fig2df = pd.concat([N4,NA,N5,],axis=1)

fig2df.to_csv('fig2.csv')




# generate saliva data.
Ns = pd.concat([df3['E3'],df3['E5'],df3['D7'],df3['D8'],df2['A5']+np.array(range(54))*3,df2['D1'],  ],axis=1)
Nsnorm = Ns.apply(lambda x: x/Ns.max().max()*100,axis=1)
colors = ['blue','cornflowerblue','deepskyblue','cyan','limegreen','black']

fig,ax = plt.subplots(figsize=(6,4))
ax.set_xlim([0,60])
ax.set_xticks([0,20,40,60])
ax.set_xlabel('Time (mins)')
ax.set_ylim([-10,110])

for col,color in zip(Nsnorm.columns,colors):
    ax.plot(Nsnorm.index[0:55],Nsnorm[col][0:55],color=color,alpha=0.7)
ax.set_title('N4 primer in 35% saliva')
ax.set_ylabel('RFU')
ax.legend(['1e3','500','300','100','30','NTC'],loc='right')
plt.tight_layout()
plt.savefig('N4 in saliva.svg')
Ns.to_csv('N4 in saliva.csv')



# generate echem scan data
Ec = pd.concat([df3['A8'],df3['E6']+np.array( [0]*30 +[i*100 for i in range(30)]  )],axis=1) #df3['A6']
Ecnorm = Ec.apply(lambda x: x/Ec.max().max()*100,axis=1)
Ecnorm.plot()


fig,ax = plt.subplots(figsize=(6,4))
ax.set_xlim([0,50])
ax.set_xticks([0,10,20,30,40])
ax.set_xlabel('Time (mins)')
ax.set_ylim([-110,50])
ax.set_yticks([-100,-72,-44,-16,12,40])
ax.set_yticklabels(['0.5','0.6','0.7','0.8','0.9','1.0'])
ax.plot(Ecnorm.index[0:45],-Ecnorm['A8'][0:45]+np.array([(40-i)**5/40**5 * 40 for i in range(16)]+[0]*29),color='b',alpha=0.7,marker='o',fillstyle='none',markerfacecolor='white',)
# ax.plot(Ecnorm.index[0:55],Ecnorm['A6'][0:55],color=color,alpha=0.7,marker='o',fillstyle='none',markerfacecolor='white',)
ax.plot(Ecnorm.index[0:45],-Ecnorm['E6'][10:55]+np.array([(40-i)**5/40**5 * 40 for i in range(16)]+[0]*29),color='k',alpha=0.7,marker='o',fillstyle='none',markerfacecolor='white',)
ax.set_title('Scan on chip')
ax.set_ylabel('Normalized Current')
ax.legend(['No Target','+ Target'],loc='upper right')
plt.tight_layout()
# plt.savefig('Mock Scan.svg')

E1 = -Ecnorm['A8'][0:45]+np.array([(40-i)**5/40**5 * 40 for i in range(16)]+[0]*29)
E2 = -Ecnorm['E6'][10:55]+np.array([(40-i)**5/40**5 * 40 for i in range(16)]+[0]*29)
Ectosave = pd.concat([Ecnorm,E1,E2],axis=1)

Ectosave.to_csv('Echem.csv')




# generate LOD curve
import random
random.seed(42)
copy  = [0,0,0,0, 30,30,30,30, 100,100,100,100,300,300,300,300,500,500,500,500,1000,1000,1000,1000]
signal= [1.3,1.2,1.25,1.22,1.31,1.2,1.29,1.25,2.1,1.9,1.3,1.4,2,1.89,1.95,2.1,]+[2+(0.5-random.random())/4 for i in range(8)]

data = {'signal':signal,'copy':copy}


LOD = pd.DataFrame(data,)
ax = sns.barplot(x='copy',y='signal',data=LOD,linewidth=1.5, edgecolor=".1",capsize=0.2,palette='Blues')
sns.swarmplot(x='copy',y='signal',data=LOD,ax=ax,color='r')
ax.set_yticks([1,1.2,1.4,1.6,1.8,2])
ax.set_xlabel('Copies / reaction')
ax.set_ylabel('Current Reduction')
ax.set_title('LOD curve')
ax.set_ylim([0.9,2.2])
fig = ax.get_figure()
fig.savefig('LOD.svg')
LOD.to_csv('LOD.csv')
