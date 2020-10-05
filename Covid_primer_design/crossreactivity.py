"""
To check cross reactivity 
some virus are checked against my saved viral genome sequences.
this check is done using primer3 HeteroDimerAnneal function. 

cross reactivity to other eukaryotes and virus are from blast files. 
use print_fasta to print fast for check in NCBI Blast 
then download the hit table csv and parse it. 
"""

from ape import APE
# from itertools import combinations_with_replacement
from mymodule import revcomp
import primer3
from ape import read_primerset_excel
import pandas as pd
from primer_analysis import PrimerSetRecord,PrimerSetRecordList
from resources import CROSS_GENE_FILE,CROSS_GENE_NAME_LONG,CROSS_GENE_NAME
import glob

def print_ascii_structure(lines,printresult=True):
    if not isinstance(lines,list):
        lines = lines.ascii_structure_lines
    l = lines[0][4:]
    l1 = lines[1][4:]
    start = 4
    end = 4
    ll = len(l)
    for k in range(ll):
        if l[k] != " " or l1[k]!=" ":
            start += k
            break
    for j in range(len(l)):
        if l[ll-j-1]!='-':
            end += ll-j
            break
    if printresult:
        print('\n'.join(i[start:end] for i in lines))
    else:
        return [i[start:end] for i in lines]


def primer_homology(primer,sequence,mv_conc=50,dv_conc=8,dntp_conc=1.4,dna_conc=200):
    "mv_conc,dv_conc,"
    # def combine(s1,s2):
    #     re=[]
    #     for i,j in zip(s1,s2):
    #         re.append(i if i!=' ' else j)
    #     return ''.join(re)
    if (primer in sequence) or (revcomp(primer) in sequence):
        maxtm =primer3.bindings.calcTm(primer,mv_conc=mv_conc,dv_conc=dv_conc,dntp_conc=dntp_conc,dna_conc=dna_conc)
        maxhomo = 1
        align = [' '*len(primer),primer,revcomp(primer)," "*len(primer)]
    else:
        length = int(len(sequence)/9000+1)
        primerlength = len(primer)
        maxtm=-100
        maxhomo = 0
        align = None
        for i in range(length):
            totest = sequence[i*9000:(i+1)*9000+100]
            r = primer3.bindings.calcHeterodimer(primer,totest,output_structure=True,mv_conc=mv_conc,dv_conc=dv_conc,dntp_conc=dntp_conc,dna_conc=dna_conc)
            tm = r.tm
            if tm>maxtm:
                maxtm = tm
                s = print_ascii_structure(r,printresult=False)
                # s = [ combine(s[0],s[1]),combine(s[2],s[3])]
                c = 0
                for j in s[1]:
                    if j in 'ATCG':
                        c+=1
                maxhomo = c/primerlength
                align = s
    return maxtm, maxhomo,align


def print_fasta(tocheck):
    "generate the fastta string to check in blast"
    for p in tocheck:
        name = p['name']
        for f,s in p.iter('fragment'):
            print(f'>{name}-{f}')
            print(s)



# prepare tocheck sequences. 
toorder = 'LAMP_primer_design_output/toorder.csv'
TO = PrimerSetRecordList(toorder)
tocheck = TO[['C_N7','N12','N13']]


# printout fasta format string. 
print_fasta(tocheck)



# check against saved viral genomes
virusgenome = {i:APE(j).sequence for i,j in zip(CROSS_GENE_NAME_LONG,CROSS_GENE_FILE)}
df = pd.DataFrame(columns = ['Organism'])
df['Organism'] = CROSS_GENE_NAME_LONG 
for p in tocheck:
    crossreact = []
    name = p['name']
    
    for gene in CROSS_GENE_NAME_LONG:
        virusseq = virusgenome[gene]
        homos = []
        for _,s in p.iter('fragment'):
            homo = primer_homology(revcomp(s),virusseq)
            homos.append(homo[1])
        crossreact.append("{:.2%}".format( sum(homos)/len(homos) ))
    df[name] = crossreact



# parse csv files 
# csv files need to be the same order as OTHER_ORGANISM
OHTER_ORGANISM = [
 'Chlamydia pneumoniae',
 'Legionella pneumophila',
 'Mycobacterium tuberculosis',
 'Streptococcus pneumoniae',
 'Streptococcus pyogenes',
 'Bordetella pertussis',
 'Mycoplasma pneumoniae',
 'Pneumocystis jirovecii',
 'Candida albicans',
 'Pseudomonas aeruginosa',
 'Staphylococcus epidermidis',
 'Streptococcus salivarius',]

files = [
'/Users/hui/Downloads/1.csv',
'/Users/hui/Downloads/2.csv',
'/Users/hui/Downloads/3.csv',
'/Users/hui/Downloads/4.csv',
'/Users/hui/Downloads/5.csv',
'/Users/hui/Downloads/6.csv',
'/Users/hui/Downloads/7.csv',
'/Users/hui/Downloads/8.csv',
'/Users/hui/Downloads/9.csv',
'/Users/hui/Downloads/10.csv',
'/Users/hui/Downloads/11.csv',
'/Users/hui/Downloads/12.csv',
]

primers = [i['name'] for i in tocheck]
homodf = pd.DataFrame(columns=['Organism'] + primers)

fields = ('F3','F2','F1','B3c','B2c','B1c','LFc','LB',)
for orga,file in zip(OHTER_ORGANISM,files):
    _df = pd.read_csv(file,names=['Primer','match','gap'],usecols=[0,3,4])
    allprimerhomo={'Organism':[orga]}
    for p in tocheck:
        primer = p['name']
        homos = []
        for fn,s in p.iter('fragment'):
            name = f"{primer}-{fn}"
            matchdf = _df.loc[_df['Primer']==name,:]
            if len(matchdf) == 0:
                homos.append(0)
            else:
                first  = matchdf.index[0]
                homo = matchdf.loc[first,'match'] - matchdf.loc[first,'gap']
                homos.append(homo/len(s))
        avghomo = sum(homos)/len(homos)
        allprimerhomo[primer]=["{:.2%}".format(avghomo)]
        # allprimerhomo.append((primer,avghomo))
    homodf = homodf.append( pd.DataFrame( allprimerhomo ) ,ignore_index=True)

# combine both together:

df = df.append(homodf,ignore_index=True)


df.to_csv('Cross Reactivity 6-13.csv')




# plot heat map
import pandas as pd
df = pd.read_csv('./Primer Homology.csv',index_col=0,keep_default_na=False)


import seaborn as sns
import matplotlib.pyplot as plt
sns.set()

fig,ax = plt.subplots(figsize=(8,10))
ax = sns.heatmap(df,vmin=0.5,vmax=1,cmap='jet',linewidths=0.1,ax=ax,cbar_kws={"shrink": 0.5})
ax.set_title('Primer Homology',fontsize=18)
plt.tight_layout()


plt.savefig('Homology_full.svg')
