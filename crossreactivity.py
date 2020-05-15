from ape import APE
# from itertools import combinations_with_replacement
from mymodule import revcomp
import primer3
from itertools import product
from mymodule import ViennaRNA
from mymodule.RNAstructure import RNA
from Levenshtein import distance
import random
from ape import read_primerset_excel
import pandas as pd


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

def determine_mutation(a,b):
    "b is the sequence in template"
    f=[]
    t=[]
    for i,j in zip(a,b):
        if i== ' ':
            continue
        if i=='-':
            f.append(i)
            t.append(j)
        else:
            _=revcomp(i)
            if _ != j:
                f.append(revcomp(i))
                t.append(j)
    return "".join(f),''.join(t)

files = glob.glob('/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/cross_reactivity/*')
genes = [
"Human Coronavirus NL63 (NC005831)",
"Human Coronavirus OC43 (KX344031.1)",
"SARS Coronavirus (NC_004718)",
"Human Parainfluenza 4a (NC_021928)",
"Human Coronavirus 229E (MN369046)",
"Human Coronavirus HKU1 (MH940245)",
"Human Coronavirus OC43 (NC_006213)",
"Human Coronavirus HKU1 (NC_006577.2)",
"Human Metapneumovirus (NC_039199)",
"Bat SARS-CoV 2015 (MG772933.1)",
"Influenza A (NC_026435)",
"MERS (NC_019843)",
"Human Coronavirus 229E (NC_002645)",
"Human rhinovirus (NC_038311)",
"Human Parainfluenza 4b (MN306032)",
"Adenovirus (AC_000017)",
"Human Coronavirus HKU1 (MK167038)",
"Bat SARS-CoV 2017 (MG772934.1)",
"Human Parainfluenza 1 (AF457102)",
"Influenza B (NC_002205)",
]

genome = {i:APE(j).sequence for i,j in zip(genes,files)}
genomelist = list(genome.values())
len(genomelist)
[set(i) for i in genomelist]



primerset = read_primerset_excel()
primerset[0]

columns=['Primer Set','Primer','Homology','Tm','Mismatch','Gene']
df = pd.DataFrame(columns=columns)
for primers in primerset[0:1]:
    setname = primers['set']
    primerseq = primers['feature']
    homos = []
    minhomo = 1
    mintm = 100
    for ID,virus in genome.items():

        for fn,_,fn_seq in primerseq:
            print('Testing {} - {}'.format(ID,fn))
            ttprimer=None
            if fn in [setname+'-F3',setname+'-F2',setname+'-B1c',setname+'-LB']:
                ttprimer = fn_seq
            elif fn in [setname+'-B3',setname+'-B2',setname+'-F1c',setname+'-LF']:
                ttprimer = revcomp(fn_seq)
            if ttprimer and (ttprimer not in virus):
                tm,homo,align=primer_homology(revcomp(ttprimer),virus)
                minhomo = min(minhomo, homo)
                mintm = min(tm,mintm)
                mismatch = '\n'.join(determine_mutation(align[0],align[3]))
                df=df.append(dict(zip(columns,[setname,fn,homo,tm,mismatch,ID])),ignore_index=True)

df.to_csv('detailed_crossreactivity.csv')





# generate cross reactivity heatmap
df = pd.read_csv('detailed_crossreactivity.csv',index_col=0,keep_default_na=False)
newdf = pd.DataFrame(columns='N2	NA	N3	N4	N5	N6	OB	O1	O2	O3	O4	E1'.split())

genes = set(df['Gene'])


for gene in genes:
    homo=[]
    for col in newdf.columns:
        data = df.loc[(df['Primer Set']==col) & (df.Gene==gene),['Homology'] ]
        value = (data.sum().item() + 8- len(data))/8
        homo.append( value )


    newdf.loc[gene,:]=homo

genes = [
"Human Coronavirus NL63 (NC005831)",
"Human Coronavirus OC43 (KX344031.1)",
"SARS Coronavirus (NC_004718)",
"Human Parainfluenza 4a (NC_021928)",
"Human Coronavirus 229E (MN369046)",
"Human Coronavirus HKU1 (MH940245)",
"Human Coronavirus OC43 (NC_006213)",
"Human Coronavirus HKU1 (NC_006577.2)",
"Human Metapneumovirus (NC_039199)",
"Bat SARS-CoV 2015 (MG772933.1)",
"Influenza A (NC_026435)",
"MERS (NC_019843)",
"Human Coronavirus 229E (NC_002645)",
"Human rhinovirus (NC_038311)",
"Human Parainfluenza 4b (MN306032)",
"Adenovirus (AC_000017)",
"Human Coronavirus HKU1 (MK167038)",
"Bat SARS-CoV 2017 (MG772934.1)",
"Human Parainfluenza 1 (AF457102)",
"Influenza B (NC_002205)",
]
new_index = []
for i in newdf.index:
    for j in genes:
        id = j[j.index('(')+1:j.index(')')]
        if id in i:
            new_index.append(j)
            break
newdf.index=new_index
newdf.to_csv('cross_reactivity_avg.csv')






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
