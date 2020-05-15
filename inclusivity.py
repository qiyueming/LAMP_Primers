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
from collections import Counter



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



primerset = read_primerset_excel()

[len(i['feature']) for i in primerset]

with open('/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/All_CoV2.fasta','rt') as f:
    data = f.read()
data = data.split('\n')
sequences={}

for line in data:
    if line.startswith('>'):
        key = line[1:].split(' ')[0]
        sequences[key]=[]
        continue
    sequences[key].append(line)
sequences = {i:''.join(j) for i,j in sequences.items()}
slist = list(sequences.values())


# for check Inclusitivity:
df = pd.DataFrame(columns=['Primer Set','Inclusitivity',"Count",'F1c','F2','F3','B1c','B2','B3','LB','LF'])
fullslist=[i for i in slist if len(i)>20000]
for primers in primerset:
    setname = primers['set']
    primerseq = primers['feature']
    count = 0
    frag = Counter()
    for virus in fullslist:
        mismatch=False
        for fn,_,fn_seq in primerseq:
            ttprimer=None
            if fn in [setname+'-F3',setname+'-F2',setname+'-B1c',setname+'-LB']:
                ttprimer = fn_seq
            elif fn in [setname+'-B3',setname+'-B2',setname+'-F1c',setname+'-LF']:
                ttprimer = revcomp(fn_seq)
            if ttprimer and (ttprimer not in virus):
                mismatch=True
                name = fn.split('-')[1]
                frag[name]+=1
        if mismatch: count +=1
    update = {'Primer Set':setname,'Inclusitivity':"{:.2%}".format(1-count/len(fullslist)),'Count':count}
    update.update(dict(zip(['F1c','F2','F3','B1c','B2','B3','LB','LF'],[0]*8)))
    update.update(frag)
    df=df.append(update,ignore_index=True)
df
df.to_csv('Inclusitivity.csv')
df = pd.read_csv('Inclusitivity.csv',index_col=0)

df.loc[:,['F1c','F2','F3','B1c','B2','B3','LB','LF']] = 1-df.loc[:,['F1c','F2','F3','B1c','B2','B3','LB','LF']]/2315

df.loc[:,['F1c','F2','F3','B1c','B2','B3','LB','LF']]=df.loc[:,['F1c','F2','F3','B1c','B2','B3','LB','LF']].applymap(lambda x: "{:.2%}".format(x))
df.to_csv('Inclusitivity update .csv')


# for check inclusivity reactivity with examine non matching homomology and Tm.
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
# for check inclusivity reactivity with examine non matching homomology and Tm.
columns=['Primer Set','Primer','Homology','Tm','Mismatch','GenID']
df = pd.DataFrame(columns=columns)
fullslist={i:j for i,j in sequences.items() if len(j)>20000}
for primers in primerset:
    setname = primers['set']
    primerseq = primers['feature']
    homos = []
    minhomo = 1
    mintm = 100
    for ID,virus in fullslist.items():
        if len(virus)<20000:
            continue
        for fn,_,fn_seq in primerseq:
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
df.to_csv('detailed_inclusivity.csv')


detail = pd.read_csv('detailed_inclusivity.csv',index_col=0,keep_default_na=False)

detail.head()
detail.loc[0,:]['Primer Set']

from collections import Counter


result = {} #dict.fromkeys(set(detail['Primer Set']))



for index in detail.index:
    row = detail.loc[index,:]
    primer = row['Primer Set']
    if primer not in result:
        result[primer] = Counter()
    mismatch = row['Mismatch']
    if not (set(mismatch) <= set('-\nATCG')):
        result[primer]['N'] +=1
    else:
        result[primer][mismatch] +=1
r=result['N4']

for k,i in result.items():
    print(k)
    mc = i.most_common(3)
    total = sum(i.values())
    tp = ["{:.2%} {}".format(b/total,a.replace('\n','->')) for a,b in mc]
    print(', '.join(tp))




primer3.bindings.calcHeterodimerTm('CTTGAACTGTTGCGACTACGTGAT','GTTCCTCATCACGTAGTCGCAACAGTTCAAGAAATTCAACTC',mv_conc=50,dv_conc=8,dntp_conc=1.4,dna_conc=200)

# N4 primer detailed analysis
primer3.bindings.calcHeterodimerTm('CTTGAACTGTTGCGACTACGTGAT',mv_conc=50,dv_conc=8,dntp_conc=1.4,dna_conc=200)
N4primerTm = {n:primer3.bindings.calcHeterodimerTm(s,revcomp(s),mv_conc=50,dv_conc=8,dntp_conc=1.4,dna_conc=200) for n,_,s in primerset[3]['feature']}
N4primerTm
N4result = detail.loc[detail['Primer Set']=='N4',:]
N4result.head()
genes = set(N4result['GenID'])


gene_dtm = {}

for gene in genes:
    tdf = N4result.loc[N4result['GenID']==gene,:]
    gene_dtm[gene] = []
    causedbymis = False
    for mis in tdf['Mismatch']:
        if not (set(mis) <= set('-\nATCG')):
            gene_dtm[gene].append('Caused By Sequencing')
            causedbymis = True
            break
    if causedbymis:continue
    for p,t,mis in zip(tdf['Primer'],tdf['Tm'],tdf['Mismatch']):
        dtm = N4primerTm[p]-t
        miscount = (len(mis) - 1) /2
        gene_dtm[gene].append((dtm,miscount))
len(gene_dtm)
gene_dtm


result_df = pd.DataFrame(columns=['GenID','dTm','MisMatch'])
for gene,mis in gene_dtm.items():
    if isinstance(mis[0],str):
        result_df.loc[gene,:] = [gene,'Sequencing','']
        continue
    # for dtm, m in mis:
    #     result_df.loc[gene,:] = [gene,dtm,m]
    maxdtm,maxmis = max(mis, key=lambda x: x[0])
    result_df.loc[gene,:] = [gene,maxdtm,maxmis]
result_df
result_df.to_csv('N4 detailed dTm.csv')

miscounter = Counter()
dTmcounter = Counter()
for i in result_df.index:
    row = result_df.loc[i,:]
    dtm = row['dTm']
    mis = row['MisMatch']
    dTmcounter[dtm] +=1
    miscounter[mis] +=1
for k,i in miscounter.items():
    print(k,i)

for k,i in :
    print(k,i)
