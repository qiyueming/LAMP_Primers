"""
from generate inclusivity data of a primer set. 
first download multialign file form GISAID;
or download from gene bank and use mafft to align. 
mafft --thread -1 --nomemsave  gisaid_cov_sequences.fasta > msa_0613.fasta
then use Reference class to read the align fasta file. 
"""

from ape import APE
# from itertools import combinations_with_replacement
from mymodule import revcomp
import primer3
from itertools import product
from mymodule import ViennaRNA
from mymodule.RNAstructure import RNA
from Levenshtein import distance
import random
from primer_analysis import PrimerSetRecord,PrimerSetRecordList
import pandas as pd
from collections import Counter
from align_sequence import Reference


# pull in the relevant align file. 
gismsa = "/Users/hui/Downloads/msa_0613/msa_0613.fasta"

REF = Reference(gismsa)


# pull in the relevant primers 
toorder = 'LAMP_primer_design_output/toorder.csv'
TO = PrimerSetRecordList(toorder)

tocheck = TO[['C_N7','N12','N13']]


# for check Inclusitivity within the REF:
df = pd.DataFrame(columns=['Primer Set','Seqs','Mismatch','100%Homology','F3','F2','F1','B3c','B2c','B1c','LFc','LB'])

for p in tocheck:
    row = {'Primer Set':p['name'],'Seqs': len(REF.aln)}
    
    for n,s in p.iter('fragment'):
        row.update({n: REF.check_inclusivity(s) * 100 })
    inclusivity =  REF.check_inclusivity(*list( i[1] for i in p.iter('fragment') if i[1]))
    row.update({'100%Homology': inclusivity * 100})
    row.update({'Mismatch': int( round( (1-inclusivity) * len(REF.aln) ,0) )})
    df=df.append(row,ignore_index=True)



def determine_mutation(s,align,reverse=True):
    """
    s is the primer sequence, align is the fragment in align. 
    reverse = True if this is for LFc or B1c ... 
    so that will check homology distance to 5' end. 
    """
    if reverse:
        s = s[::-1]
        align = align[::-1]
    f=[]
    t=[]
    pos=0
    append=False
    for k,(i,j) in enumerate(zip(s,align)):
        if i==j:
            append=False 
            continue 
        else:
            pos = k
            if append:
                f[-1]+=i 
                t[-1]+=j 
            else:
                f.append(i)
                t.append(j)
            append=True 
    return ','.join( f"{a}->{b}" for a,b in zip(f,t) ), len(s) - 1 - pos
  

# for check detailed mismatch 
mismatches = {}
for p in tocheck:
    mismatch = Counter()
    positions = []
    for n,s in p.iter('fragment'):
        _,pos = REF.find_seq(s)
        positions.append((n,_,pos))
    for virusSeq in REF.aln:
        mutation = ()
        for n,s,pos in positions:
            align = virusSeq[pos[0]:pos[1]]
            # check mutation if align here is different and align only contain AGCT-, 
            # because N or R,Y are indication of sequencing error thus not considerred. 
            if align!=s and ( set(align)  <= {'A','G','C','T','-'} ) :    
                mut,pos = determine_mutation(s,align,n.endswith('c'))
                mutation += (n,mut,pos)
        if mutation:
            mismatch[mutation]+=1
    mismatches[p['name']] = mismatch
    
# append the mismatches result to df
mostfrequent = []
singlemut = []
seqerror = []
tPrimOver5 = []

for p in tocheck:
    mis = mismatches[p['name']]
    mm = mis.most_common(5)
    total = sum(mis.values())
    
    mostfreq = []
    for i,c in mm:
        mostfreq.append( "{:.1%} {} in {}, {}n.t. to 3'".format(c/total, i[1],i[0],i[2])) 
    mostfrequent.append('; '.join(mostfreq))
    
    sc = 0 # single mutation count
    for k,value in mis.items():
        if len(k)==3:
            sc += value 
    singlemut.append("{:.2%}".format(sc/total))
    
    sr = 0  # possible sequencing error
    for k,value in mis.items():
        for i in k:
            if isinstance(i,str) and len(i.split('->'))>1 and set(i.split('->')[1]) == set('-'):
                sr += value
                break
    seqerror.append( "{:.2%}".format(sr/total) )
    
    tDo5 = 0  # three primer distance over 5
    for k,value in mis.items():
        if k[2] >=5:
            tDo5 += value 
    tPrimOver5.append("{:.2%}".format(tDo5/total))
    

df['Most Frequent'] = mostfrequent
df['Single Mutation'] = singlemut
df['Sequencing Error'] = seqerror
df["3' Distance Over 5"] = tPrimOver5


df.to_csv('inclusivit analysis 6-13.csv')