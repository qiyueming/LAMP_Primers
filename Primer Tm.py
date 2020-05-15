from ape import APE
import primer3
from itertools import product
from ape import read_primerset_excel
import pandas as pd

primerset = read_primerset_excel()
primerset[3]['feature']

primer3.bindings.calcTm('CTTGAACTGTTGCGACTACGTGAT',mv_conc=50,dv_conc=8,dntp_conc=1.4)

# calculate Tms with primer 3 and write to file.
ref = APE('/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/cov2_ref/HCoV2 NC045512.2.ape')

"Caution: The Dimer_dG calculated here doesn't make sense. "

result = {}
for p in primerset:
    name = p['set']
    tms = []
    for i,j,k in p['feature']:
        tms.append(primer3.bindings.calcTm(k,mv_conc=50,dv_conc=8,dntp_conc=1.4))
    dg = 0
    for i,j in product([i[2] for i in p['feature']],repeat=2):
        dg = min(dg,primer3.bindings.calcHeterodimer(i,j,mv_conc=50,dv_conc=8,dntp_conc=1.4).dg)
    start = ref.locate_primer(p['feature'][0][2])[0]
    end = ref.locate_primer(p['feature'][1][2])[0]
    result[name] = {'tm':[tms[i] for i in [3,2,0,5,4,1,6,7] if i<len(tms)],'dg':dg,'pos':f"{start}-{end}"}

df = pd.DataFrame(columns="F1c 	F2	F3	B1c	B2	B3	LF	LB	Dimer_dG locus".split())
for k,i in result.items():
    df.loc[k,:] = i['tm'] + [0]* (8-len(i['tm'])) + [i['dg'],i['pos']]

df.to_csv('Primer Tms.csv')
