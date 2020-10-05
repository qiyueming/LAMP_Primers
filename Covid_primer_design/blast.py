from Bio.Blast import NCBIWWW
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from ape import read_primerset_excel
import pandas as pd



result = NCBIWWW.qblast('blastn','nt','TCATCAAACGTTCGGATGCT',nucl_reward=1,
nucl_penalty=-3,megablast='on',expect=1000,hitlist_size=100,gapcosts="5 2",
filter='F',genetic_code=1,entrez_query="(txid1508227 [ORGN])")




# write fasta file of primers.
primerset = read_primerset_excel()
fields = "F1c 	F2	F3	B1c	B2	B3	LF	LB".split()
fasta = []
for primers in primerset:
    setname = primers['set']
    primerseq = primers['feature']
    pt = [ (i,j) for i,_,j in primerseq if i.replace(setname+'-','') in fields]
    for n,s in pt:
        fasta.append(f'>{n}')
        fasta.append(s)

with open('primers.fasta','wt') as f:
    f.write('\n'.join(fasta))
print('\n'.join(fasta))




# parse blast results in csv format for blast of all primers to single_organism:
# the csv file is downloaded from ncbi hit table.
organisms = """MERS
Human Coronavirus 229E
Human Coronavirus HKU1
Human Coronavirus OC43
Human Coronavirus NL63
SARS Coronavirus
Bat SARS-CoV
Adenovirus
Human Parainfluenza 1-4
Influenza A/B
Human Metapneumovirus
Human rhinovirus
Chlamydia pneumoniae
Legionella pneumophila
Mycobacterium tuberculosis
Streptococcus pneumoniae
Streptococcus pyogenes
Bordetella pertussis
Mycoplasma pneumoniae
Pneumocystis jirovecii
Candida albicans
Pseudomonas aeruginosa
Staphylococcus epidermis
Staphylococcus salivarius""".split('\n')

files = [
"./blast_hits/MERS.csv",
"./blast_hits/229E.csv",
"./blast_hits/HKU1.csv",
"./blast_hits/OC43.csv",
"./blast_hits/NL63.csv",
"./blast_hits/a.csv",
"./blast_hits/b.txt",
"./blast_hits/adeno.csv",
"./blast_hits/hpi.csv",
"./blast_hits/flu.csv",
"./blast_hits/meta.csv",
"./blast_hits/rino.csv",
"./blast_hits/1.csv",
"./blast_hits/2.csv",
"./blast_hits/3.csv",
"./blast_hits/4.csv",
"./blast_hits/5.csv",
"./blast_hits/6.csv",
"./blast_hits/7.csv",
"./blast_hits/8.csv",
"./blast_hits/9.csv",
"./blast_hits/10.csv",
"./blast_hits/11.csv",
"./blast_hits/12.csv",
]

primerlength = {i:len(j) for primer in primerset for i,_,j in primer['feature']}

fields = "F1c	F2	F3	B1c	B2	B3	LF	LB".split()
primers = "N2	NA	N3	N4	N5	N6	OB	O1	O2	O3	O4	E1".split()

homodf = pd.DataFrame(columns=primers)
for orga, file in zip(organisms,files):
    if file.endswith('.csv'):
        df = pd.read_csv(file,names=['Primer','match','gap'],usecols=[0,3,4])
    elif file.endswith('.txt'):
        with open(file,'rt') as f:
            data = f.read().split('\n')
        df = pd.DataFrame(columns=['Primer','match','gap'])
        key = ""
        match = 0
        gap = 0
        for line in data:
            if line.startswith('Query='):
                key = line.split('= ')[1]
                match,gap = 0,0
            elif line.startswith(' Identities'):
                i,g = line.split(', ')
                if match == 0:
                    match = match or int(i.split('= ')[1].split('/')[0].strip())
                    gap = gap or int(g.split('= ')[1].split('/')[0].strip())
                    df.loc[key,:]=[key,match,gap]
    allprimerhomo=[]
    for primer in primers:
        homos = []
        for field in fields:
            name = f"{primer}-{field}"
            if name in primerlength:
                matchdf = df.loc[df['Primer']==name,:]
                if len(matchdf) == 0:
                    homos.append(0)
                else:
                    first  = matchdf.index[0]
                    homo = matchdf.loc[first,'match'] - matchdf.loc[first,'gap']
                    homos.append(homo/primerlength[name])
        avghomo = sum(homos)/len(homos)
        allprimerhomo.append(avghomo)
    homodf.loc[orga,:] = allprimerhomo

homodf

homodf.to_csv('Primer Homology.csv')



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
