"""
generate analysis from primer design csv file.
"""
from ape import REFape,read_primerset_excel
from align_sequence import REF,BAT
import pandas as pd
from mymodule import revcomp,LazyProperty
import glob
import primer3
from itertools import product
from collections import Counter,OrderedDict
# from primer_para import *
mv_conc=50
dv_conc=8
dntp_conc=1.4


ps = read_primerset_excel()


def iter_primerset_excel():
    "yield:PrimerSetRecord([setname,F3,B3,FIP,BIP,LF,LB,B2c,B1c,F2,F1,gene,B3c,LFc,])"
    ps = read_primerset_excel()
    for p in ps:
        name = p['set']
        pm = p['feature']
        pd = {i.replace(name+'-',''):j for i,_,j in pm}
        F3 = pd.get('F3',None)
        B3 = pd.get('B3',None)
        F2 = pd.get('F2',None)
        F1c = pd.get('F1c',None)
        B2 = pd.get('B2',None)
        B1c = pd.get('B1c',None)
        LF = pd.get('LF',None)
        LB = pd.get('LB',None)
        gene = REFape.name_primer(revcomp(F1c))
        if all([F3,B3,F2,F1c,B2,B1c,LF,LB]):
            yield PrimerSetRecord([name,F3,B3,F1c+F2,
                  B1c+B2,LF,LB,revcomp(B2),B1c,F2,revcomp(F1c),gene,revcomp(B3),revcomp(LF)])


def iter_primerset_html(files):
    "iterate over primerdesign result from website."



files = glob.glob('./LAMP_primer_design_output/*.csv')


def iter_primerset_lamp_design(files,skiprows=0,usecols=[1,2],skipfooter=0,return_df=False):
    """
    read all csv files convert to a single DataFrame
    then iterate over primer sets,give it a name based on locus
    yield: PrimerSetRecord()
    [setname,F3,B3,FIP,BIP,LF,LB,B2c,B1c,F2,F1,gene,B3c,LFc,]
    """
    df = pd.DataFrame(columns=['name','seq'])
    dfs = [df]
    for f in files:
        _df = pd.read_csv(f,skiprows=list(range(skiprows)),usecols=usecols,skipfooter=skipfooter)
        dfs.append(_df)
    df = pd.concat(dfs,axis=0,ignore_index=True)
    if return_df:
        yield df
    else:
        name_counter = Counter()
        for i in range((len(df)//8)):
            F3, F2, F1, B1c, B2c, B3c, LFc, LB = list(df.loc[i*8:(i*8+7),'seq'])
            gene = REFape.name_primer(F1)
            name_counter[gene[0]] += 1
            setname = gene[0] + str(name_counter[gene[0]])
            yield PrimerSetRecord([setname,F3,revcomp(B3c),revcomp(F1)+F2,
                  B1c+revcomp(B2c),revcomp(LFc),LB,B2c,B1c,F2,F1,gene,B3c,LFc])



class PrimerSetRecord(OrderedDict):
    """
    Class to store primer set info.
    data stored in an ordered dict
    [setname,F3,B3,FIP,BIP,LF,LB,B2c,B1c,F2,F1,gene]
    """
    sequence_order = ('F3','B3','FIP','BIP','LF','LB','B2c','B1c','F2','F1')
    fragment_order = ('F3','F2','F1','B3c','B2c','B1c','LFc','LB',)
    primer_order = ('F3','B3','FIP','BIP','LF','LB',)
    def __init__(self,inputs):
        super().__init__(zip(
        ['name','F3','B3','FIP','BIP','LF','LB','B2c','B1c','F2','F1','gene','B3c','LFc'],inputs))

    @property
    def dg(self):
        return 0
    @property
    def tm(self):
        return 0

    def __str__(self):
        return f"PrimerSetRecord {self['name']} in {self['gene']}"

    def __repr__(self):
        return self.table.to_markdown()

    @property
    def table(self):
        "show a table for insepect"
        df=pd.DataFrame.from_dict(self,orient='index',columns=['Value'])
        df.index.name="Key"
        return df

    def iter(self,name):
        r = getattr(self,name+'_order',[])
        for i in r:
            yield i,self[i]
    def iterF(self):
        yield from zip(('F3','F2','F1','B3c','B2c','B1c','LFc','LB'),
        (self['F3'],self['F2'],self['F1'],revcomp(self['B3']),self['B2c'],self['B1c'],revcomp(self['LF']),self['LB']))

    def calc_GC_ratio(self,seq):
        return round((seq.count('G')+seq.count('C'))/len(seq),4) if seq else 0

    def GC_ratio(self,):
        for p,seq in self.iter('primer'):
            self[p+'-GC%']=self.calc_GC_ratio(seq)
        return self

    def Amplicon_pos(self):
        "locate amplicon; the positions are 0 indexed."
        F3 = self['F3']
        B3c = self['B3c']
        s,_=REFape.locate_primer(F3)
        _,e=REFape.locate_primer(B3c)
        self['A_start']=s
        self['A_end']=e
        return self

    def length(self,):
        for p,seq in self.iter('primer'):
            r = len(seq)
            self[p+'-Length']=r
        return self

    def Tm(self):
        for p,seq in self.iter('fragment'):
            tm = primer3.bindings.calcTm(seq,mv_conc,dv_conc,dntp_conc,) if seq else 0
            self[p+'-Tm']=round(tm,3)
        return self

    def Hairpin(self):
        for p,seq in self.iter('primer'):
            r=primer3.bindings.calcHairpin(seq, mv_conc,dv_conc,dntp_conc) if seq else self
            self[p+'-HpdG']=round(r.dg/1000,4)
            self[p+'-HpTm']=round(r.tm,3)
        return self

    def PrimerDimer(self):
        for (p1,s1),(p2,s2) in product(self.iter('primer'),repeat=2):
            r = primer3.bindings.calcHeterodimer(s1,s2,mv_conc,dv_conc,dntp_conc) if s1 and s2 else self
            self[p1+'-'+p2+'-PDdG']=round(r.dg/1000,4)
            self[p1+'-'+p2+'-PDTm']=round(r.tm,3)
        return self

    def Inclusivity(self):
        s = list( i[1] for i in self.iter('fragment') if i[1])
        r = REF.check_inclusivity(*s)
        self['Inclusivity']=round(r,5)
        return self

    @LazyProperty
    def gap_positions(self):
        "return a list of all 12 gap positions, 0 indexed."
        return [j for i in ['F3','F2','F1','B1c','B2c','B3c'] for j in REFape.locate_primer(self[i])]

    def Gaps(self):
        "Gap distances between fragments and amplicon length."
        ps = self.gap_positions
        r = [ i-j for i,j in zip(ps[2::2] , ps[1::2])]
        self.update(zip(('G1','G2','G3','G4','G5'), r ))
        return self

    def LoopHairpin(self):
        "check loop region of amplicon hairpin"
        lfs,_ = REFape.locate_primer(self['F2'])
        lfe,_ = REFape.locate_primer(self['F1'])
        _,lrs = REFape.locate_primer(self['B1c'])
        _,lre = REFape.locate_primer(self['B2c'])
        lf = revcomp(REFape[lfs:lfe])
        lr = REFape[lrs:lre]
        lfr=primer3.bindings.calcHairpin(lf, mv_conc,dv_conc,dntp_conc)
        lrr=primer3.bindings.calcHairpin(lr, mv_conc,dv_conc,dntp_conc)
        self['FloopdG'] = round(lfr.dg/1000,4)
        self['FloopTm'] = round(lfr.tm,3)
        self['RloopdG'] = round(lrr.dg/1000,4)
        self['RloopTm'] = round(lrr.tm,3)
        return self


    def NonTarget(self):
        "Check non target dg"
        return self

    def ExtensionStartGCratio(self,forward=8):
        "the GC ratio at extension start point"
        def f(pos):
            return self.calc_GC_ratio(REFape[pos:pos+forward])
        def r(pos):
            return self.calc_GC_ratio(REFape[pos-forward:pos])

        gs = self.gap_positions
        self.update(zip(
        ('F3extGC%','F2extGC%','F1extGC%','B1extGC%','B2extGC%','B3extGC%'),
        (f(gs[1]),f(gs[3]),f(gs[5]),r(gs[6]),r(gs[8]),r(gs[10]))
        ))
        return self



    def CrossReactivity(self):
        "check cross reactivity with other viruses; need to implement BLAST myself"
        return self

p = list(iter_primerset_excel())[1]
p
d=pd.DataFrame.from_dict(p,orient='index',columns=['Value'],)

d.index.name='new'
a = str(p)
a
p.ExtensionStartGCratio()

(p.Inclusivity()
  .Amplicon_pos()
  .Gaps()
  .GC_ratio()
  .length()
  .Tm()
  .Hairpin()
  .PrimerDimer())

for k,i in p.items():
    print(k,'=',i)
