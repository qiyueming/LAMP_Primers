from ape import REFape
from mymodule import revcomp
import primer3
# from mymodule import ViennaRNA
# from mymodule.RNAstructure import RNA
# import pandas as pd
from align_sequence import REF,BAT
from primer_para import *
from datetime import datetime
from Levenshtein import distance
# from old_primer import print_ascii_structure
from mymodule import ft,ProgressBar
from ape import read_primerset_excel
from collections import Counter



# result format:
# name:N-F1; Tm; GCratio; Inclusitivity; CrossReact;
# possibles: extension start point GCratio; Whether have ss loop structure at 65deg;
# 3' end of primer mis match in other region of target. Or total Heterodimer Tm to the entire genome.
# primer dimer formation
# check list for final product
# 0. cross react with bat CoV. check 3' end homology instead of total homology.
# 1. primer dimer √
# 2. primer hairpin structure √
# 3. 3' end of primer extension start GC content within 6 n.t. between 2 and 4
# 4. LAMP product loop region use loophairpin dG < -2 or -3.
# 5. heterodimer in other position of viral genome if 3' end stable.
# 6. check repeats of special sequence in a primer i.e. 4nt repeat, 3n.t. repeat.


def draw_hairpin(seq):
    r=primer3.bindings.calcHairpin(seq, mv_conc,dv_conc,dntp_conc,output_structure=True )
    print(r.dg)
    print(r.tm)
    print(r.ascii_structure)



def GC_ratio(seq):
    return (seq.count('G')+seq.count('C'))/len(seq)

def time_filter(F,number=1000):
    f = F()
    ft(f,args=[],number=number)

def hamming(s1,s2,k=1):
    return sum(s1[i:i+k]!=s2[i:i+k] for i in range(len(s1)-k+1))

# individual sequence filter that filter out sequence based on sequence info.
def TmFilter(Tm=(60,65), return_value=False):
    "return a filter have Tm threshold"
    def wrap(seq="GCCAAAAGGCTTCTACGCA"):
        tm = primer3.bindings.calcTm(seq,mv_conc,dv_conc,dntp_conc,)
        if return_value:return tm
        return  tm >= Tm[0] and tm <= Tm[1]
    return wrap

def GCfilter(ratio):
    "return a GC ratio filter"
    def wrap(seq):
        r = GC_ratio(seq)
        return (r >=ratio[0] and r<=ratio[1])
    return wrap

def ESfilter(count = E3):
    "return a end GC count filter"
    def wrap(seq):
        return seq[-5:].count('G')+seq[-5:].count('C') >= count
    return wrap

def ESCfilter(N=6,M=4):
    "check for 3' end self complementary, and if the last M nt are the same"
    def wrap(seq):
        return revcomp(seq[-N//2:])!=seq[-N:-N//2] and len(set(seq[-M:]))>1
    return wrap

def Hammingfilter(threshold = 0.3,k=2,return_value=False):
    lr = len(REFape.sequence)
    def wrap(seq):
        l=len(seq)
        if seq in REFape.sequence:
            seqindex = REFape.sequence.index(seq)
            minham = min( (hamming(seq, REFape.sequence[i:i+l],k ),REFape.sequence[i:i+l]) for i in range(seqindex - l) )
            minham2 = min( (hamming(seq, REFape.sequence[i:i+l],k ),REFape.sequence[i:i+l]) for i in range(seqindex +l , lr-l) )
            mindis = min(minham,minham2)
        else:
            mindis =  min( (hamming(seq, REFape.sequence[i:i+l],k ),REFape.sequence[i:i+l]) for i in range( lr-l) )
        if return_value:
            return mindis
        return mindis[0]/l >= threshold
    return wrap

def Levfilter(threshold = 0.3,k=2,return_value=False):
    lr = len(REFape.sequence)
    def wrap(seq):
        l=len(seq)
        if seq in REFape.sequence:
            seqindex = REFape.sequence.index(seq)
            md = min( (distance(seq, REFape.sequence[i:i+l] ), REFape.sequence[i:i+l]) for i in range(seqindex - l) )
            md2 = min( (distance(seq, REFape.sequence[i:i+l]),REFape.sequence[i:i+l]) for i in range(seqindex +l , lr-l) )
            mindis = min(md,md2)
        else:
             mindis = min( (distance(seq, REFape.sequence[i:i+l] ),REFape.sequence[i:i+l]) for i in range(lr - l) )
        if return_value:
            return mindis
        return mindis[0]/l >= threshold
    return wrap

def HoMofilter(homology=0.8,end=4,direction='F',return_value=False):
    "also check homology at end of the sequence"
    "This is too stringent."
    def wrap(seq="GCCAAAAGGCTTCTACGCA",):
        if return_value:
            return max(BAT.check_homology(seq,end,direction))
        return max(BAT.check_homology(seq,end,direction)) <= homology
    return wrap

def CombFilter(*filters):
    def wrap(seq):
        return all(i(seq) for i in filters)
    return wrap

def RCwrapper(func):
    "make the filter function take reversecomp sequence."
    def wrap(seq):
        return func(revcomp(seq))
    return wrap

def Hairpinfilter(dG=-4,return_value=False):
    "filter out primer with high dG or 3' end stem."
    def wrap(seq="GCCAAAAGGCTTCTACGCA"):
        r=primer3.bindings.calcHairpin(seq, mv_conc,dv_conc,dntp_conc,output_structure=True )
        if return_value: return r
        if not r.ascii_structure: return True
        if r.dg/1000>dG and r.ascii_structure_lines[0][-2:]=='--':
            return True
        return False
    return wrap

def PrimerDimerfilter(Tm=30,return_value=False):
    """
    check a new sequence with a pool of already checked sequence.
    if the Tm of the structure is higher than threshold, check if there is
    extendable 3'.
    if have dimer, return False. if not, return True.
    """
    def checkhairpin(r):
        "return True if result have extendable hairpin"
        if r.tm<Tm: return False
        lines = r.ascii_structure_lines
        for i,j in zip(lines[0][::-1],lines[1][::-1]):
            if i.isalpha(): break
            if j.isalpha(): return True
        for i,j in zip(lines[2],lines[3]):
            if i.isalpha(): break
            if j.isalpha(): return True
        return False
    def wrap(seq="GCCAAAAGGCTTCTACGCA",pool=[]):
        cH = primer3.bindings.calcHeterodimer
        sd = cH(seq,seq,mv_conc,dv_conc,dntp_conc,output_structure=True)
        if return_value: return sd
        if checkhairpin(sd): return False
        for i in pool:
            if checkhairpin(cH(seq,i,mv_conc,dv_conc,dntp_conc,output_structure=True)):
                return False
        return True
    return wrap

def PrimerComplexityfilter(trimer=2,tetramer=1,pentamer=0):
    "test if primer consist of tandom repeats"
    def wrap(seq):
        if sum(seq.count(i) for i in ['AAAAA','TTTTT','GGGGG','CCCCC']) > pentamer:
            return False
        if sum(seq.count(i) for i in ['AAAA','TTTT','GGGG','CCCC']) > tetramer:
            return False
        if sum(seq.count(i) for i in ['AAA','TTT','GGG','CCC']) > trimer:
            return False
        return True
    return wrap

class HomologyCounter(dict):
    """
    To keep track of how many remaning homology fragment needed.
    """
    def __init__(self,order):
        self.order = order
        super()
    def __call__(self, sequence,name):
        h = BAT.check_homology(sequence)
        self[name] = max(h) <= BatHomology
        for n in self.order[ self.order.index(name) + 1:]:
            self[n] = 0
        return self

    def check_homology(self,x):
        return  max(BAT.check_homology( REF.ref[x:x+25].replace('-','') ))

    def check_possibility(self,starts,ends):
        if sum(self.values()) >= BatHomologyCount: return True
        if isinstance(starts,int):
            starts = [starts]
            ends = [ends]
        return any(self.check_homology(i)<= BatHomology  for start,end in zip(starts,ends) for i in range(start,end,2))

# method to save primer set.
class PrimerSetHandler:
    def __init__(self,path,prefix,batchcount=100):
        self.count = 0
        self.data = [('set','name','seq','time')]

        # if path exist
        self.path = path.replace('.csv','')+datetime.now().strftime('%Y%m%d_%H:%M:%S')+'.csv'
        self.prefix = prefix
        self.batchcount = batchcount

        # read parameters and save.
        with open('primer_para.py','rt') as f:
            lines = f.read()
        para = [i for i in lines.split('\n') if i and i[0].isalpha()]
        with open(self.path,'wt') as f:
            f.write('\n'.join(para))
            f.write(f"\n====> Started On {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} <====")

    def __call__(self,data):
        self.count +=1
        time = f"{datetime.now().strftime('%H:%M:%S')}"
        self.data.extend((f"{self.prefix}{self.count}" ,n,s,time)
                for n,s in zip(['F3','F2','F1','B1c','B2c','B3c','LFc','LB'] , data))
        if self.count % self.batchcount == 0:
            self.write()

    def write(self):
        with open(self.path,'at') as f:
            f.write('\n')
            f.write('\n'.join(','.join(primerset) for primerset in self.data))
        self.data = []

    def done(self):
        'write a line to indicate its done.'
        with open(self.path,'at') as f:
            f.write(f"\n====> END On {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} <====")

def main(target=None,step=1,MAX_primerset=1000,savepath='./LAMP_primer.csv'):
    """
    start from F3 to B3.
    """
    pDimerfilter = PrimerDimerfilter(PrimerDimerTm)
    hPfilter = Hairpinfilter(HairpindG)
    LoopHPfilter = Hairpinfilter(LoopHairpindG)
    hOmofilterF = HoMofilter(BatHomology,BatHomologyEnd,'F')
    hOmofilterR = HoMofilter(BatHomology,BatHomologyEnd,'R')

    F3filter = CombFilter(TmFilter(P3Tm),GCfilter(GCratio),ESfilter(E3),ESCfilter(N=ESC),Hairpinfilter(HairpindG),PrimerComplexityfilter())
    F2filter = CombFilter(TmFilter(P2Tm),GCfilter(GCratio),ESfilter(E3),ESCfilter(N=ESC),PrimerComplexityfilter())
    F1filter = CombFilter(TmFilter(P1Tm),GCfilter(GCratio),ESfilter(E3),ESCfilter(N=ESC),PrimerComplexityfilter())
    B1cfilter = RCwrapper(F1filter)
    B2cfilter = RCwrapper(F2filter)
    B3cfilter = RCwrapper(F3filter)
    LBfilter = CombFilter(TmFilter(LPTm),GCfilter(GCratio),ESfilter(E3),ESCfilter(N=ESC),Hairpinfilter(HairpindG))
    LFcfilter = RCwrapper(LBfilter)

    if isinstance(target,str):
        F3_start, F3_end = REF.genes[target]
    else:
        F3_start, F3_end = target
        target = 'W'

    stop = False

    primerset_counter = 0

    SavePrimerSet=PrimerSetHandler(savepath,prefix=target[0],batchcount=5)

    progress = ProgressBar(limits=(F3_start,F3_end))
    F3iter = REF.primer_iter(F3_start,F3_end,P3L,PInclu,F3filter,step)
    for F3,(sF3,eF3) in F3iter:
        progress(sF3)
        if stop: break
        F2_start = eF3 + g1[0]
        F2_end = eF3 + g1[1]
        if not pDimerfilter(F3,[]): continue

        F2iter = REF.primer_iter(F2_start,F2_end,P2L,PInclu,F2filter,)
        for F2, (sF2,eF2) in F2iter:
            if stop:break
            F1_start = eF2 + g2[0]
            F1_end = eF2 + g2[1]
            if not pDimerfilter(F2,[F3]): continue

            F1iter = REF.primer_iter(F1_start,F1_end,P1L,PInclu,F1filter)
            for F1, (sF1,eF1) in F1iter:
                if stop:break
                B1c_start = eF1 + g3[0]
                B1c_end = eF1 + g3[1]
                FIP = revcomp(F1)+F2
                maxloopfc =  revcomp(REF.ref[sF2:sF1].replace('-','')[0:59])
                if not LoopHPfilter(maxloopfc): break # break this F1 if strong loop formed.
                if not hPfilter(FIP): continue
                if not pDimerfilter(FIP,[F3]): continue


                B1citer = REF.primer_iter(B1c_start,B1c_end,P1L,PInclu,B1cfilter)
                for B1c, (sB1c,eB1c) in B1citer:
                    if stop:break
                    B2c_start = eB1c + g4[0]
                    B2c_end = eB1c + g4[1]

                    B2citer = REF.primer_iter(B2c_start,B2c_end,P2L,PInclu,B2cfilter)
                    for B2c, (sB2c,eB2c) in B2citer:
                        if stop:break
                        if eB2c - sF2 > g6[1] or eB2c - sF2 < g6[0]: # if there is no room for B2 to satisfy g6.
                            break
                        B3c_start = eB2c + g5[0]
                        B3c_end = eB2c + g5[1]
                        BIP = B1c + revcomp(B2c)
                        maxloopr =  REF.ref[eB1c:eB2c].replace('-','')[0:59]
                        if not LoopHPfilter(maxloopr): break
                        if not hPfilter(BIP): continue
                        if not pDimerfilter(BIP,[F3,FIP]): continue

                        B3citer = REF.primer_iter(B3c_start,B3c_end,P3L,PInclu,B3cfilter)
                        for B3c,(sB3c,eB3c) in B3citer:
                            if stop: break
                            B3 = revcomp(B3c)
                            if not pDimerfilter(B3,[F3,FIP,BIP]): continue

                            LFfound = False
                            for LFc, (sLFc,eLFc) in REF.primer_iter(eF2,sF1,LPL,LoopInclu,LFcfilter):
                                LF = revcomp(LFc)
                                if not pDimerfilter(LF,[F3,B3,BIP,FIP]): continue
                                LFfound=True
                                break
                            if not LFfound: continue

                            LBfound = False
                            for LB, (sLB,eLB) in REF.primer_iter(eB1c,sB2c,LPL,LoopInclu,LBfilter):
                                if not pDimerfilter(LB,[F3,B3,BIP,FIP,LF]): continue
                                LBfound = True
                                break
                            if not LBfound: continue

                            primerset = (F3,F2,F1,B1c,B2c,B3c,LFc,LB)
                            # if sum([hOmofilterF(i) for i in primerset[0:3]] + [hOmofilterR(i) for i in primerset[3:6]])>=3: # at least 3 of the 6 are satisfied.
                            SavePrimerSet(primerset)
                            primerset_counter += 1
                            if primerset_counter > MAX_primerset:
                                stop = True
    SavePrimerSet.write()
    SavePrimerSet.done()

def main_limit_count(target=None,span=None,MAX_primerset=1000,savepath='./LAMP_primer.csv'):
    """
    start with F2 -> F1 -> B1 -> B2 -> B3 -> F3 .
    Strong check Homology to bat CoV at F2, B2 and F3, B3.
    check primer complexity.
    faster adjust on positions.
    """
    if span:
        A_start, A_end = span
        target = target or 'W'
    else:
        A_start, A_end = REF.genes[target]


    pDimerfilter = PrimerDimerfilter(PrimerDimerTm)
    hPfilter = Hairpinfilter(HairpindG)
    LoopHPfilter = Hairpinfilter(LoopHairpindG)
    hOmofilterF = HoMofilter(BatHomology,BatHomologyEnd,'F')
    hOmofilterR = HoMofilter(BatHomology,BatHomologyEnd,'R')

    F3filter = CombFilter(TmFilter(P3Tm),GCfilter(GCratio),ESfilter(E3),
                          ESCfilter(N=ESC),Hairpinfilter(HairpindG),PrimerComplexityfilter())
    F2filter = CombFilter(TmFilter(P2Tm),GCfilter(GCratio),ESfilter(E3),
                          ESCfilter(N=ESC),PrimerComplexityfilter())
    F1filter = CombFilter(TmFilter(P1Tm),GCfilter(GCratio),ESfilter(E3),ESCfilter(N=ESC),PrimerComplexityfilter())
    B1cfilter = RCwrapper(F1filter)
    B2cfilter = RCwrapper(F2filter)
    B3cfilter = RCwrapper(F3filter)

    LBfilter = CombFilter(TmFilter(LPTm),GCfilter(GCratio),ESfilter(E3),ESCfilter(N=ESC),Hairpinfilter(HairpindG))
    LFcfilter = RCwrapper(LBfilter)




    stop = False

    primerset_counter = 0

    SavePrimerSet=PrimerSetHandler(savepath,prefix=target[0],batchcount=10)

    F2iter = REF.primer_iter(A_start+g1[1]+P3L[1],A_end-g6[1]-g5[1]-P3L[1],P2L,PInclu,F2filter,)

    progress = ProgressBar(limits=(A_start+g1[1]+P3L[1],A_end-g6[1]-g5[1]-P3L[1]))
    for F2, (sF2,eF2) in F2iter:
        progress(sF2)
        F2_Count = 0
        if stop:break
        F1_start = eF2 + g2[0]
        F1_end = eF2 + g2[1]
        if not pDimerfilter(F2,[]): continue

        F1iter = REF.primer_iter(F1_start,F1_end,P1L,PInclu,F1filter)
        for F1, (sF1,eF1) in F1iter:
            F1_Count = 0
            if stop:break

            if F2_Count > F2_CountThreshold: # if toomany with samfe F2, break out of F1 loop.
                F2iter.next_pos(10) # if found many in this position, move F2 forward 10 n.t.
                break

            B1c_start = eF1 + g3[0]
            B1c_end = eF1 + g3[1]
            FIP = revcomp(F1)+F2
            maxloopfc =  revcomp(REF.ref[sF2:sF1].replace('-','')[0:59])
            if not LoopHPfilter(maxloopfc):
                F2iter.next_pos(AdjustStep) # move F2 forward a step.
                break # break this F1 if strong loop formed.
            if not hPfilter(FIP): # if FIP have hairpin move F1 forward 3 n.t.
                F1iter.next_pos(AdjustStep)
                continue

            if not pDimerfilter(FIP,[]):
                continue

            B1citer = REF.primer_iter(B1c_start,B1c_end,P1L,PInclu,B1cfilter)
            for B1c, (sB1c,eB1c) in B1citer:
                B1c_Count = 0
                if stop:break

                if F1_Count > F1_CountThreshold:
                    F1iter.next_pos(5)
                    break

                B2c_start = eB1c + g4[0]
                B2c_end = eB1c + g4[1]

                B2citer = REF.primer_iter(B2c_start,B2c_end,P2L,PInclu,B2cfilter)
                for B2c, (sB2c,eB2c) in B2citer:
                    B2c_Count = 0
                    if stop:break
                    if B1c_Count > B1c_CountThreshold:
                        B1citer.next_pos(5)
                        break

                    if eB2c - sF2 > g6[1] or eB2c - sF2 < g6[0]: # if there is no room for B2 to satisfy g6.
                        B1citer.next_pos(1000) # end B1 iteration
                        break
                    B3c_start = eB2c + g5[0]
                    B3c_end = eB2c + g5[1]
                    BIP = B1c + revcomp(B2c)
                    maxloopr =  REF.ref[eB1c:eB2c].replace('-','')[0:59]
                    if not LoopHPfilter(maxloopr):
                        B1citer.next_pos(AdjustStep) # move B1c forward a step
                        break # break this B2 if strong loop formed.
                    if not hPfilter(BIP): # if BIP have hairpin, move B2 forward.
                        B2citer.next_pos(AdjustStep)
                        continue
                    if not pDimerfilter(BIP,[FIP]):
                        continue


                    B3citer = REF.primer_iter(B3c_start,B3c_end,P3L,PInclu,B3cfilter)
                    for B3c,(sB3c,eB3c) in B3citer:
                        B3c_Count = 0
                        if stop: break
                        if B2c_Count > B2c_CountThreshold:
                            B2citer.next_pos(5)
                            break

                        B3 = revcomp(B3c)
                        if not pDimerfilter(B3,[FIP,BIP]): continue
                        F3_start = sF2-g1[1]-P3L[1]
                        F3_end = sF2 - g1[0]- P3L[0]

                        F3iter = REF.primer_iter(F3_start,F3_end,P3L,PInclu,F3filter,step=-1)
                        for F3,(sF3,eF3) in F3iter:
                            if stop: break
                            if B3c_Count >B3c_CountThreshold:
                                B3citer.next_pos(5)
                                break
                            if eF3>sF2: continue # in case overlap happened.
                            if not pDimerfilter(F3,[B3,FIP,BIP]): continue


                            LFfound = False
                            for LFc, (sLFc,eLFc) in REF.primer_iter(eF2,sF1,LPL,LoopInclu,LFcfilter):
                                LF = revcomp(LFc)
                                if not pDimerfilter(LF,[F3,B3,BIP,FIP]): continue
                                LFfound=True
                                break
                            if not LFfound: continue

                            LBfound = False
                            for LB, (sLB,eLB) in REF.primer_iter(eB1c,sB2c,LPL,LoopInclu,LBfilter):
                                if not pDimerfilter(LB,[F3,B3,BIP,FIP,LF]): continue
                                LBfound = True
                                break
                            if not LBfound: continue

                            primerset = (F3,F2,F1,B1c,B2c,B3c,LFc,LB)
                            # save if 3 of the core primers have homology specificity.
                            if sum([hOmofilterF(i) for i in primerset[0:3]] + [hOmofilterR(i) for i in primerset[3:6]])>=2:
                                SavePrimerSet(primerset)
                                primerset_counter += 1
                                F2_Count+=1
                                F1_Count+=1
                                B1c_Count+=1
                                B2c_Count+=1
                                B3c_Count+=1

                            # move F2 and everything forward.
                            if primerset_counter > MAX_primerset:
                                stop = True
    SavePrimerSet.write()
    SavePrimerSet.done()
    progress.end_bar()
    print(f'Runing {target} - {span} Finished.')

def main_Counter(target=None,span=None,MAX_primerset=1000,savepath='./LAMP_primer.csv'):
    """
    start with F2 -> F1 -> B1 -> B2 -> B3 -> F3 .
    Strong check Homology to bat CoV at F2, B2 and F3, B3.
    check primer complexity.
    faster adjust on positions.
    """
    if span:
        A_start, A_end = span
        target = target or 'W'
    else:
        A_start, A_end = REF.genes[target]

    pDimerfilter = PrimerDimerfilter(PrimerDimerTm)
    hPfilter = Hairpinfilter(HairpindG)
    LoopHPfilter = Hairpinfilter(LoopHairpindG)
    hOmofilterF = HoMofilter(BatHomology,BatHomologyEnd,'F')
    hOmofilterR = HoMofilter(BatHomology,BatHomologyEnd,'R')

    F3filter = CombFilter(TmFilter(P3Tm),GCfilter(GCratio),ESfilter(E3),
                          ESCfilter(N=ESC),Hairpinfilter(HairpindG),PrimerComplexityfilter())
    F2filter = CombFilter(TmFilter(P2Tm),GCfilter(GCratio),ESfilter(E3),
                          ESCfilter(N=ESC),PrimerComplexityfilter())
    F1filter = CombFilter(TmFilter(P1Tm),GCfilter(GCratio),ESfilter(E3),ESCfilter(N=ESC),PrimerComplexityfilter())
    B1cfilter = RCwrapper(F1filter)
    B2cfilter = RCwrapper(F2filter)
    B3cfilter = RCwrapper(F3filter)

    LBfilter = CombFilter(TmFilter(LPTm),GCfilter(GCratio),ESfilter(E3),ESCfilter(N=ESC),Hairpinfilter(HairpindG))
    LFcfilter = RCwrapper(LBfilter)

    homologyCounter = HomologyCounter(order=('F2','F3','F1','B1c','B2c','B3c'))

    Counters = [Counter() for i in range(5)]
    F1Counter,F3Counter,B1cCounter,B2cCounter,B3cCounter = Counters

    stop = False

    primerset_counter = 0

    SavePrimerSet=PrimerSetHandler(savepath,prefix=target[0],batchcount=10)

    print(f'Started {target} P({A_start}-{A_end})...\n')
    progress = ProgressBar(limits=(A_start,A_end),prefix=f'P({A_start}-{A_end})')

    F2iter = REF.primer_iter(A_start,A_end,P2L,PInclu,F2filter,)
    for F2, (sF2,eF2) in F2iter:
        progress(sF2)
        F2_Count = 0
        if stop:break
        if not homologyCounter(F2,'F2').check_possibility(sF2-g1[1]-P3L[1],sF2+g6[1]+g5[1]+P3L[1]): break
        F1_start = eF2 + g2[0]
        F1_end = eF2 + g2[1]
        if not pDimerfilter(F2,[]): continue

        F1iter = REF.primer_iter(F1_start,F1_end,P1L,PInclu,F1filter)

        for F1, (sF1,eF1) in F1iter:
            F1_Count = 0
            if stop:break
            if not homologyCounter(F1,'F1').check_possibility((sF2-g1[1]-P3L[1],eF1), (sF2, eF1+g3[1]+P1L[1]+g4[1]+P2L[1]+g5[1]+P3L[1])): break
            if F2_Count > F2_CountThreshold: # if toomany with samfe F2, break out of F1 loop.
                F2iter.next_pos(10) # if found many in this position, move F2 forward 10 n.t.
                break

            B1c_start = eF1 + g3[0]
            B1c_end = eF1 + g3[1]
            FIP = revcomp(F1)+F2
            maxloopfc =  revcomp(REF.ref[sF2:sF1].replace('-','')[0:59])
            if not LoopHPfilter(maxloopfc):
                F2iter.next_pos(AdjustStep) # move F2 forward a step.
                break # break this F1 if strong loop formed.
            if F1Counter[F1] >= SAME_Fragment_Threshold:continue
            if not hPfilter(FIP): # if FIP have hairpin move F1 forward 3 n.t.
                F1iter.next_pos(AdjustStep)
                continue
            if not pDimerfilter(FIP,[]):
                continue

            B1citer = REF.primer_iter(B1c_start,B1c_end,P1L,PInclu,B1cfilter)
            for B1c, (sB1c,eB1c) in B1citer:
                B1c_Count = 0
                if stop:break
                if not homologyCounter(B1c,'B1c').check_possibility((sF2-g1[1]-P3L[1],eB1c), (sF2, eB1c+g4[1]+P2L[1]+g5[1]+P3L[1])): break
                if F1_Count > F1_CountThreshold:
                    F1iter.next_pos(5)
                    break
                if B1cCounter[B1c] >= SAME_Fragment_Threshold:continue
                B2c_start = eB1c + g4[0]
                B2c_end = eB1c + g4[1]

                B2citer = REF.primer_iter(B2c_start,B2c_end,P2L,PInclu,B2cfilter)
                for B2c, (sB2c,eB2c) in B2citer:
                    B2c_Count = 0
                    if stop:break
                    if not homologyCounter(B2c,'B2c').check_possibility((sF2-g1[1]-P3L[1],eB2c), (sF2, eB2c+g5[1]+P3L[1])): break
                    if B1c_Count > B1c_CountThreshold:
                        B1citer.next_pos(5)
                        break
                    if eB2c - sF2 > g6[1] or eB2c - sF2 < g6[0]: # if there is no room for B2 to satisfy g6.
                        B1citer.next_pos(1000) # end B1 iteration
                        break
                    B3c_start = eB2c + g5[0]
                    B3c_end = eB2c + g5[1]
                    BIP = B1c + revcomp(B2c)
                    maxloopr =  REF.ref[eB1c:eB2c].replace('-','')[0:59]
                    if not LoopHPfilter(maxloopr):
                        B1citer.next_pos(AdjustStep) # move B1c forward a step
                        break # break this B2 if strong loop formed.
                    if B2cCounter[B2c] >= SAME_Fragment_Threshold: continue
                    if not hPfilter(BIP): # if BIP have hairpin, move B2 forward.
                        B2citer.next_pos(AdjustStep)
                        continue
                    if not pDimerfilter(BIP,[FIP]):
                        continue


                    B3citer = REF.primer_iter(B3c_start,B3c_end,P3L,PInclu,B3cfilter)
                    for B3c,(sB3c,eB3c) in B3citer:
                        B3c_Count = 0
                        if stop: break
                        if not homologyCounter(B3c,'B3c').check_possibility(sF2-g1[1]-P3L[1],sF2): break
                        if B2c_Count > B2c_CountThreshold:
                            B2citer.next_pos(5)
                            break
                        if B3cCounter[B3c] >= SAME_Fragment_Threshold:continue
                        B3 = revcomp(B3c)
                        if not pDimerfilter(B3,[FIP,BIP]): continue
                        F3_start = sF2-g1[1]-P3L[1]
                        F3_end = sF2 - g1[0]- P3L[0]

                        F3iter = REF.primer_iter(F3_start,F3_end,P3L,PInclu,F3filter,step=-1)
                        for F3,(sF3,eF3) in F3iter:
                            if stop: break
                            if B3c_Count >B3c_CountThreshold:
                                B3citer.next_pos(5)
                                break
                            if F3Counter[F3] >= SAME_Fragment_Threshold:continue
                            if eF3>sF2: continue # in case overlap happened.
                            if not pDimerfilter(F3,[B3,FIP,BIP]): continue

                            LFfound = False
                            for LFc, (sLFc,eLFc) in REF.primer_iter(eF2,sF1,LPL,LoopInclu,LFcfilter):
                                LF = revcomp(LFc)
                                if not pDimerfilter(LF,[F3,B3,BIP,FIP]): continue
                                LFfound=True
                                break
                            if not LFfound: continue

                            LBfound = False
                            for LB, (sLB,eLB) in REF.primer_iter(eB1c,sB2c,LPL,LoopInclu,LBfilter):
                                if not pDimerfilter(LB,[F3,B3,BIP,FIP,LF]): continue
                                LBfound = True
                                break
                            if not LBfound: continue

                            primerset = (F3,F2,F1,B1c,B2c,B3c,LFc,LB)
                            # save if 3 of the core primers have homology specificity.
                            if sum([hOmofilterF(i) for i in primerset[0:3]] + [hOmofilterR(i) for i in primerset[3:6]])>=BatHomologyCount:
                                SavePrimerSet(primerset)
                                primerset_counter += 1
                                F2_Count+=1
                                F1_Count+=1
                                B1c_Count+=1
                                B2c_Count+=1
                                B3c_Count+=1
                                F1Counter[F1] +=1
                                F3Counter[F3] +=1
                                B1cCounter[B1c] +=1
                                B2cCounter[B2c] +=1
                                B3cCounter[B3c] +=1
                            # move F2 and everything forward.
                            if primerset_counter > MAX_primerset:
                                stop = True
    SavePrimerSet.write()
    SavePrimerSet.done()
    progress.end_bar()
    print(f'Runing {target} ({A_start}-{A_end}) Finished.\n')
    return 0
