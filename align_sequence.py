from ape import APE,REFape
from mymodule import Alignment
import numpy as np
import pandas as pd
def read(file):
    with open(file,'rt') as f:
        return f.read().split('\n')

def lines_to_dict(lines):
    align={}
    for line in lines:
        if line.startswith('>'):
            key = line.split('|')[1]
            align[key]=[]
        else:
            align[key].append(line)
    return {i:''.join(j) for i,j in align.items()}


def bp_align(A,B,ref='NC_045512.2'):
    "Aligne two dictionary of aligned sequence according to their common ref sequence."
    def dict_append(s,d,D):
        if s == '-':
            for k in d:
                d[k].append(s)
        else:
            for k in d:
                d[k].append(D[k][s])
    a = A[ref]
    b = B[ref]
    A_align = {i:[] for i in A}
    B_align = {i:[] for i in B}
    la = len(a) -1
    lb = len(b) -1
    ia = 0
    ib = 0
    ir = 0
    while True:
        ir+=1
        if ir > 31000:
            break
        if ia > la and ib > lb:
            break
        if ia > la or ib>lb:
            if ia > la:
                dict_append('-',A_align,A)
                dict_append(ib,B_align,B)
            if ib > lb:
                dict_append('-',A_align,A)
                dict_append(ia,A_align,A)
            ia+=1
            ib+=1
        elif a[ia] == b[ib]:
            dict_append(ia,A_align,A)
            dict_append(ib,B_align,B)
            ia+=1
            ib+=1
        else:
            if a[ia]=='-':
                dict_append(ia,A_align,A)
                dict_append('-',B_align,B)
                ia+=1
            else:
                dict_append('-',A_align,A)
                dict_append(ib,B_align,B)
                ib+=1
    A_align = {i:''.join(j) for i,j in A_align.items()}
    B_align = {i:''.join(j) for i,j in B_align.items()}
    A_align.update(B_align)
    return A_align

def map_log_rate(x,low=0.999,high=1):
    res = np.log10(high)-np.log10(x) / (np.log10(high) - np.log10(low))
    return min(max(0,res),1)

def reduce_density(array,ratio=2):
    result = []
    for i,j in enumerate(array):
        if i%ratio == 0:
            result.append(j)
        elif j!=0:
            if result[-1]==0:
                result.append(j)
    return result


class AlignmentF(Alignment):
    "Augmented for deal with primer designing."
    def freq_calc(self,list_of_seq,count=False):
        """
        return frequency list with order of AGCT-;
        for position with N, doesn't count N or other letters towards frequency.
        """
        list_of_seq = [list_of_seq] if isinstance(
            list_of_seq, str) else list_of_seq
        seq_c = len(list_of_seq)
        if seq_c == 0:
            result = []
        else:
            count = (np.array(count) /
                             sum(count) if count else np.full(seq_c, 1.0/seq_c))
            count = count.reshape(-1, 1)
            A_ = np.sum(np.multiply(np.array(
                [[k == 'A' for k in list(i)] for i in list_of_seq], float), count), axis=0)
            G_ = np.sum(np.multiply(np.array(
                [[k == 'G' for k in list(i)] for i in list_of_seq], float), count), axis=0)
            C_ = np.sum(np.multiply(np.array(
                [[k == 'C' for k in list(i)] for i in list_of_seq], float), count), axis=0)
            T_ = np.sum(np.multiply(np.array(
                [[k == 'T' for k in list(i)] for i in list_of_seq], float), count), axis=0)
            gap = np.sum(np.multiply(np.array(
                [[k == '-' for k in list(i)] for i in list_of_seq], float), count), axis=0)
            N  = np.sum(np.multiply(np.array(
                [[k not in 'AGCT-' for k in list(i)] for i in list_of_seq], float), count), axis=0)
            result = [np.array(i)*seq_c/(seq_c - n*seq_c ) for n,*i in zip(N,A_, G_, C_, T_, gap)]
        return result


class Reference():
    """
    use sequence and mutation information from saved csv file.
    """
    def __init__(self,csvfile=None,alnfile=None,ref='NC_045512.2'):
        if csvfile:
            self.df = pd.read_csv(csvfile,index_col=0)
            self.ref = ''.join(self.df.nt)
        if alnfile:
            all_align = lines_to_dict(read(alnfile))
            refseq = all_align.pop(ref)
            self.aln = [refseq] + list(all_align.values())
            self.aln_count = len(self.aln)
            if not csvfile:
                self.ref = refseq
    def __getitem__(self,slice):
        return [i[slice] for i in self.aln]

    def label_gene(self,ape):
        self.genes={}
        for feature in ape.features:
            start = feature['start']
            name = feature['tag']
            end = feature['end']
            startseq = ape.sequence[start-1:start+19]
            endseq = ape.sequence[end-20:end]
            try:
                startindex = self.find_seq(startseq)[1][0]
                endindex = self.find_seq(endseq)[1][1]
                self.genes[name] = (startindex,endindex)
            except:
                print(f'feature {name} Not found')


    def check_homology(self,seq,end=0,direction='F'):
        "check homology of sequence to other sequences in the align. return a list of homology"
        f,(start,_end) = self.find_seq(seq)
        seqs = self[start:_end]
        result = []
        l = len(f)
        for i in seqs[1:]:
            result.append(sum([i==j for i,j in zip(i,f)]) / l )
        if end:
            if direction == 'F':
                _slice = slice(-end,None)
            else:
                _slice = slice(0,end)
            for i in seqs[1:]:
                result.append(sum([i==j for i,j in zip(i[_slice],f[_slice])]) / end )
        return result

    def find_seq(self,seq):
        "find a sequence based on ref sequence,return the sequence and position."
        try:
            pos = self.ref.index(seq)
        except:
            try:
                alphacount = self.ref.replace('-','').index(seq)
            except:
                print(f'<{seq}> not found')
                return 0
            pos = 0
            while alphacount:
                if self.ref[pos].isalpha(): alphacount -= 1
                pos +=1
        span = pos+len(seq)
        totest = self.ref[pos:span]
        while totest.replace('-','') != seq:
            totest += self.ref[span]
            span+=1
        return totest,(pos,span)

    def check_inclusivity(self,*seq):
        if len(seq) == 1:
            seq = seq[0]
            totest,(pos,span) = self.find_seq(seq)
            return self._check_inclusivity(totest,pos)
        else:
            # check multi sequcne match simutaneously.
            testcase = [self.find_seq(i) for i in seq]
            return sum([ all( [  i[p:q]==s or (not (set(i[p:q]) <= {'A','G','C','T','-'})) for s,(p,q) in testcase]) for i in self.aln]) /self.aln_count

    def _check_inclusivity(self,seq,pos):
        "check sequence that already found in ref"
        sl = len(seq)
        match = [ i[pos:pos+sl]==seq or (not (set(i[pos:pos+sl]) <= {'A','G','C','T','-'})) for i in self.aln]
        return sum(match)/self.aln_count


    def primer_iter_gene(self,gene,length=(19,22),inclusivity=0.99,filter_func=None,step=1):
        start,end = self.genes[gene]
        return self.primer_iter(start,end,length,inclusivity,filter_func,step)

    def primer_iterOLD(self,start,end,length=(19,22),inclusivity=0.99,filter_func=None,step=1):
        """
        iterate primers in a region, use filter_func to give additional criteria
        homology is the thereshold, to limit minimum of homology for each nt in the primer.
        both start and end are included.
        primer length will include length[0] and leng[1]
        yield the primer sequence, and position (start,end) of the primer. start end is the sequence index.
        """
        for i in range(start,end+1,step):
            for j in range(length[0],length[1]+1):
                nt = self.ref[i:i+j]
                if self._check_inclusivity(nt,i) < inclusivity:
                    break
                seq = nt.replace('-','')
                index = i+j
                fail = False
                while len(seq) < j:
                    extrant = self.ref[index]
                    if self._check_inclusivity(nt,index) < inclusivity:
                        fail = True
                        break
                    if extrant != '-': seq +=extrant
                    index +=1
                if fail: break
                if filter_func == None or filter_func(seq):
                        yield seq, (i,index)

    def primer_iter(self,start,end,length=(19,22),inclusivity=0.99,filter_func=None,step=1,iter_length=False):
        return PrimerIter(self,start,end,length,inclusivity,filter_func,step,iter_length)


class PrimerIter():
    def __init__(self,REF,start,end,length=(19,22),inclusivity=0.99,filter_func=None,step=1,iter_length=False):
        """
        if not iter_length, will move to next
        position once find a certain length at current position
        """
        self.start = start
        self.end = end
        self.step = step
        self.lengthlow,self.lengthhigh = length
        self.inclusivity = inclusivity
        self.filter_func = filter_func
        self.REF=REF
        self.posi = start if step>0 else end
        self.lengthi = self.lengthlow
        self.iter_length=iter_length

    def getnextNT(self):
        "return the next nucleotide"
        if self.lengthi > self.lengthhigh:
            self.lengthi = self.lengthlow
            self.posi += self.step
        if self.posi <= self.end and self.posi>= self.start:
            nt = self.REF.ref[self.posi:self.posi+self.lengthi]
            self.lengthi += 1
            return nt,self.posi,self.lengthi - 1
        raise StopIteration

    def next_pos(self,forward=None):
        self.posi += forward or self.step
        self.lengthi = self.lengthlow

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            nt,i,j = self.getnextNT()
            if self.REF._check_inclusivity(nt,i) < self.inclusivity:
                self.next_pos()
                continue
            seq = nt.replace('-','')
            index = i+j
            fail = False
            while len(seq) < j:
                extrant = self.REF.ref[index]
                if self.REF._check_inclusivity(nt,index) < self.inclusivity:
                    fail = True
                    break
                if extrant != '-': seq +=extrant
                index +=1
            if fail:
                self.next_pos()
                continue
            if self.filter_func == None or self.filter_func(seq):
                if not self.iter_length:
                    self.next_pos()
                return seq,(i,index)

alnfile = './viral_genome/all_align_0512.aln'

REF = Reference(alnfile=alnfile)
BAT = Reference(alnfile='./viral_genome/CoV2+Bat.aln')
REF.label_gene(REFape)

#
# REF = Reference(alnfile=alnfile)
#
# a = REF.primer_iter(1000,1002)
# for i in a :
#     print(i)
#     a.next_pos(2)
# r=list(a)
# len(r)
#
#
# b = REF.primer_iter(1000,1002,step=-1)
# r1=list(b)
# r1

if __name__ == '__main__':

    lm = APE('/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/cov2_ref/HCoV2 MN908947.3.ape')
    ref = APE('/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/cov2_ref/HCoV2 NC045512.2.ape')

    p100 = APE('P100Converve 200seq.ape')

    alignfiles = [
    '/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/world_align/A.aln',
    '/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/world_align/B.aln',
    '/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/world_align/C.aln',
    '/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/world_align/D.aln',
    '/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/world_align/E.aln',
    '/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/world_align/F.aln',
    ]

    A = lines_to_dict(read(alignfiles[0]))
    B = lines_to_dict(read(alignfiles[1]))
    C = lines_to_dict(read(alignfiles[2]))
    D = lines_to_dict(read(alignfiles[3]))
    E = lines_to_dict(read(alignfiles[4]))
    F = lines_to_dict(read(alignfiles[5]))

    # Align A BCDEF together and write to file.
    res = A
    for t in [B,C,D,E,F]:
        res = bp_align(res,t)

    with open('all_align_0512.aln','wt') as f:
        for k,i in res.items():
            f.write(f">gb|{k}\n{i}\n")

    all_align = lines_to_dict(read('/Users/hui/Desktop/WFH papers/COVID-19/Virus Genes/viral_genome/all_align_0512.aln',))
    align = AlignmentF(sequence=list(all_align.values()))

    # save align sequence and frequency as csv
    df = pd.DataFrame(columns=['nt','A','G','C','T','-'])
    seq = align.rep_seq(count=True)
    for i in range(len(align.freq)):
        row = [seq[i]] + align.freq[i].tolist()
        df.loc[i,:] = row
    df.to_csv('all_align_0512.csv')


    # plot mutation map
    import matplotlib.pyplot as plt

    freq = [i.max() for i in align.freq]
    seq = align.rep_seq(count=True)
    low = 0.989 # percentage for it to be draw as 0
    high = 0.99 # percentage for it to be draw as 1
    density_reduce_fold = 2
    cmap='YlGnBu'
    figrows =len(ref.features)


    fig,axes = plt.subplots(nrows = figrows,figsize=(16,int(1+0.5*figrows)))
    fig.subplots_adjust(top=0.9, bottom=0.01, left=0.065, right=0.99)

    gradient = [map_log_rate(i,low,high) for i in freq]
    gradient = reduce_density(gradient,int(len(gradient) / 200) + 1 )
    gradient = np.vstack((gradient,gradient))
    ax = axes[0]
    ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap('YlGnBu'))
    pos = list(ax.get_position().bounds)
    x_text = pos[0] - 0.01
    y_text = pos[1] + pos[3]/2.
    fig.text(x_text, y_text, f'SARS_CoV2\n{len(freq)}nt', va='center', ha='right', fontsize=10,)

    for ax,feature in zip(axes[1:],ref.features[0:-2]+ref.features[-1:]):
        start = feature['start']
        name = feature['tag']
        end = feature['end']
        startseq = ref.sequence[start-1:start+9]
        endseq = ref.sequence[end-10:end]
        try:
            startindex = seq.index(startseq)
            endindex = seq.index(endseq) + 10
        except:
            print(f'feature {name} Not found')
            continue
        gradient = [map_log_rate(i,low,high) for i in freq[startindex:endindex]]
        gradient = reduce_density(gradient,int(len(gradient) / 800) + 1 )
        gradient = np.vstack((gradient,gradient))
        ax.imshow(gradient,aspect='auto',cmap= plt.get_cmap(cmap))
        pos = list(ax.get_position().bounds)
        x_text = pos[0] - 0.01
        y_text = pos[1] + pos[3]/2
        if len(name)==1: name = name+' Gene'
        fig.text(x_text, y_text,name+'\n'+f'{end-start}n.t.', va='center', ha='right', fontsize=10)
        # ax.set_axis_off()

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['bottom'].set_color('#dddddd')
        ax.spines['top'].set_color('#dddddd')
        ax.spines['right'].set_color('#dddddd')
        ax.spines['left'].set_color('#dddddd')

    fig.suptitle('Mutations in SARS-CoV2',fontsize=16)
    plt.savefig('SARS_CoV2 mutation 0.989-0.99.svg')
