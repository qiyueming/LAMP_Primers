from mymodule import Alignment
import numpy as np
import re

def read(file):
    with open(file,'rt') as f:
        return f.read().split('\n')

def lines_to_dict(lines):
    align={}
    for line in lines:
        if line.startswith('>'):
            key = line
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
    build a reference genome with aligned sequences.
    """
    def __init__(self,alnfile=None,ref=None):
        all_align = lines_to_dict(read(alnfile))
        if ref:
            for k,i in all_align.items():
                if i.replace('-','') == ref.sequence:
                    break
        refseq = all_align.pop(k)
        self.aln = [refseq] + list(all_align.values())
        self.aln_count = len(self.aln)

        self.ref = refseq
        self.genes={}

    def __getitem__(self,slice):
        return [i[slice] for i in self.aln]

    def to_PE_format(self,start,end,inclusivity=0.99,save=True):
        "output a sequence and homology reprensentation for PrimerExplorer."
        seq = self.ref[start:end]
        consensus = []
        for i in range(start,end):
            inclu = self._check_inclusivity(self.ref[i],i)
            if inclu >= inclusivity:
                consensus.append('*')
            else:
                consensus.append('-')
        if save:
            if not isinstance(save,str):
                save = f"PE_Target_{start}-{end}"
            #  clean up - in sequences for primer explorer
            newseq,newcon = [],[]
            for i,j in zip(seq,consensus):
                if i!='-':
                    newseq.append(i)
                    newcon.append(j)
                else:
                    if j =='-':
                        newseq.append('A')
                        newcon.append(j)
                    else:
                        continue


            with open(save,'wt') as f:
                f.write(f'sequence={"".join(newseq)}\n')
                f.write(f'consensus={"".join(newcon)}\n')
                f.write('\n'.join(['F3_5pos=-1','F3_3pos=-1','F2_5pos=-1','F2_3pos=-1',
                'F1c_5pos=-1','F1c_3pos=-1','B3_5pos=-1','B3_3pos=-1','B2_5pos=-1',
                'B2_3pos=-1','B1c_5pos=-1','B1c_3pos=-1','target_range_type=0',
                'target_range_from=','target_range_to=']))
        else:
            return seq,''.join(consensus)


    def label_gene(self,ape):
        "label genes on self.aln from a ape file."
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
            span = pos+len(seq)
            totest = self.ref[pos:span]
            return totest,(pos,span)
        except:
            pattern = re.compile('-*'.join(list(seq)))
            match = pattern.search(self.ref)
            if match:
                totest = match.group()
                return totest, match.span()
            else:
                raise ValueError(f'{seq} not found.')

    def check_inclusivity(self,*seq):
        if len(seq) == 1:
            seq = seq[0]
            totest,(pos,span) = self.find_seq(seq)
            return self._check_inclusivity(totest,pos)
        else:
            # check multi sequcne match simutaneously.
            testcase = [self.find_seq(i) for i in seq]
            return sum([ all([i[p:q]==s or (not (set(i[p:q]) <= {'A','G','C','T','-'})) for s,(p,q) in testcase]) for i in self.aln]) /self.aln_count

    def _check_inclusivity(self,seq,pos):
        "check sequence that already found in ref"
        sl = len(seq)
        match = [ i[pos:pos+sl]==seq or (not (set(i[pos:pos+sl]) <= {'A','G','C','T','-'})) for i in self.aln]
        return sum(match)/self.aln_count

    def primer_iter_gene(self,gene,length=(19,22),inclusivity=0.99,filter_func=None,step=1):

        start,end = self.genes.get(gene,(0,len(self.ref)-1))
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

