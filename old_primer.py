from mymodule import revcomp
from mymodule import ViennaRNA
from mymodule.RNAstructure import RNA
def find_LOOP_primers(LF=None,LB=None,primers=[]):
    """
    scan for the minimum overlapping length for a primer to other primers in LAMP.
    """
    def scan(s,ps):
        for i in range(3,len(s)):
            for j in range(len(s)-i+1):
                end = s[j:j+i]
                find = True
                for p in ps:
                    if end in p:
                        find=False
                        break
                if find:

                    return i, s[j:]
        return None, None

    if LF:
        print('Finding LF primers')
        i,s = scan(LF,primers + [revcomp(LF)])
        print(f"Find LF {s}, end length {i}")
    if LB:
        print('Finding LB primers')
        i,s = scan(revcomp(LB),primers + [LB])
        print(f'Find LB {s and revcomp(s)}, end length {i}')

def predictMFE(sequence,backbone='dna'):
    if backbone=='rna':
        fc = ViennaRNA.fold_compound(sequence)
        _,mfe=fc.mfe()
        return mfe
    else:
        p = RNA.fromString(sequence,'dna')
        p.FoldSingleStrand()
        return p.GetFreeEnergy(1)

def piece_together(remainings,minimum_overlap=4,max_overlap=5,current="",pieces=[]):
    """
    piece together a list of sequences with given overlaps.
    """
    find = False
    for i in range(minimum_overlap,max_overlap+1):
        overlap = current[-i:]
        for k,s in enumerate(remainings):
            if overlap == s[0:i]:
                newremainings = remainings.copy()
                newremainings.pop(k)
                find = True
                newpieces = pieces.copy()
                newpieces.append(s)
                yield from piece_together(newremainings,minimum_overlap,max_overlap,current=current + s[i:],pieces=newpieces)
    if not find:
        yield current,pieces

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

def primer_end_scan(s1,s2):
    """
    find the maximum sequence on 3' of s1 that hybridize to s2.
    """
    s1rc = revcomp(s1)
    for i in range(1,len(s1rc)):
        totest = s1rc[0:i]
        if totest not in s2:
            break
    totest = totest[0:-1]
