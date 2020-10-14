import sys
sys.path.append('../')

from primer_design.config import config
import flu_para
# always config files first.
config(refApe="",
       refAln="./fluA_post2009_unique_aln.fasta",
       batAln="",
       para=flu_para)

from primer_design.config import REF,PARAMETER

PARAMETER.PInclu

from primer_design.design import *
from primer_design.align_sequence import *
import matplotlib.pyplot as plt
import numpy as np


all_align = lines_to_dict(read('./fluA_post2009_unique_aln.fasta'))
c = Counter(all_align.values())

aln = AlignmentF(sequence=REF.aln)



REF.find_seq('GACTTGAAGATGTCTTTGC')
len(REF)


count = 0

for seq in REF.aln:
    if 'GACTTGAAGATGTCTTTGC' in seq:
        count +=1
count

len(REF.aln)




REF.ref


REF._check_inclusivity('GACTTGAAGATGTCTTTGC')


ref = aln.rep_seq()

# main_Quality(span=(0,1215-25))
