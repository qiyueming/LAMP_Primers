import sys
sys.path.append('../')

from primer_design.config import config
import flu_para
# always config files first.
config(refApe="",
       refAln="./fluB_M_gene_aln.fast",
       batAln="",
       para=flu_para)

from primer_design.config import REF,PARAMETER

from primer_design.design import *


main_Quality(span=(0,1201-25))
