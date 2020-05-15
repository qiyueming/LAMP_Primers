"""
generate analysis from primer design csv file.
"""
from ape import REFape,read_primerset_excel
from align_sequence import REF,BAT
import pandas as pd
from primer import GC_ratio
from mymodule import revcomp

primerset = read_primerset_excel()

primerset[0]

def GC_ratio(seq):
    return (seq.count('G')+seq.count('C'))/len(seq)


def analysis_pipe(sequeces):
    """
    analysis pipe. 
    """

REF.genes
REFape[0:10]

GC_ratio(REFape[1000:5000])

primer={
'O1-F3':"TCAAGGATGCTACTCCTTCAGA",
'O1-F2':"GCTACTGCAACGATACCGAT",
'O1-F1':"TTATTGTTGGCGTTGCACTTCTTGC",
'O1-B1':"AGCGCTTCCAAAATCATAACCCTCA",
'O1-B2':"GGTGTTCACTTTGTTTGCAACT",
'O1-B3':"ACTCACACCTTTTGCTCGTT",
'O1-LF':"ACAAGCCTCACTCCCTTTCGG",
'O1-LB':"AAAAGAGATGGCAACTAGCACTCTC",}
revcomp('ACAAGCCTCACTCCCTTTCGG')

REF.check_inclusivity(*primer.values())


ape = REFape.label_from_primers(
{'F3':'TCAAGGATGCTACTCCTTCAGA',
'B3':'AACGAGCAAAAGGTGTGAGT',
'FIP':'GCAAGAAGTGCAACGCCAACAATAAGCTACTGCAACGATACCGAT',
'BIP':'AGCGCTTCCAAAATCATAACCCTCAAGTTGCAAACAAAGTGAACACC',
'LF': 'CCGAAAGGGAGTGAGGCTTGT',
'LB':'AAAAGAGATGGCAACTAGCACTCTC'
},'New',updateself=False
)

ape.save_ape('new predict.ape')
