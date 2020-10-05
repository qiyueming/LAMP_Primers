from primer_design.ape import APE
from primer_design.align_sequence import Reference
from primer_design.primer_para import *

"""
Configure the reference genes for design.
Load primer design parameters from primer_para.py
"""


global REFape,REF,BAT

def config(refApe="",refAln="",batAln=""):
    global REFape,REF,BAT
   
    REFape = APE(refApe) if refApe else None
    REF = Reference(refAln,REFape) if refAln else None
    if REFape:
        REF.label_gene(REFape)
    BAT = Reference(batAln,REFape) if batAln else None




CROSS_GENE_FILE=['./viral_genome/cross_reactivity/HCoV 229E NC_002645.gb',
 './viral_genome/cross_reactivity/HCoV 229E MN369046 .gb',
 './viral_genome/cross_reactivity/HCoV OC43 KX344031.1.gb',
 './viral_genome/cross_reactivity/HCoV OC43 NC_006213.gb',
 './viral_genome/cross_reactivity/HCoV HKU1 NC_006577.2.gb',
 './viral_genome/cross_reactivity/HCoV HKU1 MK167038.gb',
 './viral_genome/cross_reactivity/HCoV HKU1 MH940245.gb',
 './viral_genome/cross_reactivity/HCoV NL63 NC005831.gb',
 './viral_genome/cross_reactivity/SARS CoV NC_004718.gb',
 './viral_genome/cross_reactivity/MERS NC_019843.gb',
 './viral_genome/cross_reactivity/Adenovirus AC_000017.gb',
 './viral_genome/cross_reactivity/hMPV NC_039199.gb',
 './viral_genome/cross_reactivity/H parainfluenza virus 1 AF457102.gb',
 './viral_genome/cross_reactivity/H Parainfluenza 4a NC_021928.gb',
 './viral_genome/cross_reactivity/H Parainfluenza 4b MN306032.gb',
 './viral_genome/cross_reactivity/Influenza A NC_026435.gb',
 './viral_genome/cross_reactivity/Influenza B NC_002205.gb',
 './viral_genome/cross_reactivity/H rhinovirus NC_038311.gb',
 './viral_genome/cross_reactivity/Bat SARS-CoV 2015 MG772933.1 Bat.ape',
 './viral_genome/cross_reactivity/Bat SARS-CoV 2017 MG772934.1.gb']


CROSS_GENE_NAME = ['HCoV 229E',
 'HCoV 229E-1',
 'HCoV OC43',
 'HCoV OC43-1',
 'HCoV HKU1',
 'HCoV HKU1-1',
 'HCoV HKU1-2',
 'HCoV NL63',
 'SARS CoV',
 'MERS',
 'Adenovirus',
 'hMPV',
 'H parainfluenza virus 1',
 'H Parainfluenza 4a',
 'H Parainfluenza 4b',
 'Influenza A',
 'Influenza B',
 'H rhinovirus',
 'Bat SARS-CoV 2015',
 'Bat SARS-CoV 2017']



CROSS_GENE_NAME_LONG = [
"Human Coronavirus 229E (NC_002645)",
"Human Coronavirus 229E (MN369046)",
"Human Coronavirus OC43 (KX344031.1)",
"Human Coronavirus OC43 (NC_006213)",
"Human Coronavirus HKU1 (NC_006577.2)",
"Human Coronavirus HKU1 (MK167038)",
"Human Coronavirus HKU1 (MH940245)",
"Human Coronavirus NL63 (NC005831)",
"SARS Coronavirus (NC_004718)",
"MERS (NC_019843)",
"Adenovirus (AC_000017)",
"Human Metapneumovirus (NC_039199)",
"Human Parainfluenza 1 (AF457102)",
"Human Parainfluenza 4a (NC_021928)",
"Human Parainfluenza 4b (MN306032)",
"Influenza A (NC_026435)",
"Influenza B (NC_002205)",
"Human rhinovirus (NC_038311)",
"Bat SARS-CoV 2015 (MG772933.1)",
"Bat SARS-CoV 2017 (MG772934.1)",
]

OHTER_ORGANISM = [
 'Chlamydia pneumoniae',
 'Legionella pneumophila',
 'Mycobacterium tuberculosis',
 'Streptococcus pneumoniae',
 'Streptococcus pyogenes',
 'Bordetella pertussis',
 'Mycoplasma pneumoniae',
 'Pneumocystis jirovecii',
 'Candida albicans',
 'Pseudomonas aeruginosa',
 'Staphylococcus epidermis',
 'Staphylococcus salivarius',]
