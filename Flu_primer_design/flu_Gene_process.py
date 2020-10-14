import sys
sys.path.append('../')
from primer_design.align_sequence import lines_to_dict,read,AlignmentF
import matplotlib.pyplot as plt
import re
import numpy as np
from collections import Counter,namedtuple
import xml.etree.ElementTree as ET
from dateutil.parser import parse

class Record(dict):
    @property
    def year(self):
        if self.date:
            try:
                return parse(self.date).year
            except:
                return None
        return None
    def __repr__(self):
        return f"{self.accession}-{self.serotype}"
    def __getattr__(self,key):
        return self.get(key,None)
    def __setattr__(self,key,value):
        self[key] = value

class FluRecords:
    def __init__(self,fastaFile,accessionFile):
        self.records = self.parse_accession_xml(accessionFile)
        self.keys = list(self.records.keys())
        fasta = lines_to_dict(read(fastaFile))
        
        for k,seq in fasta.items():
            accession = k[1:9]
            if accession in self.records:
                self.records[accession]['sequence'] = seq
                
    def __getitem__(self,accession):
        if isinstance(accession,int):
            return self.records.get(self.keys[accession])
        return self.records.get(accession,{})
        
    def parse_accession_xml(self,file):
        tree = ET.parse(file)
        root = tree.getroot()
        result = {}
        for item in root:
            entry = Record()
            for data in item:
                entry[data.tag] = data.text
            result[entry['accession']] = entry
        return result
                
    def align(self):
        self.aln = AlignmentF(sequence=list( 
        i['sequence'] for i in self.records.values() if i.get('sequence',None)))
        return self


    def plot_year(self,hasSeq = True):
        years = []
        for record in self:
            if record.get('sequence',None) or (not hasSeq):
                date = record['date']
                if date:
                    try:
                        year = parse(date).year
                    except:
                        year = 0
                    years.append(year)

        year_counter = list(Counter(years).items())
        year_counter.sort()
        labels = [i[0] for i in year_counter]
        year_height =  [i[1] for i in year_counter]
        x_position = list(range(len(labels)))

        fig,ax = plt.subplots()
        ax.bar(x_position,year_height,)
        ax.set_xticks(x_position[::10])
        ax.set_xticklabels(labels[::10], )
        ax.set_title('Count in years')
        plt.show()
        return fig
    
    def fastA_after_year(self,year):
        fasta = {}
        for record in self:
            if record.year >= year and record.sequence:
                fasta[record.accession] = record.sequence
        return fasta
    
    def write_fasta(self,fasta,filename):
        with open(filename,'wt') as f:
            for name, seq in fasta.items():
                f.write(f'>{name}\n')
                f.write(seq+'\n')
    def plot_logo(self,ax):
        return self.aln.dna_logo(save=False,show=False,ax=ax)
    
    def plot_info(self,ax):
        en = 2.88539008/max(sum(self.aln.count),1.5)
        info = (np.log(5)/np.log(2)-en-self.aln.entropy())
        return ax.plot(info)





fluAMacessionFile = '/home/hui/Downloads/FluA_M_N.A._HxNx_accession.xml'
fluAMfasta = '/home/hui/AptitudeUsers/R&D/Users/Hui Kang/flu/FluA_genes/InfluenzaA_M1_M2.fa'


amR = FluRecords(fluAMfasta,fluAMacessionFile)
len(lines_to_dict(read(fluAMfasta)))
amR[0]
amR[1]
len([i for i in amR if i.sequence])
amR[1].sequence

len(amR.records)

types = Counter([i.serotype for i in amR])
types

years = Counter([i.year for i in amR])
years


fasta = {i.accession: i.sequence for i in amR if i.serotype == 'H1N1' and i.sequence}
len(fasta)


amR.write_fasta(fasta, 'FluA_M_N.A._H1N1.fasta')












