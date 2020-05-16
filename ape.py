"""
Add features to a plain ape file.
compare primers to virus genome, automatically generate amplicon for each primer set.
"""
import random
from mymodule import revcomp
import re
import pandas as pd
import os




def colormap(tag,rev):
    "rev=1 if fragment is reversecomp"
    color=["#fe6dbc","#c7d7c0","#6600ed","#c4d9f8","#f8e3c4","#aee5d8","#ebedff","#b8f2e6","#009468","#f8b290",
    "#fca085","#462d00","#00c290","#ffd10f","#deccad","#d3e9fe","#ed1c24","#00c38f"]

    lamps = {"F3":0 ,"F2":1, "F1":2, "LF":3,"LB":4,"B1":9,"B2":10,"B3":-2}
    for k,i in lamps.items():
        if k in tag:
            return color[i]
    return random.choice(color)



def read_primerset_excel(path="./LAMPprimers/LAMP primers.xlsx",sheet=None):
    """
    read primers for LAMP from an excel file. skipping sequence lines that contain non alphabet letter.
    """
    df = pd.read_excel(path,sheet_name=sheet, usecols=[0,1,5,9],header=None,
            names=["Set","label",'Tm',"Sequence"],na_filter=False)
    primerset = []
    if isinstance(df,dict):
        df = pd.concat(list(df.values()))
    df.index=list(range(len(df)))

    for i in range(len(df)):

        Set,label,Tm,seq = df.loc[i,:]
        if label=='label':
            primerset.append({'set':Set})
        elif seq and label and seq.replace(' ','').isalpha():
            seq = seq.replace(' ','')
            if 'feature' in primerset[-1]:
                primerset[-1]['feature'].append((label,Tm,seq))
            else:
                primerset[-1]['feature'] = [(label,Tm,seq)]
    return primerset


class APE:
    """
    a python rep of ape structure, stripped version.
    """
    def __init__(self,filepath=None):
        if filepath:
            self.path = filepath
            self.protected=True
            with open(filepath,'rt') as f:
                data = f.read()
            self.sequence = "".join(filter(lambda x:x in "ATCGN" , data[data.index('ORIGIN')+7:].upper()))
            self.raw_sequence = "".join(filter(lambda x:x in "ATCGNatgcn" , data[data.index('ORIGIN')+7:])) # without converting
            self.features = []
            if "FEATURES             Location/Qualifiers" in data:
                features = data[data.index("FEATURES             Location/Qualifiers\n")+41:data.index("ORIGIN")-1]
                for line in features.split('\n'):
                    if not line.strip():
                        continue
                    elif line[6] != " ":
                        t,p=line.split()
                        s,e = re.search(re.compile('(\d+)\.\.(\d+)'),p).groups()
                        self.features.append({'type':t,'start': int(s) ,'end': int(e),'rev': 'complement' in p})
                    elif '/locus_tag' in line:
                        self.features[-1]['tag'] = re.search(re.compile('"(.*)"'), line).groups()[0]
                    elif 'fwdcolor' in line:
                        self.features[-1]['fwdcolor'] = re.search(re.compile('"(.*)"'), line).groups()[0]
                    elif 'revcolor' in line:
                        self.features[-1]['revcolor'] = re.search(re.compile('"(.*)"'), line).groups()[0]
    def read_seq(self,file):
        "read sequence from a string or file"
        if os.path.exists(file):
            with open(file,'rt') as f:
                data = f.read()
        else:
            data = file
        self.path = 'new APE.ape'
        self.sequence = "".join(filter(lambda x:x in "ATCGN" , data.upper()))
        self.raw_sequence = "".join(filter(lambda x:x in "ATCGNatgcn" , data))
        self.features = []

    def __getitem__(self,slice):
        return self.sequence[slice]

    def __repr__(self):
        return f" APE {getattr(self,'path','Unsaved')}"

    def __len__(self):
        return len(self.sequence)

    def save_ape(self,path=None):
        if path == None: path=self.path
        with open(path,'wt') as f:
            f.write(self.dumps())

    def name_primer(self,seq):
        _,end=self.locate_primer(seq)
        for d in self.features:
            if end> d['start'] and end< d['end']:
                return d['tag']
        return 'Unknown'


    def dumps(self):
        s = ['LOCUS','FEATURES             Location/Qualifiers']
        for f in self.features:
            start = f['start']
            end = f['end']
            rev = f['rev']
            tag = f['tag']
            fwdcolor = f['fwdcolor']
            revcolor = f['revcolor']
            pos = f"complement({start}..{end})" if rev else f"{start}..{end}"
            s.append(f"     {f['type']}    {pos}")
            s.append(f'                     /locus_tag="{tag}"')
            s.append(f'                     /ApEinfo_fwdcolor="{fwdcolor}"')
            s.append(f'                     /ApEinfo_revcolor="{revcolor}"')
        s.append('ORIGIN')
        for i in range(int(len(self.raw_sequence)/60)+1):
            fragment = self.raw_sequence[60*i:60*(i+1)]
            frags = ["{:>9}".format(i*60+1)]
            frags.extend(fragment[j*10:j*10+10] for j in range(6))
            s.append(' '.join(frags))
        s.append('//')
        return '\n'.join(s)

    def locate_primer(self,sequence):
        "return start and end of a primer 5' to 3', 0 based index, "
        rev = False
        if sequence not in self.sequence:
            rev = revcomp(sequence)
            if rev in self.sequence:
                sequence = rev
            else:
                raise ValueError(f'Primer {sequence} not found')
        ind = self.sequence.index(sequence)
        if rev:
            return ind+len(sequence) , ind
        return ind, ind+len(sequence)

    def locate_primer_position(self,sequence):
        """return 1 based index, can be directly used to add features."""
        s,e = self.locate_primer(sequence)
        if s<e:
            return s+1,e
        else:
            return s,e+1

    def add_feature(self, start, end, tag,type='misc_feature',cm=colormap):
        """
        start and end are not 0 based.
        """
        if isinstance(cm,str):
            b=cm
            cm = lambda x,y:b
        self.features.append({
        'type':type,
        'start':min(start,end),
        'end':max(start,end),
        'rev': start>end,
        'tag':tag,
        'fwdcolor':cm(tag, start>end),
        'revcolor':cm(tag, start>end)
        })

    def remove_feature(self,key):
        "remove feature based on a string or function"
        if isinstance(key,str):
            self.features=list(filter(lambda x:x['tag']!=key,self.features))
        else:
            self.features=list(filter(lambda x: not key(x),self.features))

    def truncate(self,start,end,keepfeatures='gene'):
        """
        start and end are 0 based index.
        """
        features = []
        for f in self.features:
            if keepfeatures!='all' and f['type'] != keepfeatures:
                continue
            if f['start'] >= end or f['end'] <= start:
                continue
            new = f.copy()
            new['start'] = max(1, new['start'] - start )
            new['end'] = min(end-start, new['end'] - start)
            features.append(new)
        ape = APE()
        ape.sequence = self.sequence[start:end]
        ape.raw_sequence = self.raw_sequence[start:end]
        ape.features = features
        ape.protected = False
        return ape

    def label_sequence(self,sequence,name,**kwargs):
        """
        label a sequence with name, kwargs can have type, cm
        """
        start,end = self.locate_primer_position(sequence)
        self.add_feature(start,end,name,**kwargs)

    def find_LAMP_amplicon(self,primerset):
        f = primerset['feature']
        pos = []
        for tag,tm,s in f:
            s,e = self.locate_primer(s)
            pos.append(s)
            pos.append(e)
        start = min(pos)
        end = max(pos)
        ape = self.truncate(start,end,keepfeatures='gene')
        if ape.features:
            maxlengthfeature = max(ape.features,key=lambda x: x['end']-x['start'])['tag']
        else:
            maxlengthfeature=""

        # if there is only one feature, remove it.
        if len(ape.features) == 1:
            ape.features = []

        ape.add_feature(1, len(ape), f"Primer Set {primerset['set']} {maxlengthfeature} {start}-{end}",cm='#e0e0e0')
        for tag,tm,s in f:
            ape.label_sequence(s,f"{tag} {tm}")

        return ape

    def label_from_primers(self,primers,name,updateself=True,):
        """
        label lamp amplicons from its primers
        priemrs = {'F3','B3','FIP','BIP','LF','LB'}
        name is the name of the amplicon.
        """
        f3=primers.pop('F3')
        b3=primers.pop('B3')
        fip = primers.pop('FIP')
        bip = primers.pop('BIP')
        s,_ = self.locate_primer(f3)
        e,_=self.locate_primer(b3)
        ape = self.truncate(s,e,keepfeatures='gene')
        if ape.features:
            maxlengthfeature = max(ape.features,key=lambda x: x['end']-x['start'])['tag']
        else:
            maxlengthfeature=""
        if len(ape.features) == 1:
            ape.features = []
        ape.add_feature(1, len(ape), f"{name} in {maxlengthfeature} {s}-{e}",cm='#e0e0e0')
        if primers:
            for k,i in primers.items():
                ape.label_sequence(i,k)
        ape.label_sequence(f3,'F3')
        ape.label_sequence(b3,'B3')

        # label F1c
        fipl = len(fip)
        for i in range(0,fipl):
            totest = fip[fipl-1-i:fipl]
            if totest not in ape.sequence:
                break
        f2 = fip[fipl-i:]
        f1c = fip[0:fipl-i]
        ape.label_sequence(f2,'F2')
        ape.label_sequence(f1c,'F1c')

        #label Bic
        bipl = len(bip)
        for j in range(0,bipl):
            totest = bip[0:j+1]
            if totest not in ape.sequence:
                break
        b1c = bip[0:j]
        b2 = bip[j:]
        ape.label_sequence(b2,'B2')
        ape.label_sequence(b1c,'B1c')

        if updateself:
            self.add_feature(s+1,e,name)
        return ape

    def list_features(self):
        return [i['tag'] for i in self.features]

    def get_feature(self,name):
        for i in self.features:
            if i['tag'] == name:
                return self.truncate(i['start']-1,i['end'],)
    def get_feature_pos(self,name):
        " return position of a feature 0 index."
        for i in self.features:
            if i['tag'] == name:
                return (i['start']-1,i['end'],)

    def copy_features_from_ape(self,ape,type='gene'):
        "copy certain type of features from another ape"
        for feature in ape.features:
            if feature['type'] == type:
                s = ape.get_feature(feature['tag']).sequence
                try:
                    self.label_sequence(s,feature['tag'],type=type)
                except ValueError:
                    print(f'<{feature["tag"]}> not found')
                    continue

REFape = APE('./viral_genome/cov2_ref/HCoV2 NC045512.2.ape')

if __name__ == '__main__':
    print('Starting...')
    ape = input('Enter virus genome:\n')

    v = APE(ape.replace('\\',"").strip())

    print(f'Virus genome loaded. Length = {len(v)}')
    xls = input('Enter primerset excel file:\n')
    sheet = input('Enter sheet name in the excel:\n')
    primerset = read_primerset_excel(xls.replace('\\',"").strip(),sheet=sheet.strip())
    print(f"Primer sets found: {len(primerset)}")
    print([i['set'] for i in primerset])
    if input('Continue? [Y/n]\n')=='n':
        print('Aborted')
        exit(0)
    else:
        for primer in primerset:
            a = v.find_LAMP_amplicon(primer)
            a.save_ape(f'{primer["set"]} Amplicon.ape')
            print(f'{primer["set"]} Amplicon.ape saved.')
        print('All done.')



#
# v = APE(filepath='../MN908947.3.ape')
#
#
# n2 = v.find_LAMP_amplicon(primerset[0])
# n2.save_ape('n2.ape')

#ps = read_primerset_excel(sheet='N gene')
