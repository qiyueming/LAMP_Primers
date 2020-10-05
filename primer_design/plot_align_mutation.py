import pandas as pd
from ape import APE
from align_sequence import lines_to_dict,read,AlignmentF,bp_align

"""
plot mutation heat mat from alignment references.
"""


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
