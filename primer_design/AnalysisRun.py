from primer_design.config import config

# always config files first.
config(refApe="",
       refAln="",
       batAln="")


"""Run lamp primers analysis"""
from primer_design.analysis import PrimerSetRecordList,PrimerSetRecord
# from multiprocessing import Process
from mymodule import poolMap
from datetime import datetime



if __name__ == '__main__':
    print('++++++++++++'*3)
    print(f'=====Started on {datetime.now()}  ======')
    print('++++++++++++'*3)
    def task(record):
        pr = PrimerSetRecord.hyrolize(record)
        (pr.Inclusivity()
           .Amplicon_pos()
           .CrossReactivity()
           .Tm()
           .NonTarget()
           .Hairpin()
           .PrimerDimer()
           .LoopHairpin()
           .ExtensionStartGCratio(forward=8)
           .Gaps()
           .GC_ratio()
           .Length()
           .CheckFilter()
            )
        return pr.serialize()

    def callback(progress):
        print('Overall progress {:.2f}%...'.format(progress))
    pl = PrimerSetRecordList('./LAMP_primer_design_output/M_PE_design.csv')
    result = poolMap(task,[i.serialize() for i in pl],progress_callback=callback)
    plresult = PrimerSetRecordList(PrimerSetRecord.hyrolize(i) for i in result)
    plresult.save_csv('./LAMP_primer_design_output/M_PE_design_processed.csv')
    
#     print(f'=====Started N gene Processing on {datetime.now()}  ======')
#     pl = PrimerSetRecordList('./LAMP_primer_design_output/PrimerExplore_Ngene.csv')
#     result = poolMap(task,[i.serialize() for i in pl],progress_callback=callback)
#     plresult = PrimerSetRecordList(PrimerSetRecord.hyrolize(i) for i in result)
#     plresult.save_csv('./LAMP_primer_design_output/PrimerExplore_Ngene_processed.csv')

    
    print('++++++++++++'*3)
    print(f'=====Finished on {datetime.now()}  ======')
    print('++++++++++++'*3)
