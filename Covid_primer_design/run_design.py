from primer import main_Counter,main_Homology,main_Quality
# from multiprocessing import Process
import sys
from mymodule import poolMap
import os
from datetime import datetime

method = main_Quality

name = 'LAMP_100Inclusivity'

savepath = "./LAMP_primer_design_output"

if __name__ == '__main__':
    print('++++++++++++'*3)
    print(f'=====Started on {datetime.now()}  ======')
    print('++++++++++++'*3)
    def task(work):
        method(*work)
        return None
    workload = [ ('W',(i*500,(i+1)*500),10000,os.path.join(savepath,name+f'{i*0.5}-{(i+1)*0.5}K_.csv')) for i in range(2,58)] #
    def callback(progress):
        print('Overall progress {:.2f}%...'.format(progress))
    result = poolMap(task,workload,chunks=10000,progress_callback=callback)
    print('++++++++++++'*3)
    print(f'=====Finished on {datetime.now()}  ======')
    print('++++++++++++'*3)
