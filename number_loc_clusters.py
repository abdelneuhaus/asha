import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_poca_files(repertory):
    a = [os.path.join(dirpath,filename) for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('merged.txt')]
    return [x.replace("\\", "/") for x in a]

def read_poca_files(file):
    df = pd.read_csv(file)
    df.loc[df['# seq OFF'] > 5000, '# seq OFF'] = df['blinks']
    return df.iloc[:,:-1]

def read_locPALMTracer_file(file):
    data=pd.read_csv(file, sep='\t', skiprows=2)
    data['Plane'] = [int(i) for i in data['Plane']]
    data['Index'] = [int(i) for i in data['Index']]
    return data

def get_PALMTracer_files(repertory):
    a = [os.path.join(dirpath,filename) for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('locPALMTracer.txt')]
    return [x.replace("\\", "/") for x in a]


# repertory = 'D:/240925_W1_FPs/'
# # poca_files = get_poca_files(repertory)
# pt_files = get_PALMTracer_files(repertory)
# # nclusters = 0
# nlocs = 0
# # for i in poca_files:
# #     raw_file_poca = read_poca_files(i)
# #     nclusters += len(raw_file_poca)
# for i in pt_files:
#     raw_pt = read_locPALMTracer_file(i)
#     nlocs += len(raw_pt)
# # print(nclusters)
# print(nlocs)