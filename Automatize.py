import numpy as np
import math
import re
import itertools
import pandas as pd
import argparse
import pickle
import time
import os
import copy
import json


from Functions import *
from Fast_functions import *


def from_id_to_well(id):
    plate=id//96
    id=id-plate*96
    row=id//12
    column=id-row*12
    well=chr(65 + row)
    return 'Plate'+str(plate+1), str(well)+str(column)



parser = argparse.ArgumentParser(description='Parse some arguments')
parser.add_argument('--path_to_WA', type=str, help="A string argument containing the path to the pooling strategy file.")

args = parser.parse_args()

args_dict = vars(args)

dira=args_dict['path_to_WA']

WA_df=pd.read_csv(dira, index_col=0)

WA=WA_df.values

ls_to_df=[]
for i,j in zip(*np.where(WA==1)):
    d_plate, d_well= from_id_to_well(i)
    r_plate, r_well= from_id_to_well(j)
    ls_to_df.append([d_plate,d_well,r_plate,r_well])

df_out=pd.DataFrame(data=ls_to_df, columns=['Source', 'SourceWell', 'Dest', 'DestWell'])
df_out['Volume']=1

df_out.to_csv('Robot_protocol.csv', index=False)







