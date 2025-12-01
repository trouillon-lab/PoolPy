import numpy as np
import math
import re
import itertools
from itertools import combinations as comb
import pandas as pd
import argparse
import pickle
import time
import os
import copy
import json


from Functions import *
from Fast_functions import *



parser = argparse.ArgumentParser(description='Parse some arguments')
parser.add_argument('--differentiate', type=int, default=-1, help='An integer argument with default 2')
parser.add_argument('--path_to_WA', type=str, help="A string argument with default './pooling_results'")
parser.add_argument('--readout', type=str, help="A string either containing the readout or containing a path a csv of the readout (readout in the form 0,1,0,1,0,0 or 3,6,14)")

args = parser.parse_args()

args_dict = vars(args)



dira=args_dict['path_to_WA']
diff_deco=args_dict['differentiate']
readout_in=args_dict['readout']
WA_df=pd.read_csv(dira, index_col=0)

if readout_in.endswith('csv'):
    readout = np.genfromtxt('readout_in', delimiter=',', dtype=int)
else:
    readout = np.fromstring(readout_in, sep=',', dtype=int)   

WA=WA_df.values
n_compounds=WA.shape[1]

if np.max(readout)>1 or len(readout)!=n_compounds:
    readout_bin_ls = [1 if i in readout else 0 for i in range(n_compounds)]
    readout_bl=np.array(readout_bin_ls)
 

st_dir=str(dira)


mask = ~np.any((WA == 1) & (readout_bl == 0), axis=1)
                    #msg+=f'boolean readout {readout_bl}<br>'
original_indices = np.where(mask)[0]  # Get original row indices
filtered_WA = WA[mask]
n_compounds=filtered_WA.shape[0]
if diff_deco> n_compounds:
    diff_deco=n_compounds
#msg+=f'<br>we have {filtered_WA.shape} shape of the filtered WA<br>'
msg=''
if n_compounds<2:
    if n_compounds==1:
        decoded=[original_indices[0]]
        msg += '<span style="color: #2ecc40;"><b>A single positive sample was found:</b></span><br>'
        for deco in decoded:
            msg += f'<span style="color: #2ecc40;"><b>Sample: {deco}</b></span><br>'
    else:
        msg += '<span style="color: #c00;"><b>We found no matches for the given parameters, check your input or try increasing the differentiate value.</b></span>'                    
else:
    ls_combs=[comb(n_compounds,i) for i in range(diff_deco)]
    max_combs=np.sum(ls_combs)
    if max_combs>1e4:
        decoded = [int(original_indices[idx]) for idx in range(len(original_indices))]
        msg += (
            '<b>Putative positive samples were identified</b>, but the app does not have the computational power to attempt to decode the exact combination.<br>'
            'Either test all putative positive samples individually or change pooling strategy. A lower differentiate (only if it makes sense) might narrow it down.<br>'
            f'There are up to {min([diff_deco,len(decoded)])} positive samples among the following samples: <b>{decoded}</b>.'
        )   
    

    scrambler={1:np.arange(n_compounds)}
    for j in range(2,diff_deco+1):
        scrambler.update({j:np.array(list(comb(np.arange(n_compounds),j)))})

    decoded_pre=decode_precomp(well_assigner=filtered_WA, differentiate= diff_deco, scrambler=scrambler, 
            readout=readout_bl)
    # Map filtered indices back to original indices
    decoders = [combination if isinstance(combination, list) else [combination] for combination in decoded_pre]
    decoded = [[int(original_indices[idx]) for idx in combination] for combination in decoders]

    if len(decoded)==0:
        msg+= '<b>We found no matches for the given parameters, check your input or try increasing the differentiate value.'

    elif len(decoded)==1:

        if len(decoded[0])==1:
            msg += '<span style="color: #2ecc40;"><b>A single positive sample was found:</span></b><br>'
            for deco in decoded[0]:
                msg += f'<span style="color: #2ecc40;"><b>Sample: {deco}</b></span><br>' 
        else:
            msg+=f'<span style="color: #2ecc40;"><b>A single possible combination of positive samples was found. The positive samples are:</span></b><br>'
            msg += f'<span style="color: #2ecc40;"><b>Samples: {", ".join(map(str, decoded[0]))}</b></span><br>'

    elif len(decoded)>n_compounds:
        decoded_set = list(set([x for combo in decoded for x in combo]))
        decoded_set.sort()
        msg += (
            '<b>Putative positive samples were identified</b>, but the exact combination could not be pinpointed.<br>'
            'Either test all putative positive samples individually or change pooling strategy. A lower differentiate (only if it makes sense) might narrow it down.<br>'
            f'There are up to {min([diff_deco,len(decoded_set)])} positive samples among the following samples: <b>{decoded_set}</b>.'
        )
    else:
        msg += f'{len(decoded)} possible combinations of positive samples were found. The possible combinations are:<br>'
        for i,deco in enumerate(decoded):
            if i!=0:
                msg+='or<br>'
            msg += f'<b>Samples: {", ".join(map(str, deco))}</b><br>'

message = re.sub(r'</?b>', '', msg)
txt = f"Uploaded file: {dira}\nReadout: {readout_in}\n\n{message.replace('<br>', '\n')}"


# Print the message to stdout
print(txt)

# Save the message to a file under a 'decoded' folder next to the WA file
wa_dir = os.path.dirname(dira)
decoded_dir = os.path.join(wa_dir, 'decoded')
os.makedirs(decoded_dir, exist_ok=True)
base_name = os.path.splitext(os.path.basename(dira))[0]
out_path = os.path.join(decoded_dir, f"{base_name}_decoded.txt")
with open(out_path, 'w', encoding='utf-8') as f:
    f.write(txt)



