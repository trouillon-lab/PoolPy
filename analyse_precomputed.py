# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 12:03:03 2024

@author: jtrouillon
"""
import pandas as pd
import os
import pickle

# plotting
import matplotlib.pyplot as plt
import matplotlib

#Set general plot parameters
subplot_height=6
axis_label_fontsize=18
tick_label_fontsize=16
plt.rcParams['axes.linewidth'] = 1
# Set color palette
pal=['#F6DDB5','#F7CCBC','#F6B6BE','#CEA1B6','#9A9BC1','#718BBB','#5868A0']
extra_pal=['#828386','#53BBB0']
cmap_pal= matplotlib.colors.LinearSegmentedColormap.from_list("", pal)

"""define paths"""
dir_path=os.path.abspath(r'C:\Users\jtrouillon\Desktop\Python Scripts\pooling')
# dif_dir= dir_path +"\\diff_4"
# out_dir= dir_path +"\\precomputed"
plot_dir= dir_path +"\\plots"
os.chdir(dir_path)

# dif_path=os.listdir(dif_dir)
# out_path=os.listdir(out_dir)

#%% Concatenate all pickles into one df for each differentiate value
# Needed only once

''' deprecated, use file_merger '''


# dfs = []
# os.chdir(dif_dir)
# for filename in dif_path:
#     if '.pk' in filename:
#         newDF = pd.read_pickle(filename)
#         dfs.append(newDF)
        
# #concat
# full_dict={}
# for d in dfs:
#     for k, v in d.items():
#         if isinstance(k,int): #Ignores kwarg info
#             full_dict[k]=(v)

# os.chdir(dir_path)

# #save as one pickle
# file_path=os.path.join(out_dir,"diff_4.pk")

# with open(file_path, 'wb') as handle:
#     pickle.dump(full_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
#%% Get data

#Open all pickles
# os.chdir(out_dir)
# for filename in out_path:
#     if '.pk' in filename:
#         name=filename[:-3]
#         exec(name + "=pd.read_pickle(filename)")
        
# os.chdir(dir_path)

# all_precomputed=[diff_1,diff_2,diff_3,diff_4]
# diff_names=['diff_1','diff_2','diff_3','diff_4']


data_file='Final_inline_precomputed_file.pk'
#data_file='Final_precomputed_file.pk'
all_precomputed=pd.read_pickle(data_file)


all_methods = ['matrix', 'multidim', 'random', 'STD', 'Chinese trick','Binary']

#%% test plot

# Function to generate the scatter plot
def plot_specific_method(data, selected_method, to_plot='mean_experiments',plot_prefix=''):

    plt.figure(figsize=(10, 6))

    color_map = {diff: pal[i % len(pal)] for i, (diff,v) in enumerate(data.items())}

    for i, (diff, main_dict) in enumerate(data.items()):
        test_numbers = []
        mean_values = []

        for test_number, (df, _) in main_dict.items():
            for method in df.index:
                # Trim method name if it starts with 'multidim:'
                base_method = method.split(':')[0].strip()
                if selected_method==base_method:
                    test_numbers.append(test_number)
                    mean_values.append(df.loc[method, to_plot])

        plt.scatter(test_numbers, mean_values, label=diff,color=color_map[diff], alpha=0.7, s=80)


    plt.xticks(fontsize=tick_label_fontsize)
    plt.yticks(fontsize=tick_label_fontsize)
    plt.legend(edgecolor='black')
        
    plt.title(f'Method: {selected_method}',size=axis_label_fontsize)
    plt.xlabel('Test Numbers',size=axis_label_fontsize)
    plt.ylabel(f'{to_plot.replace("_", " ").capitalize()}',size=axis_label_fontsize)
    plt.legend()
    
    os.chdir(plot_dir)
    #plt.savefig(str(plot_prefix+selected_method+'_'+to_plot+'.png'), dpi=300)
    os.chdir(dir_path)

    plt.show()


def plot_specific_diff(data, diff, metrics,plot_prefix=''):

    for metric in metrics:
        plt.figure(figsize=(10, 6))
    
        main_dict = data['Differentiate '+str(diff)]
        test_numbers = []
        method_values = {}
        
        
    
        for test_number, (df, _) in main_dict.items():
            test_numbers.append(test_number)
            for method in df.index:
                # Trim method name if it starts with 'multidim:'
                base_method = method.split(':')[0].strip()
                if base_method not in method_values:
                    method_values[base_method] = []
                method_values[base_method].append(df.loc[method, metric])
    
        color_map = {method: pal[i % len(pal)] for i, method in enumerate(method_values.keys())}
    
        for method, values in method_values.items():
            plt.scatter(test_numbers, values, label=method,color=color_map[method], alpha=0.7, s=80)
            
        plt.xticks(fontsize=tick_label_fontsize)
        plt.yticks(fontsize=tick_label_fontsize)
        plt.legend(edgecolor='black')
        
        
        plt.title(f'Differentiate {diff}',size=axis_label_fontsize)
        plt.xlabel('Test Numbers',size=axis_label_fontsize)
        plt.ylabel(f'{metric.replace("_", " ").capitalize()}',size=axis_label_fontsize)
        plt.legend()
        
        os.chdir(plot_dir)
        #plt.savefig(str(plot_prefix+'diff_'+str(diff)+'_'+metric+'.png'), dpi=300)
        os.chdir(dir_path)
        
        plt.show()
    


#%% Plots

# Plot all diffs for each method
for selected_method in all_methods:
    plot_specific_method(all_precomputed, selected_method)


# Plot all methods for each diff
metrics = ['mean_experiments', 'max_compounds_per_well', 'n_wells', 'percentage_check', 'mean_extra_exp']

for diff in [1,2,3,4]:
    plot_specific_diff(all_precomputed, diff,metrics)



#single
plot_specific_method(all_precomputed, 'STD','n_wells',plot_prefix='STDonly_')
plot_specific_diff(all_precomputed, 1,metrics=['n_wells','max_compounds_per_well'],plot_prefix='STDonly_')
