# -*- coding: utf-8 -*-
"""
Created on Fri Aug 15 13:55:02 2025

@author: nanoseq
"""


import numpy as np
import pandas as pd
import os
import cmcrameri.cm as cmc
from datetime import date
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import textwrap

#%% Prepare path and data

main_path=os.path.abspath(r'C:\Users\nanoseq\Documents\GitHub\pooling')

prevalence_path=os.path.abspath(r'C:\Users\nanoseq\Documents\GitHub\pooling\OLD\prevalence')

#dir_WAs=os.path.abspath(r'D:\precomputed')
dir_WAs=os.path.abspath(r'D:\NEW_precomputed')

plt_path=os.path.join(dir_WAs,'AA_plots',str(date.today()))
os.makedirs(plt_path, exist_ok=True)

os.chdir(main_path)

save_path=False


#Selection of what to plot
N=0
start=20
max_N=101
step=1

min_diff=1
max_diff=4


#Plotting parameters

cmap=cmc.lapaz

label_fontsize=16
tick_fontsize=13
scatter_size=40
scatter_alpha=0.3
jitter=0  # Amount of jitter on x-axis
line_alpha=0.7     # Transparency for the line
line_width=5


# Helper to sample palette uniformly
def sample_palette(palette, n):
    idxs = np.linspace(0, len(palette)-1, n+1, dtype=int)
    return [palette[i] for i in idxs[:-1]]

# Convert colormap to hex colors for seaborn
cmap_hex = [mpl.colors.to_hex(c) for c in cmc.lapaz.colors]


# Get data


y_as_inf=[]
ls_metrics=[]
while N<max_N:
    Npath=os.path.join(dir_WAs,'N_'+str(N))
    if not os.path.exists(Npath):
        N+=step
        continue
    diff=int(min_diff)
    while diff<=max_diff:
        dpath=os.path.join(Npath,'diff_'+str(diff))
        if not os.path.exists(dpath):
            diff+=1
            continue
        metname=os.path.join(dpath, 'Metrics_N_'+str(N)+'_diff_'+str(diff)+'.csv')
        if not os.path.isfile(metname):
            diff +=1
            continue
        df_met=pd.read_csv(metname, header=0)
        if len(y_as_inf)==0:
            y_as_inf=list(df_met)[1:]
        df_met['N']=N
        df_met['diff']=diff

        ls_metrics.append(df_met)
        

        diff+=1


    N+=step

    
full_df_met=pd.concat(ls_metrics)

#Merge the N wells and N pools columns since some were generated with one or the others
if 'N wells' in full_df_met.columns:
    full_df_met['N wells']=full_df_met['N wells'].fillna(full_df_met['N pools'])

#Correct discrepancies between names
full_df_met=full_df_met.replace(to_replace={'Method': 'Chinese remainder'}, value='Chinese Remainder')
full_df_met=full_df_met.replace(to_replace={'Method': 'Chinese trick'}, value='Chinese Remainder')
full_df_met=full_df_met.replace(to_replace={'Method': 'Chinese reminder'}, value='Chinese Remainder')
full_df_met=full_df_met.replace(to_replace={'Method': 'Ch. Rm. Bktrk'}, value='Ch. Rem. Backtrack')
full_df_met=full_df_met.replace(to_replace={'Method': 'Ch. Rm. Special'}, value='Ch. Rem. Special')
full_df_met=full_df_met.replace(to_replace={'Method': 'STD'}, value='Shifted transversal')

#
full_df_met=full_df_met.reset_index().drop(['N pools','Unnamed: 0','index'],axis=1)
full_df_met=full_df_met[full_df_met['Method']!='multidim-2'] #remove multidim-2 since redundant with matrix


full_df_met['test_per_sample']=full_df_met['Mean experiments']/full_df_met['N']
full_df_met['groupsize_vs_N']=full_df_met['Max samples per pool']/full_df_met['N']   


#%% Calculate metric quintiles

#Get method ranks within each combination of N and diff
for N in range(start,max_N):
    if N in full_df_met['N'].values:
        for diff in range(min_diff,max_diff+1):
            temp_df=full_df_met[full_df_met['N']==N]
            temp_df_diff=temp_df[temp_df['diff']==diff]
            
            temp_df_diff['Mean_Exp_rank'] = temp_df_diff['Mean experiments'].rank(pct=True,method='min')
            temp_df_diff['Group_size_rank'] = temp_df_diff['Max samples per pool'].rank(pct=True,method='min')
            temp_df_diff['Mean_step_rank'] = temp_df_diff['Mean steps'].rank(pct=True,method='min')
            
            for idx in temp_df_diff.index:
                full_df_met.loc[idx,'Mean_Exp_rank']=temp_df_diff.loc[idx,'Mean_Exp_rank']
                full_df_met.loc[idx,'Group_size_rank']=temp_df_diff.loc[idx,'Group_size_rank']
                full_df_met.loc[idx,'Mean_step_rank']=temp_df_diff.loc[idx,'Mean_step_rank']
                
#average rank per method
methods = list(full_df_met['Method'].unique())
df_method_ranks=pd.DataFrame(index=methods)    

#check only for diff 1
for method in methods:
    temp_df=full_df_met[full_df_met['Method']==method]
    for diff in range(min_diff,max_diff+1):
        temp_df_diff=temp_df[temp_df['diff']==diff]
        df_method_ranks.loc[method,f'AVG_Mean_Exp_rank_diff_{diff}']=temp_df_diff.loc[:,'Mean_Exp_rank'].mean()
        df_method_ranks.loc[method,f'AVG_Group_size_rank_diff_{diff}']=temp_df_diff.loc[:,'Group_size_rank'].mean()
        df_method_ranks.loc[method,f'AVG_Mean_step_rank_diff_{diff}']=temp_df_diff.loc[:,'Mean_step_rank'].mean()

'''Export'''
#df_method_ranks.to_csv('method_ranking.csv')

#%% Plot metrics across methods

metrics_to_plot=['Mean experiments','Max samples per pool', 'Mean steps','test_per_sample']
ylabels=['Number of tests','Maximum group size','Number of steps','Test per sample']

#keep binary methods for diff>1 or not
keep_binary=True


methods = list(full_df_met['Method'].unique())
n_methods = len(methods)
colors = [cmap(x) for x in np.linspace(0, 0.9, n_methods)]
color_dict = dict(zip(methods, colors))
rng = np.random.default_rng(seed=42)  # For reproducibility


for i in range(len(metrics_to_plot)):
    
    metric_to_plot=metrics_to_plot[i]
    ylabel=ylabels[i]
    
    x_col='N'
    xlabel='Number of samples'
    
    title=False
    
    
    #loop across diffs
    for diff in range(min_diff,max_diff+1):
    
        #subset df
        df_filtered=full_df_met[full_df_met['diff']==diff]
        
        ylim_max=max(df_filtered[metric_to_plot])*1.05 #to keep same ylim across diffs
        
        fig, ax = plt.subplots(figsize=(8, 5))
        
        plt.text(0.05,0.9,'Differentiate '+str(diff), transform = ax.transAxes,fontsize=tick_fontsize)
        
        for method in methods:
            
            #Filter binary methods above diff 1
            if keep_binary==False:
                if diff>1 and method=='Binary': ######################################################################### Binary removed here
                    continue
            
            try:
                grp = df_filtered[df_filtered['Method'] == method]
                y_vals = grp[metric_to_plot].values
                # Jitter x values
                x_vals = grp[x_col].values + rng.uniform(-jitter, jitter, size=len(grp)) if jitter > 0 else grp[x_col].values
                # Sort by jittered x for line plotting
                sort_idx = np.argsort(x_vals)
                x_sorted = x_vals[sort_idx]
                y_sorted = y_vals[sort_idx]
                # Line: through jittered points

                ax.plot(
                    x_sorted, y_sorted,
                    color=color_dict[method],
                    alpha=line_alpha,
                    linewidth=line_width,
                    zorder=1
                )
                ax.annotate(xy=(x_sorted[-1],y_sorted[-1]), xytext=(5,0), textcoords='offset points',
                            text=method, va='center',color=color_dict[method],fontsize=tick_fontsize,alpha=0.7)
                
                # Scatter: with jitter
                # ax.scatter(
                #     x_vals, y_vals,
                #     label=method,
                #     color=color_dict[method],
                #     s=scatter_size,
                #     alpha=scatter_alpha,
                #     zorder=2
                # )
            except:
                continue
            
                          
        if metric_to_plot=='test_per_sample':
            ax.hlines(y=1, xmin=0, xmax=max_N, linewidth=1, color='lightgrey',linestyles='--')
            ax.axhspan(1, 1000, alpha=0.75,color='white')
        
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        
        ax.set_xlim(xmin=0)
        ax.set_ylim(ymin=0,ymax=ylim_max)
        
        if title:
            ax.set_title(title, fontsize=label_fontsize + 2)
        ax.tick_params(axis='both', labelsize=tick_fontsize)
        
        
        #ax.spines[['top', 'right']].set_visible(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
            
        plt.tight_layout()
        
        # Place legend outside the plot
        # ax.legend(
        #     title='Method',
        #     fontsize=label_fontsize,
        #     bbox_to_anchor=(1.05, 1),
        #     loc='upper left',
        #     borderaxespad=0.
        # )
        
        
        suffix=''
        
        if diff>1 and keep_binary==False:
            suffix='_ExclBinary'
        
        save_path = os.path.join(
            plt_path, f"{metric_to_plot.replace(' ','_')}_N_{max_N}_diff_{diff}{suffix}"
        )
        
        
        if plt_path:
            plt.savefig(save_path+'.png', bbox_inches='tight',format='png',dpi=150)
            plt.savefig(save_path+'.svg', bbox_inches='tight',format='svg',dpi=150)
            
        plt.show()
        
        
    
#%% Plot test/sample vs group size 
    
# #metric_to_plot='Max samples per pool'
# metric_to_plot='groupsize_vs_N'

# x_col='test_per_sample'

# ylabel='Relative group size'
# xlabel='Test per sample'


# title=False


# methods = list(full_df_met['Method'].unique())
# n_methods = len(methods)
# colors = [cmap(x) for x in np.linspace(0, 0.9, n_methods)]
# color_dict = dict(zip(methods, colors))
# rng = np.random.default_rng(seed=42)  # For reproducibility

# #loop across diffs
# for diff in range(min_diff,max_diff+1):

#     #subset df
#     df_filtered=full_df_met[full_df_met['diff']==diff]
#     df_filtered=df_filtered[df_filtered['N'].isin([500])]
    
#     fig, ax = plt.subplots(figsize=(4, 5))
    
#     plt.text(0.05,0.9,'Differentiate '+str(diff), transform = ax.transAxes,fontsize=tick_fontsize)
    
#     for method in methods:
        
#         #Filter binary methods above diff 1
#         if diff>1 and method=='Binary': ######################################################################### Binary removed here
#             continue
        
#         try:
#             grp = df_filtered[df_filtered['Method'] == method]
#             y_vals = grp[metric_to_plot].values
#             # Jitter x values
#             x_vals = grp[x_col].values + rng.uniform(-jitter, jitter, size=len(grp)) if jitter > 0 else grp[x_col].values
#             # Sort by jittered x for line plotting
#             sort_idx = np.argsort(x_vals)
#             x_sorted = x_vals[sort_idx]
#             y_sorted = y_vals[sort_idx]
#             # Line: through jittered points

#             ax.scatter(
#                 x_sorted, y_sorted,s=200,
#                 color=color_dict[method],
#                 alpha=0.5,
#                 zorder=1
#             )
#             # ax.plot(
#             #     x_sorted, y_sorted,
#             #     color=color_dict[method],
#             #     alpha=0.3,
#             #     linewidth=line_width,
#             #     zorder=1
#             # )
#             ax.annotate(xy=(x_sorted[-1],y_sorted[-1]), xytext=(15,0), textcoords='offset points',
#                         text=method, va='center',color=color_dict[method],fontsize=tick_fontsize,alpha=0.7)
#             # Scatter: with jitter
#             # ax.scatter(
#             #     x_vals, y_vals,
#             #     label=method,
#             #     color=color_dict[method],
#             #     s=scatter_size,
#             #     alpha=scatter_alpha,
#             #     zorder=2
#             # )
#         except:
#             continue
    
#     ax.set_xlabel(xlabel, fontsize=label_fontsize)
#     ax.set_ylabel(ylabel, fontsize=label_fontsize)
    
#     ax.set_xlim(xmin=0,xmax=0.3)
#     ax.set_ylim(ymin=0,ymax=0.55)
    
#     if title:
#         ax.set_title(title, fontsize=label_fontsize + 2)
#     ax.tick_params(axis='both', labelsize=tick_fontsize)
    
    
#     #ax.spines[['top', 'right']].set_visible(False)
#     for spine in ax.spines.values():
#         spine.set_visible(False)
        
#     plt.tight_layout()
    
#     # Place legend outside the plot
#     # ax.legend(
#     #     title='Method',
#     #     fontsize=label_fontsize,
#     #     bbox_to_anchor=(1.05, 1),
#     #     loc='upper left',
#     #     borderaxespad=0.
#     # )
    
    
#     save_path = os.path.join(
#         plt_path, f"{metric_to_plot.replace(' ','_')}_N_{max_N}_diff_{diff}"
#     )
    
    
#     if plt_path:
#         plt.savefig(save_path+'_dualplot_v2.png', bbox_inches='tight',format='png',dpi=150)
#         plt.savefig(save_path+'_dualplot.svg', bbox_inches='tight',format='svg',dpi=150)
        
#     plt.show()
    
    
    
    
    
#%% Plot by diff

metrics_to_plot=['Mean experiments','Max samples per pool', 'Mean steps','test_per_sample']
ylabels=['Number of tests','Maximum group size','Number of steps','Test per sample']

#keep binary methods for diff>1 or not
keep_binary=True


methods = list(full_df_met['Method'].unique())
n_methods = len(methods)

diffs=[x for x in range(min_diff,max_diff+1)]
n_diff=len(diffs)

colors = [cmap(x) for x in np.linspace(0, 0.6, n_diff)]
color_dict = dict(zip(diffs, colors))
rng = np.random.default_rng(seed=42)  # For reproducibility


for method in methods:
    
    #Filter binary methods above diff 1
    if keep_binary==False:
        if diff>1 and method=='Binary': ######################################################################### Binary removed here
            continue
    
    for i in range(len(metrics_to_plot)):
        
        metric_to_plot=metrics_to_plot[i]
        ylabel=ylabels[i]
        
        x_col='N'
        xlabel='Number of samples'
        
        title=False
    
        
        #subset df
        df_filtered=full_df_met[full_df_met['Method']==method]
        
        ylim_max=max(df_filtered[metric_to_plot])*1.05
        
        fig, ax = plt.subplots(figsize=(8, 5))
        
        plt.text(0.05,0.9,str(method)+' design', transform = ax.transAxes,fontsize=tick_fontsize)
        
        #loop across diffs
        for diff in diffs:


            try:
                grp = df_filtered[df_filtered['diff'] == diff]
                y_vals = grp[metric_to_plot].values
                # Jitter x values
                x_vals = grp[x_col].values + rng.uniform(-jitter, jitter, size=len(grp)) if jitter > 0 else grp[x_col].values
                # Sort by jittered x for line plotting
                sort_idx = np.argsort(x_vals)
                x_sorted = x_vals[sort_idx]
                y_sorted = y_vals[sort_idx]
                # Line: through jittered points

                ax.plot(
                    x_sorted, y_sorted,
                    color=color_dict[diff],
                    alpha=line_alpha,
                    linewidth=line_width,
                    label=diff,
                    zorder=1
                )
                # ax.annotate(xy=(x_sorted[-1],y_sorted[-1]), xytext=(5,0), textcoords='offset points',
                #             text=method, va='center',color=color_dict[method],fontsize=tick_fontsize,alpha=0.7)
                
                # Scatter: with jitter
                # ax.scatter(
                #     x_vals, y_vals,
                #     label=method,
                #     color=color_dict[method],
                #     s=scatter_size,
                #     alpha=scatter_alpha,
                #     zorder=2
                # )
            except:
                continue
                
                          
        if metric_to_plot=='test_per_sample':
            ax.hlines(y=1, xmin=0, xmax=max_N, linewidth=1, color='lightgrey',linestyles='--')
            ax.axhspan(1, 1000, alpha=0.75,color='white')
        
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        
        ax.set_xlim(xmin=0)
        ax.set_ylim(ymin=0,ymax=ylim_max)
        
        if title:
            ax.set_title(title, fontsize=label_fontsize + 2)
        ax.tick_params(axis='both', labelsize=tick_fontsize)
        
        
        #ax.spines[['top', 'right']].set_visible(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
            
        plt.tight_layout()
        
        # Place legend outside the plot
        ax.legend(
            title='Differentiate',
            fontsize=8,
        )
        
        
        suffix=''
        
        if diff>1 and keep_binary==False:
            suffix='_ExclBinary'
        
        save_path = os.path.join(
            plt_path, f"{method.replace(' ','_').replace('.','')}_method_{metric_to_plot.replace(' ','_')}_N_{max_N}{suffix}"
        )
        
        
        if plt_path:
            plt.savefig(save_path+'.png', bbox_inches='tight',format='png',dpi=150)
            plt.savefig(save_path+'.svg', bbox_inches='tight',format='svg',dpi=150)
            
        plt.show()







#%% Plotting prevalence

# Get prevalence df
os.chdir(prevalence_path)
df_prevalence=pd.read_csv('prevalence_results.csv')
os.chdir(main_path)


N_min=10
N_max=100
diff_min=1
diff_max=4


df_prevalence_subset=df_prevalence[df_prevalence['N'].between(N_min,N_max)][df_prevalence['Diff'].between(diff_min,diff_max)]


prevalences=df_prevalence_subset['prevalence'].unique()
n_prevalence=len(prevalences)

colors_prevalence = [cmap(x) for x in np.linspace(0, 0.8, n_prevalence)]
color_dict_prevalence = dict(zip(prevalences, colors_prevalence))

diffs=[x for x in range(min_diff,max_diff+1)]
n_diff=len(diffs)

colors_diff = [cmap(x) for x in np.linspace(0, 0.8, n_diff)]
color_dict_diff = dict(zip(diffs, colors_diff))


# Plot 1: One plot per prevalence
for i, prevalence in enumerate(prevalences):
    
    
    plt.figure(figsize=(12,6))  # Wider figure
    subset = df_prevalence_subset[df_prevalence_subset['prevalence'] == prevalence]
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    for diff in diffs:
        subset_diff=subset[subset['Diff']==diff]
        ax.plot(
            subset_diff['N'], subset_diff['p_error'],
            color=color_dict_diff[diff],
            alpha=line_alpha,
            linewidth=line_width,
            label=diff,
            zorder=1
        )
    
    plt.text(0.05,0.9,f'Prevalence {prevalence}', transform = ax.transAxes,fontsize=tick_fontsize)

    ax.set_xlabel("Number of samples", fontsize=label_fontsize)
    ax.set_ylabel("Probability of error", fontsize=label_fontsize)
    
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0,ymax=1)
    
    if title:
        ax.set_title(title, fontsize=label_fontsize + 2)
    ax.tick_params(axis='both', labelsize=tick_fontsize)
    
    
    #ax.spines[['top', 'right']].set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
        
    plt.tight_layout()
    
    # Place legend outside the plot
    ax.legend(
        title='Differentiate',
        fontsize=8,
    )
    plot_file = os.path.join(plt_path, f'Prevalence_{N_max}_N_vs_p_error_prevalence_{str(prevalence).replace('.','-')}')
    plt.savefig(plot_file+'.svg', dpi=150, format='svg')
    plt.savefig(plot_file+'.png', dpi=150, format='png')
    plt.show()
    
  
    
  
# Plot 2: One plot per diff
for i, diff in enumerate(diffs):
    
    plt.figure(figsize=(12,6))  # Wider figure
    subset = df_prevalence[df_prevalence['Diff'] == diff]
    
    fig, ax = plt.subplots(figsize=(8, 5))
    
    for prevalence in prevalences:
        subset_prevalence=subset[subset['prevalence']==prevalence]
        ax.plot(
            subset_prevalence['N'], subset_prevalence['p_error'],
            color=color_dict_prevalence[prevalence],
            alpha=line_alpha,
            linewidth=line_width,
            label=prevalence,
            zorder=1
        )
    
    plt.text(0.05,0.9,f'Differentiate {diff}', transform = ax.transAxes,fontsize=tick_fontsize)

    ax.set_xlabel("Number of samples", fontsize=label_fontsize)
    ax.set_ylabel("Probability of error", fontsize=label_fontsize)
    
    ax.set_xlim(xmin=0)
    ax.set_ylim(ymin=0,ymax=1)
    
    if title:
        ax.set_title(title, fontsize=label_fontsize + 2)
    ax.tick_params(axis='both', labelsize=tick_fontsize)
    
    
    #ax.spines[['top', 'right']].set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
        
    plt.tight_layout()
    
    # Place legend outside the plot
    ax.legend(
        title='Prevalence',
        fontsize=8,
    )
    plot_file = os.path.join(plt_path, f'Prevalence_{N_max}_N_vs_p_error_Diff_{str(diff).replace('.','-')}')
    plt.savefig(plot_file+'.svg', dpi=150, format='svg')
    plt.savefig(plot_file+'.png', dpi=150, format='png')
    plt.show()




#%% Plot line plots for metrics across diffs for specific Ns

metrics_to_plot=['test_per_sample','groupsize_vs_N']
ylabels=['Test per sample','Relative group size']

N_to_plot=[20,40,60,80,99]

methods = list(full_df_met['Method'].unique())
n_methods = len(methods)
colors = [cmap(x) for x in np.linspace(0, 0.9, n_methods)]
color_dict = dict(zip(methods, colors))
rng = np.random.default_rng(seed=42)  # For reproducibility


# IF WANT to plot specific subsets of methods independently
methods_nonadaptive=['Shifted transversal','Chinese Remainder','Ch. Rem. Backtrack','Ch. Rem. Special']
methods_adaptive=['Hierarchical','Binary','Random','Matrix','multidim-3','multidim-4']
methods_split=[methods_nonadaptive,methods_adaptive]

for j,methods in enumerate(methods_split):
    
    n_methods = len(methods)
    if j==0:
        colors = [cmap(x) for x in np.linspace(0.1, 0.65, n_methods)]
        color_dict = dict(zip(methods, colors))
    elif j==1:
        colors = [cmap(x) for x in np.linspace(0.35, 0.9, n_methods)]
        color_dict = dict(zip(methods, colors))
    else: print('Need to re-work color split')
    
    for i in range(len(metrics_to_plot)):
        
        metric_to_plot=metrics_to_plot[i]
        ylabel=ylabels[i]
        
        x_col='N'
        xlabel='Number of positive samples'
        
        title=False
        
        
        #loop across diffs
        for N in N_to_plot:
        
            #subset df
            df_filtered=full_df_met[full_df_met['N']==N]
            
            fig, ax = plt.subplots(figsize=(4, 5))
            
            plt.text(0.05,0.9,'n = '+str(N), transform = ax.transAxes,fontsize=tick_fontsize)
            
            for method in methods:
                
                try:
                    grp = df_filtered[df_filtered['Method'] == method]
                    y_vals = grp[metric_to_plot].values
                    # Jitter x values
                    x_vals = grp['diff'].values
                    # Sort by jittered x for line plotting
                    sort_idx = np.argsort(x_vals)
                    x_sorted = x_vals[sort_idx]
                    y_sorted = y_vals[sort_idx]
                    # Line: through jittered points
    
                    ax.plot(
                        x_sorted, y_sorted,
                        color=color_dict[method],
                        alpha=line_alpha,
                        linewidth=line_width,
                        zorder=1,
                        label=method
                    )

                except:
                    continue
                
                              
            if metric_to_plot=='test_per_sample':
                ax.hlines(y=1, xmin=0, xmax=max_N, linewidth=1, color='lightgrey',linestyles='--')
                ax.axhspan(1, 1000, alpha=0.75,color='white')
                
                ax.set_ylim(ymin=0,ymax=2)
            
            else:
                ax.set_ylim(ymin=0,ymax=0.55)
            
            ax.set_xlabel(xlabel, fontsize=label_fontsize)
            ax.set_ylabel(ylabel, fontsize=label_fontsize)
            
            ax.set_xlim(xmin=min_diff-0.3,xmax=max_diff+0.3)

                
            ax.tick_params(axis='both', labelsize=tick_fontsize)
            
            
            #ax.spines[['top', 'right']].set_visible(False)
            for spine in ax.spines.values():
                spine.set_visible(False)
                
            plt.tight_layout()
            
            ax.legend(
                title='Method',
                fontsize=8,
            )
            
            
            suffix=''
            
            if diff>1 and keep_binary==False:
                suffix='_ExclBinary'
            
            save_path = os.path.join(
                plt_path, f"Diff_lines_{metric_to_plot.replace(' ','_')}_N_{N}_MethodGroup_{j}{suffix}"
            )
            
            
            if plt_path:
                plt.savefig(save_path+'.png', bbox_inches='tight',format='png',dpi=150)
                plt.savefig(save_path+'.svg', bbox_inches='tight',format='svg',dpi=150)
                
            plt.show()




#%% Plot line plots for prevalence across diffs for specific Ns

metrics_to_plot=['p_error']
ylabels=['Probability of error']


prevalences=df_prevalence['prevalence'].unique()
n_prevalence=len(prevalences)

colors_prevalence = [cmap(x) for x in np.linspace(0, 0.8, n_prevalence)]
color_dict_prevalence = dict(zip(prevalences, colors_prevalence))


N_to_plot=[20,40,60,80,99]


for i in range(len(metrics_to_plot)):
    
    metric_to_plot=metrics_to_plot[i]
    ylabel=ylabels[i]
    
    x_col='N'
    xlabel='Number of positive samples'
    
    title=False
    
    
    #loop across diffs
    for N in N_to_plot:
    
        #subset df
        df_filtered=df_prevalence[df_prevalence['N']==N]
        
        fig, ax = plt.subplots(figsize=(4, 5))
        
        plt.text(0.05,0.9,'n = '+str(N), transform = ax.transAxes,fontsize=tick_fontsize)
        
        for prevalence in prevalences:

            try:
                grp = df_filtered[df_filtered['prevalence'] == prevalence]
                y_vals = grp[metric_to_plot].values
                # Jitter x values
                x_vals = grp['Diff'].values
                # Sort by jittered x for line plotting
                sort_idx = np.argsort(x_vals)
                x_sorted = x_vals[sort_idx]
                y_sorted = y_vals[sort_idx]
                # Line: through jittered points

                ax.plot(
                    x_sorted, y_sorted,
                    color=color_dict_prevalence[prevalence],
                    alpha=line_alpha,
                    linewidth=line_width,
                    zorder=1,
                    label=prevalence
                )

            except:
                continue
            
                          
        # if metric_to_plot=='test_per_sample':
        #     ax.hlines(y=1, xmin=0, xmax=max_N, linewidth=1, color='lightgrey',linestyles='--')
        #     ax.axhspan(1, 1000, alpha=0.75,color='white')
            
        #     ax.set_ylim(ymin=0,ymax=2)
        
        # else:
            ax.set_ylim(ymin=0,ymax=1)
        
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        
        ax.set_xlim(xmin=min_diff-0.3,xmax=max_diff+0.3)

            
        ax.tick_params(axis='both', labelsize=tick_fontsize)
        
        
        #ax.spines[['top', 'right']].set_visible(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
            
        plt.tight_layout()
        
        ax.legend(
            title='Prevalence',
            fontsize=8,
        )
        
        
        suffix=''
        
        if diff>1 and keep_binary==False:
            suffix='_ExclBinary'
        
        save_path = os.path.join(
            plt_path, f"Diff_lines_{metric_to_plot.replace(' ','_')}_N_{N}_MethodGroup_{j}{suffix}"
        )
        
        
        if plt_path:
            plt.savefig(save_path+'.png', bbox_inches='tight',format='png',dpi=150)
            plt.savefig(save_path+'.svg', bbox_inches='tight',format='svg',dpi=150)
            
        plt.show()




#%%





