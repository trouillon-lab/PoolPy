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

from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The Axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)
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


# Get datadraws


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
        
        ylim_max=max(df_filtered[metric_to_plot])*1.05
        
        # fig, ax = plt.subplots(figsize=(8, 5))
        fig, ax = plt.subplots(figsize=(7, 5))
        
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
                # ax.annotate(xy=(x_sorted[-1],y_sorted[-1]), xytext=(5,0), textcoords='offset points',
                #             text=method, va='center',color=color_dict[method],fontsize=tick_fontsize,alpha=0.7)
                ax.annotate(xy=(x_sorted[-1],y_sorted[-1]), xytext=(20,0), textcoords='offset points',
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
            ylim_max=1.1
        
        ax.set_xlabel(xlabel, fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
        
        ax.set_xlim(xmin=0)
        ax.set_ylim(ymin=0,ymax=ylim_max)
        
        if title:
            ax.set_title(title, fontsize=label_fontsize + 2)
        ax.tick_params(axis='both', labelsize=tick_fontsize)
        
        
        #ax.spines[['top', 'right']].set_visible(False)
        # for spine in ax.spines.values():
        #     spine.set_visible(False)
            
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
    
#metric_to_plot='Max samples per pool'
metric_to_plot='groupsize_vs_N'

x_col='test_per_sample'

ylabel='Relative group size'
xlabel='Test per sample'


title=False


methods = list(full_df_met['Method'].unique())
n_methods = len(methods)
colors = [cmap(x) for x in np.linspace(0, 0.9, n_methods)]
color_dict = dict(zip(methods, colors))
rng = np.random.default_rng(seed=42)  # For reproducibility

#loop across diffs
for diff in range(min_diff,max_diff+1):

    #subset df
    df_filtered=full_df_met[full_df_met['diff']==diff]
    df_filtered=df_filtered[df_filtered['N'].isin([100])]
    
    fig, ax = plt.subplots(figsize=(5.5, 5))
    
    plt.text(0.05,0.9,'Differentiate '+str(diff), transform = ax.transAxes,fontsize=tick_fontsize)
    
    for method in methods:
        
        #Filter binary methods above diff 1
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

            ax.scatter(
                x_sorted, y_sorted,s=200,
                color=color_dict[method],
                edgecolor='black',
                alpha=0.8,
                zorder=1
            )
            # ax.plot(
            #     x_sorted, y_sorted,
            #     color=color_dict[method],
            #     alpha=0.3,
            #     linewidth=line_width,
            #     zorder=1
            # )
            
            # ax.annotate(xy=(x_sorted[-1],y_sorted[-1]), xytext=(15,0), textcoords='offset points',
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
    
    if metric_to_plot=='groupsize_vs_N':
        ax.vlines(x=1, ymin=0, ymax=1, linewidth=1, color='lightgrey',linestyles='--')
        ax.axvspan(1, 1000, alpha=0.75,color='white')
    
        
    ax.set_xlabel(xlabel, fontsize=label_fontsize)
    ax.set_ylabel(ylabel, fontsize=label_fontsize)
    
    ax.set_xlim(xmin=0,xmax=1.2)
    ax.set_ylim(ymin=0,ymax=0.55)
    
    if title:
        ax.set_title(title, fontsize=label_fontsize + 2)
    ax.tick_params(axis='both', labelsize=tick_fontsize)
    
    
    #ax.spines[['top', 'right']].set_visible(False)
    # for spine in ax.spines.values():
    #     spine.set_visible(False)
        
    plt.tight_layout()
    
    # Place legend outside the plot
    # ax.legend(
    #     title='Method',
    #     fontsize=label_fontsize,
    #     bbox_to_anchor=(1.05, 1),
    #     loc='upper left',
    #     borderaxespad=0.
    # )
    
    
    save_path = os.path.join(
        plt_path, f"{metric_to_plot.replace(' ','_')}_N_{max_N}_diff_{diff}"
    )
    
    
    if plt_path:
        plt.savefig(save_path+'_dualplot.png', bbox_inches='tight',format='png',dpi=150)
        plt.savefig(save_path+'_dualplot.svg', bbox_inches='tight',format='svg',dpi=150)
        
    plt.show()
    

#%% Plot test/sample vs group size  OVER ALL DIFFs   


#metric_to_plot='Mean steps'
#metric_to_plot='groupsize_vs_N'
metrics_to_plot_=['groupsize_vs_N','Mean steps']

for metric_to_plot in metrics_to_plot_:

    N_to_plot=100
    
    x_col='test_per_sample'
    
    ylabel='Relative group size'
    xlabel='Test per sample'
    
    marker_size_circle=200
    marker_size_square=120
    
    
    if metric_to_plot=='groupsize_vs_N':
        ylim_max=0.55
        ylim_min=0
        ylabel='Relative group size'
        arrow_head_size=0.02
    elif metric_to_plot=='Mean steps':
        ylim_max=5
        ylim_min=0
        ylabel='Number of steps'
        arrow_head_size=0.04
    
    title=False
    
    #Plot only specific methods ?
    methods=list(full_df_met['Method'].unique())
    methods_list=methods.copy()
    #methods=['Shifted transversal']
    
    #subset df to take only N to plot
    df_filtered_alldiff=full_df_met[full_df_met['N'].isin([N_to_plot])].copy()
    
    #get markers
    df_filtered_alldiff['marker_style']="o" #default
    df_filtered_alldiff['marker_style']=df_filtered_alldiff['marker_style'].astype('str')
    df_filtered_alldiff['marker_size']=marker_size_circle #default
    
    #Add marker info for number of steps
    # for idx in df_filtered_alldiff.index:
    #     if df_filtered_alldiff.loc[idx,'Mean steps']>1:
    #         df_filtered_alldiff.loc[idx,'marker_style']='s'
    #         df_filtered_alldiff.loc[idx,'marker_size']=marker_size_square
            
    #Alternatively, use marker for diff
    for idx in df_filtered_alldiff.index:
        if df_filtered_alldiff.loc[idx,'diff']>1:
            if df_filtered_alldiff.loc[idx,'diff']==2:
                df_filtered_alldiff.loc[idx,'marker_style']='s'
                df_filtered_alldiff.loc[idx,'marker_size']=marker_size_square
            elif df_filtered_alldiff.loc[idx,'diff']==3:
                df_filtered_alldiff.loc[idx,'marker_style']='^'
            elif df_filtered_alldiff.loc[idx,'diff']==4:
                df_filtered_alldiff.loc[idx,'marker_style']='p'
    
    n_methods = len(methods)
    colors = [cmap(x) for x in np.linspace(0, 0.9, n_methods)]
    color_dict = dict(zip(methods, colors))
    
    # for method in methods_list:
    #     methods=[method]
    
    fig, ax = plt.subplots(figsize=(5.5, 5))
    
    #Plot arrows if more than one diff
    if max_diff>1:
        for method in methods:
        
            #subset df
            df_filtered=df_filtered_alldiff[df_filtered_alldiff['Method']==method]
        
            for diff in range(min_diff,max_diff+1):
                
                try:
                    x0=df_filtered[df_filtered['diff']==diff][x_col].values[0]
                    y0=df_filtered[df_filtered['diff']==diff][metric_to_plot].values[0]
                    x1=df_filtered[df_filtered['diff']==diff+1][x_col].values[0]
                    y1=df_filtered[df_filtered['diff']==diff+1][metric_to_plot].values[0]
                    dx=(x1-x0)*0.9
                    dy=(y1-y0)*0.9
                    ax.arrow(x0,y0,dx,dy,
                              shape='full', color=color_dict[method], lw=1, length_includes_head=True, 
                              zorder=0, head_length=arrow_head_size, head_width=arrow_head_size, alpha=0.3)
                    
                except:
                    continue
    
    
    #Plot markers
    for diff in range(min_diff,max_diff+1):
    
        #subset df
        df_filtered=df_filtered_alldiff[df_filtered_alldiff['diff']==diff]
    
        for method in methods:
            
            # #Filter binary methods above diff 1
            # if diff>1 and method=='Binary': ######################################################################### Binary removed here
            #     continue
            
            try:
                grp = df_filtered[df_filtered['Method'] == method]
                markers_to_plot=grp['marker_style'].astype('string').values[0]
    
                ax.scatter(
                    grp[x_col], grp[metric_to_plot],s=grp['marker_size'],
                    color=color_dict[method],
                    edgecolor='black',
                    marker=markers_to_plot,
                    alpha=(0.8/diff)+0.2,
                    zorder=1
                )
                
                # ax.annotate(xy=(grp[x_col].values[-1],grp[metric_to_plot].values[-1]), xytext=(15,0), textcoords='offset points',
                #             text=method, va='center',color=color_dict[method],fontsize=tick_fontsize,alpha=0.7)
                
            except:
                continue
        
    if x_col=='test_per_sample':
        ax.axvspan(1, 1000, alpha=0.75,color='white')
        ax.vlines(x=1, ymin=0, ymax=10, linewidth=1, color='lightgrey',linestyles='--',zorder=1)
    
        
    ax.set_xlabel(xlabel, fontsize=label_fontsize)
    ax.set_ylabel(ylabel, fontsize=label_fontsize)
    
    ax.set_xlim(xmin=0,xmax=1.2)
    ax.set_ylim(ymin=ylim_min,ymax=ylim_max)
    
    
    if title:
        ax.set_title(title, fontsize=label_fontsize + 2)
    ax.tick_params(axis='both', labelsize=tick_fontsize)
    
    
    #ax.spines[['top', 'right']].set_visible(False)
    # for spine in ax.spines.values():
    #     spine.set_visible(False)
        
    plt.tight_layout()
    
    save_path = os.path.join(
        plt_path, f"{metric_to_plot.replace(' ','_')}_N_{max_N}_diff_{diff}"
    )
    
    
    if plt_path:
        plt.savefig(save_path+'_dualplot_allDiffs.png', bbox_inches='tight',format='png',dpi=150)
        plt.savefig(save_path+'_dualplot_allDiffs.svg', bbox_inches='tight',format='svg',dpi=150)
        
    plt.show()

#%% Plot FITS of test/sample vs group size  OVER ALL DIFFs   

N_to_plot=100

#metric_to_plot='Max samples per pool'
metric_to_plot='groupsize_vs_N'

x_col='test_per_sample'

ylabel='Relative group size'
xlabel='Test per sample'

marker_size_circle=200
marker_size_square=120

title=False

diffs=[x for x in range(min_diff,max_diff+1)]
n_diff=len(diffs)

colors = [cmap(x) for x in np.linspace(0, 0.6, n_diff)]
color_dict = dict(zip(diffs, colors))



#subset df to take only N to plot
df_filtered_alldiff=full_df_met[full_df_met['N'].isin([N_to_plot])].copy()

fig, ax = plt.subplots(figsize=(5.5, 5))


#Plot fits
for diff in range(min_diff,max_diff+1):
    
    #subset df
    df_filtered=df_filtered_alldiff[df_filtered_alldiff['diff']==diff]
    
    x=df_filtered[x_col]
    y=df_filtered[metric_to_plot]
    
    if diff==1: 
        markers_to_plot='o'
        marker_size=marker_size_circle
    elif diff==2: 
        markers_to_plot='s'
        marker_size=marker_size_square
    elif diff==3: 
        markers_to_plot='^'
        marker_size=marker_size_circle
    elif diff==4: 
        markers_to_plot='p'
        marker_size=marker_size_circle
    
    # ax.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),
    #         color=color_dict[diff],
    #         linewidth=3,
    #         zorder=1
    #         )
    
    ax.scatter(
        x, y,
        edgecolor='black',
        color=color_dict[diff],
        marker=markers_to_plot,
        s=marker_size/2,
        alpha=0.3,
    )
    
    #Plot confidence ellipse over n_std standard deviations
    confidence_ellipse(x, y, ax, n_std=1 , edgecolor=color_dict[diff],linewidth=2)
    
if metric_to_plot=='groupsize_vs_N':
    ax.vlines(x=1, ymin=0, ymax=1, linewidth=1, color='lightgrey',linestyles='--')
    ax.axvspan(1, 1000, alpha=0.75,color='white')

    
ax.set_xlabel(xlabel, fontsize=label_fontsize)
ax.set_ylabel(ylabel, fontsize=label_fontsize)

ax.set_xlim(xmin=0,xmax=1.2)
ax.set_ylim(ymin=0,ymax=0.55)


if title:
    ax.set_title(title, fontsize=label_fontsize + 2)
ax.tick_params(axis='both', labelsize=tick_fontsize)


#ax.spines[['top', 'right']].set_visible(False)
# for spine in ax.spines.values():
#     spine.set_visible(False)
    
plt.tight_layout()

save_path = os.path.join(
    plt_path, f"{metric_to_plot.replace(' ','_')}_N_{max_N}_diff_{diff}"
)


if plt_path:
    plt.savefig(save_path+'_dualplot_FIT_allDiffs.png', bbox_inches='tight',format='png',dpi=150)
    plt.savefig(save_path+'_dualplot_FIT_allDiffs.svg', bbox_inches='tight',format='svg',dpi=150)
    
plt.show()
    
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
        
        
        # if plt_path:
        #     plt.savefig(save_path+'.png', bbox_inches='tight',format='png',dpi=150)
        #     plt.savefig(save_path+'.svg', bbox_inches='tight',format='svg',dpi=150)
            
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
    # plt.savefig(plot_file+'.svg', dpi=150, format='svg')
    # plt.savefig(plot_file+'.png', dpi=150, format='png')
    plt.show()




#%% Plot line plots for metrics across diffs for specific Ns

metrics_to_plot=['test_per_sample','groupsize_vs_N']
ylabels=['Test per sample','Relative group size']

N_to_plot=[20,40,60,80,100]

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
            
            fig, ax = plt.subplots(figsize=(4, 4))
            
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
                
                ax.set_ylim(ymin=0,ymax=1.1)
            
            else:
                ax.set_ylim(ymin=0,ymax=0.55)
            
            ax.set_xlabel(xlabel, fontsize=label_fontsize)
            ax.set_ylabel(ylabel, fontsize=label_fontsize)
            
            ax.set_xlim(xmin=min_diff-0.3,xmax=max_diff+0.3)

                
            ax.tick_params(axis='both', labelsize=tick_fontsize)
            
            
            #ax.spines[['top', 'right']].set_visible(False)
            # for spine in ax.spines.values():
            #     spine.set_visible(False)
                
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


N_to_plot=[20,40,60,80,100]


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
        
        fig, ax = plt.subplots(figsize=(4, 4))
        
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
        # for spine in ax.spines.values():
        #     spine.set_visible(False)
            
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





