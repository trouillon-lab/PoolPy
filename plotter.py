import sys
print(sys.executable)

import numpy as np
import pandas as pd
import os
import cmcrameri.cm as cmc
from datetime import date
import matplotlib.pyplot as plt
import argparse


def plot_with_custom_labels(
    full_df_met,
    y_col,
    cmap,
    diff_val=None,
    plot_diff_on_x=False,
    selected_N=None,
    xlabel='N',
    ylabel=None,
    title=None,
    label_fontsize=12,
    tick_fontsize=10,
    save_path=None,
    scatter_size=40,
    scatter_alpha=0.8
):
    """
    Plots a scatter plot grouped by 'Method' from a DataFrame, using a cmcrameri colormap,
    with options for x-axis as 'N' or 'diff', custom labels/font sizes, and saves the plot.
    """
    import numpy as np

    # Filter data based on x-axis mode
    if plot_diff_on_x:
        df_filtered = full_df_met[full_df_met['N'] == selected_N]
        x_col = 'diff'
        x_label = xlabel
    else:
        df_filtered = full_df_met[full_df_met['diff'] == diff_val]
        x_col = 'N'
        x_label = xlabel

    fig, ax = plt.subplots(figsize=(8, 6))
    methods = list(df_filtered['Method'].unique())
    n_methods = len(methods)
    colors = [cmap(x) for x in np.linspace(0, 1, n_methods)]
    color_dict = dict(zip(methods, colors))

    for method, grp in df_filtered.groupby('Method'):
        ax.scatter(
            grp[x_col], grp[y_col],
            label=method,
            color=color_dict[method],
            s=scatter_size,
            alpha=scatter_alpha
        )

    ax.set_xlabel(x_label, fontsize=label_fontsize)
    ax.set_ylabel(ylabel if ylabel else y_col, fontsize=label_fontsize)
    if title:
        ax.set_title(title, fontsize=label_fontsize + 2)
    ax.tick_params(axis='both', labelsize=tick_fontsize)
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.tight_layout()

    # Place legend outside the plot
    ax.legend(
        title='Method',
        fontsize=label_fontsize,
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        borderaxespad=0.
    )

    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def plot_all_combinations(
    full_df_met,
    y_as,
    x_as,
    cmap,
    xlabel_dict=None,
    base_title=None,
    label_fontsize=12,
    tick_fontsize=10,
    save_dir='plots',
    scatter_size=40,
    scatter_alpha=0.8
):
    """
    Plots all combinations of y_as vs. x_as ('N' or 'diff'), for all unique values of diff/N.
    Saves plots in subfolders 'N/' and 'diff/' under save_dir, with clear filenames.

    Parameters:
        full_df_met (pd.DataFrame): The data source.
        y_as (list of str): List of column names for y-axis.
        x_as (list of str): List containing 'N', 'diff', or both, specifying x-axes to use.
        cmap (matplotlib colormap): Colormap to use for methods.
        xlabel_dict (dict, optional): Mapping of 'N' and 'diff' to custom x-axis labels.
        base_title (str, optional): Base plot title.
        label_fontsize (int): Font size for axis labels and legend.
        tick_fontsize (int): Font size for tick labels.
        save_dir (str): Directory to save plots (creates 'N' and 'diff' subfolders).
        scatter_size (int): Size of scatter points.
        scatter_alpha (float): Alpha (transparency) of scatter points (0.0 to 1.0).
    """
    if xlabel_dict is None:
        xlabel_dict = {'N': 'N', 'diff': 'diff'}

    # Create subdirectories for 'N' and 'diff'
    n_dir = os.path.join(save_dir, 'N')
    diff_dir = os.path.join(save_dir, 'diff')
    os.makedirs(n_dir, exist_ok=True)
    os.makedirs(diff_dir, exist_ok=True)

    for x_mode in x_as:
        if x_mode == 'N':
            unique_diffs = full_df_met['diff'].unique()
            for diff_val in unique_diffs:
                for y_col in y_as:
                    save_path = os.path.join(
                        n_dir, f"{y_col.replace(' ', '_')}_diff_{diff_val}.png"
                    )
                    title = (f"{y_col} vs N (diff={diff_val})"
                             if base_title is None else f"{base_title} (diff={diff_val})")
                    plot_with_custom_labels(
                        full_df_met=full_df_met,
                        y_col=y_col,
                        cmap=cmap,
                        diff_val=diff_val,
                        plot_diff_on_x=False,
                        selected_N=None,
                        xlabel=xlabel_dict.get('N', 'N'),
                        ylabel=y_col,
                        title=title,
                        label_fontsize=label_fontsize,
                        tick_fontsize=tick_fontsize,
                        save_path=save_path,
                        scatter_size=scatter_size,
                        scatter_alpha=scatter_alpha
                    )
        elif x_mode == 'diff':
            unique_Ns = full_df_met['N'].unique()
            for n_val in unique_Ns:
                for y_col in y_as:
                    save_path = os.path.join(
                        diff_dir, f"{y_col.replace(' ', '_')}_N_{n_val}.png"
                    )
                    title = (f"{y_col} vs diff (N={n_val})"
                             if base_title is None else f"{base_title} (N={n_val})")
                    plot_with_custom_labels(
                        full_df_met=full_df_met,
                        y_col=y_col,
                        cmap=cmap,
                        diff_val=None,
                        plot_diff_on_x=True,
                        selected_N=n_val,
                        xlabel=xlabel_dict.get('diff', 'diff'),
                        ylabel=y_col,
                        title=title,
                        label_fontsize=label_fontsize,
                        tick_fontsize=tick_fontsize,
                        save_path=save_path,
                        scatter_size=scatter_size,
                        scatter_alpha=scatter_alpha
                    )
        else:
            raise ValueError(f"Unknown x axis: {x_mode}. Use 'N' and/or 'diff'.")


def plotter(dir_WAs, max_diff, min_diff, start=0, stop=np.inf, step=1, x_axis='both', y_axis='all', save_df=False, **kwargs):
    N=start
    #ls_names_met=['Method', 'Mean experiments', 'Max compunds per well', 'N wells', 'Percentage check', 'Mean extra experiments', 'Mean steps', 'N', 'diff' ]
    
    #full_df_met=pd.DataFrame(columns=ls_names_met)

    ls_metrics=[]
    print(stop)
    while N<stop:
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
            df_met=pd.read_csv(metname, header=0)
            df_met['N']=N
            df_met['diff']=diff

            ls_metrics.append(df_met)
            

            diff+=1


        N+=step
    
        
    full_df_met=pd.concat(ls_metrics)
    
    print(full_df_met.shape)

    if save_df:
        full_df_met.to_csv(metname)

    if x_axis=='both':
        x_as=['N', 'diff']
    else:
        x_as=[x_axis]

    if y_axis=='all':
        y_as=['Mean experiments', 'Max compunds per well', 'N wells', 'Percentage check', 'Mean extra experiments', 'Mean steps']
    else:
        y_as=[y_axis]

    today = date.today()

    plt_path=os.path.join(dir_WAs,'plots',str(today))

    plot_all_combinations(
        full_df_met=full_df_met,
        y_as=y_as,
        x_as=x_as,
        cmap=cmc.batlow,
        xlabel_dict={'N': 'Sample Size', 'diff': 'Maximun number of positives'},
        base_title="Method Comparison",
        label_fontsize=14,
        tick_fontsize=12,
        save_dir=plt_path,
        scatter_size=80,      # Larger points
        scatter_alpha=0.5     # More transparent
    )
        


parser = argparse.ArgumentParser(description="Run the plotter utility for experiment metrics.")

parser.add_argument("--dir_WAs", type=str, default=' ./',  help="Root directory containing N_x/diff_y folders.")
parser.add_argument("--max_diff", type=float,  default=3, help="Maximum diff value to process.")
parser.add_argument("--min_diff", type=float, default=1, help="Minimum diff value to process.")

parser.add_argument("--start", type=int, default=0, help="Start value for N (default: 0)")
parser.add_argument("--stop", type=float, default=np.inf, help="Stop value for N (default: infinity)")
parser.add_argument("--step", type=int, default=1, help="Step for N (default: 1)")
parser.add_argument("--x_axis", type=str, default='both', choices=['N', 'diff', 'both'], help="Which x axis to use (default: both)")
parser.add_argument("--y_axis", type=str, default='all', help="Which y axis to use (default: all)")
parser.add_argument("--save_df", action='store_true', help="Whether to save the combined DataFrame (default: False)")

# Parse arguments and convert to dictionary
args = parser.parse_args()
args_dict = vars(args)

# Call your function using the dictionary of arguments
plotter(**args_dict)



    
                

