import numpy as np
import pandas as pd
from . import data_utils

import seaborn as sns; sns.set()
import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib import rcParams

plt.rcParams['figure.figsize'] = [15, 10]

font = {'size'   : 15, }
axes = {'labelsize': 'medium', 'titlesize': 'medium'}

sns.set_context("talk")
sns.set_style("white")
mpl.rc('font', **font)
mpl.rc('axes', **axes)
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
# plt.rcParams['text.usetex'] = True
plt.rcParams['axes.unicode_minus'] = False

import pandas as pd
from . import posterior_analysis
import matplotlib.style as style
import math

import subprocess
from . import data_plotting

import os
from glob import glob
from pathlib import Path
from . import motif_counting
from . import nearest_neighbours
from . import distance_analysis
from . import NMF_analysis
from . import system_summary

def merge_model_space_report_df_list(df_list):
    merge_func = lambda x, y, suff: pd.merge(x, y, on='model_idx', suffixes=suff)
    for i in range(1, len(df_list)):
        df_list[0] = merge_func(df_list[0], df_list[i], [str("_" + str(i-1)), str("_" + str(i))])
    model_space_report_df = df_list[0]

    return model_space_report_df


def generate_model_space_statistics(df, target_column_name):
    # Get appropriate columns to generate stdev
    column_names = list(df)
    target_cols = [col for col in column_names if target_column_name in col]

    model_stds = []
    model_means = []
    for idx, row in df.iterrows():
        model_stds.append(np.std(row[target_cols].values))
        model_means.append(np.mean(row[target_cols].values))

    df[target_column_name + "_std"] = model_stds
    df[target_column_name + "_mean"] = model_means

def generate_marginal_probability_boxplot(pop_dir_list, output_dir, hide_x_ticks=True, show_median=True, show_BF=False, drop_eqless=-1):
    print("Generating marginal probability distribution... ")

    model_space_report_list = []
    num_pops = len(pop_dir_list)

    model_space_report_path = pop_dir_list + "model_space_report_long.csv"
    model_space_report_df = pd.read_csv(model_space_report_path, index_col=0)
    model_space_report_df['model_marginal_median'] = np.nan

    # Calculate median for each model
    model_idxs = model_space_report_df['model_idx'].unique()
    median_marginals = []
    for m_idx in model_idxs:
        sub_df = model_space_report_df.loc[model_space_report_df['model_idx'] == m_idx]
        m_med_marg = np.median(sub_df['model_marginal'].values)
        model_space_report_df.loc[(model_space_report_df['model_idx'] == m_idx), 'model_marginal_median'] = m_med_marg
        median_marginals.append(m_med_marg)
    

    model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_median'] <= drop_eqless].index, inplace=True)


    # Sort model indexes by median marginal descending order
    median_marginals, model_idxs = [list(x) for x in zip(*sorted(zip(median_marginals, model_idxs), key=lambda pair: pair[0], reverse=True))]
    unique_model_idxs = model_space_report_df['model_idx'].unique()
    model_idxs = [x for x in model_idxs if x in unique_model_idxs]

    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    sns.boxplot(x='model_idx', y='model_marginal', order=model_idxs, data=model_space_report_df, ax=ax, dodge=False, showfliers=False)

    ax.unicode_minus = True

    if hide_x_ticks:
        ax.set(xticklabels=[])
        ax.set(xlabel='')
        ax.legend().remove()    
    
    else:
        ax.set(xticklabels=model_space_report_df['model_idx'])
        ax.set(xlabel='Model')
        ax.legend()


    ax.set_ylabel('')
    ax.set(xlim=(-0.5,None))
    ax.set(ylim=(-0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.tick_params(labelsize=15)
    ax.margins(x=0)
    ax.margins(y=0)
    fig.tight_layout()
    plt.show()
    # model_space_report_df.drop(model_space_report_df.filter(regex="Unname"),axis=1, inplace=True)
    # model_space_report_df = model_space_report_df.sort_values(by='model_idx', ascending=False).reset_index(drop=True)
    # print(model_space_report_df.columns)
    # model_space_report_list.append(model_space_report_df)


def generate_marginal_probability_distribution(pop_dir_list, output_dir, hide_x_ticks=True, show_median=True, show_BF=False, drop_eqless=-1):
    print("Generating marginal probability distribution... ")

    model_space_report_list = []
    num_pops = len(pop_dir_list)

    for data_dir in pop_dir_list:
        # Load model space report
        model_space_report_path = data_dir + "model_space_report.csv"
        model_space_report_df = pd.read_csv(model_space_report_path, index_col=0)
        model_space_report_df.drop(model_space_report_df.filter(regex="Unname"),axis=1, inplace=True)
        model_space_report_df = model_space_report_df.sort_values(by='model_idx', ascending=False).reset_index(drop=True)
        print(model_space_report_df.columns)
        model_space_report_list.append(model_space_report_df)

    model_space_report_df = merge_model_space_report_df_list(model_space_report_list)
    print(model_space_report_df.columns)
    generate_model_space_statistics(model_space_report_df, "model_marginal")

    # model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] == 0].index, inplace=True)
    model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] <= drop_eqless].index, inplace=True)
    print(model_space_report_df)

    model_space_report_df = model_space_report_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)

    output_path = output_dir + "model_marginal_probability.pdf"
    
    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))

    bar_color = sns.xkcd_palette(["amber"])[0]
    bar_color = "#333333"

    sns.barplot(model_space_report_df.index, model_space_report_df.model_marginal_mean, linewidth=1.0,
                     data=model_space_report_df, alpha=0.9, ax=ax, dodge=False, color=bar_color) #, palette="BuGn_r"
    # sns.barplot(model_space_report_df.index, model_space_report_df.model_marginal_mean, scale=0.1, join=False, linewidth=1.0,
    # data=model_space_report_df, alpha=0.9, ax=ax, dodge=False, color=bar_color) #, palette="BuGn_r"

    # ax.stackplot(model_space_report_df.index, model_space_report_df.model_marginal_mean)
    if num_pops > 1:
        ax.errorbar(model_space_report_df.index, 
                    model_space_report_df['model_marginal_mean'], 
                    yerr=model_space_report_df['model_marginal_std'], fmt='None', color='black', alpha=1,
                    label=None, elinewidth=0.1)

    ax.unicode_minus = True

    if hide_x_ticks:
        ax.set(xticklabels=[])
        ax.set(xlabel='')
        ax.legend().remove()    
    
    else:
        ax.set(xticklabels=model_space_report_df['model_idx'])
        ax.set(xlabel='Model')
        ax.legend()

    if show_median:
        median = np.median(model_space_report_df['model_marginal_mean'].values)
        ax.axhline(median, ls='--', label='median', linewidth=1.0)
        ax.legend()

    if show_BF:
        BF = model_space_report_df.iloc[0]['model_marginal_mean'] / 3.0
        ax.axhline(BF, ls='--', label='Bayes\' Factor = 3.0', linewidth=1.0)
        ax.legend().remove()

    ax.set_ylabel('')
    ax.set(xlim=(-0.5,None))
    ax.set(ylim=(-0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.tick_params(labelsize=15)
    ax.margins(x=0)
    ax.margins(y=0)

    fig.tight_layout()
    plt.savefig(output_path, dpi=500, bbox_inches='tight')

    output_path = output_dir + "model_marginal_probability_log_scale.eps"
    ax.set(yscale="log")
    ax.set(ylim=(None, None))
    plt.savefig(output_path, dpi=500)

    print(output_path)
    print("\n")

def get_total_simulations(rep_dirs):
    total_sims = 0

    for rep in rep_dirs:
        sub_dirs = glob(rep + "*/")
        pop_dirs = [f for f in sub_dirs if "Population" in f.split('/')[-2]]
        pop_dirs = sorted(pop_dirs, key=lambda a: a[-2])

        for pop in pop_dirs:
            model_space_report_path = pop + "model_space_report.csv"
            df = pd.read_csv(model_space_report_path)
            total_sims += sum(df['simulated_count'].values)

    return total_sims

def plot_all_model_param_distributions(pop_dir, inputs_dir, figure_output_dir):
    data_utils.make_folder(figure_output_dir)
    model_space_report_path = pop_dir + "combined_model_space_report.csv"
    model_space_report_df = pd.read_csv(model_space_report_path, index_col=0)
    # generate_model_space_statistics(model_space_report_df, "model_marginal")
    print(list(model_space_report_df))

    model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] == 0].index, inplace=True)

    for model_idx in model_space_report_df['model_idx'].values:
        model_output_parameter_path = pop_dir + "combined_model_params/" + "model_" + str(model_idx) + "_population_all_params"
        model_input_parameter_path = inputs_dir + "input_files/params_" + str(model_idx) + ".csv"
        model_input_species_path = inputs_dir + "input_files/species_" + str(model_idx) + ".csv"

        R_script = os.path.dirname(os.path.realpath(__file__)) + "/dens_plot_2D.R"

        subprocess.call(['Rscript', R_script, model_output_parameter_path, model_input_parameter_path, model_input_species_path, str(model_idx), figure_output_dir, 'TRUE', 'TRUE'])



def population_analysis(pop_dir, inputs_dir):
    analysis_dir = pop_dir + "analysis/"
    data_utils.make_folder(analysis_dir)
    print(pop_dir)
    generate_marginal_probability_distribution([pop_dir], analysis_dir, hide_x_ticks=True, drop_unnacepted=True)
    write_model_order([pop_dir], analysis_dir)
    
    param_dists_dir = analysis_dir + "param_dists/"
    data_utils.make_folder(param_dists_dir)
    plot_all_model_param_distributions(pop_dir, inputs_dir, param_dists_dir)

def spock_vs_others(data_dir, adj_mat_dir, output_dir, drop_unnacepted=False, show_median=False):
    hide_x_ticks = False

    model_space_report_path = data_dir + "combined_model_space_report.csv"
    model_space_report_df = pd.read_csv(model_space_report_path)
        
    model_idxs = [79, 65, 129, 130, 131]
    model_labels = ['Spock_v1', 'Spock_v2', 'Balagadde', 'Scott', 'Mccardell']

    # Make subset of only the best models

    sub_model_space_df = model_space_report_df.loc[model_space_report_df['model_idx'].isin(model_idxs)]
    sub_model_space_df = sub_model_space_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)

    ordered_labels = []
    for x in sub_model_space_df['model_idx'].values:
        for y, label in zip(model_idxs, model_labels):
            if y == x:
                ordered_labels.append(label)

    print(ordered_labels)
    sub_model_space_df['names'] = ordered_labels

    # Sort data frame in order of highest acceptance ratio to lowest

    # Plot 
    output_path = output_dir + "spock_vs_others.pdf"
    
    fig, ax = plt.subplots()

    ax.errorbar(sub_model_space_df.index, 
                sub_model_space_df['model_marginal_mean'], 
                yerr=sub_model_space_df['model_marginal_std'], fmt=',', color='black', alpha=1,
                label=None, elinewidth=0.5)

    sns.barplot(x='names', y='model_marginal_mean', 
                     data=sub_model_space_df, alpha=0.9, ax=ax)
    ax.unicode_minus = True

    if hide_x_ticks:
        ax.set(xticklabels=[])
        ax.set(xlabel='')
        ax.legend().remove()

    else:
        ax.set(xticklabels=sub_model_space_df['names'])
        ax.set(xlabel='Model')
        ax.legend()

    if show_median:
        median = np.median(sub_model_space_df['model_marginal_mean'].values)
        ax.axhline(median, ls='--', label='Median', linewidth=1.0)
        ax.legend()


    ax.set(ylabel='Model marginal probability')
    ax.set(xlim=(-0.5,None))
    ax.set(ylim=(-0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    fig.tight_layout()
    plt.savefig(output_path, dpi=500)

    # output_path = output_dir + "model_marginal_probability_log_scale.eps"
    # ax.set(yscale="log")
    # ax.set(ylim=(None, None))
    # plt.savefig(output_path, dpi=500)

    print(output_path)
    print("\n")
    exit()


def compare_top_models_by_parts(data_dir, adj_mat_dir, output_dir, drop_eqless=-1, show_median=False, hide_x_ticks=False):
    model_space_report_path = data_dir + "combined_model_space_report.csv"
    model_space_report_df = pd.read_csv(model_space_report_path)

        # model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] == 0].index, inplace=True)
    model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] <= drop_eqless].index, inplace=True)


    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    all_num_parts = data_utils.make_num_expressed_parts(model_space_report_df, adj_matrix_path_template)
    # all_num_parts = data_utils.make_num_species(model_space_report_df, adj_matrix_path_template)

    # all_num_parts = data_utils.make_num_parts_alt(model_space_report_df, adj_matrix_path_template)
    model_space_report_df['num_parts'] = all_num_parts

    unique_num_parts = model_space_report_df['num_parts'].unique()
    unique_num_parts.sort()
    best_model_idxs = []

    # Get best model for each part class
    for num_parts in unique_num_parts:
        sub_model_space_df = model_space_report_df.loc[model_space_report_df['num_parts'] == num_parts]

        # Sort data frame in order of highest acceptance ratio to lowest
        sub_model_space_df = sub_model_space_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)
        best_model_idxs.append(sub_model_space_df.iloc[0]['model_idx'])

        # Make subset of only the best models
        sub_model_space_df = model_space_report_df.loc[model_space_report_df['model_idx'].isin(best_model_idxs)]

        # Sort data frame in order of highest acceptance ratio to lowest
        sub_model_space_df = sub_model_space_df.sort_values(by='num_parts', ascending=True).reset_index(drop=True)

    # Plot 
    output_path = output_dir + "best_models_marginal_by_parts.pdf"

    
    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))

    ax.errorbar(sub_model_space_df.index, 
                sub_model_space_df['model_marginal_mean'], 
                yerr=sub_model_space_df['model_marginal_std'], fmt=',', color='black', alpha=1,
                label=None, elinewidth=0.5)


    sub_model_space_df['model_idx'] = sub_model_space_df['model_idx'].astype('category')
    sub_model_space_df['index'] = sub_model_space_df.index

    bar_color = sns.xkcd_palette(["amber"])[0]

    bar_color = "#333333"

    sns.barplot(x=sub_model_space_df.index, y='model_marginal_mean', color=bar_color, linewidth=0,
                     data=sub_model_space_df, alpha=0.9, ax=ax)
    print(sub_model_space_df)
    ax.unicode_minus = True

    if hide_x_ticks:
        ax.set(xticklabels=[])
        ax.set(xlabel='')
        ax.legend().remove()

    else:
        print(hide_x_ticks)
        ax.set(xticklabels=sub_model_space_df['model_idx'])
        ax.set(xlabel='Model')
        ax.legend()

    if show_median:
        median = np.median(sub_model_space_df['model_marginal_mean'].values)
        ax.axhline(median, ls='--', label='Median', linewidth=1.0)
        ax.legend()

    ax.set(ylabel='')
    # ax.set(yticklabels=[])
    ax.set(xlim=(-0.5,None))
    ax.set(ylim=(-0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    # ax.spines["left"].set_visible(False)

    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    fig.tight_layout()
    plt.savefig(output_path, dpi=500, bbox_inches='tight')

    # output_path = output_dir + "model_marginal_probability_log_scale.eps"
    # ax.set(yscale="log")
    # ax.set(ylim=(None, None))
    # plt.savefig(output_path, dpi=500)

    print(output_path)
    print("\n")

    # best_model_idxs = [62, 48, 66]
    # Generate bayes factors
    model_posterior_probs = [sub_model_space_df.loc[sub_model_space_df['model_idx']==idx]['model_marginal_mean'].values[0] for idx in best_model_idxs]
    bayes_factor_mat = np.zeros([len(model_posterior_probs), len(model_posterior_probs)])

    B_ij = lambda p_m1, p_m2: p_m1/p_m2
    sig_func = lambda x: 1/ np.exp(-x)

    for i in range(len(best_model_idxs)):
        for j in range(len(best_model_idxs)):
            p_m1 = model_posterior_probs[i]
            p_m2 = model_posterior_probs[j]

            bayes_factor_mat[i, j] = (p_m1/p_m2) 

    print(bayes_factor_mat)
    print(bayes_factor_mat[0, 1])


def ammensal_vs_cooperative_systems(pop_dir_list, output_dir, adj_mat_dir, hide_x_ticks=True, drop_unnacepted=True, show_median=True):
    show_median = True

    model_space_report_list = []
    num_pops = len(pop_dir_list)

    for data_dir in pop_dir_list:
        # Load model space report
        model_space_report_path = data_dir + "model_space_report.csv"
        model_space_report_df = pd.read_csv(model_space_report_path, index_col=0)
        model_space_report_list.append(model_space_report_df)

    model_space_report_df = merge_model_space_report_df_list(model_space_report_list)
    generate_model_space_statistics(model_space_report_df, "model_marginal")

    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    sys_category = []
    sys_category_dummy = []
    for idx in model_space_report_df['model_idx'].values:
        adj_mat = pd.read_csv(adj_matrix_path_template.replace("#REF#", str(idx)))
        total_coop = sum(adj_mat['S_aa1']) + sum(adj_mat['S_aa2'])
        total_amen = abs(sum(adj_mat['B_1'])) + abs(sum(adj_mat['B_2']))

        print(adj_mat)
        print(list(adj_mat.columns).index('B_1'))
        N_1_idx = list(adj_mat.columns).index('N_1') - 1
        N_2_idx = list(adj_mat.columns).index('N_2') - 1
        S_aa1_idx =  list(adj_mat.columns).index('S_aa1') - 1
        S_aa2_idx =  list(adj_mat.columns).index('S_aa2') - 1
        B_1_idx = list(adj_mat.columns).index('B_1') - 1
        B_2_idx = list(adj_mat.columns).index('B_2') - 1

        adj_mat = adj_mat.drop(adj_mat.columns[0], axis=1).as_matrix()

        ammensal = False
        self_regulation = False
        cross_feed = False

        print(adj_mat[:, N_1_idx][B_1_idx:B_2_idx+1])
        exit()

        if total_coop > 0 and total_amen > 0:
            sys_category.append('Both')
            sys_category_dummy.append(1)

        elif total_coop > 0:
            sys_category.append('Cooperative only')
            sys_category_dummy.append(2)

        elif total_amen > 0:
            sys_category.append('Ammensal only')
            sys_category_dummy.append(3)


    model_space_report_df['sys_cat'] = sys_category
    model_space_report_df['sys_category_dummy'] = sys_category_dummy

    if drop_unnacepted:
        model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] == 0].index, inplace=True)

    model_space_report_df = model_space_report_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)

    output_path = output_dir + "coop_amen_model_marginal_probability.pdf"

    norm = plt.Normalize(min(model_space_report_df.sys_category_dummy.values), max(model_space_report_df.sys_category_dummy.values))
    cmap = plt.get_cmap("magma")

    fig, ax = plt.subplots()
    if num_pops > 1:
        ax.errorbar(model_space_report_df.index, 
                    model_space_report_df['model_marginal_mean'], 
                    yerr=model_space_report_df['model_marginal_std'], fmt=',', color='black', alpha=1,
                    label=None, elinewidth=0.5)

    ax.bar(model_space_report_df.index, model_space_report_df.model_marginal_mean, color=cmap(norm(model_space_report_df.sys_category_dummy.values)), 
        data=model_space_report_df, alpha=0.9)

    ax.unicode_minus = True
    if hide_x_ticks:
        ax.set(xticklabels=[])
        ax.set(xlabel='')
        ax.legend().remove()    
    
    else:
        ax.set(xticklabels=model_space_report_df['model_idx'])
        ax.set(xlabel='Model')
        ax.legend()

    if show_median:
        median = np.median(model_space_report_df['model_marginal_mean'].values)
        ax.axhline(median, ls='--', label='median', linewidth=1.0)
        ax.legend()


    ax.set(ylabel='Model marginal probability')
    ax.set(xlim=(-0.5,None))
    ax.set(ylim=(-0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    fig.tight_layout()
    plt.savefig(output_path, dpi=500)

    output_path = output_dir + "model_marginal_probability_log_scale.eps"
    ax.set(yscale="log")
    ax.set(ylim=(None, None))
    plt.savefig(output_path, dpi=500)

    print(output_path)
    print("\n")

def split_by_num_parts(data_dir, adj_mat_dir, output_dir, drop_unnacepted=False):
    hide_x_ticks = False
    show_median = True

    model_space_report_path = data_dir + "combined_model_space_report.csv"
    model_space_report_df = pd.read_csv(model_space_report_path)

    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    # all_num_parts, all_AHL_num_parts, all_microcin_num_parts = data_utils.make_num_parts(model_space_report_df, adj_matrix_path_template)
    all_num_parts, all_AHL_num_parts, all_microcin_num_parts = data_utils.make_num_parts(model_space_report_df, adj_matrix_path_template)

    model_space_report_df['num_parts'] = all_num_parts

    unique_num_parts = model_space_report_df['num_parts'].unique()

    # Acceptance rate
    for num_parts in unique_num_parts:
        sub_model_space_df = model_space_report_df.loc[model_space_report_df['num_parts'] == num_parts]
        valid_model_refs = sub_model_space_df['model_idx'].unique()

        if drop_unnacepted:
            sub_model_space_df.drop(sub_model_space_df[sub_model_space_df['model_marginal_mean'] == 0].index, inplace=True)
        
        if len(sub_model_space_df) == 0:
            continue

        # Sort data frame in order of highest acceptance ratio to lowest
        sub_model_space_df = sub_model_space_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)


        # Generate standard deviation
        output_path = output_dir + "model_marginal_probability_" + str(int(num_parts)) + "_parts.pdf"

        fig, ax = plt.subplots()

        ax.errorbar(sub_model_space_df.index, 
                    sub_model_space_df['model_marginal_mean'], 
                    yerr=sub_model_space_df['model_marginal_std'], fmt=',', color='black', alpha=1,
                    label=None, elinewidth=0.5)
        print(sub_model_space_df)

        sns.barplot(sub_model_space_df.index, sub_model_space_df.model_marginal_mean, 
                         data=sub_model_space_df, alpha=0.9, ax=ax)
        ax.unicode_minus = True

        if hide_x_ticks:
            ax.set(xticklabels=[])
            ax.set(xlabel='')
            ax.legend().remove()    
        
        else:
            ax.set(xticklabels=sub_model_space_df['model_idx'])
            ax.set(xlabel='Model')
            ax.legend()

        if show_median:
            median = np.median(sub_model_space_df['model_marginal_mean'].values)
            ax.axhline(median, ls='--', label='Median', linewidth=1.0)
            ax.legend()


        ax.set(ylabel='Model marginal probability')
        ax.set(xlim=(-0.5,None))
        ax.set(ylim=(-0))
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_alpha(0.5)
        ax.spines["bottom"].set_alpha(0.5)
        fig.tight_layout()
        plt.savefig(output_path, dpi=500)

        # output_path = output_dir + "model_marginal_probability_log_scale.eps"
        # ax.set(yscale="log")
        # ax.set(ylim=(None, None))
        # plt.savefig(output_path, dpi=500)

        print(output_path)
        print("\n")


        print("\n")

def self_regulation_bar_plot(pop_dir_list, adj_mat_dir, output_dir, hide_x_ticks=True, drop_unnacepted=False, show_median=False):
    model_space_report_list = []
    num_pops = len(pop_dir_list)

    for data_dir in pop_dir_list:
        # Load model space report
        model_space_report_path = data_dir + "model_space_report.csv"
        model_space_report_df = pd.read_csv(model_space_report_path, index_col=0)
        model_space_report_list.append(model_space_report_df)


    model_space_report_df = merge_model_space_report_df_list(model_space_report_list)
    generate_model_space_statistics(model_space_report_df, "model_marginal")
    
    # total_accepted = sum(model_space_report_df['accepted_count'].values)
    # new_model_marginals = [x/total_accepted for x in model_space_report_df['accepted_count'].values]
    # model_space_report_df['model_marginal_mean'] = new_model_marginals

    if drop_unnacepted:
        model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] == 0.0].index, inplace=True)

    print(len(model_space_report_df))
    model_space_report_df = model_space_report_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)
    output_path = output_dir + "model_marginal_probability_self_reg.pdf"

    model_idxs = model_space_report_df['model_idx'].values



    target_species_list = []
    # model_space_report_df['self_reg_loops'] = data_utils.get_self_regulators(model_idxs, target_species_list, adj_mat_dir)
    strain_loop_bal, sum_pos_loops, sum_neg_loops = data_utils.make_strain_feedback_loop_balance(model_idxs, adj_mat_dir)
    model_space_report_df['strain_loop_bal'] = strain_loop_bal
    model_space_report_df['sum_pos_loops'] = sum_pos_loops
    model_space_report_df['sum_neg_loops'] = sum_neg_loops

    model_space_report_df['mat_distance'] = data_utils.get_adj_mat_distances(model_idxs[0], model_idxs, adj_mat_dir)
    model_space_report_df['sum_interactions'] = data_utils.get_num_interactions(model_idxs, adj_mat_dir)

    model_space_report_df['symmetrical'] = [4 if x == True else 0 for x in  data_utils.get_two_strain_symmetric_adj_mats(model_idxs, adj_mat_dir)]
    symm_models = model_space_report_df.loc[model_space_report_df['symmetrical'] == 4]['model_idx'].values


    sym_col_names = []
    for m in symm_models:
        col_name = str(m) + "_dist"
        sym_col_names.append(col_name)
        model_space_report_df[col_name] = data_utils.get_adj_mat_distances(m, model_idxs, adj_mat_dir)
    
    min_symm_neighbour = []
    min_dist_val = []
    for m_idx in model_idxs:
        min_symm_vals = model_space_report_df.loc[model_space_report_df['model_idx'] == m_idx][sym_col_names].values[0]
        min_symm_neighbour.append(int(sym_col_names[np.argmin(min_symm_vals)].split('_')[0]))
        min_dist_val.append(min(model_space_report_df.loc[model_space_report_df['model_idx'] == m_idx][sym_col_names].values[0]))

    model_space_report_df['min_symm_neighbour'] = min_symm_neighbour
    model_space_report_df['min_symm_dist'] = min_dist_val
    # model_space_report_df.drop(model_space_report_df[model_space_report_df['143_dist'] > model_space_report_df['min_symm_dist']].index, inplace=True)
    # model_space_report_df.drop(model_space_report_df[model_space_report_df['143_dist'] > model_space_report_df['min_symm_dist']].index, inplace=True)

    adjusted_vals = []
    for idx, row in model_space_report_df.iterrows():
        if row['224_dist'] > row['min_symm_dist']:
            adjusted_vals.append(0)

        else:
            adjusted_vals.append(row['model_marginal_mean'])
    model_space_report_df['model_marginal_mean'] = adjusted_vals
    print(model_space_report_df)
    # if drop_unnacepted:
    # model_space_report_df.drop(model_space_report_df[model_space_report_df['symmetrical'] == False].index, inplace=True)
    # model_idxs = model_space_report_df['model_idx'].values
    # model_space_report_df = model_space_report_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)
    model_space_report_df.to_csv(output_dir + "dists_test.csv")

    fig, ax = plt.subplots()

    print(model_space_report_df)
    if num_pops > 1:
        pass
        # ax.errorbar(model_space_report_df.index, 
        #             model_space_report_df['model_marginal_mean'], 
        #             yerr=model_space_report_df['model_marginal_std'], fmt=',', color='black', alpha=1,
        #             label=None, elinewidth=0.5)

    sns.barplot(model_space_report_df.index, model_space_report_df.model_marginal_mean, 
                     data=model_space_report_df, alpha=0.9, ax=ax)
    ax.unicode_minus = True

    if hide_x_ticks:
        ax.set(xticklabels=[])
        ax.set(xlabel='')
        ax.legend().remove()    
    
    else:
        ax.set(xticklabels=model_space_report_df['model_idx'])
        ax.set(xlabel='Model')
        ax.legend()

    if show_median:
        median = np.median(model_space_report_df['model_marginal_mean'].values)
        ax.axhline(median, ls='--', label='median', linewidth=1.0)
        ax.legend()


    ax2 = ax.twinx()
    ax2.set(ylim=(0, 5))
    # sns.scatterplot(x=model_space_report_df.index, y=model_space_report_df['194_dist'].values, ax=ax2)

    ax.set(ylabel='Model marginal probability')
    ax.set(xlim=(-0.5,None))
    ax.set(ylim=(-0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)


    fig.tight_layout()
    plt.savefig(output_path, dpi=500)

    output_path = output_dir + "model_marginal_probability_self_reg_log_scale.eps"
    ax.set(yscale="log")
    ax.set(ylim=(None, None))
    plt.savefig(output_path, dpi=500)

    print(output_path)
    print("\n")


def generate_critical_parameter_bar_plot(data_dir, KS_data_dir, output_dir, num_params):
    print("Generating critical parameter bar plot... ")

    # Load model space report 
    model_space_report_path = data_dir + "model_space_report.csv"
    model_space_report_df = pd.read_csv(model_space_report_path)
    

    model_indexes = model_space_report_df['model_idx'].values
    KS_data_path_template = KS_data_dir + "model_NUM_KS.csv"
    
    model_space_report_df['crit_param_1'] = np.nan

    for model_idx in (model_indexes):
        ks_data_path = KS_data_path_template.replace('NUM', str(int(model_idx)))

        try:
            df = pd.read_csv(ks_data_path, sep=',')
    
        except(FileNotFoundError, ValueError):
            continue

        param_names = df.columns[3:]
        row = df.iloc[0]
        param_KS_values = row[3:].values

        
        # Reorder parameters according to their KS value
        params_order = [param for _,param in sorted(zip(param_KS_values, param_names), reverse=True)]
        critical_parameters = [param for param in params_order]
        model_space_report_df.loc[model_idx, 'crit_param_1'] = critical_parameters[0]

    # Make plot!
    style.use('seaborn-muted')
    current_palette = sns.color_palette()

    fig, axes = plt.subplots(ncols=1)
    crit_params = [1]

    for idx, col_num in enumerate(crit_params):
        print(col_num)
        crit_param_col = 'crit_param_NUM'.replace('NUM', str(col_num))
        names = data_utils.translate_param_names(model_space_report_df[crit_param_col].value_counts().index)
        value_counts = model_space_report_df[crit_param_col].value_counts()

        value_counts = value_counts[:num_params]
        names = names[:num_params]
        
        sns.barplot(value_counts, names, ax=axes)
        
        axes.spines["right"].set_visible(False)
        axes.spines["top"].set_visible(False)
        axes.set_title('Rank ' + str(col_num) + ' critical parameter')
        axes.tick_params(labelsize=40)
        axes.set_xlabel('Count')

    plt.tight_layout()
    plt.savefig(output_dir+'rank_one_params.pdf', dpi=500)

    print("\n")


def write_model_order(pop_dir_list, output_dir):
    model_space_report_list = []
    for data_dir in pop_dir_list:
        # Load model space report
        model_space_report_path = data_dir + "model_space_report.csv"
        model_space_report_df = pd.read_csv(model_space_report_path, index_col=0)
        model_space_report_list.append(model_space_report_df)

    model_space_report_df = merge_model_space_report_df_list(model_space_report_list)
    generate_model_space_statistics(model_space_report_df, "model_marginal")

    model_space_report_df = model_space_report_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)
    
    file = open(output_dir + "model_order.txt", "w") 
    for item in list(model_space_report_df['model_idx'].values):
        file.write("%s\n" % item)
    file.close()

    print("\n")


def write_combined_model_space_report(pop_dir_list, output_dir):
    output_path = output_dir + "combined_model_space_report.csv"
    model_space_report_list = []
    for data_dir in pop_dir_list:
        # Load model space report
        model_space_report_path = data_dir + "model_space_report.csv"
        model_space_report_df = pd.read_csv(model_space_report_path, index_col=0)
        model_space_report_df.drop(model_space_report_df.filter(regex="Unname"),axis=1, inplace=True)
        model_space_report_df = model_space_report_df.sort_values(by='model_idx', ascending=False).reset_index(drop=True)
        print(model_space_report_df.columns)
        model_space_report_list.append(model_space_report_df)

    model_space_report_df = merge_model_space_report_df_list(model_space_report_list)
    generate_model_space_statistics(model_space_report_df, "model_marginal")

    model_space_report_df = model_space_report_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)
    model_space_report_df.to_csv(output_path)

    return model_space_report_df


def write_combined_model_space_with_motif_counts(pop_dir_list, output_dir, adj_mat_dir, window=0, normalise=True, plot=False):
    output_path = output_dir + "combined_model_space_report_with_motifs.csv"
    model_space_report_list = []

    for data_dir in pop_dir_list:
        # Load model space report
        model_space_report_path = data_dir + "model_space_report.csv"
        model_space_report_df = pd.read_csv(model_space_report_path, index_col=0)
        model_space_report_df.drop(model_space_report_df.filter(regex="Unname"),axis=1, inplace=True)
        model_space_report_df = model_space_report_df.sort_values(by='model_idx', ascending=False).reset_index(drop=True)
        model_space_report_list.append(model_space_report_df)

    model_space_report_df = merge_model_space_report_df_list(model_space_report_list)
    generate_model_space_statistics(model_space_report_df, "model_marginal")
    model_space_report_df = model_space_report_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)

    motif_columns = ['permissive_counts', 'dependent_counts',
       'submissive_counts', 'hedonistic_counts',
       'defensive_counts', 'logistic_counts', 'opportunistic_counts',
       'exponential_counts']

    motif_columns = ['SL1', 'SL2', 'SL3', 'SL4', 
    'OL1', 'OL2', 'OL3', 'OL4']

    # model_space_report_df = motif_counting.count_direct_self_limiting(model_space_report_df, adj_mat_dir, window=window, normalise=normalise)
    model_space_report_df = motif_counting.count_permissive(model_space_report_df, adj_mat_dir, window=window, normalise=normalise)
    model_space_report_df = motif_counting.count_dependent(model_space_report_df, adj_mat_dir, window=window, normalise=normalise)
    model_space_report_df = motif_counting.count_submissive(model_space_report_df, adj_mat_dir, window=window, normalise=normalise)
    model_space_report_df = motif_counting.count_hedonistic(model_space_report_df, adj_mat_dir, window=window, normalise=normalise)

    # model_space_report_df = motif_counting.count_direct_competitive(model_space_report_df, adj_mat_dir, window=window, normalise=normalise)
    model_space_report_df = motif_counting.count_defensive(model_space_report_df, adj_mat_dir, window=window, normalise=normalise)
    model_space_report_df = motif_counting.count_logistic(model_space_report_df, adj_mat_dir, window=window, normalise=normalise)
    model_space_report_df = motif_counting.count_opportunistic(model_space_report_df, adj_mat_dir, window=window, normalise=normalise)
    model_space_report_df = motif_counting.count_exponential(model_space_report_df, adj_mat_dir, window=window, normalise=normalise)

    model_marginal_means = model_space_report_df['model_marginal_mean'].values
    normalised_marginal_means = motif_counting.normalise_list(model_marginal_means)
    model_space_report_df['norm_marginal_means'] = normalised_marginal_means
    
    model_space_report_df.to_csv(output_path)

    if plot:
        output_path = output_dir + "motif_counts_plot.pdf"
        fig, ax = plt.subplots()

        # sns.barplot(x='name', y='mean', 
        #                  data=analysis_df, alpha=0.9, ax=ax, palette=diverging_colours, vert=False)

        # analysis_df.plot("name", "mean", kind="barh", color=diverging_colours, ax=ax, title='', width=1)
        x = range(len(model_space_report_df['defensive_rolling'].values))
        sns.lineplot(x=x, y='dependent_rolling', data=model_space_report_df, label='dependent', ax=ax, size=2)
        # sns.scatterplot(x=x, y='permissive_rolling', data=model_space_report_df,label='permissive', ax=ax, size=2)
        sns.lineplot(x=x, y='defensive_rolling', data=model_space_report_df, label='defensive', ax=ax, size=2)
        sns.lineplot(x=x, y='exponential_rolling', data=model_space_report_df, label='exponential', ax=ax, size=2)

        # sns.scatterplot(x=x, y='dependent_rolling', data=model_space_report_df, ax=ax)
        # sns.scatterplot(x=x, y='dependent_rolling', data=model_space_report_df, ax=ax)

        ax.set_xlabel('index')
        ax.set_ylabel('rolling motif count')


        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_visible(True)

        ax.spines["bottom"].set_alpha(0.5)
        ax.spines["left"].set_alpha(0.5)

        ax.legend()

        fig.tight_layout()
        plt.savefig(output_path, dpi=500)



    return model_space_report_df


def write_experiment_summary(population_size, n_repeats, final_pop_dirs, inputs_dir, output_dir):
    # Number of repeats
    # Number of models
    # Population size
    # Total number of simulations

    exp_rep_dirs = []

    for d in final_pop_dirs:
        split_dir = d.split('/')

        # Take step back in directory
        exp_rep_dirs.append("/".join(split_dir[0:-2]) + "/")


    total_sims = get_total_simulations(exp_rep_dirs)

    model_space_report_df = pd.read_csv(output_dir + "combined_model_space_report.csv")
    n_models = len(model_space_report_df['model_idx'].values)
    n_dead_models = len(model_space_report_df.loc[model_space_report_df['model_marginal_mean'] == 0])

    file = open(output_dir + "experiment_summary.txt", "w") 
    file.write("Number of repeats: " + str(n_repeats) + "\n")
    file.write("Number of models: " + str(n_models) + "\n")
    file.write("Number of dead models: " + str(n_dead_models) + "\n")
    file.write("Population size: " + str(population_size) + "\n")
    file.write("Total simulations: " + str(total_sims) + "\n")
    file.close()


def find_finished_experiments(data_dir):
    repeats_exp_dirs = list(glob(data_dir + "**/"))

    final_pop_dirs = []


    for exp_dir in repeats_exp_dirs:
        sub_dirs = glob(exp_dir + "*/")
        pop_dirs = [f for f in sub_dirs if "Population" in f.split('/')[-2]]

        if len(pop_dirs) == 0:
            continue
        pop_dirs = sorted(pop_dirs, key=lambda a: a[-2])

        final_pop_dirs.append(pop_dirs[-1])

    finished_pop_dirs = []
    for pop in final_pop_dirs:
        model_space_report = Path(pop + "model_space_report.csv")
        if model_space_report.is_file():
            finished_pop_dirs.append(pop)

        else:
            print(pop)


    return finished_pop_dirs

def find_latest_pop_dirs(data_dir):
    repeats_exp_dirs = list(glob(data_dir + "**/"))
    final_pop_dirs = []


    for exp_dir in repeats_exp_dirs:
        sub_dirs = glob(exp_dir + "*/")
        pop_dirs = [f for f in sub_dirs if "Population" in f.split('/')[-2]]

        if len(pop_dirs) == 0:
            continue
        pop_dirs = sorted(pop_dirs, key=lambda a: a[-2])

        final_pop_dirs.append(pop_dirs[-1])
    return final_pop_dirs

def combine_model_params(pop_dirs, output_dir):
    data_utils.make_folder(output_dir + "combined_model_params")

    # Get all model numbers
    model_space_report_df = pd.read_csv(output_dir + "combined_model_space_report.csv")
    model_idxs = model_space_report_df['model_idx'].values
    model_idxs = sorted(model_idxs)
    model_params_template = "model_#IDX#_population_all_params"

    for idx in model_idxs:
        first_instance = True
        master_df = None
        model_params_file_name = model_params_template.replace('#IDX#', str(idx))

        for pop in pop_dirs:
            model_params_file_path = pop + "model_sim_params/" + model_params_file_name

            try:
                model_params_df = pd.read_csv(model_params_file_path)

                if first_instance:
                    master_df = model_params_df
                    first_instance = False

                else:
                    master_df = pd.concat([master_df, model_params_df])

            except(FileNotFoundError):
                continue

        print("Saving model params idx: ", idx)

        # Normalise particle weights
        if len(master_df) > 1:
            max_weight = max(master_df['particle_weight'].values)
            min_weight = min(master_df['particle_weight'].values)

            if max_weight != min_weight:
                min_max_scale = lambda x: (x - min_weight) / (max_weight - min_weight)
                new_weights = [min_max_scale(x) for x in master_df['particle_weight'].values]
                master_df['particle_weight'] = new_weights


        master_df.to_csv(output_dir + "combined_model_params/" + model_params_file_name)

def combine_model_distances(pop_dirs, output_dir):
    data_utils.make_folder(output_dir + "combined_model_distances/")

    model_space_report_df = pd.read_csv(output_dir + "combined_model_space_report.csv")
    model_idxs = model_space_report_df['model_idx'].values
    model_idxs = sorted(model_idxs)
    model_distances_template = "model_#IDX#_population_distances.csv"

    # Open all distance files
    distance_dfs = []
    print("Loading model distance dfs")
    for pop in pop_dirs:
        pop_distances_file_path = pop + "distances.csv"
        dist_df = pd.read_csv(pop_distances_file_path)
        print(dist_df.columns)
        distance_dfs.append(dist_df)


    for idx in model_idxs:
        df_list = []
        model_distance_file_name = model_distances_template.replace('#IDX#', str(idx))

        for df in distance_dfs:
            sub_df = df.loc[df['model_ref'] == idx]
            sub_df = sub_df.loc[sub_df['Accepted'] == True]

            df_list.append(sub_df)

        model_distances_df = pd.concat(df_list, ignore_index=True) 
        model_distances_df.to_csv(output_dir + "combined_model_distances/" + model_distance_file_name)


def combine_model_distances_individually(pop_dirs, output_dir):
    data_utils.make_folder(output_dir + "combined_model_distances/")

    model_space_report_df = pd.read_csv(output_dir + "combined_model_space_report.csv")
    model_idxs = model_space_report_df['model_idx'].values
    model_idxs = sorted(model_idxs)
    model_distances_template = "model_#IDX#_distances.csv"

    # Open all distance files

    for idx in model_idxs:
        first_instance = True
        master_df = None
        model_distance_file_name = model_distances_template.replace('#IDX#', str(idx))

        for pop in pop_dirs:
            model_distance_path = pop + "model_sim_distances/" + model_distance_file_name

            try:
                model_dist_df = pd.read_csv(model_distance_path)

                if first_instance:
                    master_df = model_dist_df
                    first_instance = False

                else:
                    master_df = pd.concat([master_df, model_dist_df])

            except(FileNotFoundError):
                continue

        print("Saving model distance idx: ", idx)

        if not first_instance:
            master_df.to_csv(output_dir + "combined_model_distances/" + model_distance_file_name)

