import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

import matplotlib as mpl
from . import network_vis

font = {'size'   : 5, }
axes = {'labelsize': 'medium', 'titlesize': 'medium'}

sns.set_context("talk")
sns.set_style("white")
mpl.rc('font', **font)
mpl.rc('axes', **axes)
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
# plt.rcParams['text.usetex'] = True
plt.rcParams['axes.unicode_minus'] = False

def count_bac(adj_mat_df):
    col_names = adj_mat_df.columns

    microcin_indexes = [i for i in col_names if 'B_' in i]

    num_bac = 0

    for bac in microcin_indexes:
        x = adj_mat_df[bac].values
        x = sum(abs(x))

        if x > 0:
            num_bac += 1

    return num_bac

def count_QS(adj_mat_df):
    col_names = adj_mat_df.columns

    AHL_indexes = [i for idx, i in enumerate(col_names) if 'A_' in i]
    microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

    num_QS = 0

    for QS in AHL_indexes:
        x = adj_mat_df[QS].values
        x = sum(abs(x))

        if x > 0:
            num_QS += 1

    return num_QS

def regulation_categorisation(adj_mat_df):
    col_names = adj_mat_df.columns

    AHL_indexes = [i for idx, i in enumerate(col_names) if 'A_' in i]
    microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

    num_QS = 0
    ind_count = 0
    rep_count = 0

    for QS in AHL_indexes:
        QS_interactions = adj_mat_df[QS].values

        if np.sum(abs(QS_interactions)) > 0:
            for x in QS_interactions:
                if x > 0:
                    ind_count += 1

                elif x < 0:
                    rep_count += 1

    return ind_count, rep_count

def make_bayes_factors(data_list):
    BF_arr = np.zeros((len(data_list), len(data_list)))

    for idx_x, x in enumerate(data_list):
        for idx_y, y in enumerate(data_list):
            BF_arr[idx_x, idx_y] = x/y

    return BF_arr


def model_space_summary(combined_analysis_output_dir, adj_mat_dir, OL_only=False):
    model_space_report_path = combined_analysis_output_dir + "combined_model_space_report_with_motifs.csv"
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"

    if OL_only:
        QS_output_path = combined_analysis_output_dir + "OL_only_QS_subset_means.pdf"
        B_output_path = combined_analysis_output_dir + "OL_only_B_subset_means.pdf"
        reg_output_path = combined_analysis_output_dir + "OL_only_reg_subset_means.pdf"

    else:
        QS_output_path = combined_analysis_output_dir + "QS_subset_means.pdf"
        B_output_path = combined_analysis_output_dir + "B_subset_means.pdf"
        reg_output_path = combined_analysis_output_dir + "reg_subset_means.pdf"

    model_space_report_df = pd.read_csv(model_space_report_path)
    if OL_only:
        model_space_report_df['node_type'] = model_space_report_df.apply(network_vis.model_type_conditional, axis=1)
        model_space_report_df = model_space_report_df.loc[model_space_report_df['node_type'] == 'red']

    model_space_report_df = model_space_report_df.sort_values('model_idx')

    model_idxs = model_space_report_df['model_idx'].values
    
    num_QS_counts = []
    num_bac_counts = []
    num_ind_counts = []
    num_rep_counts = []

    for m_idx in model_idxs:

        model_self_limiting_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)

        model_QS_count = count_QS(adj_mat_df)
        model_bac_count = count_bac(adj_mat_df)
        ind_count, rep_count = regulation_categorisation(adj_mat_df)

        num_ind_counts.append(ind_count)
        num_rep_counts.append(rep_count)

        num_QS_counts.append(model_QS_count)
        num_bac_counts.append(model_bac_count)

    model_space_report_df['num_QS'] = num_QS_counts
    model_space_report_df['num_bac'] = num_bac_counts
    model_space_report_df['num_ind'] = num_ind_counts
    model_space_report_df['num_rep'] = num_rep_counts

    # models with 1 QS
    sub_1_QS_models = model_space_report_df.loc[model_space_report_df['num_QS'] == 1.0]
    QS_1_mean_pprob = np.mean(sub_1_QS_models['model_marginal_mean'].values)
    QS_1_std = np.std(sub_1_QS_models['model_marginal_mean'].values)

    # models with 2 QS
    sub_2_QS_models = model_space_report_df.loc[model_space_report_df['num_QS'] == 2.0]
    QS_2_mean_pprob = np.mean(sub_2_QS_models['model_marginal_mean'].values)
    QS_2_std = np.std(sub_2_QS_models['model_marginal_mean'].values)

    # models with 1 B
    sub_1_B_models = model_space_report_df.loc[model_space_report_df['num_bac'] == 1.0]
    B_1_mean_pprob = np.mean(sub_1_B_models['model_marginal_mean'].values)
    B_1_std = np.std(sub_1_B_models['model_marginal_mean'].values)

    # models with 2 B
    sub_2_B_models = model_space_report_df.loc[model_space_report_df['num_bac'] == 2.0]
    B_2_mean_pprob = np.mean(sub_2_B_models['model_marginal_mean'].values)
    B_2_std = np.std(sub_2_B_models['model_marginal_mean'].values)

    # models with 2 B
    sub_3_B_models = model_space_report_df.loc[model_space_report_df['num_bac'] == 3.0]
    B_3_mean_pprob = np.mean(sub_3_B_models['model_marginal_mean'].values)
    B_3_std = np.std(sub_2_B_models['model_marginal_mean'].values)

    # Models with induction only
    sub_ind_models = model_space_report_df.loc[model_space_report_df['num_ind'] > 0.0]
    sub_ind_models = sub_ind_models.loc[sub_ind_models['num_rep'] == 0.0]
    ind_pprob = np.mean(sub_ind_models['model_marginal_mean'].values)


    # Models with rep only
    sub_rep_models = model_space_report_df.loc[model_space_report_df['num_rep'] > 0.0]
    sub_rep_models = sub_rep_models.loc[sub_rep_models['num_ind'] == 0.0]
    rep_pprob = np.mean(sub_rep_models['model_marginal_mean'].values)

    # Models with mixed rep ind
    sub_mixed_models = model_space_report_df.loc[model_space_report_df['num_ind'] > 0.0]
    sub_mixed_models = sub_mixed_models.loc[sub_mixed_models['num_rep'] > 0.0]
    mixed_pprob = np.mean(sub_mixed_models['model_marginal_mean'].values)

    print(len(sub_ind_models))
    print(len(sub_rep_models))
    print(len(sub_mixed_models))

    QS_summary_data = [QS_1_mean_pprob, QS_2_mean_pprob]
    reg_summary_data = [ind_pprob, rep_pprob, mixed_pprob]
    B_summary_data = [B_1_mean_pprob, B_2_mean_pprob, B_3_mean_pprob]

    QS_BF = make_bayes_factors(QS_summary_data)
    reg_BFs = make_bayes_factors(reg_summary_data)
    B_BFs = make_bayes_factors(B_summary_data)

    print("QS BFs: \n", QS_BF)
    print("")

    print("reg BFs: \n", reg_BFs)
    print("")

    print("B BFs: \n", B_BFs)
    print("")

    width_inches = 110/2*2 / 25.4
    height_inches = 45*2 / 25.4

    bar_color = "#333333"

    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    
    
    sns.barplot(x=[0, 1], y=QS_summary_data, color=bar_color, linewidth=0, alpha=0.9, ax=ax)

    ax.set(xticklabels=[])
    ax.set(yticklabels=[])

    ax.set(xlabel='')
    ax.legend().remove()    

    ax.set(xlim=(-0.5,None))
    ax.set(ylim=(None, 10**-3))
    ax.set_yscale('log')

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    fig.tight_layout()

    plt.savefig(QS_output_path, dpi=500)

    fig, ax = plt.subplots(figsize=(width_inches, height_inches))    
    sns.barplot(x=[0, 1, 2], y=B_summary_data, color=bar_color, linewidth=0, alpha=0.9, ax=ax)

    ax.set(xticklabels=[])
    ax.set(yticklabels=[])

    ax.set(xlabel='')
    ax.legend().remove()    

    ax.set(xlim=(-0.5,None))
    ax.set(ylim=(None, 10**-3))
    # ax.set_yscale('log')

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    fig.tight_layout()

    plt.savefig(B_output_path, dpi=500)


    fig, ax = plt.subplots(figsize=(width_inches, height_inches))

    sns.barplot(x=[0, 1, 2], y=reg_summary_data, color=bar_color, linewidth=0, alpha=0.9, ax=ax)

    ax.set(xticklabels=[])
    ax.set(yticklabels=[])

    ax.set(xlabel='')
    ax.legend().remove()    

    ax.set(xlim=(-0.5,None))
    ax.set(ylim=(None, 10**-3))
    # ax.set_yscale('log')

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    fig.tight_layout()

    plt.savefig(reg_output_path, dpi=500)



