import numpy as np
from sklearn.decomposition import NMF

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import matplotlib as mpl
from matplotlib import rcParams
import itertools

import csv
from scipy import stats

from sklearn.decomposition.nmf import _beta_divergence 
from . import data_utils

def convert_QS_column(adj_mat):
    col_names = adj_mat.columns
    A_columns = [x for x in col_names if 'A_' in x]

    # Iterate AHLs and split into inducing and repressing
    for x in A_columns:
        col_vals = adj_mat[x]

        ind_vals = [0.0 if x != 1 else 1.0 for x in col_vals ]
        rpr_vals = [0.0 if x != -1 else 1.0 for x in col_vals ]

        ind_col = x + "_ind"
        rpr_col = x + "_rpr"

        adj_mat[ind_col] = ind_vals
        adj_mat[rpr_col] = rpr_vals

        adj_mat.drop(x, inplace=True, axis=1)

    return adj_mat


def write_H_to_csv(H, out_path, row_names, column_names, normalise=True, W_max_list=None):
    with open(out_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        
        for base_element in H:
            # median_base_element = np.median(base_element)
            # base_element = base_element - median_base_element

            # Normalise columns
            if normalise:
                max_base_element = np.max(base_element)
                min_base_element = np.min(base_element)

                min_max_normalise = lambda x: (x - min_base_element) / (max_base_element - min_base_element)

                base_element = min_max_normalise(base_element)
                base_element = np.around(base_element, decimals=3)

                for idx, _ in enumerate(base_element):
                    col_mode = stats.mode(base_element[:, idx])
                    base_element[:, idx] = base_element[:, idx] - col_mode[0]

            writer.writerow( ['-']  + column_names.tolist())
            for idx, i in enumerate(base_element):
                i = i.tolist()
                i = [row_names[idx]] + i
                writer.writerow(i)

            writer.writerows("\n")
            writer.writerows("\n")

def write_W_to_csv(W_array, model_idxs, out_path):
    W_element_col_names = ['V_' + str(i) for i in range(np.shape(W_array)[1])]

    df = pd.DataFrame(data=W_array, columns=W_element_col_names)

    df['model_idx'] = model_idxs
    
    df.to_csv(out_path)

def nmf_motif_count_decomposition(output_dir):
    model_space_report_df = pd.read_csv(output_dir + "combined_model_space_report_with_motifs.csv")

    # model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] <= 0.000].index, inplace=True)
    model_space_report_df = model_space_report_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)

    model_idxs = model_space_report_df['model_idx'].values

    all_motif_data = [0 for x in range(len(model_idxs))]

    motif_columns = ['permissive_counts', 'dependent_counts',
       'submissive_counts', 'hedonistic_counts',
       'defensive_counts', 'logistic_counts', 'opportunistic_counts',
       'exponential_counts']

    for i, m_idx in enumerate(model_idxs):
        model_sub_df = model_space_report_df.loc[model_space_report_df['model_idx'] == m_idx]
        motif_data = model_sub_df[motif_columns].values
        motif_data = motif_data.reshape(1, -1)
        print(motif_data)

        all_motif_data[i] = motif_data

    for i in range(2, 9):
        n_features = all_motif_data[0].shape[1]
        n_samples = len(model_idxs)

        motif_data = np.array(all_motif_data).reshape(n_samples, n_features)

        model = NMF(n_components=i, init='random', verbose=1, max_iter=1000000, tol=4e-16)

        W = np.array(model.fit_transform(motif_data))

        for idx in range(i):
            W_col_max = np.max(W[:, idx])
            W[:, idx] = W[:, idx] / W_col_max

        H = model.components_
        H = np.around(H, decimals=3)

        H = H.reshape(i, 1, n_features)

        write_H_to_csv(H, "H_motif_" + str(i) + ".csv", normalise=False)

        data_dict = {}
        data_dict['model_idx'] = model_space_report_df['model_idx']
        data_dict['model_marginal'] = model_space_report_df['model_marginal_mean'].values

        df = pd.DataFrame(data_dict)

        width_inches = 200 / 25.4
        height_inches = 150 / 25.4
        fig, ax = plt.subplots(figsize=(width_inches, height_inches))

        # sns.scatterplot(x='X', y='Y', hue='model_marginal', cmap=cmap, data=data_dict, ax=ax, size=0.5)
        # print(np.argmax(W[:, 0]))
        ax2 = ax.twinx()
        sns.heatmap(W.T, ax=ax, cbar=False)
        sns.scatterplot(y='model_marginal', x=df.index, data=df, ax=ax2, s=1, alpha=0.5, color='white')
        ax2.legend().remove()

        ax.set_xlabel('Models ordered high to low Pprob')
        # ax.set_xlabel('mean marginal probability change')
        ax.set_xticklabels('')
        # ax.set_yticklabels([])
        ax.set_ylabel('Components')

        plt.savefig('NMF_motif_' + str(i) + '_comp.pdf')
        # plt.show()
        plt.close()


        # # H = np.array(H).reshape(i, -1)
        # print(np.shape(H))
        # print(np.shape(W))
        # print(H[0].reshape(1, -1), np.argmax(H[0].reshape(1, -1)))
        # print(H[1].reshape(1, -1), np.argmax(H[1].reshape(1, -1)))
        # print(H[2].reshape(1, -1), np.argmax(H[2].reshape(1, -1)))
        # print(H[3].reshape(1, -1), np.argmax(H[3].reshape(1, -1)))
        # print(H[4].reshape(1, -1), np.argmax(H[4].reshape(1, -1)))
        # print(H[5].reshape(1, -1), np.argmax(H[5].reshape(1, -1)))
        # print(H[6].reshape(1, -1), np.argmax(H[6].reshape(1, -1)))
        # print(H[7].reshape(1, -1), np.argmax(H[7].reshape(1, -1)))

        # # np.multiply(W[0], H[:, 0])

        # exit()


def measure_sparseness(H_mat):
    # Flatten H
    H_mat = np.copy(H_mat).reshape(1, -1)
    H_mat = np.around(H_mat, decimals=3)

    # Count non zero elements
    num_elements = np.shape(H_mat)[1]
    non_zero_elements = np.count_nonzero(H_mat)

    sparsity = non_zero_elements / num_elements

    return sparsity



def optimize_NMF_rank_fuv(data, n_samples, plot_output_dir, train_size=0.8, k_min_max=[2, 30]):
    k_range = range(k_min_max[0], k_min_max[1])
    k_fuv_dict = {}
    k_fuv_dict['rep'] = []
    k_fuv_dict['k'] = []
    k_fuv_dict['fuv_vals'] = []
    k_fuv_dict['error_variance'] = []
    k_fuv_dict['non_zero_ratio'] = []
    # k_fuv_dict['sparsity_var_ratio'] = []
    k_fuv_dict['SS_err'] = []

    group_dict = {}
    group_dict['rep_err_var'] = []
    group_dict['reconstruct_X_test'] = []
    group_dict['X_test_flat'] = []
    group_dict['k'] = []


    n_repeats = 15

    for k in k_range:
        # group_dict['X_test_flat'] = []
        # group_dict['reconstruct_X_test'] = []
        for rep in range(n_repeats):

            # Generate test and train data
            model_indexes = list(range(n_samples))
            train_indexes = np.random.choice(model_indexes, size=int(n_samples*train_size), replace=False)
            test_indexs = [i for i in model_indexes if i not in train_indexes]

            X_test = np.copy(data[test_indexs])
            X_train = np.copy(data[train_indexes])

            # perturb_mat = np.random.normal(0.0, scale=10, size=np.shape(X_train))

            # X_train = np.random.normal(0.0, scale=10, size=np.shape(X_train))
            X_train = abs(X_train)

            # Apply speckled mask
            mask = np.random.choice([0, 1], size=X_train.shape, p=[0.2, 0.8]).astype(np.bool)
            # mask = np.random.randint(0,2,size=X_train.shape, weights=[0.2, 0.8]).astype(np.bool)
            # print(mask)
            r = np.zeros(X_train.shape)
            
            X_train[mask] = r[mask]

            model = NMF(n_components=k, init='nndsvda', verbose=0, max_iter=100, tol=4e-18, l1_ratio=1).fit(X_train)

            # Transform test set and reconstruct
            W_test = model.transform(X_test)
            reconstruct_X_test = model.inverse_transform(W_test).reshape(1, -1)

            # reconstruct_X_test = np.round(reconstruct_X_test, decimals=0)

            # Flatten elements
            X_test_flat = np.copy(X_test).reshape(1, -1)
            X_test_mean = np.mean(X_test_flat)

            SS_err = np.sum((X_test_flat - reconstruct_X_test)**2)
            SS_tot = np.sum((X_test_flat - X_test_mean)**2)
            fuv = SS_err / SS_tot

            error_variance = np.mean(SS_err)

            sparsity = measure_sparseness(model.components_)

            k_fuv_dict['rep'].append(rep)
            k_fuv_dict['k'].append(k)
            k_fuv_dict['fuv_vals'].append(fuv)
            k_fuv_dict['error_variance'].append(error_variance)
            k_fuv_dict['SS_err'].append(SS_err)
            k_fuv_dict['non_zero_ratio'].append(sparsity)
            # k_fuv_dict['sparsity_var_ratio'].append( error_variance / sparsity)

            # group_dict['X_test_flat'].extend(X_test_flat)
            # group_dict['reconstruct_X_test'].extend(reconstruct_X_test)

        # X_test_mean = np.mean(group_dict['X_test_flat'])
        # X_test_flat = np.array(group_dict['X_test_flat']).reshape(1, -1)
        # reconstruct_X_test = np.array(group_dict['reconstruct_X_test']).reshape(1, -1)

        # SS_err = np.sum((X_test_flat - reconstruct_X_test)**2)
        # group_dict['rep_err_var'].append(SS_err)
        # group_dict['k'].append(k)
        print(k)


    df = pd.DataFrame(k_fuv_dict)


    group_dict.pop('X_test_flat', None)
    group_dict.pop('reconstruct_X_test', None)
    group_df = pd.DataFrame(group_dict)

    df.to_csv('fuv_vals.csv')

    width_inches = 200 / 25.4
    height_inches = 150 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    sns.lineplot(x='k', y='SS_err', data=df)
    # sns.lineplot(x='k', y='error_variance', data=df)
    plt.savefig(plot_output_dir + 'NMF_optim_SS_err.pdf', dpi=500, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    sns.lineplot(x='k', y='error_variance', data=df)
    # sns.lineplot(x='k', y='error_variance', data=df)
    plt.savefig(plot_output_dir + 'NMF_optim_err_var.pdf', dpi=500, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    sns.lineplot(x='k', y='non_zero_ratio', data=df)
    # sns.lineplot(x='k', y='error_variance', data=df)
    plt.savefig(plot_output_dir + 'NMF_sparsity.pdf', dpi=500, bbox_inches='tight')
    plt.close()



    # fig, ax = plt.subplots(figsize=(width_inches, height_inches))

    # sns.lineplot(x='k', y='rep_err_var', data=group_df)
    # # sns.lineplot(x='k', y='error_variance', data=df)
    # plt.savefig(plot_output_dir + 'NMF_optim_err_var.pdf', dpi=500, bbox_inches='tight')
    # plt.close()


def nmf_decomposition(output_dir, adj_mat_dir, n_components=4):
    figure_output_dir = output_dir + 'nmf_analysis/'
    data_utils.make_folder(figure_output_dir)


    model_space_report_df = pd.read_csv(output_dir + "combined_model_space_report_with_motifs.csv")

    # model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] <= 0.000].index, inplace=True)
    model_space_report_df = model_space_report_df.sort_values(by='model_marginal_mean', ascending=False).reset_index(drop=True)

    model_idxs = model_space_report_df['model_idx'].values
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"

    all_adj_mats = [0 for x in range(len(model_idxs))]

    for i, m_idx in enumerate(model_idxs):
        model_self_limiting_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat = pd.read_csv(adj_mat_path, index_col=0)
        adj_mat = convert_QS_column(adj_mat)
        row_names = adj_mat.index.values
        column_names = adj_mat.columns.values

        adj_mat = abs(adj_mat.values)
        adj_mat = adj_mat.reshape(1, -1)
        all_adj_mats[i] = adj_mat


    n_features = all_adj_mats[0].shape[1]
    print(n_features)
    n_samples = len(model_idxs)
    model_adj_mats = np.array(all_adj_mats).reshape(n_samples, -1)

    model = NMF(n_components=n_components, init='random', random_state=2, verbose=0, max_iter=100000, tol=4e-18)


    W = np.array(model.fit_transform(model_adj_mats))

    for idx in range(n_components):
        W_col_max = np.max(W[:, idx])
        W[:, idx] = W[:, idx] / W_col_max
    
    W_csv_path = figure_output_dir +  "W_" + str(n_components) + ".csv"
    write_W_to_csv(W, model_idxs, W_csv_path)

    H = model.components_
    # H = H.reshape(n_components, 9, 11)
    H = H.reshape(n_components, 7, 9)

    H_csv_path = figure_output_dir +  "H_" + str(n_components) + ".csv"
    write_H_to_csv(H, H_csv_path, row_names, column_names, normalise=True)

    data_dict = {}
    data_dict['model_idx'] = model_space_report_df['model_idx']
    data_dict['model_marginal'] = model_space_report_df['model_marginal_mean'].values
    df = pd.DataFrame(data_dict)

    vmin, vmax = 0, 1
    cmap = sns.dark_palette("palegreen", as_cmap=True)

    width_inches = 175 / 25.4
    height_inches = 150 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))

    # sns.scatterplot(x='X', y='Y', hue='model_marginal', cmap=cmap, data=data_dict, ax=ax, size=0.5)
    # print(np.argmax(W[:, 0]))
    sns.heatmap(W.T, ax=ax, cbar=False)
    # sns.scatterplot(y='model_marginal', x=df.index, data=df, ax=ax2, s=1, alpha=0.5, color='white', orient='V')

    ax.set_xlabel('')
    # ax.set_xlabel('mean marginal probability change')
    ax.set_xticklabels('')
    ax.set_yticklabels('')

    # ax.set_yticklabels([])
    ax.set_ylabel('')
    ax.margins(x=0)
    ax.margins(y=0)

    ax.legend().remove()

    fig.tight_layout()

    plt.savefig(figure_output_dir + 'NMF_' + str(n_components) + '_comp.pdf', dpi=500, bbox_inches='tight')

    # plt.show()
    plt.close()
