import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
from sklearn.tree import export_graphviz
from itertools import combinations
from sklearn.ensemble import RandomForestRegressor
# # Importing the Keras libraries and packages
# from keras.models import Sequential
# from keras.layers import Convolution2D
# from keras.layers import MaxPooling2D
# from keras.layers import Flatten
# from keras.layers import Dense
import os
import subprocess
from collections import Counter

import pyunicorn
from sklearn import tree

import graphviz
import pickle

def process_data(df):
    df = df[df['integ_error'] != 'species_decayed']
    df.reset_index(inplace=True)

    return df

def make_interaction_strings():
    n_species = 4
    interaction_template = 'm_TO_FROM'
    interaction_strings = []

    for strain_TO_idx in range(1, n_species + 1):
        for strain_FROM_idx in range(1, n_species + 1):
            interaction = interaction_template.replace('TO', str(strain_TO_idx))
            interaction = interaction.replace('FROM', str(strain_FROM_idx))
            interaction_strings.append(interaction)

    return interaction_strings

def permute_params(df):
    param_columns = [x for x in df.columns if 'param_' in x]
    species_labels = ['1', '2', '3', '4']

    interaction_columns = [x for x in df.columns if 'param_m_' in x]
    growth_param_cols = [x for x in df.columns if 'param_mu' in x]
    init_N_cols = [x for x in df.columns if 'N_' in x]

    interaction_matrix = np.reshape(interaction_columns, (4, 4))
    # print(interaction_columns[:, 0])
    interaction_columns = np.reshape(interaction_matrix, (-1))

    permuted_param_columns = []

    for i, _ in enumerate(species_labels):
        for j, _ in enumerate(species_labels):
            if i == j :
                continue
            permuted_matrix = np.copy(interaction_matrix)

            # Swap rows
            permuted_matrix[i] = np.copy(interaction_matrix[j])
            permuted_matrix[j] = np.copy(interaction_matrix[i])

            # Swap columns
            col_i = np.copy(permuted_matrix[:, i])
            col_j = np.copy(permuted_matrix[:, j])
            permuted_matrix[:, i] = col_j
            permuted_matrix[:, j] = col_i
            permuted_matrix = np.reshape(permuted_matrix, (-1))


            # Swap growth params
            premuted_growth = np.copy(growth_param_cols)
            mu_i = growth_param_cols[i]
            mu_j = growth_param_cols[j]

            premuted_growth[i] = mu_i
            premuted_growth[j] = mu_j

            # Swap init N
            permuted_N = np.copy(init_N_cols)
            N_i = init_N_cols[i]
            N_j = init_N_cols[j]

            permuted_N[i] = N_i
            permuted_N[j] = N_j

            concat_permuted_params = np.concatenate([permuted_matrix, premuted_growth, permuted_N], axis=0) 
            permuted_param_columns.append(list(concat_permuted_params))

    param_idx_cols = ['param_sim_idx', 'param_batch_idx']

    non_param_columns = [x for x in df.columns if not 'param_' in x]
    new_rows = []
    for idx, row in df.iterrows():
        for p_cols in permuted_param_columns:
            row_n = row[non_param_columns + param_idx_cols + p_cols]
            new_rows.append(row_n.values)

    permuted_df = pd.DataFrame(data=new_rows, columns=df.columns)
    return permuted_df


def plot_params(df, model_idx, data_dir, output_dir):
    output_path = output_dir + "params_model_" + str(model_idx) + ".pdf" 
    model_params_dir = data_dir + "model_sim_params/"
    model_params_path = model_params_dir + "model_#IDX#_all_params.csv".replace('#IDX#', str(model_idx))
    params_df = pd.read_csv(model_params_path)
    # species_labels = ['1', '2', '3', '4']
    # interaction_strings = make_interaction_strings()

    params_list = []
    for idx, row in df.iterrows():
        model_ref = int(row['model_ref'])
        sim_idx = row['sim_idx']
        batch_idx = row['batch_idx']

        sub_df = params_df.loc[params_df['sim_idx'] == sim_idx]
        sub_df = sub_df.loc[sub_df['batch_idx'] == batch_idx]

        interaction_columns = [x for x in sub_df.columns if 'm_' in x]
        param_columns = [x for x in sub_df.columns]

        params = sub_df[param_columns].values[0]
        params_list.append(params)

    params_df = pd.DataFrame(params_list, columns=param_columns)
    params_df.drop('sim_idx', axis=1, inplace=True)
    params_df.drop('batch_idx', axis=1, inplace=True)
       
    interaction_columns = [x for x in params_df.columns if 'm_' in x]
    plot_params = interaction_columns + ['D']

    master_df = pd.concat([df, params_df], axis=1)
    subplot_cols = 3
    subplot_rows = int(np.ceil(len(plot_params) / subplot_cols))

    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches), nrows=subplot_rows, ncols=subplot_cols, sharex=True, sharey=True)
    ax = ax.reshape(-1)
    
    for idx, interaction in enumerate(plot_params):
        # sns.scatterplot(x=interaction, y='max_LE', ax=ax[idx], data=master_df, hue='label')
        sns.distplot(master_df[interaction], ax=ax[idx], bins=50)

        ax[idx].set_xlim(0.1, 2.5)
        ax[idx].legend().remove()
        ax[idx].spines["right"].set_visible(False)
        ax[idx].spines["top"].set_visible(False)
        ax[idx].spines["left"].set_alpha(0.5)
        ax[idx].spines["bottom"].set_alpha(0.5)

    plt.tight_layout()
    plt.savefig(output_path)

def plot_2d_param_dens(model_distance_df, model_idx, data_dir, priors_dir, output_dir):
    output_path = output_dir + "params_model_" + str(model_idx) + ".pdf" 
    model_params_dir = data_dir + "model_sim_params/"
    model_params_path = model_params_dir + "model_#IDX#_all_params.csv".replace('#IDX#', str(model_idx))    
    params_df = pd.read_csv(model_params_path)

    temp_model_params_path = model_params_dir + "temp_model_#IDX#_all_params.csv".replace('#IDX#', str(model_idx))    

    model_input_parameter_path = priors_dir + '/params_#IDX#.csv'.replace('#IDX#', str(model_idx))
    model_input_species_path = priors_dir + '/species_#IDX#.csv'.replace('#IDX#', str(model_idx))
    
    params_list = []
    for idx, row in model_distance_df.iterrows():
        model_ref = int(row['model_ref'])
        sim_idx = row['sim_idx']
        batch_idx = row['batch_idx']

        sub_df = params_df.loc[params_df['sim_idx'] == sim_idx]
        sub_df = sub_df.loc[sub_df['batch_idx'] == batch_idx]

        interaction_columns = [x for x in sub_df.columns if 'm_' in x]
        param_columns = [x for x in sub_df.columns]

        params = sub_df[param_columns].values[0]
        params_list.append(params)

    params_df = pd.DataFrame(params_list, columns=param_columns)
    # params_df.drop('sim_idx', axis=1, inplace=True)
    # params_df.drop('batch_idx', axis=1, inplace=True)

    params_df.to_csv(temp_model_params_path)

    R_script = os.path.dirname(os.path.realpath(__file__)) + "/data_analysis/dens_plot_2D_LE.R"
    subprocess.call(['Rscript', R_script, temp_model_params_path, model_input_parameter_path, model_input_species_path, str(model_idx), output_dir, 'TRUE', 'TRUE'])


def plot_max_LE_distribution(df, output_dir):
    bins = 100
    plt.hist(df.max_LE, bins=bins)
    plt.savefig(output_dir + "max_LE_dist.pdf", dpi=500)
    # plt.show()

def PCA_analysis(df):
    param_columns = [x for x in df.columns if 'param_' in x]
    print(param_columns)

    n_components = 6
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(df[param_columns])

    component_columns = ['pca' + str(i) for i in range(n_components)]

    pca_data = pd.DataFrame(X_pca, columns=component_columns)
    pca_data['label'] = df['label'].values

    col_combinations = list(combinations(component_columns, r=2))


    n_cols = 2
    nrow = int(np.ceil(len(col_combinations) / n_cols))

    fig, axes = plt.subplots(nrows=nrow, ncols=n_cols)
    axes = np.array(axes)

    for idx, ax in enumerate(axes.reshape(-1)):
        if idx >= len(col_combinations):
            break

        x = col_combinations[idx][0]
        y = col_combinations[idx][1]
        sns.scatterplot(data=pca_data, x=x, y=y, hue='label', ax=ax, s =4)
    
    plt.show()


def multi_linear_regression(output_dir):
    merge_df = pd.read_csv(output_dir + 'net_analysis_df.csv')
    merge_df.fillna(0, inplace=True)

    X_columns = ['total_mat_weight', 'transitivity', 'transitivity_4', 'link_density', 
    'diameter', 'edge_betweenness', 'spreading', 'assortativity']


    # X_columns = ['transitivity','edge_betweenness', 'spreading', 'assortativity']
    X_columns = ['transitivity', 'edge_betweenness', 'spreading', 'assortativity', 'diameter']


    Y_column = ['surv_pprob']
    Y_column = ['strange_pprob']

    X_data = merge_df[X_columns].values
    Y_data = merge_df[Y_column].values
    print(np.shape(X_data))
    print(np.shape(Y_data))

    reg = LinearRegression().fit(X_data, Y_data)
    print(reg.score(X_data, Y_data))
    Y_pred = reg.predict(X_data)

    Y_data = Y_data.reshape(-1)
    Y_pred = Y_pred.reshape(-1)

    sns.scatterplot(x=Y_data, y=Y_pred)

    plt.savefig(output_dir + "mlr.pdf")
    plt.show()



def k_means(output_dir):
    merge_df = pd.read_csv(output_dir + 'net_analysis_df.csv')
    merge_df.fillna(0, inplace=True)

    X_columns = ['total_mat_weight', 'transitivity', 'transitivity_4', 'link_density', 
    'diameter', 'edge_betweenness', 'spreading', 'assortativity']


    # X_columns = ['transitivity','edge_betweenness', 'spreading', 'assortativity']
    X_columns = ['transitivity', 'edge_betweenness', 'spreading', 'assortativity', 'diameter']


    Y_column = ['surv_pprob']
    Y_column = ['strange_pprob']

    X_data = merge_df[X_columns].values
    Y_data = merge_df[Y_column].values

    random_state = 42
    y_pred = KMeans(n_clusters=2, random_state=random_state).fit_predict(X_data)
    print(y_pred[0])
    Y_data = Y_data.reshape(-1)
    # y_pred = y_pred.reshape(-1)

    print(np.shape(y_pred))

    sns.scatterplot(x=y_pred[:,0], y=y_pred[:,1], s=4, hue=Y_data)
    plt.show()

def random_forest(output_dir):
    merge_df = pd.read_csv(output_dir + 'net_analysis_df.csv')
    merge_df.fillna(-1, inplace=True)

    merge_df.drop(merge_df.loc[merge_df['strange_pprob']==0.0].index, inplace=True)

    X_columns = ['total_mat_weight', 'transitivity', 'transitivity_4', 'link_density', 
    'diameter', 'edge_betweenness', 'spreading', 'assortativity']


    # X_columns = ['transitivity','edge_betweenness', 'spreading', 'assortativity']
    X_columns = ['total_mat_weight', 'transitivity', 'edge_betweenness', 'spreading', 'assortativity', 'diameter']

    X_columns = ['total_mat_weight', 'transitivity', 'spreading', 'synchronizability']#, 'edge_betweenness', 'spreading', 'assortativity', 'diameter']
    # X_columns = ['transitivity', 'spreading'] #, 'transitivity', 'spreading']#, 'edge_betweenness', 'spreading', 'assortativity', 'diameter']
    X_columns = ['norm_total_mat_weight', 'transitivity', 'norm_spreading', 'norm_synchronizability']#, 'edge_betweenness', 'spreading', 'assortativity', 'diameter']

    Y_column = ['surv_pprob']
    Y_column = ['strange_pprob']

    Y_column = ['norm_strange_pprob']


    mask = np.random.rand(len(merge_df)) < 1.0
    train_df = merge_df[mask]
    test_df = merge_df[~mask]
    test_df = merge_df[mask]


    X_train = train_df[X_columns].values
    Y_train = train_df[Y_column].values.reshape(-1)
    
    X_test = test_df[X_columns].values
    Y_test = test_df[Y_column].values

    regr = RandomForestRegressor(max_depth=25, random_state=0)
    regr.fit(X_train, Y_train)

    print(regr.feature_importances_)

    Y_train_pred = regr.predict(X_train)
    Y_test_pred = regr.predict(X_test)


    Y_test = Y_test.reshape(-1)
    Y_train = Y_train.reshape(-1)
    Y_test_pred = Y_test_pred.reshape(-1)
    Y_train_pred = Y_train_pred.reshape(-1)


    x_line = [0, 0.1, 1.0]
    y_line = x_line

    # sns.scatterplot(x=Y_train, y=Y_train_pred, s=6, label='train_set')
    sns.scatterplot(x=Y_test, y=Y_test_pred, s=2, label='test_set', linewidth=0)

    sns.lineplot(x=x_line, y=y_line)

    plt.savefig(output_dir + "rdn_frst.pdf")
    plt.close()

    with open(output_dir + 'strange_rdn_frst.pkl', 'wb') as fid:
        pickle.dump(regr, fid)    




def CART(output_dir):
    merge_df = pd.read_csv(output_dir + 'net_analysis_df.csv')
    merge_df.fillna(-1, inplace=True)

    # merge_df.drop(merge_df.loc[merge_df['strange_pprob']==0.0].index, inplace=True)
    X_columns = ['total_mat_weight', 'transitivity', 'transitivity_4', 'link_density', 
    'diameter', 'edge_betweenness', 'spreading', 'assortativity']


    # X_columns = ['transitivity','edge_betweenness', 'spreading', 'assortativity']
    X_columns = ['transitivity', 'edge_betweenness', 'spreading', 'assortativity', 'diameter']
    X_columns = ['norm_total_mat_weight', 'transitivity', 'norm_spreading', 'norm_synchronizability']#, 'edge_betweenness', 'spreading', 'assortativity', 'diameter']


    Y_column = ['surv_pprob']
    Y_column = ['strange_pprob']

    Y_column = ['norm_surv_pprob']
    Y_column = ['norm_strange_pprob']

    mask = np.random.rand(len(merge_df)) < 1.0
    train_df = merge_df[mask]
    test_df = merge_df[~mask]
    test_df = merge_df[mask]


    X_train = train_df[X_columns].values
    Y_train = train_df[Y_column].values


    clf = tree.DecisionTreeRegressor()
    clf = clf.fit(X_train, Y_train)
    print(clf.feature_importances_)

    X_test = test_df[X_columns].values
    Y_test = test_df[Y_column].values


    Y_train_pred = clf.predict(X_train)
    Y_test_pred = clf.predict(X_test)


    dot_data = export_graphviz(clf, out_file =None, feature_names = X_columns)
    graph = graphviz.Source(dot_data)
    # graph.view()


    tree.plot_tree(clf)
    plt.savefig(output_dir + "tree.pdf")
    plt.close()

    Y_test = Y_test.reshape(-1)
    Y_train = Y_train.reshape(-1)
    Y_test_pred = Y_test_pred.reshape(-1)
    Y_train_pred = Y_train_pred.reshape(-1)

    x_line = [0, 0.01,1.0]
    y_line = x_line

    # sns.scatterplot(x=Y_train, y=Y_train_pred, s=6, label='train_set')
    sns.scatterplot(x=Y_test, y=Y_test_pred, s=3, label='test_set', linewidth=0)

    sns.lineplot(x=x_line, y=y_line, alpha=0.5, color='black', linewidth=1)

    plt.savefig(output_dir + "tree_fit.pdf")
    plt.close()
    
    with open(output_dir + 'strange_dectree.pkl', 'wb') as fid:
        pickle.dump(clf, fid)    

    # cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
    # sns.scatterplot(x='synchronizability', y='strange_pprob', data=merge_df, hue='strange_pprob', s=5, linewidth=0, palette=cmap)
    # plt.show()

    cmap = sns.cubehelix_palette(rot=-.2, as_cmap=True)
    sns.scatterplot(x='synchronizability', y='total_mat_weight', data=merge_df, hue='strange_pprob', palette=cmap)
    # plt.show()

def min_max_normalise_list(data_list):
    max_data = np.max(data_list)
    min_data = np.min(data_list)
    norm_list = [(x - min_data) / (max_data - min_data) for x in data_list]

    return norm_list


def make_network_analysis_data(counts_df, adj_mat_dir, output_dir):
    adj_mat_path_templ = adj_mat_dir + "model_#IDX#_adj_mat.csv"
    
    make_csv = True
    if make_csv:
        data = {
        'model_idx': [],
        'total_mat_weight': [],
        'transitivity': [],
        'link_density': [], 
        'diameter': [], 
        'edge_betweenness': [],
        'spreading': [],
        'synchronizability': [],
        'assortativity': []
        }

        for model_idx in counts_df.model_idx.values:
            print(model_idx)

            adj_mat_path = adj_mat_path_templ.replace('#IDX#', str(model_idx))
            adj_mat = pd.read_csv(adj_mat_path)
            adj_mat = adj_mat.loc[:, ~adj_mat.columns.str.contains('^Unnamed')]
            adj_mat = adj_mat.values

            if np.sum(adj_mat) == 0:
                data['model_idx'].append(model_idx)
                data['total_mat_weight'].append(0)
                data['transitivity'].append(0)
                data['link_density'].append(0)
                data['diameter'].append(0)
                data['edge_betweenness'].append(0)
                data['spreading'].append(0)
                data['synchronizability'].append(0)
                data['assortativity'].append(0)

            else:
                net = pyunicorn.core.network.Network(adjacency=adj_mat)

                data['model_idx'].append(model_idx)
                data['total_mat_weight'].append(np.sum(adj_mat))
                data['transitivity'].append(net.transitivity())
                data['link_density'].append(net.link_density)
                data['diameter'].append(net.diameter())
                data['edge_betweenness'].append(np.sum(net.edge_betweenness()))
                
                data['spreading'].append(np.sum(net.spreading()))

                # print(net.local_cliquishness(3))
                data['synchronizability'].append(net.msf_synchronizability())
                # print(net.msf_synchronizability())

                try:
                    data['assortativity'].append(net.assortativity())

                except ZeroDivisionError as e:
                    data['assortativity'].append(0)

        sum_strange = np.sum(counts_df['strange_counts'].values)
        strange_pprobs = [x /sum_strange for x in counts_df['strange_counts'].values]

        sum_surv = np.sum(counts_df['all_counts'].values)
        surv_pprobs = [x / sum_surv for x in counts_df['all_counts'].values]

        data['strange_pprob'] = strange_pprobs
        data['surv_pprob'] = surv_pprobs

        data['norm_surv_pprob'] = min_max_normalise_list(surv_pprobs)
        data['norm_strange_pprob'] = min_max_normalise_list(strange_pprobs)
        data['norm_total_mat_weight'] = min_max_normalise_list(data['total_mat_weight'])
        data['norm_spreading'] = min_max_normalise_list(data['spreading'])
        data['norm_synchronizability'] = min_max_normalise_list(data['synchronizability'])

        net_analysis_df = pd.DataFrame(data)

        merge_df = pd.merge(counts_df, net_analysis_df, on='model_idx')
        merge_df.sort_values(by=['model_idx'], ascending=True, inplace=True)
        merge_df.to_csv(output_dir + 'net_analysis_df.csv')

    merge_df = pd.read_csv(output_dir + 'net_analysis_df.csv')
    merge_df.sort_values(by=['all_counts'], ascending=False, inplace=True)

    merge_df = merge_df.loc[merge_df['total_mat_weight'] != 0]

    # width_inches = 95*4 / 25.4
    # height_inches = 51*4 / 25.4
    # fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    # sns.scatterplot(data=merge_df, y='surv_pprob', x='transitivity', ax=ax)
    # ax.legend().remove()
    # ax.spines["right"].set_visible(False)
    # ax.spines["top"].set_visible(False)
    # ax.spines["left"].set_alpha(0.5)
    # ax.spines["bottom"].set_alpha(0.5)
    # ax.set(xticklabels=[])
    # # ax.set(xlabel='')

    # plt.show()
    # plt.close()


    # width_inches = 95*4 / 25.4
    # height_inches = 51*4 / 25.4
    # fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    # sns.scatterplot(data=merge_df, y='surv_pprob', x='transitivity_4', ax=ax)
    # ax.legend().remove()
    # ax.spines["right"].set_visible(False)
    # ax.spines["top"].set_visible(False)
    # ax.spines["left"].set_alpha(0.5)
    # ax.spines["bottom"].set_alpha(0.5)
    # ax.set(xticklabels=[])
    # # ax.set(xlabel='')

    # plt.show()
    # plt.close()

    # width_inches = 95*4 / 25.4
    # height_inches = 51*4 / 25.4
    # fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    # sns.scatterplot(data=merge_df, y='strange_pprob', x='edge_betweenness', ax=ax)
    # ax.legend().remove()
    # ax.spines["right"].set_visible(False)
    # ax.spines["top"].set_visible(False)
    # ax.spines["left"].set_alpha(0.5)
    # ax.spines["bottom"].set_alpha(0.5)
    # ax.set(xticklabels=[])
    # # ax.set(xlabel='')

    # plt.show()
    # plt.close()


    # width_inches = 95*4 / 25.4
    # height_inches = 51*4 / 25.4
    # fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    # sns.barplot(data=merge_df, x='model_idx', y='strange_pprob', ax=ax)
    # ax.legend().remove()
    # ax.spines["right"].set_visible(False)
    # ax.spines["top"].set_visible(False)
    # ax.spines["left"].set_alpha(0.5)
    # ax.spines["bottom"].set_alpha(0.5)
    # ax.set(xticklabels=[])
    # # ax.set(xlabel='')

    # plt.show()
    # plt.close()

    # width_inches = 95*4 / 25.4
    # height_inches = 51*4 / 25.4
    # fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    # sns.scatterplot(data=merge_df, x='transitivity', y='strange_pprob', ax=ax)
    # ax.legend().remove()
    # ax.spines["right"].set_visible(False)
    # ax.spines["top"].set_visible(False)
    # ax.spines["left"].set_alpha(0.5)
    # ax.spines["bottom"].set_alpha(0.5)
    # ax.set(xticklabels=[])
    # # ax.set(xlabel='')

    # plt.show()
    # plt.close()


    # width_inches = 95*4 / 25.4
    # height_inches = 51*4 / 25.4
    # fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    # sns.scatterplot(data=merge_df, y='strange_pprob', x='total_mat_weight', ax=ax)
    # ax.legend().remove()
    # ax.spines["right"].set_visible(False)
    # ax.spines["top"].set_visible(False)
    # ax.spines["left"].set_alpha(0.5)
    # ax.spines["bottom"].set_alpha(0.5)
    # ax.set(xticklabels=[])
    # # ax.set(xlabel='')

    # plt.show()
    # plt.close()


def combine_distances(data_dir, model_idxs, output_dir):
    distances_dir = data_dir + "model_sim_distances/"
    model_distance_template = distances_dir + "model_#IDX#_distances.csv"

    df_list = []
    for m_idx in model_idxs:
        m_distance_path = model_distance_template.replace('#IDX#', str(m_idx))
        try:
            df = pd.read_csv(m_distance_path)

            df_list.append(df)

        except FileNotFoundError:
            continue

    master_df = pd.concat(df_list)

    master_df.to_csv(output_dir + 'master_distances.csv')

    return master_df

def label_attractor(row):
    max_LE = row['max_LE']
    # max_LE = LE_list[0]

    if max_LE > 0.05:
        label= 'STRANGE'

    elif max_LE > -0.05 and max_LE <=0.05:
        label= 'LIMIT'

    else:
        label= 'STABLE'

    return label

def make_max_LE(df):
    LE_list = ['d1', 'd2', 'd3', 'd4']
    LE_arr = df[LE_list].values
    max_LE = [np.max(x) for x in LE_arr]

    return max_LE

def make_counts_df(df, model_idxs): 
    strange_sub_df = df.loc[df['label'] == 'STRANGE']
    strange_counts = Counter(strange_sub_df['model_ref'].values)
    all_counts = Counter(df['model_ref'].values)

    all_model_idxs = []
    strange_counts_list = []
    all_counts_list = []
    ratio = []

    data = {
    'model_idx': [],
    'strange_counts': [],
    'all_counts': [],
    'ratio': []
    }

    threshold = 0.05
    for x in model_idxs:
        if all_counts[x] > 0:
            data['model_idx'].append(x)
            data['strange_counts'].append(strange_counts[x])
            data['all_counts'].append(all_counts[x])
            data['ratio'].append(strange_counts[x] / all_counts[x])

        else:
            data['model_idx'].append(x)
            data['strange_counts'].append(0)
            data['all_counts'].append(0)
            data['ratio'].append(0)

    counts_df = pd.DataFrame(data)

    return counts_df

def compare_prediction_and_real(data_dir):
    real_df = pd.read_csv(data_dir + "net_analysis_df.csv")
    pred_df = pd.read_csv(data_dir + "sampled_models.csv")

    pred_df = pred_df[['model_idx', 'rdn_frst_predictions']]


    sum_pred = np.sum(pred_df['rdn_frst_predictions'].values)

    rescale_func = lambda x: x / sum_pred

    pred_df['norm_pred'] = np.vectorize(rescale_func)(pred_df['rdn_frst_predictions'])

    print(pred_df.columns)
    print(real_df.columns)

    merge_df = real_df.merge(pred_df, how='outer')

    sns.scatterplot(y='rdn_frst_predictions', x='strange_pprob', data=merge_df)
    plt.show()


def main():
    data_dir = "./output/gen_LV_four_chem_0_rej/chunk_0/Population_end/"
    data_dir = "./output/gen_LV_four_chem_0_rej/gen_LV_four_chem_0_rej_10/Population_0/"
    data_dir = '~/Documents/AutoCD/output/gen_LV_four_chem_0_rej_2/gen_LV_four_chem_0_rej_1/Population_0/'
    data_dir = '/home/behzad/Documents/AutoCD/output/gen_LV_four_chem_0_rej_4/chunk_0/Population_end/'

    data_dir = '/home/behzad/Documents/AutoCD/output/gen_LV_four_chem_3_rej_0/chunk_0/Population_end/'

    priors_dir = './ABC_input_files/input_files_gen_LV_four_chemostat_3/input_files/'
    adj_mat_dir = './ABC_input_files/input_files_gen_LV_four_chemostat_3/adj_matricies/'



    # data_dir = '/home/behzad/Documents/AutoCD/output/gen_LV_five_chem_0_rej_0/chunk_0/Population_end/'

    # priors_dir = './ABC_input_files/input_files_gen_LV_five_chemostat_0/input_files/'
    # adj_mat_dir = './ABC_input_files/input_files_gen_LV_five_chemostat_0/adj_matricies/'


    output_dir = data_dir

    # compare_prediction_and_real(data_dir)
    # exit()

    model_idxs = list(range(1017))
    model_idxs = list(range(1))
    df = combine_distances(data_dir, model_idxs, output_dir)
    df = pd.read_csv(output_dir + 'master_distances.csv')
    # exit()
    # # k_means(output_dir)
    # # multi_linear_regression(output_dir)
    # random_forest(output_dir)
    # CART(output_dir)
    # exit()
    # CART(output_dir)
    # exit()


    LE_list = ['d1', 'd2', 'd3', 'd4', 'd5']
    max_LE_func = lambda row: np.max(row[LE_list].values)
   
    print("Making max LE column")
    df['max_LE'] = make_max_LE(df)

    # print(df.loc[df['max_LE'].idxmax()])
    # exit()
    print("Making attractor labels")
    df['label'] = df[['max_LE']].apply(label_attractor, axis=1)

    # plot_max_LE_distribution(df, output_dir)

    counts_df = make_counts_df(df, model_idxs)
    counts_df.sort_values(by=['strange_counts'], ascending=False, inplace=True)
    make_network_analysis_data(counts_df, adj_mat_dir, output_dir)

    random_forest(output_dir)

    # multi_linear_regression(output_dir)
    # k_means(output_dir)
    df = df.loc[df['model_ref'] == 0]
    strange_sub_df = strange_sub_df.loc[strange_sub_df['model_ref'] == 0]

    plot_2d_param_dens(df, 0, data_dir, priors_dir, output_dir)

    exit()
    # exit()
    # print(counts_df)
    # strange_sub_df = strange_sub_df.loc[strange_sub_df['model_ref'] == 307]
    # plot_params(strange_sub_df, 243, data_dir, output_dir)

    exit()

    print(x)
    exit()
    for idx, row in strange_sub_df.iterrows():
        print(row['model_ref'])
        # print(row[param_columns])

        # print(row[param_columns].values)
        print("")
    # exit()
    plot_params(df, data_dir, output_dir)
    plot_max_LE_distribution(df, output_dir)
    # CNN(df)
    # exit()
    # k_means(df)
    # exit()
    # PCA_analysis(df)
    exit()

if __name__ == "__main__":
    main()