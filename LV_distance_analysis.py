import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

from itertools import combinations

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


def k_means(df):
    param_columns = [x for x in df.columns if 'param_' in x]
    print(param_columns)

    X_data = df[param_columns].values

    random_state = 42
    y_pred = KMeans(n_clusters=2, random_state=random_state).fit_predict(X_data)

    print(y_pred)

    sns.scatterplot(x=df.max_LE, y=y_pred, s=3)
    plt.show()
    exit(0)


def CART(output_dir):
    merge_df = pd.read_csv(output_dir + 'net_analysis_df.csv')
    merge_df.fillna(0, inplace=True)

    X_columns = ['total_node_weight', 'transitivity', 'transitivity_4', 'link_density', 
    'diameter', 'edge_betweenness', 'spreading', 'assortativity']


    # X_columns = ['transitivity','edge_betweenness', 'spreading', 'assortativity']
    X_columns = ['transitivity', 'edge_betweenness', 'spreading', 'assortativity', 'diameter']


    Y_column = ['surv_pprob']
    Y_column = ['surv_pprob']

    X_data = merge_df[X_columns].values
    Y_data = merge_df[Y_column].values

    print(X_columns)

    clf = tree.DecisionTreeRegressor()
    clf = clf.fit(X_data, Y_data)
    print(clf.feature_importances_)

    Y_pred = clf.predict(X_data)
    
    # tree.plot_tree(clf)
    # plt.show()

    Y_pred = Y_pred.reshape(-1)
    Y_data = Y_data.reshape(-1)

    x_line = [0, 0.01, 0.0175]
    y_line = x_line

    print(np.shape(Y_pred))
    print(np.shape(Y_data))
    sns.scatterplot(x=Y_data, y=Y_pred)
    sns.lineplot(x=x_line, y=y_line)

    # plt.yscale('symlog')
    plt.show()

    print(diff)

def make_network_analysis_data(counts_df, adj_mat_dir, output_dir):
    adj_mat_path_templ = adj_mat_dir + "model_#IDX#_adj_mat.csv"

    make_csv = True
    if make_csv:
        data = {
        'model_idx': [],
        'total_node_weight': [],
        'transitivity': [],
        'transitivity_4': [],
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
                data['total_node_weight'].append(0)
                data['transitivity'].append(0)
                data['transitivity_4'].append(0)
                data['link_density'].append(0)
                data['diameter'].append(0)
                data['edge_betweenness'].append(0)
                data['spreading'].append(0)
                data['synchronizability'].append(0)
                data['assortativity'].append(0)

            else:
                net = pyunicorn.core.network.Network(adjacency=adj_mat)

                data['model_idx'].append(model_idx)
                data['total_node_weight'].append(net.total_node_weight)
                data['transitivity'].append(net.transitivity())
                data['transitivity_4'].append(net.higher_order_transitivity(4))
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
        net_analysis_df = pd.DataFrame(data)

        merge_df = pd.merge(counts_df, net_analysis_df, on='model_idx')
        merge_df.to_csv(output_dir + 'net_analysis_df.csv')

    merge_df = pd.read_csv(output_dir + 'net_analysis_df.csv')
    merge_df.sort_values(by=['all_counts'], ascending=False, inplace=True)


    merge_df = merge_df.loc[merge_df['model_idx'] != 1739]

    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    sns.scatterplot(data=merge_df, y='surv_pprob', x='transitivity', ax=ax)
    ax.legend().remove()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.set(xticklabels=[])
    # ax.set(xlabel='')

    plt.show()
    plt.close()


    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    sns.scatterplot(data=merge_df, y='surv_pprob', x='transitivity_4', ax=ax)
    ax.legend().remove()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.set(xticklabels=[])
    # ax.set(xlabel='')

    plt.show()
    plt.close()

    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    sns.scatterplot(data=merge_df, y='strange_pprob', x='edge_betweenness', ax=ax)
    ax.legend().remove()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.set(xticklabels=[])
    # ax.set(xlabel='')

    plt.show()
    plt.close()


    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    sns.barplot(data=merge_df, x='model_idx', y='strange_pprob', ax=ax)
    ax.legend().remove()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.set(xticklabels=[])
    # ax.set(xlabel='')

    plt.show()
    plt.close()

    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    sns.scatterplot(data=merge_df, x='transitivity', y='strange_pprob', ax=ax)
    ax.legend().remove()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.set(xticklabels=[])
    # ax.set(xlabel='')

    plt.show()
    plt.close()


    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    sns.scatterplot(data=merge_df, y='strange_pprob', x='total_node_weight', ax=ax)
    ax.legend().remove()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.set(xticklabels=[])
    # ax.set(xlabel='')

    plt.show()
    plt.close()


    exit()



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

    if max_LE >=0.001:
        label= 'STRANGE'

    elif max_LE > -0.001 and max_LE <=0.001:
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


def main():
    data_dir = "./output/gen_LV_four_chem_0_rej/chunk_0/Population_end/"
    data_dir = "./output/gen_LV_four_chem_0_rej/gen_LV_four_chem_0_rej_10/Population_0/"
    data_dir = '~/Documents/AutoCD/output/gen_LV_four_chem_0_rej_2/gen_LV_four_chem_0_rej_1/Population_0/'
    data_dir = '/home/behzad/Documents/AutoCD/output/gen_LV_four_chem_0_rej_4/chunk_0/Population_end/'
    
    priors_dir = './ABC_input_files/input_files_gen_LV_four_chemostat/input_files/'
    adj_mat_dir = './ABC_input_files/input_files_gen_LV_four_chemostat/adj_matricies/'

    distance_path = data_dir + "distances.csv"
    output_dir = data_dir


    model_idxs = list(range(1740))
    # df = combine_distances(data_dir, model_idxs, output_dir)
    df = pd.read_csv(output_dir + 'master_distances.csv')

    # df = df.loc[df['model_ref'] == 307]
    LE_list = ['d1', 'd2', 'd3', 'd4']
    max_LE_func = lambda row: np.max(row[LE_list].values)
   
    print("Making max LE column")
    # df['max_LE'] = df[LE_list].apply(max_LE_func, axis=1)
    df['max_LE'] = make_max_LE(df)

    print("Making attractor labels")
    df['label'] = df[['max_LE']].apply(label_attractor, axis=1)

    # plot_max_LE_distribution(df, output_dir)

    counts_df = make_counts_df(df, model_idxs)
    counts_df.sort_values(by=['strange_counts'], ascending=False, inplace=True)
    
    CART(output_dir)
    exit()

    make_network_analysis_data(counts_df, adj_mat_dir, output_dir)
    exit()

    df = df.loc[df['model_ref'] == 307]
    strange_sub_df = strange_sub_df.loc[strange_sub_df['model_ref'] == 307]

    plot_2d_param_dens(df, 307, data_dir, priors_dir, output_dir)

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