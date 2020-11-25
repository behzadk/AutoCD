import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import LogNorm


def normalise_list(x):
    x_min = np.min(x)
    x_max = np.max(x)

    norm_lambd = lambda v: (v - x_min) / (x_max - x_min)
    norm_x = [norm_lambd(v) for v in x ]
    return norm_x

def set_state(sep_coeff, threshold):
    states = []
    for x in sep_coeff:
        if x < 0:
            states.append('STABLE')

        elif x > 0 and x < threshold:
            states.append('OSC')

        else:
            states.append('CHAOS')

    return states

def label_segment(m_values):
    states = []
    for m in m_values:
        if m > 0.48 and m < 0.55:
            states.append('B_d')

        elif m < 0.48 and m > 0.2:
            states.append('B')

        else:
            states.append('A')

    return states



def plot_figure_2():
    # https://www.pnas.org/content/pnas/86/1/142.full.pdf
    data_dir = "./output/autocat_0_rej/autocat_0_rej_2/Population_0/"
    distances_df = pd.read_csv(data_dir + "distances.csv")
    params_df = pd.read_csv(data_dir + "model_sim_params/model_0_population_all_params")
    print(params_df)
    print(distances_df)

    plot_data = {'m': params_df['m'].values, 'K_1': params_df['K_1'].values, 'd1': distances_df['d1'].values}
    plot_df = pd.DataFrame(plot_data)
    # plot_df['d1'] = normalise_list(plot_df['d1'].values)
    plot_df['state'] = set_state(plot_df['d1'].values, threshold=0.0004)
    print(plot_df['state'])

    width_inches = 95*4 / 25.4
    height_inches = 95*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches,height_inches))
    sns.scatterplot(data=plot_df, x='K_1', y='m', hue='state', ax=ax)
    ax.set_xlim(1, 20)
    ax.set_ylim(0.1, 1)
    plt.show()
    exit()

def m_change():
    # https://www.pnas.org/content/pnas/86/1/142.full.pdf
    data_dir = "./output/autocat_2_rej/autocat_2_rej_10/Population_0/"
    distances_df = pd.read_csv(data_dir + "distances.csv")
    params_df = pd.read_csv(data_dir + "model_sim_params/model_0_population_all_params")

    plot_data = {'m': params_df['m'].values, 'K_1': params_df['K_1'].values, 'd1': distances_df['d1'].values}
    plot_df = pd.DataFrame(plot_data)
    # plot_df['d1'] = normalise_list(plot_df['d1'].values)
    plot_df['state'] = label_segment(plot_df['m'].values)
    print(plot_df['state'])

    width_inches = 95*4 / 25.4
    height_inches = 95*4 / 25.4
    fig, ax = plt.subplots(figsize=(width_inches,height_inches))
    sns.scatterplot(data=plot_df, x='m', y='d1', ax=ax)
    ax.set_xlim(0.5, 0.54)
    ax.set_xlabel('m')
    ax.set_ylabel('separation coefficient')

    # ax.set_ylim(0.1, 1)
    # plt.show()
    plt.savefig(data_dir + 'autocatalytic_test.pdf', dpi=500)
    exit()

def entropy_heatmap():
    data_dir = "./output/autocat_3_rej/autocat_3_rej_3/Population_0/"
    distances_df = pd.read_csv(data_dir + "distances.csv")
    params_df = pd.read_csv(data_dir + "model_sim_params/model_0_population_all_params")

    inf_idxs = distances_df.index[np.isinf(distances_df[['d1', 'd2', 'd3']]).any(1)]
    # distances_df.drop(index=inf_idxs, inplace=True)
    # params_df.drop(index=inf_idxs, inplace=True)

    distances_df.reset_index()
    params_df.reset_index()

    plot_data = {
    'm': params_df['m'].values, 'K_1': params_df['K_1'].values, 'd1': distances_df['d1'].values, 
    'd2': distances_df['d2'].values, 'd3': distances_df['d3'].values, 'd4': distances_df['d4'].values,
    'sep_coeff': distances_df['d5'].values, 'DET': distances_df['d6'].values, 'ENT': distances_df['d7'].values, 
    'LAM': distances_df['d8'].values, 'max_LE': np.max(distances_df[['d9', 'd10', 'd11', 'd12']].values)
    }
    plot_df = pd.DataFrame(plot_data)

    n_bins = 30

    m_bins = np.linspace(0.1, 1.0, n_bins)
    K_1_bins = np.linspace(1.0, 20.0, n_bins)

    m_bins = list(np.around(m_bins,2))
    K_1_bins = list(np.around(K_1_bins,2))

    plot_metrics = ['sep_coeff', 'DET', 'ENT', 'LAM', 'max_LE']

    for metric in plot_metrics:
        fig, ax = plt.subplots(figsize=(13,10)) 

        binned_grid = stats.binned_statistic_2d(x=plot_df['K_1'], y=plot_df['m'], values=plot_df[metric], statistic='mean', bins=[K_1_bins, m_bins], range=None, expand_binnumbers=False)
        sns.heatmap(binned_grid.statistic.T, ax=ax)
        # ax.set_xticklabels(K_1_bins, rotation=90)
        ax.set_facecolor("yellow")
        ax.set_title(metric)

        ax.invert_yaxis()
        # ax.set_yticklabels(m_bins, rotation=0)

        output_path = data_dir + metric + '_autocat.pdf'

        plt.savefig(output_path, dpi=500)

def make_LE_data():
    data_dir = "./output/autocat_7_sprott_bene/chunk_0/Population_end/"
    distances_df = pd.read_csv(data_dir + "model_sim_distances/model_0_distances.csv")
    params_df = pd.read_csv(data_dir + "model_sim_params/model_0_all_params.csv")

    inf_idxs = distances_df.index[np.isinf(distances_df[['d1', 'd2', 'd3']]).any(1)]
    # distances_df.drop(index=inf_idxs, inplace=True)
    # params_df.drop(index=inf_idxs, inplace=True)
    # distances_df = distances_df.round({'d1': 3, 'd2': 3, 'd3': 3, 'd4': 3})

    distances_df.reset_index()
    params_df.reset_index()

    LE_classification = []
    # plt.hist(distances_df['d1'].values, 100)
    # plt.show()
    # exit(0)

    max_LE = []

    min_species = []
    m_values = []

    for idx, row in distances_df.iterrows():
        X1_min = 1/row['d6']
        X2_min = 1/row['d7']
        X3_min = 1/row['d8']
        X4_min = 1/row['d9']

        print(np.min([X1_min, X2_min, X3_min, X4_min]))
        min_species.append(np.min([X1_min, X2_min, X3_min, X4_min]))

    plot_data = {
    'm': params_df['m'].values, 'K_1': params_df['K_1'].values, 'sprott_LE': distances_df['d1'].values, 
    'bene_LE': distances_df['d2'].values, 'min_species': min_species
    }

    plot_df = pd.DataFrame(plot_data)

    plot_df.to_csv(data_dir + 'LE_data.csv')


def RQA_plots():
    # https://www.pnas.org/content/pnas/86/1/142.full.pdf
    data_dir = "./output/autocat_3_rej/chunk_0/Population_end/"
    distances_df = pd.read_csv(data_dir + "distances.csv")
    params_df = pd.read_csv(data_dir + "model_sim_params/model_0_population_all_params")

    plot_data = {
    'm': params_df['m'].values, 'K_1': params_df['K_1'].values, 'd1': distances_df['d1'].values, 
    'd2': distances_df['d2'].values, 'd3': distances_df['d3'].values, 'd4': distances_df['d4'].values
    }

    plot_df = pd.DataFrame(plot_data)
    # plot_df['d1'] = normalise_list(plot_df['d1'].values)
    plot_df['state'] = label_segment(plot_df['m'].values)
    # print(plot_df['state'])

    width_inches = 95*4 / 25.4
    height_inches = 95*4 / 25.4
    fig, axes = plt.subplots(figsize=(width_inches,height_inches), nrows=4, ncols=1)
    sns.scatterplot(data=plot_df, x='m', y='d1', ax=axes[0])
    axes[0].set_xlim(0.5, 0.54)
    axes[0].set_xlabel('m')
    axes[0].set_ylabel('separation coefficient')
    # axes[0].set_yscale('symlog')

    sns.scatterplot(data=plot_df, x='m', y='d2', ax=axes[1])
    axes[1].set_xlim(0.5, 0.54)
    axes[1].set_xlabel('m')
    axes[1].set_ylabel('DET')
    axes[1].set_ylim(0.99, 1.01)

    # axes[1].set_yscale('log')

    sns.scatterplot(data=plot_df, x='m', y='d3', ax=axes[2])
    axes[2].set_xlim(0.5, 0.54)
    axes[2].set_xlabel('m')
    axes[2].set_ylabel('ENT')
    # axes[2].set_ylim(4, 5)

    sns.scatterplot(data=plot_df, x='m', y='d4', ax=axes[3])
    axes[3].set_xlim(0.5, 0.54)
    axes[3].set_xlabel('m')
    axes[3].set_ylabel('LAM')
    # axes[3].set_yscale('log')



    axes[3].set_ylim(0.97, 1.056)
    # plt.show()
    plt.savefig(data_dir + 'autocatalytic_test.pdf', dpi=500)
    exit()


def LE_heatmaps():
    data_dir = "./output/autocat_7_sprott_bene/chunk_0/Population_end/"
    distances_df = pd.read_csv(data_dir + "model_sim_distances/model_0_distances.csv")
    params_df = pd.read_csv(data_dir + "model_sim_params/model_0_all_params.csv")

    plot_df = pd.read_csv(data_dir + "LE_data.csv")


    extinct_func = lambda x: True if (x < 1e-20) else False
    # extinct_idxs = plot_df.index[plot_df['min_species'].map(extinct_func)]
    plot_df.drop(plot_df[plot_df['min_species'] < 1e-4].index, inplace=True)

    print(plot_df.loc[plot_df['sprott_LE'].idxmax()])
    print(plot_df.loc[plot_df['bene_LE'].idxmax()])

    n_bins = 50
    m_bins = np.linspace(0.1, 1.0, n_bins)
    K_1_bins = np.linspace(1.0, 20.0, n_bins)

    m_bins = list(np.around(m_bins,3))
    K_1_bins = list(np.around(K_1_bins,3))

    plot_metrics = ['sprott_LE', 'bene_LE', 'min_species']
    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4

    for metric in plot_metrics:
        fig, ax = plt.subplots(figsize=(13,10))

        binned_grid = stats.binned_statistic_2d(x=plot_df['K_1'], y=plot_df['m'], values=plot_df[metric], statistic='median', 
            bins=[K_1_bins, m_bins], range=None, expand_binnumbers=False)



        xbins = binned_grid.x_edge[:-1]
        ybins = binned_grid.y_edge[:-1]

        print(np.shape(xbins))
        print(np.shape(ybins))
        print(xbins)
        print(ybins)
        # exit()

        if metric == "xag":
            sns.heatmap(binned_grid.statistic, ax=ax, color='mako')

        else:
            sns.heatmap(binned_grid.statistic.T, ax=ax, cmap='mako',vmin=-0.01, vmax=0.04)
        
        for _, spine in ax.spines.items():
            spine.set_visible(True)

        ax.set_xticklabels(xbins, rotation=90)
        ax.set_yticklabels(ybins, rotation=0)

        ax.set_facecolor("#e30b17")
        ax.patch.set_alpha(0.1)
        ax.set_title(metric)
        ax.invert_yaxis()
        # ax.xaxis.set_major_locator(ticker.MultipleLocator(5))

        output_path = data_dir + metric + '_autocat.pdf'

        fig.tight_layout()

        plt.savefig(output_path, dpi=500)
        plt.close()


    # width_inches = 95*4 / 25.4
    # height_inches = 51*4 / 25.4

    # fig, ax = plt.subplots(figsize=(width_inches,height_inches))
    # sns.scatterplot(x="K_1", y="m", data=plot_df, hue = 'LE_classification', edgecolors=None, alpha = 0.7)
    # output_path = data_dir + 'LE_class_autocat.pdf'
    # plt.savefig(output_path, dpi=500)
    # plt.close()

def max_LE_hist():
    data_dir = "./output/autocat_6_rej/chunk_0/Population_end/"
    distances_df = pd.read_csv(data_dir + "model_sim_distances/model_0_distances.csv")
    params_df = pd.read_csv(data_dir + "model_sim_params/model_0_all_params.csv")

    plt.hist(distances_df['d1'], bins=100)
    plt.show()


def main():
    # max_LE_hist()
    # exit()
    make_LE_data()
    LE_heatmaps()
    exit()
    entropy_heatmap()
    RQA_plots()
    m_change()
    plot_figure_2()


if __name__ == "__main__":
    main()
