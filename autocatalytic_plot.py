import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from matplotlib.ticker import FormatStrFormatter

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
    data_dir = "./output/autocat_2_rej/autocat_2_rej_4/Population_0/"
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
    'LAM': distances_df['d8'].values
    }
    plot_df = pd.DataFrame(plot_data)

    n_bins = 30

    m_bins = np.linspace(0.1, 1.0, n_bins)
    K_1_bins = np.linspace(1.0, 20.0, n_bins)

    m_bins = list(np.around(m_bins,2))
    K_1_bins = list(np.around(K_1_bins,2))

    plot_metrics = ['sep_coeff', 'DET', 'ENT', 'LAM']

    for metric in plot_metrics:
        fig, ax = plt.subplots(figsize=(13,10)) 

        binned_grid = stats.binned_statistic_2d(x=plot_df['K_1'], y=plot_df['m'] ,  values=plot_df[metric], statistic='mean', bins=[K_1_bins, m_bins], range=None, expand_binnumbers=False)
        sns.heatmap(binned_grid.statistic.T, ax=ax)
        ax.set_xticklabels(K_1_bins, rotation=90)
        ax.set_facecolor("yellow")
        ax.set_title(metric)


        ax.invert_yaxis()
        ax.set_yticklabels(m_bins, rotation=0)

        output_path = data_dir + metric + '_autocat.pdf'

        plt.savefig(output_path, dpi=500)

def RQA_plots():
    # https://www.pnas.org/content/pnas/86/1/142.full.pdf
    data_dir = "./output/autocat_1_rej/autocat_1_rej_1/Population_0/"
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


def main():
    entropy_heatmap()
    RQA_plots()
    m_change()
    plot_figure_2()


if __name__ == "__main__":
    main()
