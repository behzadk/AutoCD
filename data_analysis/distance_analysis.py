import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib as mpl
from matplotlib import rcParams
from . import data_utils

import ternary
from scipy import stats
import math
from sklearn.neighbors import KernelDensity
from sklearn import mixture
from itertools import combinations
# from sklearn.linear_model import LinearRegression

import statsmodels.regression as sm_reg

from sklearn import svm
from scipy.stats.stats import pearsonr
import subprocess
from glob import glob
import os

plt.rcParams['figure.figsize'] = [15, 15]

font = {'size'   : 30, }
axes = {'labelsize': 'medium', 'titlesize': 'medium'}

sns.set_context("talk")
sns.set_style("white")
mpl.rc('font', **font)
mpl.rc('axes', **axes)
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
# plt.rcParams['text.usetex'] = True
plt.rcParams['axes.unicode_minus'] = False



def shannon_entropy(p):
    """Computes the Shannon Entropy at a distribution in the simplex."""
    s = 0.
    for i in range(len(p)):
        try:
            s += p[i] * math.log(p[i])
        except ValueError:
            continue
    return -1. * s

# def my_heatmapf(func, step=0.001, scale=10, boundary=True, cmap=None, ax=None,
#              style='triangular', permutation=None, vmin=None, vmax=None):
    

#     for i in np.arange(start=0, stop=scale, step=step)


# def simplex_iterator(scale=1.0, boundary=True, num_steps=50):
#     start = 0

#     for i in np.arange(start=start, stop=scale, step=step)
#         for j in np.arange(start=start, stop=scale + (1 - start) - i, step=step):
#             k = scale - i - j
#             yield (i, j, k)

class kde_func:
    def __init__(self, data, scale, step):
        # self.kde = kde = stats.gaussian_kde(data.T)        # 
        # self.sklearn_model = KernelDensity(bandwidth=0.1, kernel='gaussian', leaf_size=100)
        self.sklearn_model = KernelDensity(bandwidth=0.02, kernel='gaussian', leaf_size=200)

        # self.sklearn_model = mixture.GaussianMixture(n_components=3, covariance_type='diag')
        self.sklearn_model.fit(data)

        # self.sklearn_model = svm.LinearSVR(tol=1e-5)
        # print(data)
        # exit()
        self.sklearn_model.fit(data)

        self.scale = scale
        self.step = step
        self.v_max = 0

    def run_kde(self, x):
        x = np.array(x)
        x = x.reshape(1, -1)
        # ret = self.kde(x) * self.scale

        ret = self.sklearn_model.score_samples(x)
        # ret = self.sklearn_model.predict(x)

        # print(np.exp(ret))
        if ret[0] > self.v_max:
            self.v_max = ret[0]
            print(x)
            print(self.v_max)
            print("")

        return ret[0]



def plot_steady_state_ratio(output_dir, distance_columns, figure_output_dir):
    model_space_report_df = pd.read_csv(output_dir + "combined_model_space_report.csv")
    data_utils.make_folder(figure_output_dir)

    model_idxs = model_space_report_df['model_idx'].values
    model_idxs = sorted(model_idxs)

    model_distances_path_template = output_dir +  "combined_model_distances/model_#IDX#_distances.csv"
    scatter_output_path_template = figure_output_dir + "distance_scatter_tern_#IDX#.pdf"
    heatmap_output_path_template = figure_output_dir + "distance_heatmap_tern_#IDX#.pdf"

    print(model_distances_path_template)
    for model_idx in model_idxs:
        # if model_idx not in [4128, 226, 872, 2939, 1416]:
        if model_idx not in [4128, 4125, 4119, 3938]:
            continue

        try:
            model_distances_path =  model_distances_path_template.replace("#IDX#", str(model_idx))
            model_dist_df = pd.read_csv(model_distances_path)
        
        except FileNotFoundError:
            continue
        
        # Convert distances back to OD values

        if len(distance_columns) == 3:
            # print(model_dist_df[distance_columns])
            t_data = model_dist_df[distance_columns]
            t_data = t_data
            if len(t_data) < 100:
                continue

            t_data = t_data.applymap(lambda x: (1/x))
            # Convert to fractional populations
            t_data = t_data.div(t_data.sum(axis=1), axis=0)
            # t_data = t_data.applymap(lambda x: (x))

            model_dist_df = t_data


            t_data = t_data.values
            t_data = t_data.astype(np.float64)

            fontsize = 40
            offset = 0.04


            scale = int(100)
            func = kde_func(t_data, scale, step=0.1)

            figure, t_ax = ternary.figure(scale=scale)
            t_ax.boundary(linewidth=2.0)
            t_ax.gridlines(color="black", multiple=scale/10)
            # t_ax.gridlines(color="red", multiple=scale/4, linewidth=0.5)
            t_ax.ticks(axis='lbr', linewidth=1, multiple=scale/5, tick_formats="%d", fontsize=fontsize, offset=offset)
            t_ax.clear_matplotlib_ticks()
            t_ax.get_axes().axis('off')

            # Set Axis labels and Title
            # t_ax.set_title("Model: " + str(model_idx), fontsize=fontsize)

            t_ax.heatmapf(func.run_kde, scale=scale, boundary=True, style='hexagonal', colorbar=False, cmap='hot')

            figure_output_path = heatmap_output_path_template.replace("#IDX#", str(model_idx))
            figure.tight_layout()

            plt.savefig(figure_output_path, dpi=500)
            plt.close()

            ### Do scatter ###
            scale = int(100)
            t_data = model_dist_df[distance_columns]
            t_data = t_data.applymap(lambda x: (1/x))
            t_data = t_data.div(t_data.sum(axis=1), axis=0)
            print(t_data)


            model_dist_df.loc[  (model_dist_df['d3'] >  0.1)  &  
            (model_dist_df['d6'] >  0.1)  &  
            (model_dist_df['d9'] >  0.1) , 'target'] = '#F8766D'


            model_dist_df.loc[  (model_dist_df['d3'] <  0.1) |  
            (model_dist_df['d6'] <  0.1)  |
            (model_dist_df['d9'] <  0.1), 'target'] = '#00BFC4'

            plot_data = np.multiply(t_data.values, 100)
            figure, t_ax = ternary.figure(scale=scale)
            t_ax.boundary(linewidth=2.0)
            t_ax.gridlines(color="black", multiple=scale/10)
            # t_ax.gridlines(color="red", multiple=scale/4, linewidth=0.5)
            t_ax.ticks(axis='lbr', linewidth=1, multiple=scale/5, tick_formats="%d", fontsize=fontsize, offset=offset)
            t_ax.clear_matplotlib_ticks()
            t_ax.get_axes().axis('off')
            # Set Axis labels and Title
            fontsize = 20
            # t_ax.set_title("Model: " + str(model_idx), fontsize=fontsize)


            t_ax.scatter(points=plot_data, s=0.5, c=model_dist_df['target'])

            figure_output_path = scatter_output_path_template.replace("#IDX#", str(model_idx))
            # plt.show()
            figure.tight_layout()

            plt.savefig(figure_output_path, dpi=500)
            plt.close()


        else:
            fig, ax = plt.subplots()

            pop_distances_df = model_dist_df[distance_columns]

            pop_distances_df = pop_distances_df.applymap(lambda x: 1/x)
            pop_distances_df = pop_distances_df.assign(pop_ratio=lambda x: np.log(x[distance_columns[0]] /  x[distance_columns[1]] ))
            print(pop_distances_df)
            # sns.scatterplot(x=distance_columns[0], y=distance_columns[1], ax=ax, data=pop_distances_df)

            sns.distplot(pop_distances_df['pop_ratio'], hist=False, rug=False)

            plt.show()
        # exit()


def get_two_species_ratio(mat_row):
    return mat_row[0]/mat_row[1]

def get_abs_diff_populations(mat_row):
    row_combos = combinations(mat_row, r=2)
    
    diff = []
    for i in row_combos:
        diff.append((i[0] - i[1])**2)


    return sum(diff)


def ratio_dimensionality_reduction(output_dir, distance_columns, figure_output_dir):
    model_space_report_df = pd.read_csv(output_dir + "combined_model_space_report.csv")
    data_utils.make_folder(figure_output_dir)

    model_idxs = model_space_report_df['model_idx'].values
    model_idxs = sorted(model_idxs)

    model_distances_path_template = output_dir +  "combined_model_distances/model_#IDX#_distances.csv"
    scatter_output_path_template = figure_output_dir + "distance_scatter_tern_#IDX#.pdf"
    heatmap_output_path_template = figure_output_dir + "distance_heatmap_tern_#IDX#.pdf"

    df_list = []
    
    min_pop_bal = 2.0

    for model_idx in model_idxs:
        if model_idx not in [226, 872, 1294, 1416]:
            continue

        try:
            model_distances_path =  model_distances_path_template.replace("#IDX#", str(model_idx))
            model_dist_df = pd.read_csv(model_distances_path)
        except FileNotFoundError:
            continue
        
        if len(distance_columns) == 3:
            # print(model_dist_df[distance_columns])
            t_data = model_dist_df[distance_columns]
            if len(t_data) < 100:
                continue

            t_data = t_data.applymap(lambda x: (1/x))
            # Convert to fractional populations
            t_data = t_data.div(t_data.sum(axis=1), axis=0)
            t_data['pop_balance'] = np.apply_along_axis(get_abs_diff_populations, axis=1, arr=t_data.values)
            # print(t_data.iloc[0])
            # exit()
            t_data['model_idx'] = [model_idx for i in range(len(t_data))]

            q75, q25 = np.percentile(t_data['pop_balance'].values, [75 ,25])

            # print(q75, q25)

            if q25 < min_pop_bal:
                min_pop_bal = q25
                print(model_idx)
                print(q25)
                print("")

            df_list.append(t_data)


    master_df = pd.concat(df_list)
    output_path = output_dir + "motif_comparison_boxplot.pdf"
    fig, ax = plt.subplots()

    # sns.stripplot(x="pop_balance", y="model_idx", data=master_df, ax=ax, size=4,
    #     orient="h", zorder=10)
    sns.boxplot(x="pop_balance",  y="model_idx", data=master_df, ax=ax, orient="h",
        boxprops={'facecolor':'None'}, showfliers=False, linewidth=2, width=0.9)
    ax.set_xlabel('normalised marginal change')
    # ax.set_xscale('log')
    # ax.set_yticklabels('')
    ax.set_ylabel('')

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(True)

    ax.spines["bottom"].set_alpha(0.5)
    ax.spines["left"].set_alpha(0.5)

    ax.legend().remove()

    fig.tight_layout()
    plt.show()
            # plt.savefig(output_path, dpi=500)


def ratio_parameter_correlations(pop_dirs, output_dir, distance_columns, figure_output_dir):
    model_space_report_df = pd.read_csv(output_dir + "combined_model_space_report.csv")
    data_utils.make_folder(figure_output_dir)

    model_idxs = model_space_report_df['model_idx'].values
    model_idxs = sorted(model_idxs)

    model_distances_path_template = output_dir +  "combined_model_distances/model_#IDX#_distances.csv"
    model_params_path_template = output_dir +  "combined_model_params/model_#IDX#_population_all_params"

    for model_idx in model_idxs:
        if model_idx not in [4128, 4125, 4119, 3938]:
            continue

        model_distances_path =  model_distances_path_template.replace("#IDX#", str(model_idx))
        model_params_path =  model_params_path_template.replace("#IDX#", str(model_idx))

        try:
            model_dist_df = pd.read_csv(model_distances_path, index_col=0)
            model_params_df = pd.read_csv(model_params_path, index_col=0)

        except FileNotFoundError:
            continue

        # Merge distances and parameters on population number, batch_idx and sim_idx
        drop_columns = ['Accepted', 'particle_weight', 'integ_error', 'd1', 'd2', 'd4', 'd5', 'd7', 'd8']
        merge_columns = ['sim_idx', 'batch_idx', 'population_num', 'exp_num']
        merge_df = model_dist_df.merge(model_params_df, how='inner', on=merge_columns, left_index=False)
        merge_df = merge_df[merge_df.columns.difference(drop_columns)]

        params_df = merge_df[merge_df.columns.difference(merge_columns)]
        params_df = params_df[params_df.columns.difference(distance_columns)]
        clean_param_cols = []
        for c in params_df.columns:
            if min(params_df[c]) != max(params_df[c]):
                clean_param_cols.append(c)



        # Add ones colmn to dataset to represent b_0 constant of multiple linear regression
        # params_df['b_0'] = [1 for i in range(len(params_df))]
        param_values_array = params_df.values

        distances_df = merge_df[distance_columns]
        distances_df = distances_df.applymap(lambda x: (1/x))
        distances_df['pop_balance'] = np.apply_along_axis(get_two_species_ratio, axis=1, arr=distances_df.values)

        distance_values_array = distances_df['pop_balance'].values

        merge_df['pop_balance'] = distance_values_array
        fig, ax = plt.subplots(ncols=3, nrows=int(np.ceil(len(clean_param_cols)/3)))
        ax = np.concatenate(ax).ravel()
        for idx, x in enumerate(clean_param_cols):
            print(pearsonr(merge_df[x].values, merge_df['pop_balance'].values))
            sns.scatterplot(x=x, y='pop_balance', data=merge_df, ax=ax[idx])


            ax[idx].set_ylabel('')
            ax[idx].set_yscale('log')
        # plt.show()
        plt.savefig(figure_output_dir+"pop_bal_correlation.pdf")
        plt.close()
        # print(param_values_array.shape)
        # print(distance_values_array.shape)

        # optimal_parameters = list(range(len(params_df.columns)))
        # ommited_parameters = []
        # bad_p_values_exist = True

        # while bad_p_values_exist:
        #     X_opt = param_values_array[:, optimal_parameters]
        #     regress_OLS = sm_reg.linear_model.OLS(endog=distance_values_array, exog=X_opt).fit()
            
        #     print(regress_OLS.pvalues)

        #     max_p_val = 0 
        #     max_idx = 0
        #     for idx, p_val in zip(optimal_parameters, regress_OLS.pvalues):
        #         if p_val > max_p_val:
        #             max_p_val = p_val
        #             max_idx = idx

        #     if max_p_val <= 0.05:
        #         bad_p_values_exist = False
        #         break

        #     else:
        #         print(max_idx)
        #         ommited_parameters.append(max_idx)
        #         optimal_parameters = [i for i in optimal_parameters if i != max_idx]

        # print(regress_OLS.pvalues)
        # print(optimal_parameters)



def KS_test_distance_subset(pop_dirs, output_dir, inputs_dir, distance_columns, R_script, figure_output_dir):
    model_space_report_df = pd.read_csv(output_dir + "combined_model_space_report.csv")
    data_utils.make_folder(figure_output_dir)

    model_idxs = model_space_report_df['model_idx'].values
    model_idxs = sorted(model_idxs)


    model_distances_path_template = output_dir +  "combined_model_distances/model_#IDX#_distances.csv"
    model_params_path_template = output_dir +  "combined_model_params/model_#IDX#_population_all_params"

    figure_output_name_template = "cum_dist_KS_model_#IDX#.pdf"

    for model_idx in model_idxs:
        if model_idx not in [4119]:
            continue

        output_name = figure_output_dir + "model_" + str(model_idx) + "_2D_dual.pdf"
        model_distances_path =  model_distances_path_template.replace("#IDX#", str(model_idx))
        model_params_path =  model_params_path_template.replace("#IDX#", str(model_idx))

        model_input_parameter_path = inputs_dir + "input_files/params_" + str(model_idx) + ".csv"
        model_input_species_path = inputs_dir + "input_files/species_" + str(model_idx) + ".csv"

        try:
            model_dist_df = pd.read_csv(model_distances_path, index_col=0)

        except FileNotFoundError:
            continue

        model_params_df = pd.read_csv(model_params_path, index_col=0)

        # Merge distances and parameters on population number, batch_idx and sim_idx
        drop_columns = ['Accepted', 'integ_error', 'd1', 'd2', 'd4', 'd5', 'd7', 'd8']
        merge_columns = ['sim_idx', 'batch_idx', 'population_num', 'exp_num']
        
        merge_df = model_dist_df.merge(model_params_df, how='inner', on=merge_columns, left_index=False)
        merge_df = merge_df[merge_df.columns.difference(drop_columns)]

        params_df = merge_df[merge_df.columns.difference(merge_columns)]
        params_df = params_df[params_df.columns.difference(distance_columns)]
        param_columns = params_df.columns
        
        clean_param_cols = []
        for c in params_df.columns:
            if min(params_df[c]) != max(params_df[c]):
                clean_param_cols.append(c)

        distances_df = merge_df[distance_columns]
        distances_df = distances_df.applymap(lambda x: (1/x))
        distances_df['pop_balance'] = np.apply_along_axis(get_abs_diff_populations, axis=1, arr=distances_df.values)

        merge_df['pop_balance'] = distances_df['pop_balance'].values
        distance_values_array = distances_df['pop_balance'].values

        merge_df[distance_columns] = merge_df[distance_columns].applymap(lambda x: (1/x))

        # target_pop_balance = merge_df.loc[ (merge_df['pop_balance'] < 0.5)]
        target_pop_balance = merge_df.loc[ (merge_df['d3'] >  0.1)]
        target_pop_balance = target_pop_balance.loc[ (target_pop_balance['d6'] > 0.1)]
        target_pop_balance = target_pop_balance.loc[ (target_pop_balance['d9']  >  0.1)]

        anti_target_pop_balance = merge_df.loc[ (merge_df['d3'] <  0.1) | (merge_df['d6'] <  0.1) | (merge_df['d9'] <  0.1) ]


        # print(anti_target_pop_balance)
        # exit()
        # anti_target_pop_balance = merge_df

        target_pop_balance = target_pop_balance[model_params_df.columns]
        anti_target_pop_balance = anti_target_pop_balance[model_params_df.columns]

        if len(target_pop_balance) <= 5:
            continue


        plot_param_cols = ["D", "kB_max_1", "kB_max_2", "kB_max_3"]

        fig, axes = plt.subplots(ncols=2, nrows=int(np.ceil(len(plot_param_cols)/2)))#, sharey='row')
        axes = axes.reshape(-1)

        D_crit_vals = []
        KS_vals = []
        target_pop_data = []
        all_pop_data = []


        # Save as tmp csv
        anti_target_pop_balance.to_csv(output_dir + "params_1_tmp.csv")
        target_pop_balance.to_csv(output_dir + "params_2_tmp.csv")

        output_pdf = os.path.abspath(output_name)

        subprocess.call(['Rscript', R_script, os.path.abspath(output_dir + "params_1_tmp.csv"), os.path.abspath(output_dir + "params_2_tmp.csv"), 
            model_input_parameter_path, model_input_species_path, 
            str(model_idx), output_pdf])

        # Delete tmp csv

        max_param = ["NONE", 0]
        for idx, param_name in enumerate(plot_param_cols):
            target_pop_param_vals = target_pop_balance[param_name].values
            all_pop_param_vals = anti_target_pop_balance[param_name].values

            # log_scale = False

            # if min(target_pop_param_vals) < 0:
            #     log_scale = True
            #     target_pop_param_vals = np.exp(target_pop_param_vals)
            #     all_pop_param_vals = np.exp(all_pop_param_vals)

            D_crit = 1.36 * math.sqrt(1 / len(all_pop_param_vals) + 1 / len(target_pop_param_vals))
            D_crit_vals.append(D_crit)

            sp_ks = stats.ks_2samp(all_pop_param_vals, target_pop_param_vals)
            KS_vals.append(sp_ks[0])

            target_pop_data.append(target_pop_param_vals)
            all_pop_data.append(all_pop_data)

            # sns.kdeplot(target_pop_param_vals, cumulative=True, ax=axes[idx])
            # sns.kdeplot(all_pop_param_vals, cumulative=True, ax=axes[idx])

            if sp_ks[0] > D_crit:
                print(param_name, sp_ks)

            colour_blue = '#00BFC4'
            colour_pink = '#F8766D'

            sns.distplot(target_pop_param_vals, norm_hist=True, kde=True, hist=False, ax=axes[idx], label='balanced_pop', color=colour_pink)
            sns.distplot(all_pop_param_vals, norm_hist=True, kde=True, hist=False, ax=axes[idx], label='all_data', color=colour_blue)

            # axes[idx].text(0.02, 0.98, param_name + " - " + str(np.format_float_positional(sp_ks[0], precision=3, fractional=False, trim='k')), fontsize=9, ha="left", va="top", transform=axes[idx].transAxes)
        
            # axes[idx].set(xticklabels=[])

            # if log_scale:
            #     axes[idx].set_xscale('log')

            if min(target_pop_param_vals) < 0:
                axes[idx].set_ylim(0, 0.25)

            # axes[idx].set(xlabel='')
            axes[idx].spines["right"].set_visible(False)
            axes[idx].spines["top"].set_visible(False)
            axes[idx].spines["left"].set_alpha(0.5)
            axes[idx].spines["bottom"].set_alpha(0.5)
            axes[idx].legend().remove()


        if max_param[0] == "D":
            print("DILUTION IS MAX: ", model_idx)

        # handles, labels = axes[0].get_legend_handles_labels()
        # fig.legend(handles, labels, loc='lower center')
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0.35)
        plt.savefig(figure_output_dir + figure_output_name_template.replace('#IDX#', str(model_idx)))
        plt.close()



def test():
    mu=np.array([1,10,20])
    # Let's change this so that the points won't all lie in a plane...
    sigma=np.matrix([[20,10,10],
                     [10,25,1],
                     [10,1,50]])

    data=np.random.multivariate_normal(mu,sigma,1000)
    values = data.T

    kde = stats.gaussian_kde(values)


if __name__ == "__main__":
    test_dict = {"A": [0, 5, 7],
                "B": [9, 9, 9], 
                "C": [1000, 1000, 1000]}

    df = pd.DataFrame(test_dict)
    print(df)
    column_names = ["A", "B", "C"]
    x = df.div(df.sum(axis=1), axis=0)
    # x = 

    print(x)
