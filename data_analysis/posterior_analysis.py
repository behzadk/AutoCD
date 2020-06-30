import csv
from sklearn.preprocessing import MinMaxScaler
import sklearn
import math
import matplotlib.style as style
import os
import glob
from . import data_utils
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import seaborn as sns

plt.rcParams['figure.figsize'] = [30, 30]

# font = {'size': 25, }
# axes = {'labelsize': 'medium', 'titlesize': 'medium'}

# mpl.rc('font', **font)
# mpl.rc('axes', **axes)


def generate_file_paths(posterior_dir, priors_dir):
    prior_param_paths = [file_path for file_path in glob.iglob(priors_dir + "*.csv")]
    posterior_param_paths = [file_path for file_path in glob.iglob(posterior_dir + "*_all_params")]
    num_models = len(posterior_param_paths)

    ordered_prior_species_paths = [0 for f in range(num_models)]
    ordered_prior_param_paths = [0 for f in range(num_models)]
    ordered_posterior_paths = [0 for f in range(num_models)]

    # Split prior paths into species and parameter priors in model number order
    for f in prior_param_paths:
        file_name = os.path.basename(f)
        split_name = file_name.split('_')
        model_num = int(''.join(list(filter(str.isdigit, split_name[1]))))

        if split_name[0] == "species":
            ordered_prior_species_paths[model_num] = f

        if split_name[0] == "params":
            ordered_prior_param_paths[model_num] = f

    for f in posterior_param_paths:
        file_name = os.path.basename(f)
        split_name = file_name.split('_')
        model_num = int(split_name[1])
        ordered_posterior_paths[model_num] = f

    return ordered_posterior_paths, ordered_prior_param_paths, ordered_prior_species_paths


def import_input_file(input_path):
    data_dict = {}

    with open(input_path) as fin:
        reader = csv.reader(fin, skipinitialspace=True)

        for row in reader:
            data_dict[row[0]] = [float(i) for i in row[1:]]

    return data_dict


def generate_posterior_KS_csv(data_dir, priors_dir, output_dir):
    print("Generating posterior KS csv files ...")
    posterior_dir = data_dir + "model_sim_params/"
    ordered_posterior_paths, \
        ordered_prior_param_paths, \
        ordered_prior_species_paths = generate_file_paths(posterior_dir, priors_dir)

    out_path = output_dir + "model_NUM_KS.csv"

    for model_idx, f in enumerate(ordered_posterior_paths):
        model_posterior_df = pd.read_csv(f, sep=',')
        KS_df = data_utils.make_KS_df(model_idx, model_posterior_df)
        if KS_df is None:
            continue

        KS_df.to_csv(out_path.replace('NUM', str(model_idx)))

    print("\n")

def generate_posterior_KS_csv_ABCSMC(final_pop_dir, first_pop_dir, priors_dir, output_dir):
    print("Generating posterior KS csv files ...")
    first_posterior_dir = first_pop_dir + "model_sim_params/"
    final_posterior_dir = final_pop_dir + "model_sim_params/"

    final_pop_ordered_posterior_paths, _, _ = generate_file_paths(final_posterior_dir, priors_dir)
    first_pop_ordered_posterior_paths, _, _ = generate_file_paths(first_posterior_dir, priors_dir)

    out_path = output_dir + "model_NUM_KS.csv"

    for model_idx, (f_1, f_2) in enumerate(zip(final_pop_ordered_posterior_paths, first_pop_ordered_posterior_paths)):
        model_posterior_df = pd.read_csv(f_1, sep=',')
        model_prior_df = pd.read_csv(f_2, sep=',')

        KS_df = data_utils.make_KS_df_alt(model_idx, model_posterior_df, model_prior_df)
        if KS_df is None:
            continue

        KS_df.to_csv(out_path.replace('NUM', str(model_idx)))

    print("\n")



def main():
    # Two species

    wd = "/home/behzad/Documents/barnes_lab/cplusplus_software/speed_test/repressilator/cpp/"
    data_dir = wd + "output/two_species_stable_6/Population_0/"

    posterior_params_dir = data_dir + "model_sim_params/"
    priors_dir = wd + "input_files_two_species/priors/"
    output_dir = wd + "data_analysis_notebook/posterior_analysis/"
    distributions_dir = output_dir + "distributions/"
    KS_out_dir = output_dir + "KS_data/"

    # generate_posterior_KS_csv(posterior_params_dir, priors_dir, KS_out_dir)
    visualise_posterior_distributions(posterior_params_dir, priors_dir, distributions_dir)
    exit()

    # Three species
    wd = "/home/behzad/Documents/barnes_lab/cplusplus_software/speed_test/repressilator/cpp/"
    data_dir = wd + "output/three_species_stable_comb/Population_0/"

    posterior_params_dir = data_dir + "model_sim_params/"
    priors_dir = wd + "input_files_three_species/priors/"
    output_dir = wd + "data_analysis_notebook/posterior_analysis/three_species/"
    distributions_dir = output_dir + "distributions/"
    KS_out_dir = output_dir + "KS_data/"

    # generate_posterior_KS_csv(posterior_params_dir, priors_dir, KS_out_dir)
    visualise_posterior_distributions(posterior_params_dir, priors_dir, distributions_dir)

if __name__ == "__main__":
    print("hello world")
    main()
