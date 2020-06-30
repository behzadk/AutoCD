import math
from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set()

import scipy
import os
from . import symmetry_analysis
from itertools import combinations

def make_folder(folder_dir):
    try:
        os.mkdir(folder_dir)

    except FileExistsError:
        pass


def KL_divergence(sample_one, sample_two):
    n_1 = len(sample_one)
    n_2 = len(sample_two)

    sample_one.sort()
    sample_two.sort()

    sample_one_cdf = [i*1 / n_1 for i in range(1, n_1+1)]
    sample_two_cdf = [i*1 / n_2 for i in range(1, n_2+1)]

    entropy = scipy.stats.entropy(sample_one_cdf, sample_two_cdf)
    return entropy


def kolmogorov_smirnov_test(sample_one, sample_two):
    n_1 = len(sample_one)
    n_2 = len(sample_two)

    sample_one.sort()
    sample_two.sort()

    sample_one_cdf = [i*1 / n_1 for i in range(1, n_1+1)]
    sample_two_cdf = [i*1 / n_2 for i in range(1, n_2+1)]

    s_one_idx = 0
    s_two_idx = 0

    combined_sample = []

    max_diff = 0


    while s_one_idx < n_1 and s_two_idx < n_2:
        v1 = sample_one[s_one_idx]
        v2 = sample_two[s_two_idx]

        print(v1)
        print(v2)
        print(sample_one_cdf[s_one_idx])
        print(sample_two_cdf[s_two_idx])
        print("")

        diff = abs(sample_one_cdf[s_one_idx] - sample_two_cdf[s_two_idx])
        if diff > max_diff:
            max_diff = diff

        if v1 < v2:
            combined_sample.append(v1)
            s_one_idx += 1

        elif v2 < v1:
            combined_sample.append(v2)
            s_two_idx += 1

        else:
            combined_sample.append(v1)
            s_one_idx += 1
            s_two_idx += 1


    return max_diff


def normalise_parameters(model_posterior_df):
    # Get param names
    param_names = model_posterior_df.columns[3:]
    free_params = []

    # Extract parameter columns that are not constants
    for param in param_names:
        if model_posterior_df[param].nunique() == 1:
            continue

        else:
            free_params.append(param)

    # Normalise free parameter columns
    for param in free_params:
        scaler = MinMaxScaler()
        model_posterior_df[param] = scaler.fit_transform(model_posterior_df[[param]])

    return model_posterior_df, free_params


def recur_get_connected_nodes(adj_mat, path_length, path_sign, visited_paths, from_idx, target_idx, loop_found, path_log, output, path_log_list):
    n_species = np.shape(adj_mat)[0]
    neighbours = [i for i in range(n_species) if adj_mat[i, from_idx] != 0]

    if loop_found:
        return 0

    if from_idx is target_idx:
        print("path log: \t", path_log, "current sign: \t", path_sign)
        loop_found = True
        path_log_list.append(path_log)
        # path_log.clear()
        # path_log.append(from_idx)
        # path_sign = 1
        output.append([path_sign, path_length])
        return 0

    for n in neighbours:
        if [n, from_idx] not in visited_paths:
            branched_loop_found = False
            branch_visited_nodes = visited_paths.copy()
            branch_path_length = path_length + 1
            branch_path_sign = path_sign * adj_mat[n, from_idx]
            branch_path_log = path_log.copy()
            branch_path_log.append(n * path_sign)

            # if n != target_idx:
            branch_visited_nodes.append([n, from_idx])

            recur_get_connected_nodes(adj_mat, branch_path_length, branch_path_sign, branch_visited_nodes, n, target_idx, branched_loop_found, branch_path_log, output, path_log_list)

    return path_log_list


def get_two_strain_symmetric_adj_mats(model_idxs, adj_mat_dir):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    
    symmetrical_col = []

    for m_idx in model_idxs:
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)
        adj_mat_df = adj_mat_df.drop("S_glu", axis=1).drop("S_glu", axis=0)
        new_adj_mat = symmetry_analysis.rearrange_matrix_two_strain(adj_mat_df)
        symmetrical = symmetry_analysis.check_matrix_symmetrical(new_adj_mat)

        symmetrical_col.append(symmetrical)

    return symmetrical_col


def get_n_strain_symmetric_adj_mats(model_idxs, adj_mat_dir):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    
    symmetrical_col = []

    for m_idx in model_idxs:
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)
        adj_mat_df = adj_mat_df.drop("S_glu", axis=1).drop("S_glu", axis=0)

        strain_col = [col for col in adj_mat_df if col.startswith('N_')]
        # strain_pairs = [x for x in combinations(strain_col, 2)]

        symmerical = True
        for strain_pair in combinations(strain_col, 2):
            tmp_adj_mat = adj_mat_df.copy(deep=True)
            # print(strain_pair)
            for x in strain_col:
                if x not in strain_pair:
                    tmp_adj_mat = tmp_adj_mat.drop(x, axis=1).drop(x, axis=0)

            new_adj_mat = symmetry_analysis.rearrange_matrix_three_strain(tmp_adj_mat)
            symmetrical = symmetry_analysis.check_matrix_symmetrical(new_adj_mat)
            print(symmetrical)
            if symmetrical is False:
                break

        if symmetrical:
            print(adj_mat_path)


def translate_param_names(param_names):
    translated_param_names = []

    for param in param_names:
        split_name = param.split('_')
        num_splits = len(split_name)

        if num_splits > 1:
            new_name = split_name[0]

            for split in split_name[1:]:
                new_name = new_name + "_{" + split


            new_name = new_name + "}" * (num_splits - 1)

        else:
            new_name = split_name[0]

        new_name = new_name.replace("}{", "}_{")
        new_name = "$" + new_name + "$"
        new_name = new_name.replace("omega", "\omega")
        new_name = new_name.replace("mu", "\mu")
        # new_name = new_name.replace("B", "\beta")

        translated_param_names.append(new_name)

    return translated_param_names

# Get model idxs that both produce a species and are influenced by the species directly 
def get_self_regulators(model_idxs, target_species_list, adj_mat_dir):
    # Get all species beginning with N
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"

    self_reg_loops = []

    for m_idx in model_idxs:
        self_reg_loop_count = 0
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)
        get_feedback_loops(adj_mat_df)
        # Get columns beginning with N_
        strain_col = [col for col in adj_mat_df if col.startswith('N_')]

        for strain in strain_col:
            for col in adj_mat_df.columns:
                if adj_mat_df[strain][col] == 1 and adj_mat_df[col][strain] == -1:
                    self_reg_loop_count += 1

        self_reg_loops.append(self_reg_loop_count)

    return self_reg_loops

def get_adj_mat_distances(ref_model_idx, model_idxs, adj_mat_dir):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"

    ref_adj_mat = pd.read_csv(adj_matrix_path_template.replace("#REF#", str(ref_model_idx)), index_col=0).as_matrix()
    mat_distances_list = []

    for m_idx in model_idxs:
        candidate_adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        candidate_adj_mat_df = pd.read_csv(candidate_adj_mat_path, index_col=0)

        candidate_adj_mat = candidate_adj_mat_df.as_matrix()
        distance = np.sum(abs(candidate_adj_mat - ref_adj_mat))
        mat_distances_list.append(distance)

    return mat_distances_list


# Counts feedback existing between all.
def get_feedback_loops(adj_mat_df):
    # Drop row names column
    adj_mat_df = adj_mat_df.drop("S_glu", axis=1).drop("S_glu", axis=0)
    print(adj_mat_df)
    col_names = adj_mat_df.columns
    strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
    adj_mat = adj_mat_df.values

    # Remove column 0, which contains row names
    # adj_mat = adj_mat[:, 1:]
    n_species = np.shape(adj_mat)[0]
    path_log_list = []
    output = []

    for strain_idx in strain_indexes:
        for row_idx in range(n_species):

            # Find strain to species interaction
            if adj_mat[row_idx, strain_idx] != 0:

                loop_closed = False
                path_length = 1

                from_idx = row_idx
                path_log = [strain_idx, row_idx]
                path_sign = adj_mat[row_idx, strain_idx]
                strain_feedback_loops = []
                recur_get_connected_nodes(adj_mat, path_length, path_sign, [strain_idx, row_idx], from_idx, strain_idx, False, path_log, strain_feedback_loops, path_log_list)
                
                print(strain_feedback_loops)
                for l in strain_feedback_loops:
                    l.append(strain_idx)
                    output.append(l)

    return output, path_log_list

def get_num_interactions(model_idxs, adj_mat_dir):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"

    sum_interactions = []
    for m_idx in model_idxs:
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)

        adj_mat = adj_mat_df.as_matrix()
        adj_mat = abs(adj_mat)
        sum_interactions.append(np.sum(adj_mat))

    return sum_interactions

def get_num_parts(adj_mat_df):
    # Drop row names column
    adj_mat_df.drop([adj_mat_df.columns[0]], axis=1, inplace=True)

    col_names = adj_mat_df.columns

    strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
    AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
    microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

    adj_mat = adj_mat_df.values

    n_species = np.shape(adj_mat)[0]
    total_num_parts = 0

    for strain_idx in strain_indexes:
        total_num_parts += abs(sum(adj_mat[:, strain_idx]))

    num_AHL_parts = 0
    for idx in AHL_indexes:
        if sum(adj_mat[:, idx]) >= 1:
            num_AHL_parts += 1

    num_microcin_parts = 0
    for idx in microcin_indexes:
        if abs(sum(adj_mat[:, idx])) >= 1:
            num_microcin_parts += 1

    return total_num_parts, num_AHL_parts, num_microcin_parts


def make_KS_df(model_idx, model_posterior_df):
    accepted_sims = model_posterior_df.loc[model_posterior_df['Accepted'] == True]

    if len(accepted_sims) <= 1:
        return None
        
    model_posterior_df, free_params = normalise_parameters(model_posterior_df)

    D_crit = 1.36 * math.sqrt(1 / len(accepted_sims) + 1 / len(model_posterior_df))

    KS_data_df = pd.DataFrame(columns=['model_idx', 'D_crit'])
    KS_data_df['model_idx'] = [model_idx]
    KS_data_df['D_crit'] = [D_crit]

    for param in free_params:
        D_n = kolmogorov_smirnov_test(accepted_sims[param].values, model_posterior_df[param].values)

        KS_data_df[param] = D_n

    return KS_data_df


def make_KS_df_alt(model_idx, model_posterior_df, model_prior_df):
    accepted_sims = model_posterior_df.loc[model_posterior_df['Accepted'] == True]

    if len(accepted_sims) <= 1:
        return None
        
    # model_posterior_df, free_params = normalise_parameters(model_posterior_df)
    # model_prior_df, free_params = normalise_parameters(model_prior_df)

    accepted_sims = model_posterior_df.loc[model_posterior_df['Accepted'] == True]

    D_crit = 1.36 * math.sqrt(1 / len(accepted_sims) + 1 / len(model_prior_df))

    KS_data_df = pd.DataFrame(columns=['model_idx', 'D_crit'])
    KS_data_df['model_idx'] = [model_idx]
    KS_data_df['D_crit'] = [D_crit]

    param_names = model_posterior_df.columns[3:]

    for param in param_names:

        if param == "particle_weight":
            continue

        min_val = min(model_posterior_df[param].values)
        max_val = max(accepted_sims[param].values)

        if min_val == max_val:
            continue
        # print(accepted_sims[param].values)
        # print(model_posterior_df[param].values)
        # exit()

        D_n = kolmogorov_smirnov_test(accepted_sims[param].values, model_posterior_df[param].values)
        KS_data_df[param] = D_n

    return KS_data_df


def make_num_parts(model_space_report_df, adj_matrix_path_template):
    models = model_space_report_df.model_idx.values
    all_num_parts = []
    all_AHL_num_parts = []
    all_microcin_num_parts = []

    for m in models:
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m))

        adj_mat_df = pd.read_csv(adj_mat_path)
        total_num_parts, num_AHL_parts, num_microcin_parts = get_num_parts(adj_mat_df)
        all_num_parts.append(total_num_parts)
        all_AHL_num_parts.append(num_AHL_parts)
        all_microcin_num_parts.append(num_microcin_parts)

    return all_num_parts, all_AHL_num_parts, all_microcin_num_parts


def make_num_species(model_space_report_df, adj_matrix_path_template):
    models = model_space_report_df.model_idx.values
    total_species = []

    for m in models:
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m))

        adj_mat_df = pd.read_csv(adj_mat_path)
        adj_mat_df.drop([adj_mat_df.columns[0]], axis=1, inplace=True)

        col_names = adj_mat_df.columns

        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        num_AHL = 0
        num_microcin = 0

        adj_mat_vals = adj_mat_df.values

        for a in AHL_indexes:
            for s in strain_indexes:
                if adj_mat_vals[a, s] == 1:
                    num_AHL += 1
                    break

        for m in microcin_indexes:
            for s in strain_indexes:
                if adj_mat_vals[m, s] == 1:
                    num_microcin += 1
                    break
        
        total_species.append(num_AHL + num_microcin)

    return total_species


def make_num_parts_alt(model_space_report_df, adj_matrix_path_template):
    models = model_space_report_df.model_idx.values
    adj_mat_sum = []

    for m in models:
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m))

        adj_mat_df = pd.read_csv(adj_mat_path)
        adj_mat_df.drop([adj_mat_df.columns[0]], axis=1, inplace=True)

        col_names = adj_mat_df.columns

        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat_sum.append(np.sum(np.abs(adj_mat_df.values)))

    return adj_mat_sum


def make_num_expressed_parts(model_space_report_df, adj_matrix_path_template):
    models = model_space_report_df.model_idx.values
    all_num_parts = []
    all_AHL_num_parts = []
    all_microcin_num_parts = []

    min_x = 100

    for m in models:
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m))


        adj_mat_df = pd.read_csv(adj_mat_path)
        col_names = adj_mat_df.columns


        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]

        strain_expression = [adj_mat_df.values[:, x] for x in strain_indexes]

        strain_num_parts = [sum(x) for x in strain_expression]
        min_x = min(min(strain_num_parts), min_x)
        model_parts = np.sum(strain_expression)

        all_num_parts.append(model_parts)

    return all_num_parts


def make_feedback_loop_counts(data_dir, input_files_dir):
    model_space_report_df = pd.read_csv(data_dir + "model_space_report.csv")
    adj_mat_name_template = "model_#REF#_adj_mat.csv"
    adj_matrix_path_template = input_files_dir + "adj_matricies/" + adj_mat_name_template

    models = model_space_report_df.model_idx.values
    positive_loops = []
    negative_loops = []

    for m in models:
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m))
        adj_mat_df = pd.read_csv(adj_mat_path)
        if m != 33:
            continue

        feedback_loops = get_feedback_loops(adj_mat_df)
        positive_count = 0
        negative_count = 0

        for loop in feedback_loops:
            loop_sign = loop[0]
            loop_length = loop[1]

            if loop_sign == 1:
                positive_count += 1#/loop_length

            if loop_sign == -1:
                negative_count += 1#/loop_length

            if loop_sign != 1 and loop_sign != -1:
                print("Loop sign is wrong: ", loop_sign)

        positive_loops.append(positive_count)
        negative_loops.append(negative_count)

    return positive_loops, negative_loops


def make_strain_feedback_loop_balance(model_idxs, adj_mat_dir):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"

    loop_balances = []
    sum_positive = []
    sum_negative = []
    model_path_log_lists = []
    for m_idx in model_idxs:

        total_pos = 0
        total_neg = 0

        strain_0_loops = []
        strain_1_loops = []

        self_reg_loop_count = 0
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)


        model_path_log_lists
        loops_out, paths_log_out = get_feedback_loops(adj_mat_df)
        model_path_log_lists.append(paths_log_out)
        for loop in loops_out:
            if loop[0] == 1:
                total_pos += 1

            if loop[0] == -1:
                total_neg += 1

            if loop[2] == 0:
                strain_0_loops.append(loop)

            else:
                strain_1_loops.append(loop)

        strain_0_bal = sum([l[0] for l in strain_0_loops])
        strain_1_bal = sum([l[0] for l in strain_1_loops])

        loop_balances.append(abs(strain_0_bal - strain_1_bal))
        sum_positive.append(total_pos)
        sum_negative.append(total_neg)

    return loop_balances, sum_positive, sum_negative, model_path_log_lists


def generate_replicates_and_std(all_sims_df, model_space_report_df, num_replicates):
    models = model_space_report_df.model_idx.values
    batches = all_sims_df.batch_num.unique()
    batch_sim_indexes = all_sims_df.sim_idx.unique()

    total_simulations = len(all_sims_df)

    all_sims_df['replicate'] = pd.qcut(
        range(total_simulations), num_replicates, labels=np.arange(num_replicates))

    # Create replicate subset
    for rep_num in range(num_replicates):
        rep_subset_df = all_sims_df.loc[all_sims_df['replicate'] == rep_num]
        total_accepted = len(rep_subset_df.loc[rep_subset_df['Accepted'] == True])
        print("Replicate: ", rep_num, " ",
              " Simulations: ", len(rep_subset_df.index))


        replicate_acceptance_ratios = []

        # Calculate acceptance ratio for each model within a subset
        for m in models:
            models_subset_df = rep_subset_df.loc[rep_subset_df['model_ref'] == m]
            num_accepted = len(
                models_subset_df.loc[models_subset_df['Accepted'] == True].index)
            num_rejected = len(
                models_subset_df.loc[models_subset_df['Accepted'] == False].index)

            try:
                acceptance_ratio = num_accepted / total_accepted
            except(ZeroDivisionError):
                acceptance_ratio = 0

            replicate_acceptance_ratios.append(acceptance_ratio)

        replicate_name = 'replicate_NUM'.replace('NUM', str(rep_num))
        model_space_report_df[replicate_name] = replicate_acceptance_ratios

    # Calculate standard deviation between replicates
    model_standard_deviations = []
    mean_acceptance_ratios = []
    for m in models:
        model_subset_df = model_space_report_df.loc[model_space_report_df['model_idx'] == m]
        model_replicate_acceptance_ratios = []

        for rep_num in range(num_replicates):
            replicate_name = 'replicate_NUM'.replace('NUM', str(rep_num))
            replicate_ratio = model_subset_df[replicate_name].values[0]
            model_replicate_acceptance_ratios.append(replicate_ratio)

        model_standard_deviations.append(
            np.std(model_replicate_acceptance_ratios))

    model_space_report_df['stdev'] = model_standard_deviations
    sims_per_replicate = len(all_sims_df.loc[all_sims_df['replicate'] == 0])

    return model_space_report_df, sims_per_replicate


##
# Adds new bool column for species 1 and 2 indicating whether the species was susitaned or not
# by the end of the simulation. Uses a threshold where above that number, the species is sustained.
#
##
def species_sustained(df, threshold=1e3):

    N_1_final = df.d3
    N_1_sustained = []

    for i in N_1_final:
        if i >= threshold and i < 1e100:
            N_1_sustained.append(True)

        else:
            N_1_sustained.append(False)

    N_2_final = df.d6
    N_2_sustained = []

    for i in N_2_final:
        if float(i) >= threshold and i < 1e100:
            N_2_sustained.append(True)

        else:
            N_2_sustained.append(False)

    df['N_1_sustained'] = N_1_sustained
    df['N_2_sustained'] = N_2_sustained

    return df


def make_max_eig(df):
    # Make subset containing only eigenvalues
    col_names = df.columns
    eign_cols = [x for x in col_names if 'eig' in x]
    eign_cols = [x for x in eign_cols if 'real' in x]

    eig_df = df[eign_cols]
    max_real_eigs = []

    for row in eig_df.values:
        abs_vals = [abs(x) for x in row]
        try:
            max_val_idx = np.nanargmax(abs_vals)
            max_real_eigs.append(row[max_val_idx])

        except(ValueError):
            max_real_eigs.append(np.nan)

    # print(max_real_eigs)
    max_real_eigs

    return max_real_eigs


def make_sum_eig(df):
    # Make subset containing only eigenvalues
    col_names = df.columns
    eign_cols = [x for x in col_names if 'eig' in x]
    eign_cols = [x for x in eign_cols if 'real' in x]

    eig_df = df[eign_cols]
    all_sum_eigs = []

    for row in eig_df.values:
        sum_eig = 0
        for val in row:
            if pd.isnull(val):
                continue

            else:
                sum_eig = sum_eig + val

        all_sum_eigs.append(sum_eig)

    # print(max_real_eigs)
    all_sum_eigs

    return all_sum_eigs


##
# Adds a column noting whether all the real part eigenvalues are negative
#
##
def all_negative_eigs(df):
    # Make subset containing only eigenvalues
    col_names = df.columns
    eign_cols = [x for x in col_names if 'eig' in x]
    eign_cols = [x for x in eign_cols if 'real' in x]
    eig_df = df[eign_cols]
    all_negative = []

    for row in eig_df.values:
        row_negative = True
        for val in row:
            if pd.isnull(val):
                continue

            elif val < 0:
                continue

            else:
                row_negative = False

        all_negative.append(row_negative)

    # print(max_real_eigs)
    all_negative

    return all_negative

##
# Adds a column noting whether all the real part eigenvalues are positive
#
##


def all_positive_eigs(df):
    # Make subset containing only eigenvalues
    col_names = df.columns
    eign_cols = [x for x in col_names if 'eig' in x]
    eign_cols = [x for x in eign_cols if 'real' in x]
    eig_df = df[eign_cols]
    all_positive = []

    for row in eig_df.values:
        row_positive = True
        for val in row:
            if pd.isnull(val):
                continue

            elif val > 0:
                continue

            else:
                row_positive = False

        all_positive.append(row_positive)

    # print(max_real_eigs)
    all_positive

    return all_positive


def all_real_eigs(df):
    # Make subset containing only eigenvalues
    col_names = df.columns
    eign_cols = [x for x in col_names if 'eig' in x]
    eign_cols = [x for x in eign_cols if 'imag' in x]

    eig_df = df[eign_cols]
    all_real = []

    for row in eig_df.values:
        is_real = True

        for val in row:
            if val == 0:
                continue

            elif pd.isnull(val):
                continue

            else:
                is_real = False

        # print(is_real, ": ", row)
        all_real.append(is_real)

    # print(max_real_eigs)
    all_real

    return all_real


def all_zero_eigs(df):
    # Make subset containing only eigenvalues
    col_names = df.columns
    eign_cols = [x for x in col_names if 'eig' in x]
    real_part_cols = [x for x in eign_cols if 'real' in x]

    eig_df = df[real_part_cols]
    all_zero = []

    for row in eig_df.values:
        is_zero = True

        for val in row:
            if val == 0:
                continue

            elif pd.isnull(val):
                continue

            else:
                is_zero = False

        # print(is_real, ": ", row)
        all_zero.append(is_zero)

    # print(max_real_eigs)
    all_zero

    return all_zero


##
#   Finds conjugate pairs of eigenvalues from a list format of [ [real, imaginary], [real, imaginary] ]
#
#   Returns indexes of eigenvalues which are pairs
##
def get_conjugate_pairs(df):
    eign_cols = [x for x in df.columns if 'eig' in x]
    real_part_cols = [x for x in eign_cols if 'real' in x]
    imag_part_cols = [x for x in eign_cols if 'imag' in x]

    df_num_conj_pairs = []
    for idx, row in df.iterrows():
        real_parts = row[real_part_cols]
        imag_parts = row[imag_part_cols].values

        # Find matching real parts
        real_part_pairs = []
        for idx, i in enumerate(real_parts):
            if np.isnan(i):
                continue

            matches = np.equal(real_parts, i)
            match_indexes = [i for i, x in enumerate(matches) if x]

            # Add pairs to array
            if len(match_indexes) == 2 and match_indexes not in real_part_pairs:
                real_part_pairs.append(match_indexes)

        conjugate_pairs = []
        for p in real_part_pairs:
            i_0 = imag_parts[p[0]]
            i_1 = imag_parts[p[1]]

            if i_0 == 0 and i_1 == 0:
                continue

            if i_0 == i_1 * -1:
                conjugate_pairs.append(p)

        df_num_conj_pairs.append(len(conjugate_pairs))

    return df_num_conj_pairs


def set_number_imaginary(df):
    eign_cols = [x for x in df.columns if 'eig' in x]
    real_part_cols = [x for x in eign_cols if 'real' in x]
    imag_part_cols = [x for x in eign_cols if 'imag' in x]

    df_num_imag = []
    for idx, row in df.iterrows():
        imag_parts = row[imag_part_cols].values
        imag_count = 0
        for idx, i in enumerate(imag_parts):
            if np.isnan(i):
                continue

            if i == 0.0:
                continue

            else:
                imag_count += 1

        df_num_imag.append(imag_count)

    return df_num_imag


def distances_pre_processing(df):
    col_names = df.columns
    dist_cols = ['d1', 'd2', 'd3', 'd4', 'd5', 'd6']

    for c in dist_cols:
        df[c] = [np.float128(i) for i in df[c].values]

    return df


def set_accepted_column(df):
    col_names = df.columns
    dist_cols = ['d1', 'd2', 'd3', 'd4', 'd5', 'd6']

    df.loc[df['d1']]


def main():
    output_dir = "./output/two_species_big_3/Population_0"
    distances_path = output_dir + "/distances.csv"
    eigenvalues_path = output_dir + "/eigenvalues.csv"
    model_space_report_path = output_dir + "/model_space_report.csv"

    distaces_df = pd.read_csv(distances_path)
    distaces_df = distances_pre_processing(distaces_df)

    eigenvalues_df = pd.read_csv(eigenvalues_path)

    # Ignore failed simulations
    eigenvalues_df = eigenvalues_df.loc[eigenvalues_df['integ_error'].isnull()]

    joint_df = pd.merge(left=eigenvalues_df, right=distaces_df, how='inner', on=['sim_idx', 'batch_num', 'model_ref'])
    joint_df.reset_index()

    joint_df = species_sustained(joint_df)
    joint_df = make_max_eig(joint_df)
    joint_df = all_negative_eigs(joint_df)
    joint_df = all_real_eigs(joint_df)

    print(len(joint_df))
    print(joint_df.columns)

    loc_all_negative
    # loc_all_negative = joint_df.loc[joint_df['all_negative_eigs'] == True]
    print(len(loc_all_negative))
    loc_all_negative['sum_std'] = loc_all_negative['d2'] + loc_all_negative['d3']

    # loc_all_negative = loc_all_negative.loc[loc_all_negative['N_1_sustained'] == True]
    # loc_all_negative = loc_all_negative.loc[loc_all_negative['N_2_sustained'] == True]

    all_real = loc_all_negative.loc[loc_all_negative['all_real_eigs'] == True]
    not_real = loc_all_negative.loc[loc_all_negative['all_real_eigs'] == False]
    print(len(all_real))
    exit()
    ax = sns.scatterplot(x='sum_std', y='max_eig', data=all_real)
    ax = sns.scatterplot(x='sum_std', y='max_eig', data=not_real)

    ax.set(xscale="symlog", yscale='symlog')

    plt.plot()
    plt.show()


if __name__ == "__main__":
    adj_mat_dir = "/home/behzad/Documents/barnes_lab/sympy_consortium_framework/output/two_species_no_symm/adj_matricies/model_0_adj_mat.csv"
    adj_mat_dir = "/home/behzad/Documents/barnes_lab/cplusplus_software/speed_test/repressilator/cpp/input_files/input_files_two_species_0/adj_matricies/"
    adj_mat_dir = "/home/behzad/Documents/barnes_lab/cplusplus_software/speed_test/repressilator/cpp/input_files/input_files_three_species_0/adj_matricies/"
    model_idxs = range(510, 512)
    target_species_list = ['B_1', 'B_2']
    # get_adj_mat_distances(42, model_idxs, adj_mat_dir)
    get_n_strain_symmetric_adj_mats(model_idxs, adj_mat_dir)
    # get_symmetric_adj_mats(model_idxs, adj_mat_dir)
    # get_num_interactions(model_idxs, adj_mat_dir)
    # get_self_regulators(model_idxs, target_species_list, adj_mat_dir)