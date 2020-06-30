import pandas as pd
import numpy as np
import os
import glob


def combine_model_space_reports(exp_1_dir, exp_2_dir, output_dir):
    print("Starting combine model space report")

    exp_1_df = pd.read_csv(exp_1_dir + 'model_space_report.csv')
    exp_2_df = pd.read_csv(exp_2_dir + 'model_space_report.csv')

    # Check model refs match
    exp_1_models = np.array(exp_1_df['model_idx'].values)
    exp_2_models = np.array(exp_2_df['model_idx'].values)

    if not np.array_equal(exp_1_models, exp_2_models):
        print("model references do not match")
        exit()

    row_list = []
    # Iterate both dataframes
    for (row_idx, s1), (_, s2) in zip(exp_1_df.iterrows(), exp_2_df.iterrows()):
        model_idx = s1['model_idx']
        s1_accepted = s1['accepted_count']
        s2_accepted = s2['accepted_count']
        s1_simulated = s1['simulated_count']
        s2_simulated = s2['simulated_count']

        data_dict = {'model_idx': model_idx,
                     'accepted_count': s1_accepted + s2_accepted,
                     'simulated_count': s1_simulated + s2_simulated}

        row_list.append(data_dict)

    new_df = pd.DataFrame(row_list, columns=['model_idx', 'accepted_count', 'simulated_count'])
    new_df.to_csv(output_dir + 'model_space_report.csv')
    print("Combine model space report finished")
    print("")


def combine_distances(exp_1_dir, exp_2_dir, output_dir):
    print("Combining distances")

    print("Reading in dataframes")
    exp_1_df = pd.read_csv(exp_1_dir + 'distances.csv')
    exp_2_df = pd.read_csv(exp_2_dir + 'distances.csv')

    exp_1_df = exp_1_df.loc[exp_1_df['Accepted'] == True]
    exp_2_df = exp_1_df.loc[exp_1_df['Accepted'] == True]

    # start_batches = exp_1_df.iloc[-1]['batch_num'] + 1

    # print("changing batch numbers for exp 2")
    # exp_2_df['batch_num'] = exp_2_df['batch_num'].apply(lambda x: x + start_batches)

    print("Concatenating dataframes")
    concat_df = pd.concat([exp_1_df, exp_2_df], ignore_index=True)

    print("Saving distances.csv")
    concat_df.to_csv(output_dir + 'distances.csv', index=False)
    print("Combine distances finished")
    print("")


def combine_distances_individually(exp_1_dir, output_dir):
    master_dir = output_dir + "model_sim_distances/"
    try:
        os.mkdir(master_dir)
    except FileExistsError:
        pass

    model_distance_template_path = master_dir + "model_#IDX#_distances.csv"
    # Read exp_1
    exp_1_df = pd.read_csv(exp_1_dir + 'distances.csv')
    model_idxs = exp_1_df['model_ref'].unique()

    for idx in model_idxs:
        model_data = exp_1_df.loc[exp_1_df['model_ref'] == idx]

        master_data_path = model_distance_template_path.replace('#IDX#', str(idx))

        # Check if csv file exists
        try:
            master_model_data = pd.read_csv(master_data_path)

        # Master data doesn't exist, so save current df
        except FileNotFoundError:
            model_data.to_csv(master_data_path, index=False)
            continue

        concat_df = pd.concat([master_model_data, model_data], ignore_index=True)
        concat_df.to_csv(master_data_path, index=False)


def combine_eigenvalues(exp_1_dir, exp_2_dir, output_dir):
    print("Combining eigenvalues")

    print("Reading in dataframes")
    exp_1_df = pd.read_csv(exp_1_dir + 'eigenvalues_do_fsolve_state.csv')
    exp_2_df = pd.read_csv(exp_2_dir + 'eigenvalues_do_fsolve_state.csv')

    start_batches = exp_1_df.iloc[-1]['batch_num'] + 1

    print("changing batch numbers for exp 2")
    exp_2_df['batch_num'] = exp_2_df['batch_num'].apply(lambda x: x + start_batches)

    print("Concatenating dataframes")
    concat_df = pd.concat([exp_1_df, exp_2_df], ignore_index=True)

    print("Saving eigenvalues_do_fsolve_state.csv")
    concat_df.to_csv(output_dir + 'eigenvalues_do_fsolve_state.csv')

    print("Combine eigenvalues finished")
    print("")


def combine_model_sim_params(exp_1_dir, exp_2_dir, output_dir):
    output_dir = output_dir + "model_sim_params/"
    try:
        os.mkdir(output_dir)
    except FileExistsError:
        pass

    exp_1_params_dir = exp_1_dir + "model_sim_params/"
    exp_2_params_dir = exp_2_dir + "model_sim_params/"

    exp_1_params_path = [file_path for file_path in glob.iglob(exp_1_params_dir + "*_all_params") if
                         "population" in file_path]
    exp_2_params_path = [file_path for file_path in glob.iglob(exp_2_params_dir + "*_all_params") if
                         "population" in file_path]

    if len(exp_1_params_path) != len(exp_2_params_path):
        print("model sims missing")
        exit()
    else:
        num_models = len(exp_1_params_path)

    # Get paths in order of model idx
    exp_1_ordered_paths = [0 for f in range(num_models)]
    exp_2_ordered_paths = [0 for f in range(num_models)]

    for idx, f1 in enumerate(exp_1_params_path):
        f1_base_name = os.path.basename(f1)
        f1_model_num = int(os.path.basename(f1).split('_')[1])
        match = False
        for f2 in exp_2_params_path:
            f2_base_name = os.path.basename(f2)
            f2_model_num = int(os.path.basename(f2).split('_')[1])

            if f1_model_num == f2_model_num:
                model_num = int(os.path.basename(f1).split('_')[1])

                exp_1_ordered_paths[idx] = f1
                exp_2_ordered_paths[idx] = f2
                idx += 1
                break

    count = 0
    for (exp_1_f, exp_2_f) in zip(exp_1_ordered_paths, exp_2_ordered_paths):
        if count % 100 == 0:
            print("combining model " + str(count))
        if os.path.basename(exp_1_f) != os.path.basename(exp_2_f):
            print("file mismatch, exiting")
            print(os.path.basename(exp_1_f))
            print(os.path.basename(exp_2_f))
            exit()

        exp_1_df = pd.read_csv(exp_1_f)
        exp_2_df = pd.read_csv(exp_2_f)

        outfile_name = os.path.basename(exp_1_f)

        if len(exp_2_df) == 0:
            exp_1_df.to_csv(output_dir + outfile_name, index=False)
            count += 1
            continue

        # start_batches = max(exp_2_df['batch_num'].values) + 1
        # exp_2_df['batch_num'] = exp_2_df['batch_num'].apply(lambda x: x + start_batches)
        concat_df = pd.concat([exp_1_df, exp_2_df], ignore_index=True)
        concat_df.to_csv(output_dir + outfile_name, index=False)
        count += 1


def main():
    data_dir = "./output/"
    exp_1_dir = data_dir + "three_species_stable_0_comb/Population_0/"
    exp_2_dir = data_dir + "three_species_stable_0d/Population_0/"

    new_output_dir = data_dir + "three_species_stable_0_comb/Population_0/"

    try:
        os.makedirs(new_output_dir)
    except FileExistsError:
        pass

    # New df has batch number added
    exp_1_df = pd.read_csv(exp_1_dir + 'distances.csv')
    start_batches = exp_1_df.iloc[-1]['batch_num'] + 1
    del exp_1_df

    combine_model_space_reports(exp_1_dir, exp_2_dir, new_output_dir)
    combine_model_sim_params(exp_1_dir, exp_2_dir, new_output_dir)
    combine_eigenvalues(exp_1_dir, exp_2_dir, new_output_dir)
    combine_distances(exp_1_dir, exp_2_dir, new_output_dir)


if __name__ == "__main__":
    main()
