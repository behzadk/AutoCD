import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def normalise_list(input_list):
    max_list = np.max(input_list)
    min_list = np.min(input_list)
    output_list = [(x - min_list)/(max_list - min_list) for x in input_list]

    return output_list


def motif_heatmap(output_dir):
    figure_output_path = output_dir + 'motif_heatmap.pdf'

    model_space_report_path = output_dir + "combined_model_space_report_with_motifs.csv"
    model_space_report_df = pd.read_csv(model_space_report_path)

    # Columns 
    motif_columns = ['SL1', 'SL2', 'SL3', 'SL4', 
    'OL1', 'OL2', 'OL3', 'OL4']

    data = model_space_report_df[motif_columns].values
    print(data)



    width_inches = 175 / 25.4
    height_inches = 150 / 25.4
    # fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    fig, ax = plt.subplots()

    # sns.scatterplot(x='X', y='Y', hue='model_marginal', cmap=cmap, data=data_dict, ax=ax, size=0.5)
    # print(np.argmax(W[:, 0]))
    sns.heatmap(data.T, ax=ax, cbar=False)
    # sns.scatterplot(y='model_marginal', x=df.index, data=df, ax=ax2, s=1, alpha=0.5, color='white', orient='V')

    ax.set_xlabel('')
    # ax.set_xlabel('mean marginal probability change')
    ax.set_xticklabels('')
    ax.set_yticklabels('')

    # ax.set_yticklabels([])
    ax.set_ylabel('')
    ax.margins(x=0)
    ax.margins(y=0)
    # ax.set_aspect("equal")

    ax.legend().remove()

    fig.tight_layout()

    plt.savefig(figure_output_path, dpi=500, bbox_inches='tight')

    # plt.show()
    plt.close()


#### Self limiting ####
def count_direct_self_limiting(model_space_report_df, adj_mat_dir, window, normalise=False):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    model_idxs = model_space_report_df['model_idx'].values
    direct_self_limiting_counts = []

    for m_idx in model_idxs:
        model_self_limiting_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)
        # print(adj_mat_df)
        col_names = adj_mat_df.columns


        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat = adj_mat_df.values

        # adj_mat[to][from]
        for strain in strain_indexes:

            # Get microcins expressed by the strain
            strain_microcins = [i for i in microcin_indexes if adj_mat[i][strain] == 1]

            # Filter microcins for only those the strain is sensitive to
            strain_microcins = [i for i in strain_microcins if adj_mat[strain][i] == -1]

            direct_self_limiting_microcins = []

            # Filter microcins that are constitutive
            for m in strain_microcins:
                is_direct_self_limiting = True

                for a in AHL_indexes:
                    if adj_mat[m][a] != 0:
                        is_direct_self_limiting = False

                if is_direct_self_limiting:
                    direct_self_limiting_microcins.append(m)

            if len(direct_self_limiting_microcins) > 0:
                model_self_limiting_count += 1
        
        direct_self_limiting_counts.append(model_self_limiting_count)

    if normalise:
        direct_self_limiting_counts = normalise_list(direct_self_limiting_counts)


    print("finished making list")
    model_space_report_df['direct_self_limiting_counts'] = direct_self_limiting_counts

    if window > 0:
        model_space_report_df['direct_self_limiting_rolling'] = model_space_report_df['direct_self_limiting_counts'].rolling(window=window).mean()

    return model_space_report_df

##
# Dependent relationships arise from a self-killing strain whose rate of self-killing is regulated by
# a quorum molecule that is expressed by another strain
##
def count_dependent(model_space_report_df, adj_mat_dir, window, normalise=False):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    model_idxs = model_space_report_df['model_idx'].values

    dependent_counts = []

    for m_idx in model_idxs:
        model_dependent_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)
        # print(adj_mat_df)
        col_names = adj_mat_df.columns

        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat = adj_mat_df.values

        # adj_mat[to][from]
        for strain in strain_indexes:

            # Get microcins expressed by the strain
            strain_microcins = [i for i in microcin_indexes if adj_mat[i][strain] == 1]

            # Filter microcins for only those the strain is sensitive to
            strain_microcins = [i for i in strain_microcins if adj_mat[strain][i] == -1]

            for other_strain in strain_indexes:
                if other_strain == strain:
                    continue


                other_strain_AHLs = [i for i in AHL_indexes if adj_mat[i][other_strain] == 1]
                
                is_repressed = False
                
                # Filter microcins for those induced by the AHL
                for m in strain_microcins:
                    for a in other_strain_AHLs:
                        if adj_mat[m][a] == -1:
                            is_repressed = True

                if is_repressed:
                    model_dependent_count += 1

        dependent_counts.append(model_dependent_count)


    if normalise:
        dependent_counts = normalise_list(dependent_counts)

    model_space_report_df['SL2'] = dependent_counts

    if window > 0:
        model_space_report_df['dependent_rolling'] = model_space_report_df['SL2'].rolling(window=window).mean()

    # sns.scatterplot(x=model_space_report_df.index, y='predator_prey_counts', data=model_space_report_df)
    # plt.show()

    return model_space_report_df


##
# Strain which would propel their current state. Represses self-killing microcin with self-expressed AHL
##
def count_hedonistic(model_space_report_df, adj_mat_dir, window, normalise=False):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    model_idxs = model_space_report_df['model_idx'].values

    hedonistic_counts = []

    for m_idx in model_idxs:
        model_hedonistic_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)

        # print(adj_mat_df)
        col_names = adj_mat_df.columns


        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat = adj_mat_df.values

        # adj_mat[to][from]
        for strain in strain_indexes:

            # Get microcins expressed by the strain
            strain_microcins = [i for i in microcin_indexes if adj_mat[i][strain] == 1]

            # Filter microcins for only those the strain is sensitive to
            strain_microcins = [i for i in strain_microcins if adj_mat[strain][i] == -1]

            # Get AHLs expressed by the strain
            strain_AHLs = [i for i in AHL_indexes if adj_mat[i][strain] == 1]


            # Filter microcins for those repressed by the AHL
            repressed_microcins = []
            for m in strain_microcins:
                is_repressed = False
                for a in strain_AHLs:
                    if adj_mat[m][a] == -1:
                        is_repressed = True
                        break

                if is_repressed:
                    repressed_microcins.append(m)
                    break

            strain_microcins = repressed_microcins

            if len(strain_microcins) > 0:
                model_hedonistic_count += 1

        hedonistic_counts.append(model_hedonistic_count)

    if normalise:
        hedonistic_counts = normalise_list(hedonistic_counts)

    print("finished making list")
    model_space_report_df['SL4'] = hedonistic_counts

    if window > 0:
        model_space_report_df['hedonistic_rolling'] = model_space_report_df['SL4'].rolling(window=window).mean()

    # sns.lineplot(x=model_space_report_df.index, y='hedonistic_rolling', data=model_space_report_df)
    # # sns.tsplot(data=model_space_report_df.hedonistic_count, estimator=rolling_mean)
    # # sns.scatterplot(x=model_space_report_df.index, y='hedonistic_count', data=model_space_report_df)
    # plt.show()

    return model_space_report_df


##
# Permissive motif described asa  strain that expresses a microcin that it is sensitive to.
# The microcin is induced by the QS molecule the strain also expresses
##
def count_permissive(model_space_report_df, adj_mat_dir, window, normalise=False):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    model_idxs = model_space_report_df['model_idx'].values

    permissive_counts = []

    for m_idx in model_idxs:
        model_permissive_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)
        # print(adj_mat_df)
        col_names = adj_mat_df.columns


        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat = adj_mat_df.values

        # adj_mat[to][from]
        for strain in strain_indexes:

            # Get microcins expressed by the strain
            strain_microcins = [i for i in microcin_indexes if adj_mat[i][strain] == 1]

            # Filter microcins for only those the strain is sensitive to
            strain_microcins = [i for i in strain_microcins if adj_mat[strain][i] == -1]

            # Get AHLs expressed by the strain
            strain_AHLs = [i for i in AHL_indexes if adj_mat[i][strain] == 1]


            # Filter microcins for those induced by the AHL
            induced_microcins = []
            for m in strain_microcins:
                is_induced = False
                for a in strain_AHLs:
                    if adj_mat[m][a] == 1:
                        is_induced = True

                if is_induced:
                    induced_microcins.append(m)

            strain_microcins = induced_microcins

            if len(strain_microcins) > 0:
                model_permissive_count += 1

        permissive_counts.append(model_permissive_count)

    if normalise:
        permissive_counts = normalise_list(permissive_counts)

    model_space_report_df['SL1'] = permissive_counts

    if window > 0:
        model_space_report_df['permissive_rolling'] = model_space_report_df['SL1'].rolling(window=window).mean()

    return model_space_report_df


def count_submissive(model_space_report_df, adj_mat_dir, window, normalise=False):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    model_idxs = model_space_report_df['model_idx'].values

    submissive_counts = []

    for m_idx in model_idxs:
        model_submissive_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)
        # print(adj_mat_df)
        col_names = adj_mat_df.columns

        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat = adj_mat_df.values

        # adj_mat[to][from]
        for strain in strain_indexes:

            # Get microcins expressed by the strain
            strain_microcins = [i for i in microcin_indexes if adj_mat[i][strain] == 1]

            # Filter microcins for only those the strain is sensitive to
            strain_microcins = [i for i in strain_microcins if adj_mat[strain][i] == -1]

            for other_strain in strain_indexes:
                if other_strain == strain:
                    continue


                other_strain_AHLs = [i for i in AHL_indexes if adj_mat[i][other_strain] == 1]
                
                is_induced = False
                
                # Filter microcins for those induced by the AHL
                for m in strain_microcins:
                    for a in other_strain_AHLs:
                        if adj_mat[m][a] == 1:
                            is_induced = True

                if is_induced:
                    model_submissive_count += 1

        submissive_counts.append(model_submissive_count)


    if normalise:
        submissive_counts = normalise_list(submissive_counts)

    model_space_report_df['SL3'] = submissive_counts

    if window > 0:
        model_space_report_df['dependent_rolling'] = model_space_report_df['SL3'].rolling(window=window).mean()

    # sns.scatterplot(x=model_space_report_df.index, y='predator_prey_counts', data=model_space_report_df)
    # plt.show()

    return model_space_report_df


#### Competitive ####
##
# Competitive counted as when a strain is sensitive to a microcin it does not itself produce
##
def count_direct_competitive(model_space_report_df, adj_mat_dir, window, normalise=False):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    model_idxs = model_space_report_df['model_idx'].values

    direct_competitive_count = []

    for m_idx in model_idxs:
        model_direct_competitive_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)

        # print(adj_mat_df)
        col_names = adj_mat_df.columns


        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat = adj_mat_df.values

        # adj_mat[to][from]
        for strain in strain_indexes:
            # Get microcins expressed by the strain
            strain_microcins = [i for i in microcin_indexes if adj_mat[i][strain] == 1]


            for other_strain in strain_indexes:
                is_competitive = False

                if other_strain == strain:
                    continue

                for m in strain_microcins:
                    if adj_mat[other_strain][m] == -1:
                        is_competitive = True
                

                for a in AHL_indexes:
                    if adj_mat[m][a] != 0:
                        is_competitive = False

                if is_competitive:
                    model_direct_competitive_count += 1

        direct_competitive_count.append(model_direct_competitive_count)

    if normalise:
        direct_competitive_count = normalise_list(direct_competitive_count)


    print("finished making list")
    model_space_report_df['direct_competitive_counts'] = direct_competitive_count

    if window > 0:
        model_space_report_df['direct_competitive_rolling'] = model_space_report_df['direct_competitive_counts'].rolling(window=window).mean()

    return model_space_report_df


def count_defensive(model_space_report_df, adj_mat_dir, window, normalise=False):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    model_idxs = model_space_report_df['model_idx'].values

    defensive_counts = []

    for m_idx in model_idxs:
        model_defensive_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)

        # print(adj_mat_df)
        col_names = adj_mat_df.columns


        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat = adj_mat_df.values

        # adj_mat[to][from]
        for strain in strain_indexes:

            # Get microcins expressed by the strain
            strain_microcins = [i for i in microcin_indexes if adj_mat[i][strain] == 1]

            for other_strain in strain_indexes:
                is_defensive = False
                if other_strain == strain:
                    continue

                other_strain_AHLs = [i for i in AHL_indexes if adj_mat[i][other_strain] == 1]
                
                # Filter microcins for those induced by the AHL
                induced_microcins = []
                for m in strain_microcins:
                    # If strain is sensitive 
                    if adj_mat[other_strain][m] == -1:

                        for a in other_strain_AHLs:

                            # If other strain is inducing microcin expression
                            if adj_mat[m][a] == 1:
                                is_defensive = True

                if is_defensive:
                    model_defensive_count += 1

        defensive_counts.append(model_defensive_count)

    if normalise:
        defensive_counts = normalise_list(defensive_counts)


    print("finished making list")
    model_space_report_df['OL1'] = defensive_counts

    if window > 0:
        model_space_report_df['defensive_rolling'] = model_space_report_df['OL1'].rolling(window=window).mean()

    return model_space_report_df

def count_exponential(model_space_report_df, adj_mat_dir, window, normalise=False):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    model_idxs = model_space_report_df['model_idx'].values

    exponential_counts = []

    for m_idx in model_idxs:
        model_exponential_counts = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)

        # print(adj_mat_df)
        col_names = adj_mat_df.columns


        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat = adj_mat_df.values

        # adj_mat[to][from]
        for strain in strain_indexes:

            # Get microcins expressed by the strain
            strain_microcins = [i for i in microcin_indexes if adj_mat[i][strain] == 1]

            # Get AHLS expressed by the strain
            strain_AHLs = [a for a in AHL_indexes if adj_mat[a][strain] == 1]


            # Get microcins induced by strain AHLs
            induced_microcins = []
            for m in strain_microcins:
                is_induced = False
                for a in strain_AHLs:
                    if adj_mat[m][a] == 1:
                        is_induced = True

                if is_induced:
                    induced_microcins.append(m)

            for other_strain in strain_indexes:
                is_sensitive = False
                if other_strain == strain:
                    continue

                for m in induced_microcins:
                    if adj_mat[other_strain][m] == -1:
                        is_sensitive = True

                if is_sensitive:
                    model_exponential_counts += 1


        exponential_counts.append(model_exponential_counts)
    

    print("finished making list")

    if normalise:
        exponential_counts = normalise_list(exponential_counts)

    model_space_report_df['OL4'] = exponential_counts

    if window > 0:
        model_space_report_df['exponential_rolling'] = model_space_report_df['OL4'].rolling(window=window).mean()

    return model_space_report_df

def count_logistic(model_space_report_df, adj_mat_dir, window, normalise=False):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    model_idxs = model_space_report_df['model_idx'].values

    logistic_counts = []

    for m_idx in model_idxs:
        model_logistic_counts = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)

        # print(adj_mat_df)
        col_names = adj_mat_df.columns


        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat = adj_mat_df.values

        # adj_mat[to][from]
        for strain in strain_indexes:

            # Get microcins expressed by the strain
            strain_microcins = [i for i in microcin_indexes if adj_mat[i][strain] == 1]

            # Get AHLS expressed by the strain
            strain_AHLs = [a for a in AHL_indexes if adj_mat[a][strain] == 1]


            # Get microcins induced by strain AHLs
            repressed_microcins = []
            for m in strain_microcins:
                is_repressed = False
                for a in strain_AHLs:
                    if adj_mat[m][a] == -1:
                        is_repressed = True

                if is_repressed:
                    repressed_microcins.append(m)

            for other_strain in strain_indexes:
                is_sensitive = False
                if other_strain == strain:
                    continue

                for m in repressed_microcins:
                    if adj_mat[other_strain][m] == -1:
                        is_sensitive = True

                if is_sensitive:
                    model_logistic_counts += 1


        logistic_counts.append(model_logistic_counts)
    

    if normalise:
        logistic_counts = normalise_list(logistic_counts)


    print("finished making list")
    model_space_report_df['OL2'] = logistic_counts

    if window > 0:
        model_space_report_df['logistic_rolling'] = model_space_report_df['OL2'].rolling(window=window).mean()

    return model_space_report_df

def count_opportunistic(model_space_report_df, adj_mat_dir, window, normalise=False):
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"
    model_idxs = model_space_report_df['model_idx'].values

    achillies_count = []

    for m_idx in model_idxs:
        model_achillies_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)

        # print(adj_mat_df)
        col_names = adj_mat_df.columns


        strain_indexes = [idx for idx, i in enumerate(col_names) if 'N_' in i]
        AHL_indexes = [idx for idx, i in enumerate(col_names) if 'A_' in i]
        microcin_indexes = [idx for idx, i in enumerate(col_names) if 'B_' in i]

        adj_mat = adj_mat_df.values

        # adj_mat[to][from]
        for strain in strain_indexes:

            # Get microcins expressed by the strain
            strain_microcins = [i for i in microcin_indexes if adj_mat[i][strain] == 1]

            for other_strain in strain_indexes:
                is_achillies = False
                if other_strain == strain:
                    continue

                other_strain_AHLs = [i for i in AHL_indexes if adj_mat[i][other_strain] == 1]
                
                # Filter microcins for those induced by the AHL
                for m in strain_microcins:
                    # If strain is sensitive 
                    if adj_mat[other_strain][m] == -1:

                        for a in other_strain_AHLs:
                            # If other strain is inducing microcin expression
                            if adj_mat[m][a] == -1:
                                is_achillies = True

                if is_achillies:
                    model_achillies_count += 1

        achillies_count.append(model_achillies_count)

    if normalise:
        achillies_count = normalise_list(achillies_count)


    print("finished making list")
    model_space_report_df['OL3'] = achillies_count

    if window > 0:
        model_space_report_df['opportunistic_rolling'] = model_space_report_df['OL3'].rolling(window=window).mean()

    return model_space_report_df

