import numpy as np
import pandas as pd


##
#   Returns bool describing whether matrix is symmetric or not.
#   matrix must already be ordered e.g N_1, B_1, A_1, A_2, B_2, N_2
#
##
def check_matrix_symmetrical(adj_mat):
    n_species = np.shape(adj_mat)[0]

    mid_index = (n_species - 1) / 2
    sect_one_start = 0
    sect_one_end = int(np.floor(mid_index) + 1)

    sect_two_start = int(np.ceil(mid_index))
    sect_two_end = np.shape(adj_mat)[0]


    # Split across midpoint at both axes
    sect_one_adj_mat = adj_mat[:, sect_one_start:sect_one_end]
    sect_two_adj_mat =  adj_mat[:, sect_two_start:sect_two_end]

    split_adj_mat = np.zeros([2, 2], dtype=object)
    split_adj_mat[0, 0] = sect_one_adj_mat[sect_one_start:sect_one_end]
    split_adj_mat[1, 0] = sect_one_adj_mat[sect_two_start:sect_two_end]
    split_adj_mat[0, 1] = sect_two_adj_mat[sect_one_start:sect_one_end]
    split_adj_mat[1, 1] = sect_two_adj_mat[sect_two_start:sect_two_end]

    symmetrical = False
    x = np.equal(np.flip(split_adj_mat[0, 0]), split_adj_mat[1, 1])

    if np.all(np.equal(np.flip(split_adj_mat[0, 0]), split_adj_mat[1, 1])) and np.all(np.equal(np.flip(split_adj_mat[1, 0]), split_adj_mat[0, 1])):
        symmetrical = True

    return symmetrical




##
#   Rearranges matrix for symmetry analysis
#   e.g input N_1, N_2, B_1, B_2, A_1, A_2
#       output N_1, B_1, A_1, A_2, B_2, N_2
#
##
def rearrange_matrix_two_strain(adj_mat_df):
    strain_col = [col for col in adj_mat_df if col.startswith('N_')]
    AHL_col = [col for col in adj_mat_df if col.startswith('A_')]
    microcin_col = [col for col in adj_mat_df if col.startswith('B_')]
    new_col_order = [0] * (len(adj_mat_df.columns) )

    # Swap columns
    rev_idx = -1
    fwd_idx = 0
    for x, _ in enumerate(strain_col):
        if x % 2 == 0:
            new_col_order[fwd_idx] = strain_col[x]
            fwd_idx += 1

        else:
            new_col_order[rev_idx] = strain_col[x]
            rev_idx += -1

    for x, _ in enumerate(microcin_col):
        if x % 2 == 0:
            new_col_order[fwd_idx] = microcin_col[x]
            fwd_idx += 1

        else:
            new_col_order[rev_idx] = microcin_col[x]
            rev_idx += -1

    for x, _ in enumerate(AHL_col):
        if x % 2 == 0:
            new_col_order[fwd_idx] = AHL_col[x]
            fwd_idx += 1

        else:
            new_col_order[rev_idx] = AHL_col[x]
            rev_idx += -1

    adj_mat = adj_mat_df.as_matrix()
    new_adj_mat = np.zeros([len(new_col_order), len(new_col_order)])

    for new_idx, col in enumerate(new_col_order):

        old_idx = list(adj_mat_df.columns).index(col)
        # Get column and rearrange
        col_vals = np.copy(adj_mat)[:, old_idx]
        sorted_col_vals = np.copy(col_vals)
        # Rearrange column
        for new_row_idx, row in enumerate(new_col_order):
            old_row_idx = list(adj_mat_df.columns).index(row)
            sorted_col_vals[new_row_idx] = np.copy(adj_mat)[:, old_idx][old_row_idx]

        new_adj_mat[:, new_idx] = sorted_col_vals

    return new_adj_mat


def rearrange_matrix_three_strain(adj_mat_df):
    strain_col = [col for col in adj_mat_df if col.startswith('N_')]
    AHL_col = [col for col in adj_mat_df if col.startswith('A_')]
    microcin_col = [col for col in adj_mat_df if col.startswith('B_')]
    new_col_order = [0] * (len(adj_mat_df.columns) )

    # Swap columns
    rev_idx = -1
    fwd_idx = 0
    for x, _ in enumerate(strain_col):
        if x % 2 == 0:
            new_col_order[fwd_idx] = strain_col[x]
            fwd_idx += 1

        else:
            new_col_order[rev_idx] = strain_col[x]
            rev_idx += -1

    for x, _ in enumerate(AHL_col):
        if x % 2 == 0:
            new_col_order[fwd_idx] = AHL_col[x]
            fwd_idx += 1

        else:
            new_col_order[rev_idx] = AHL_col[x]
            rev_idx += -1

    for x, _ in enumerate(microcin_col):
        if x % 2 == 0:
            new_col_order[fwd_idx] = microcin_col[x]
            fwd_idx += 1

        else:
            new_col_order[rev_idx] = microcin_col[x]
            rev_idx += -1


    print(new_col_order)
    adj_mat = adj_mat_df.as_matrix()
    new_adj_mat = np.zeros([len(new_col_order), len(new_col_order)])

    for new_idx, col in enumerate(new_col_order):

        old_idx = list(adj_mat_df.columns).index(col)
        # Get column and rearrange
        col_vals = np.copy(adj_mat)[:, old_idx]
        sorted_col_vals = np.copy(col_vals)
        # Rearrange column
        for new_row_idx, row in enumerate(new_col_order):
            old_row_idx = list(adj_mat_df.columns).index(row)
            sorted_col_vals[new_row_idx] = np.copy(adj_mat)[:, old_idx][old_row_idx]

        new_adj_mat[:, new_idx] = sorted_col_vals

    return new_adj_mat
