import os
import numpy as np


def make_folder(folder_dir):
    try:
        os.mkdir(folder_dir)

    except FileExistsError:
        pass


def check_model_list_repeats(model_list, array_a, array_b):
    # Get number of array rows

    keep_idx = []

    for idx_a, model_a in enumerate(model_list):
        model_unique = True

        array_a = model_a.adjacency_matrix
        rows = np.shape(array_a)[0]

        for idx_b, model_b in enumerate(model_list):
            array_b = model_b.adjacency_matrix

            if idx_a == idx_b:
                continue

            model_match = True
            for row_idx in range(rows):
                if list(array_a[row_idx]) != list(array_b[row_idx]):
                    model_match = False
                    break

            if model_match == True:
                model_unique = False
                break

        if model_unique:
            keep_idx.append(idx_a)
