import numpy as np
from itertools import chain
from decimal import *


##
#   Finds conjugate pairs of eigenvalues from a list format of [ [real, imaginary], [real, imaginary] ]
#
#   Returns indexes of eigenvalues which are pairs
##
def get_conjugate_pairs(eigen_values):
    real_parts = [i[0] for i in eigen_values]
    imag_parts = [i[1] for i in eigen_values]

    # Find matching real parts
    real_part_pairs = []
    for idx, i in enumerate(real_parts):
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

    return conjugate_pairs


##
#   Calculates the eigenvalue product. If conjugate pairs exist, it will compute the product of these first, followed
#   by multiplication of the remaining real only eigenvalues.
#
#   Returns product of eigenvalues
##
def eigenvalue_product(eigenvalues):
    real_parts = [i[0] for i in eigenvalues]
    imag_parts = [i[1] for i in eigenvalues]

    # Check all imaginary parts are zero
    if all(i == 0.0 for i in imag_parts):
        output = np.float64(1)
        for i in real_parts:
            output = output * i

        return output

    # Get indexes of conjugate pairs and indexes of eigenvalues that are not paired
    conjugate_pairs = get_conjugate_pairs(eigenvalues)
    not_paired = [i for i in range(len(real_parts)) if i not in chain.from_iterable(conjugate_pairs)]

    print("Paired: ", conjugate_pairs)
    print("Not paired: ", not_paired)

    # Calculate products of conjugate pairs
    conjugate_pair_products = []
    for p in conjugate_pairs:
        prod = np.float64(real_parts[p[0]] * real_parts[p[1]] + abs(imag_parts[p[0]] * imag_parts[p[1]]))
        conjugate_pair_products.append(prod)

    output = np.float64(np.product(conjugate_pair_products))

    # Check if non paired parts have no imaginary part
    for p in not_paired:
        if imag_parts[p] != 0:
            print("Non eigenvalue imaginary part does not equal 0")

        else:
            output = output * real_parts[p]

    return output


##
#   Calculates the sum of real part eigenvalues
#
##
def sum_eigenvalues(eigen_values):
    real_parts = [i[0] for i in eigen_values]
    return sum(real_parts)


def delete_matrix_row_column(matrix, del_column):
    # Delete first row
    matrix = np.delete(matrix, np.s_[0], axis=0)

    # Delete specified column
    matrix = np.delete(matrix, np.s_[del_column], axis=1)

    matrix = np.reshape(matrix, (len(matrix), len(matrix)))
    return matrix


def get_determinant(matrix):
    if (len(matrix) == 2):
        return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0])

    matrix_size = len(matrix)

    answer = 0
    for j in range(matrix_size):
        answer += (-1) ** j * matrix[0][j] * get_determinant(delete_matrix_row_column(matrix, j))

    return answer


##
#   Takes a set of eigenvalues to classify a system.
#   Eigen values come in pairs of real and imaginary parts
##
def classify_eigensystem(eigenvalues):
    real_parts = [i[0] for i in eigenvalues]
    imag_parts = [i[1] for i in eigenvalues]

    only_real_parts = all(i == 0.0 for i in imag_parts)
    conjugate_paairs = get_conjugate_pairs(eigenvalues)

    all_neg_real_parts = all(i < 0 for i in real_parts)

    centers = []

    print("")
    print("Only real parts: \t", only_real_parts)
    print("All -ve real parts: \t", all_neg_real_parts)
    print("Number of conjugate pairs: \t", len(conjugate_paairs))
    print("")

    for eig in eigenvalues:
        print(eig)
