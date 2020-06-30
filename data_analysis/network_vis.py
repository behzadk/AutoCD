import numpy as np
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import netwulf
import seaborn as sns

import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import scipy
from .NMF_analysis import convert_QS_column

import os
import glob
from pathlib import Path
# font = {'size'   : 15, }
# axes = {'labelsize': 'medium', 'titlesize': 'medium'}
# import matplotlib as mpl
# sns.set_context("talk")
# sns.set_style("white")
# mpl.rc('font', **font)
# mpl.rc('axes', **axes)
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
# plt.rcParams['text.usetex'] = True
plt.rcParams['axes.unicode_minus'] = False
sns.set_style("white")



min_max_scaling = lambda a, b, min_x, max_x, x: a + ((x - min_x) * (b - a)) / (max_x - min_x)

def two_species_vis_parameters():
    preset = {'zoom': 3.3366265659632246, 'node_charge': -45, 'node_gravity': 0.4087170329446852, 'link_distance': 15, 'link_distance_variation': 0, 
    'node_collision': True, 'wiggle_nodes': False, 'freeze_nodes': False, 'node_fill_color': '#79aaa0', 
    'node_stroke_color': '#000000', 'node_label_color': '#000000', 'display_node_labels': False, 'scale_node_size_by_strength': False, 
    'node_size': 5, 'node_stroke_width': 0, 'node_size_variation': 1.0, 'link_color': '#000000', 'link_width': 0.5, 
    'link_alpha': 0.5, 'link_width_variation': 0.5, 'display_singleton_nodes': True, 'min_link_weight_percentile': 0, 'max_link_weight_percentile': 1}

    return preset


def find_nearest_neighbours(current_combo, all_combos):
    neighbours = []

    min_nighbour_distance = 1e10

    for idx, potential_neighbour in enumerate(all_combos):
        if np.array_equal(potential_neighbour, current_combo):
            continue

        diff_combos = np.sum(np.absolute(np.subtract(potential_neighbour, current_combo)))


        if abs(diff_combos) == min_nighbour_distance:
            neighbours.append(idx)


        elif abs(diff_combos) < min_nighbour_distance:
            min_nighbour_distance = abs(diff_combos)
            neighbours = [idx]

        else:
            continue


    return neighbours


def find_adjacent_neighbours(current_combo, all_combos):
    neighbours = []

    for idx, potential_neighbour in enumerate(all_combos):
        if np.array_equal(potential_neighbour, current_combo):
            continue

        diff_combos = np.sum(np.absolute(np.subtract(potential_neighbour, current_combo)))


        if abs(diff_combos) == 1:
            neighbours.append(idx)

        else:
            continue


    return neighbours
    

def model_type_conditional(row, col_SL_only='blue', col_OL_only='red', col_mixed='grey'):
    SL_columns = ['SL1', 'SL2', 'SL3', 'SL4']
    OL_columns = ['OL1', 'OL2', 'OL3', 'OL4']

    has_SL_motifs = False
    has_OL_motifs = False

    if np.sum(row[SL_columns]) > 0:
        has_SL_motifs = True
    
    if np.sum(row[OL_columns]) > 0:
        has_OL_motifs = True


    if has_SL_motifs and has_OL_motifs:
        return col_mixed

    if has_SL_motifs:
        return col_SL_only

    if has_OL_motifs:
        return col_OL_only


def make_model_motif_colours(model_space_df, col_SL_only='blue', col_OL_only='red', col_mixed='grey'):
    model_space_df['node_type'] = model_space_df.apply(model_type_conditional, args=(col_SL_only, col_OL_only, col_mixed), axis=1)
    return model_space_df


def make_model_network(combined_analysis_output_dir, adj_mat_dir, drop_eqless=0.00, colour_by_motif=False, use_nearest_neighbour=False, use_adjacent_neighbour=False):
    if use_nearest_neighbour:
        output_name = "nearest_neighbour_vis.pdf"
    
    elif use_adjacent_neighbour:
        output_name = "adjacent_neighbour_vis.pdf"

    if use_nearest_neighbour and use_adjacent_neighbour:
        print("use only one of nearest or adjacent neighbour")
        exit()

    model_space_report_path = combined_analysis_output_dir + "combined_model_space_report_with_motifs.csv"
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"

    model_space_report_df = pd.read_csv(model_space_report_path)
    model_space_report_df = model_space_report_df.sort_values('model_idx')

    model_idxs = model_space_report_df['model_idx'].values
    model_marginals = model_space_report_df['model_marginal_mean'].values

    flat_adj_mats = []

    # Load and flatten all adj mats
    for m_idx in model_idxs:
        model_self_limiting_count = 0

        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)
        flat_adj_mats.append(adj_mat_df.values.flatten())

    model_connectivity = np.zeros(shape=(len(flat_adj_mats), len(flat_adj_mats)))

    for m_idx in model_idxs:
        if use_nearest_neighbour:
            neighbours = find_nearest_neighbours(flat_adj_mats[m_idx], flat_adj_mats)

        elif use_adjacent_neighbour:
            neighbours = find_adjacent_neighbours(flat_adj_mats[m_idx], flat_adj_mats)

        for n in neighbours:
            model_connectivity[m_idx, n] = 1

    graph = nx.convert_matrix.from_numpy_matrix(model_connectivity)

    max_marginal = np.max(model_marginals)
    min_marginal = 0

    scale_max = 50
    scale_min = 10
    node_sizes = [min_max_scaling(scale_min, scale_max, min_marginal, max_marginal, x) for x in model_marginals]

    if colour_by_motif:
        col_blue = '#1d71b8' # SL blue
        col_red = '#ac182b' # OL red
        col_grey = '#706f6f' # Mixed grey 

        model_space_report_df = make_model_motif_colours(model_space_report_df, col_SL_only=col_blue, col_OL_only=col_red, col_mixed=col_grey)
        node_colours = model_space_report_df['node_type'].values

    else:
        node_colours = ['#73ba65' for _ in model_idxs]

    # #1d71b8 
    # #ac182b - OL red
    # #706f6f - grey

    for n, data in graph.nodes(data=True):
        data['size'] = node_sizes[n]
        data['color'] = node_colours[n]
        print(data)

    vis_config = two_species_vis_parameters()

    graph, config = netwulf.visualize(graph, config=vis_config, plot_in_cell_below=False, is_test=False)

    if colour_by_motif:
        output_name = output_name + "_col"

    output_name = output_name + ".pdf"
    output_path = combined_analysis_output_dir + output_name

    print("plotting network... ")

    fig, ax = netwulf.draw_netwulf(graph, figsize=15.0)
    fig.tight_layout()

    plt.savefig(output_path)


def get_flat_adjacency_matricies_df(model_indexes, adj_matrix_path_template):
    flat_adj_mats = []

    # Unpack all adjacency matricies
    for m_idx in model_indexes:
        adj_mat_path = adj_matrix_path_template.replace("#REF#", str(m_idx))
        adj_mat_df = pd.read_csv(adj_mat_path, index_col=0)

        print(adj_mat_df)
        adj_mat_df = convert_QS_column(adj_mat_df)
        print(adj_mat_df)
        print(np.shape(adj_mat_df))
        exit()
        flat_adj_mats.append(abs(adj_mat_df.values).flatten())
    
    interaction_indexes = ["a_" + str(i) for i in range(np.shape(flat_adj_mats)[1])]

    # Make adjacency matrix dataframe
    columns = ["model_idx"] + interaction_indexes
    adj_mat_df = pd.DataFrame(columns=columns)
    adj_mat_df['model_idx'] = model_indexes
    adj_mat_df.loc[:, adj_mat_df.columns != 'model_idx'] = flat_adj_mats

    adj_mat_df.set_index('model_idx')
    adj_mat_df = adj_mat_df.loc[:, adj_mat_df.columns != 'model_idx']
    

    # Remove columns that are only zeros or all the same
    for col in adj_mat_df.columns:
        if np.max(adj_mat_df[col]) == np.min(adj_mat_df[col]):
            adj_mat_df = adj_mat_df.loc[:, adj_mat_df.columns != col]


    print(np.shape(adj_mat_df))
    exit()

    return adj_mat_df


def find_next_children(start_path_array, path_arrays):
    line_end_x_1 = start_path_array.vertices[0][0]
    line_end_x_2 = start_path_array.vertices[3][0]

    line_min_y = start_path_array.vertices[0][1]
    line_max_y = start_path_array.vertices[2][1]

    next_children = []

    for candidate in path_arrays:
        candidate_start = candidate.vertices[1][0]
        candidate_min_y = candidate.vertices[0][1]
        candidate_max_y = candidate.vertices[2][1]
        candidate_end = candidate.vertices[0][0]

        # If start of line equals end of previous line
        if candidate_start == line_end_x_1 or candidate_start == line_end_x_2:
            is_connected = False

            # if candidate_end < line_end_x_1: # and candidate_min_y < line_min_y or candidate_max_y > line_max_y:
            #     is_connected = False

            # Lower branch
            if candidate_min_y < line_min_y and candidate_max_y > line_min_y:
                is_connected = True

            elif candidate_max_y > line_max_y and candidate_min_y < line_max_y:
                is_connected = True

            # if candidate_centre == line_min_y or candidate_centre == line_max_y:
            if is_connected:
                next_children.append(candidate)

    return next_children

def find_last_children(start_path_array, path_arrays, last_children, keep_all_descendents=True):
    line_end_x_1 = start_path_array.vertices[0][0]
    line_end_x_2 = start_path_array.vertices[3][0]

    line_min_y = start_path_array.vertices[0][1]
    line_max_y = start_path_array.vertices[2][1]

    is_last_child = True

    for candidate in path_arrays:
        candidate_start = candidate.vertices[1][0]
        candidate_min_y = candidate.vertices[0][1]
        candidate_max_y = candidate.vertices[2][1]
        candidate_end = candidate.vertices[0][0]


        # If start of line equals end of previous line
        if candidate_start == line_end_x_1 or candidate_start == line_end_x_2:
            is_connected = False

            # if candidate_end < line_end_x_1: # and candidate_min_y < line_min_y or candidate_max_y > line_max_y:
            #     is_connected = False

            # Lower branch
            if candidate_min_y < line_min_y and candidate_max_y > line_min_y:
                is_connected = True

            elif candidate_max_y > line_max_y and candidate_min_y < line_max_y:
                is_connected = True

            # if candidate_centre == line_min_y or candidate_centre == line_max_y:
            if is_connected:
                is_last_child = False

                if keep_all_descendents:
                    is_last_child = True

                # Find children from end line
                find_last_children(candidate, path_arrays, last_children, keep_all_descendents)

    if is_last_child:
        last_children.append(start_path_array)


def find_last_children_levels(start_path_array, path_arrays, last_children, max_level=5, current_level=0, keep_all_descendents=True):

    if current_level == max_level:
        if not keep_all_descendents:
            last_children.append(start_path_array)
        return 0

    line_end_x_1 = start_path_array.vertices[0][0]
    line_end_x_2 = start_path_array.vertices[3][0]

    line_min_y = start_path_array.vertices[0][1]
    line_max_y = start_path_array.vertices[2][1]

    is_last_child = True

    remove_candidates = []

    for candidate in path_arrays:
        candidate_start = candidate.vertices[1][0]
        candidate_min_y = candidate.vertices[0][1]
        candidate_max_y = candidate.vertices[2][1]
        candidate_end = candidate.vertices[0][0]

        candidate_centre = candidate_min_y + (candidate_max_y - candidate_min_y)
        
        is_connected = False

        # If start of line equals end of previous line
        if candidate_start == line_end_x_1 or candidate_start == line_end_x_2:

            # if candidate_end < line_end_x_1: # and candidate_min_y < line_min_y or candidate_max_y > line_max_y:
            #     is_connected = False

            # Lower branch
            if candidate_min_y < line_min_y and candidate_max_y > line_min_y:
                is_connected = True

            elif candidate_max_y > line_max_y and candidate_min_y < line_max_y:
                is_connected = True

            # if candidate_centre == line_min_y or candidate_centre == line_max_y:
            if is_connected:
                is_last_child = False

                if keep_all_descendents:
                    is_last_child = True

                    last_children.append(candidate)


                # Clean path arrays so we only include children (start level is not higher 
                # than the start of new candidate
                clean_path_arrays = [p for p in path_arrays if p.vertices[1][0] <= candidate_start]

                # Find children from end line
                find_last_children_levels(candidate, clean_path_arrays, last_children, keep_all_descendents=keep_all_descendents, 
                    current_level=current_level+1, max_level=max_level)


    if is_last_child:
        last_children.append(start_path_array)


def find_path_root(path_arrays):
    sorted(path_arrays, key=lambda x: x.vertices[0][0])

def assign_clusters(verticies, final_clusters):


    for f in final_clusters:
        print(f)
    # exit()
    min_f = 1000

    unique_ends = []

    # Get unique ends
    for p in verticies:
        if p[0][0] == 0:
            unique_ends.append(p[0][1])

        if p[3][0] == 0:
            unique_ends.append(p[3][1])

    unique_ends = list(set(unique_ends))
    unique_ends = sorted(unique_ends)
    cluster_assignments = [-1 for x in range(len(unique_ends))]

    # Build singlet clusters
    for e_idx, end in enumerate(unique_ends):
        for f_idx, f in enumerate(final_clusters):
            if end >= f[0] and end <= f[1]:
                cluster_assignments[e_idx] = f_idx

            elif end == f[0] or end == f[1]:
                cluster_assignments[e_idx] = f_idx

        if cluster_assignments[e_idx] == -1:
            final_clusters.append([end, end])


    # Redorder clusters including the singlets
    final_clusters = sorted(final_clusters, key= lambda x: x[0])
    # print(final_clusters)
    # print(len(final_clusters))
    # final_ends = []
    # for x in final_clusters:
    #     final_ends.append(x[0])
    #     final_ends.append(x[1])

    # Repeat cluster assignment
    cluster_assignments = [-1 for x in range(len(unique_ends))]
    
    # Build singlet clusters
    for e_idx, end in enumerate(unique_ends):
        for f_idx, f in enumerate(final_clusters):
            if end >= f[0] and end <= f[1]:
                cluster_assignments[e_idx] = f_idx

            elif end == f[0] or end == f[1]:
                cluster_assignments[e_idx] = f_idx

    return cluster_assignments

def cut_dendrogram_by_level(line_collection, max_level, set_min_x=0.0):
    path_arrays = (line_collection.get_paths())
    root_path = sorted(path_arrays, key=lambda x: x.vertices[0][0])[-1]

    colour_mask = []
    for path in path_arrays:
        colour_mask.append('white')

    level_legal_paths = []
    last_children = []

    print("finding last children levels with descenents")

    # find_last_children(p, legal_path_arrays, last_children, keep_all_descendents=False)
    find_last_children_levels(root_path, path_arrays, last_children, max_level=max_level, 
        current_level=0, keep_all_descendents=True)

    level_legal_paths.extend(last_children)
    
    print("Cleaning...")

    clean_last_paths = []
    for l in level_legal_paths:
        if l not in clean_last_paths:
            clean_last_paths.append(l)

    print("Making colours")
    for idx, p in enumerate(path_arrays):
        if p in clean_last_paths:
            colour_mask[idx] = 'black'

        if p == root_path:
            colour_mask[idx] = 'black'

    print("Finding last children without descendents...")
    last_children = []
    find_last_children_levels(root_path, path_arrays, last_children, max_level=max_level, 
        current_level=0, keep_all_descendents=False)

    clean_last_paths = last_children

    print("Number of clusters: ", len(clean_last_paths))

    cluster_min_maxes = []
    cluster_final_children = []

    # Find last children of last legal paths
    for idx, p in enumerate(clean_last_paths):
        print(idx)
        last_children = []
        print("finding first children... ")
        new_children = find_next_children(p, path_arrays)

        # Remove children from path_arrays

        final_children = []
        next_children_batch = new_children
        i = 0

        p_min_y = p.vertices[0][1]
        p_max_y = p.vertices[2][1]
        p_min_x = p.vertices[0][0]
        child_iteration_count = 0

        while len(next_children_batch) != 0:
            child_iteration_count += 1

            path_arrays = [x for x in path_arrays if x not in new_children]

            new_children = next_children_batch[:]
            next_children_batch = []

            for child in new_children:
                batch_children = find_next_children(child, path_arrays)
                child_min_y = child.vertices[0][1]
                child_max_y = child.vertices[2][1]
                child_x_end_1 = child.vertices[0][0]
                child_x_end_2 = child.vertices[3][0]

                p_min_y = min([child_min_y, p_min_y])
                p_max_y = max([child_max_y, p_max_y])
                p_min_x = min([child_x_end_1, child_x_end_2, p_min_x])

                if len(batch_children) == 0:
                    final_children.append(child)

                else:
                    next_children_batch.extend(batch_children)

        p.vertices[0][1] = p_min_y
        p.vertices[1][1] = p_min_y
        p.vertices[2][1] = p_max_y
        p.vertices[3][1] = p_max_y
        p.vertices[0][0] = set_min_x
        p.vertices[3][0] = set_min_x
        cluster_final_children.append(final_children)

        cluster_min_maxes.append([p_min_y, p_max_y])

    return colour_mask, cluster_min_maxes, clean_last_paths

# Trims linkages
def cut_dendrogram(line_collection, max_x):
    path_arrays = (line_collection.get_paths())

    colour_mask = []
    line_starts = []

    for path in path_arrays:
        # print(path.vertices)
        # exit()
        # Start is legal, end is not legal
        if path.vertices[1][0] > max_x and path.vertices[0][0] < max_x:
            colour_mask.append('white')
        
        elif path.vertices[1][0] > max_x and path.vertices[3][0] < max_x:
            colour_mask.append('white')

        # end is not legal
        elif path.vertices[0][0] < max_x:
            colour_mask.append('white')

        else:
            colour_mask.append('black')

    legal_path_arrays = [p for p, m in zip(path_arrays, colour_mask) if m == "black"]
    masked_path_arrays = [p for p, m in zip(path_arrays, colour_mask) if m == "white"]

    last_legal_paths = []
    # Get final lines that are below threshold
    for p in legal_path_arrays:
        last_children = []
        # find_last_children(p, legal_path_arrays, last_children, keep_all_descendents=False)
        find_last_children_levels(legal_path_arrays[-1], path_arrays, last_children, max_level=max_x, 
            current_level=0, keep_all_descendents=True)
        # last_legal_paths.extend(last_children[0])
        last_legal_paths.extend(last_children)


    clean_last_paths = []
    for l in last_legal_paths:
        if l not in clean_last_paths:
            clean_last_paths.append(l)

    for idx, p in enumerate(path_arrays):
        if p in clean_last_paths:
            colour_mask[idx] = 'red'

    print("Last path")


    # Find last children of last legal paths
    for p in clean_last_paths:
        last_children = []

        find_last_children(p, path_arrays, last_children)
        last_children = last_children[:-1]

        p_min_y = p.vertices[0][1]
        p_max_y = p.vertices[2][1]

        for child_path in last_children:
            if child_path == p:
                continue

            child_min_y = child_path.vertices[0][1]
            child_max_y = child_path.vertices[2][1]

            p_min_y = min([child_min_y, p_min_y])
            p_max_y = max([child_max_y, p_max_y])



        p.vertices[0][1] = p_min_y
        p.vertices[1][1] = p_min_y
        p.vertices[2][1] = p_max_y
        p.vertices[3][1] = p_max_y
        p.vertices[0][0] = 0.0
        p.vertices[3][0] = 0.0

    return colour_mask


def plot_tree(P, pos=None):
    plt.clf()
    icoord = scipy.array(P['icoord']) # Branches
    dcoord = scipy.array(P['dcoord']) # Tree
    print(dcoord)
    color_list = scipy.array(P['color_list'])
    xmin, xmax = icoord.min(), icoord.max()
    ymin, ymax = dcoord.min(), dcoord.max()
    if pos:
        icoord = icoord[pos]
        dcoord = dcoord[pos]
        color_list = color_list[pos]
    for xs, ys, color in zip(icoord, dcoord, color_list):
        plt.plot(xs, ys, color)
    plt.xlim(xmin-10, xmax + 0.1*abs(xmax))
    plt.ylim(ymin, ymax + 0.1*abs(ymax))
    plt.show()

def make_hierarchical_clustering(combined_analysis_output_dir, adj_mat_dir, drop_eqless=-1, hide_x_ticks=True, max_level=1000, use_bar=False, plot_error=True, average_clusters=False, log_scale=False):
    model_space_report_path = combined_analysis_output_dir + "combined_model_space_report_with_motifs.csv"
    adj_matrix_path_template = adj_mat_dir + "model_#REF#_adj_mat.csv"

    print("Reading df...")
    model_space_report_df = pd.read_csv(model_space_report_path)
    model_space_report_df = model_space_report_df.loc[:, ~model_space_report_df.columns.str.contains('^Unnamed')]
    model_space_report_df = model_space_report_df.sort_values('model_idx').reset_index(drop=True)

    model_idxs = model_space_report_df['model_idx'].values
    model_marginals = model_space_report_df['model_marginal_mean'].values

    # Generate df of adjacency matricies
    print("Getting flat adj mats...")
    adj_mat_df = get_flat_adjacency_matricies_df(model_idxs, adj_matrix_path_template)

    # Make hierarchical linkage 
    print("Generating linkage... ")
    linkage = hc.linkage(adj_mat_df.values, method='average', optimal_ordering=True)

    # dmin = 1.2
    # dmax = 10
    # pos = scipy.all( (linkage[:,2] >= dmin, linkage[:,2] <= dmax), axis=0 ).nonzero()
    # print(pos)
    print("Making interactions dendrogram and heatmap")
    colors = ["black", "faded green"]
    print(sns.xkcd_palette(colors))
    pal = sns.xkcd_palette(colors)
    pal[0] = "#333333"

    # model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] <= drop_eqless].index, inplace=True)

    # Plot clustering for using linkage
    print("Making cluster map... ")
    cm = sns.clustermap(
        data=adj_mat_df, row_linkage=linkage, col_cluster=False, 
        cmap=pal, yticklabels=1, xticklabels=1,
        cbar_kws={"ticks":[]}, linewidths=0.0, linecolor='white'
        )

    print("Pruning dendrogram... ")
    original_verticies = [p.vertices for p in cm.ax_row_dendrogram.collections[0].get_paths()]

    for line in cm.ax_row_dendrogram.collections:
        colour_mask, cluster_min_maxes, cluster_final_children = cut_dendrogram_by_level(line, max_level=max_level, set_min_x=0.0)
        cluster_assignments = assign_clusters(original_verticies, cluster_min_maxes)

        # print(mask)
        # print(len(mask))
        # colours = ['black' if m else 'white' for m in mask]
        line.set_color(colour_mask)

        # line.set_visible

    # exit()

    # Assign cluster ranks
    model_index_order = [model_space_report_df.ix[int(t.get_text())]['model_idx'] for t in cm.ax_heatmap.get_yticklabels()]
    # print(model_index_order)
    # exit()
    sorterIndex = dict(zip(model_index_order, range(len(model_index_order))))
    model_space_report_df['cluster_rank'] = model_space_report_df['model_idx'].map(sorterIndex)
    model_space_report_df.sort_values(['cluster_rank'], inplace=True)
    model_space_report_df['cluster_assignment'] = cluster_assignments
    model_space_report_df.to_csv(combined_analysis_output_dir + 'cluster_model_space_report.csv')

    hm = cm.ax_heatmap.get_position()

    cm.ax_heatmap.set_position([hm.x0 * 0.9, hm.y0, hm.width*0.25, hm.height])
    cm.ax_heatmap.set(yticklabels=[])
    cm.ax_heatmap.set(xticklabels=[])
    # cm.ax_heatmap.set(linewidths=0.5, rasterized=True)

    cm.ax_heatmap.set(xlabel='')
    cm.ax_heatmap.set(ylabel='')
    cm.ax_heatmap.tick_params(left=False, bottom=False, right=False)
    cm.cax.set_visible(False)

    col = cm.ax_col_dendrogram.get_position()
    cm.ax_col_dendrogram.set_position([col.x0 * 0.9, col.y0, col.width*0.25, col.height])
    # cm.ax_row_dendrogram.axvline(cut_off)

    col = cm.ax_col_dendrogram.get_position()
    hm = cm.ax_heatmap.get_position()

    # row = cm.ax_row_dendrogram.get_position()
    # cm.ax_row_dendrogram.set_position([row.x0, row.y0, row.width*0.25, row.height])

    int_link_x = hm.x0
    int_link_y0 = hm.y0
    int_link_width = hm.width
    int_link_height = hm.height

    output_path = combined_analysis_output_dir + "interactions_dendrogram_" + str(max_level) + ".pdf"
    plt.savefig(output_path, dpi=500, bbox_inches='tight', transparent=False)
    plt.close()

    # print("Making posterior probability heatmap")
    # # Plot with posterior probability as heatmap
    # adj_mat_df['p_prob'] = model_space_report_df['model_marginal_mean'].values
    # linkage_data = adj_mat_df.loc[:, adj_mat_df.columns != 'p_prob'].values
    # linkage = hc.linkage(sp.distance.pdist(linkage_data), method='average', metric='euclidean')
    # cm = sns.clustermap(data=adj_mat_df['p_prob'], row_linkage=linkage, col_cluster=False, yticklabels=1)
    # pprob_row_colours = cm.row_colors
    # hm = cm.ax_heatmap.get_position()
    # cm.ax_heatmap.set_position([int_link_x, int_link_y0, int_link_width*0.1, int_link_height])
    
    # # cm.ax_heatmap.set(yticklabels=[])
    # # cm.ax_heatmap.set(xticklabels=[])

    # # cm.ax_heatmap.set(xlabel='')
    # # cm.ax_heatmap.set(ylabel='')
    # cm.ax_heatmap.tick_params(left=False, bottom=False, right=False)

    # col = cm.ax_col_dendrogram.get_position()
    # cm.ax_heatmap.tick_params(left=False, bottom=False, right=False)

    # cm.ax_row_dendrogram.set_visible(False)
    # cm.cax.set_visible(False)

    # plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

    # output_path = combined_analysis_output_dir + "pprob_dendrogram_" + str(max_level) + ".pdf"
    # plt.savefig(output_path, dpi=500, bbox_inches='tight', transparent=True)
    # plt.close()


    print("Making posterior probability barchart")

    # Sort model space report df
    sorterIndex = dict(zip(model_index_order, range(len(model_index_order))))
    model_space_report_df['cluster_rank'] = model_space_report_df['model_idx'].map(sorterIndex)
    model_space_report_df.sort_values(['cluster_rank'], inplace=True)

    # print(model_space_report_df.cluster_rank)

    bar_color = sns.xkcd_palette(["amber"])[0]
    bar_color = "#333333"
    width_inches = 77.5*1.5 / 25.4
    height_inches = 19*1.5 / 25.4

    fig, ax = plt.subplots(figsize=(width_inches, height_inches))

    unique_cluster_assignments = model_space_report_df['cluster_assignment'].unique()

    if average_clusters:
        for c in unique_cluster_assignments:
            sub_df = model_space_report_df.loc[model_space_report_df['cluster_assignment'] == c]
            average_marg = np.sum(sub_df['model_marginal_mean'].values) / len(sub_df)
            model_space_report_df.loc[model_space_report_df['cluster_assignment'] == c, 'model_marginal_mean'] = average_marg


    if use_bar:
        model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] <= drop_eqless].index, inplace=True)
        sns.barplot(x='cluster_rank', y='model_marginal_mean', data=model_space_report_df, ax=ax, dodge=False, color=bar_color, linewidth=0)

    else:
        model_space_report_df.drop(model_space_report_df[model_space_report_df['model_marginal_mean'] <= drop_eqless].index, inplace=True)
        sns.scatterplot(x='cluster_rank', y='model_marginal_mean', s=1, marker='s',
        data=model_space_report_df, alpha=1.0, ax=ax, color=bar_color, linewidth=0) #, palette="BuGn_r"
        # print(dropped_m_space.cluster_rank)
    despine_offfset = {'top': 0, 'bottom': 0, 'right': 0, 'left': 1.0}
    sns.despine(offset=despine_offfset, trim=False)

    ax.tick_params(left=True, bottom=False, right=False)

    if plot_error:
        ax.errorbar(x=model_space_report_df['cluster_rank'], 
                    y=model_space_report_df['model_marginal_mean'], 
                    yerr=model_space_report_df['model_marginal_std'], fmt='None', color='black', alpha=1,  join=False,
                    label=None, elinewidth=0.1)

    print(model_space_report_df.cluster_rank)

    ax.unicode_minus = True

    if hide_x_ticks:
        ax.set(xticklabels=[])
        ax.set(xlabel='')
        ax.legend().remove()    
    
    else:
        # ax.set(xticklabels=model_space_report_df['model_idx'])
        ax.set(xlabel='Model')
        ax.set_xticklabels(model_space_report_df['model_idx'], fontsize = 1)
        ax.legend()

    if log_scale:
        ax.set(yscale="log")

    ax.set_ylabel('')
    ax.set(xlim=(-0,  len(model_index_order)))
    ax.set(ylim=(-0))
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.tick_params(labelsize=5)

    ax.margins(x=0)
    ax.margins(y=0)

    if use_bar:
        out_name = "dend_order_posterior_bar_" + str(max_level)
        if log_scale:
            out_name = out_name + "_log"
        out_path_name = out_name + ".pdf"

    else:
        out_name = "dend_order_posterior_scatter_" + str(max_level)

        if log_scale:
            out_name = out_name + "_log"

        out_path_name =  out_name + ".pdf"

    output_path = combined_analysis_output_dir + out_path_name

    plt.tight_layout()
    plt.savefig(output_path, dpi=500, bbox_inches='tight', transparent=True)


def make_hierarchical_cluster_with_NMF_weights(combined_analysis_output_dir, NMF_dir, drop_eqless=-1, hide_x_ticks=True):
    model_space_report_path = combined_analysis_output_dir + "combined_model_space_report_with_motifs.csv"

    model_space_report_df = pd.read_csv(model_space_report_path)
    model_space_report_df = model_space_report_df.loc[:, ~model_space_report_df.columns.str.contains('^Unnamed')]
    model_space_report_df = model_space_report_df.sort_values('model_idx').reset_index(drop=True)

    model_idxs = model_space_report_df['model_idx'].values
    model_marginals = model_space_report_df['model_marginal_mean'].values

    for file in glob.glob(NMF_dir + "W_*.csv"):
        # Generate df of adjacency matricies
        w_df = pd.read_csv(file)
        w_df = w_df.loc[:, ~w_df.columns.str.contains('^Unnamed')]

        # reorder to match modelspace df
        w_df = w_df.sort_values('model_idx').reset_index(drop=True)
        w_df = w_df.loc[:, w_df.columns != 'model_idx']
        w_values = w_df.values

        # Make hierarchical linkage 
        linkage = hc.linkage(w_values, method='average', optimal_ordering=True)

        print("Making interactions dendrogram and heatmap")
        colors = ["faded green", "white"]


        # Plot clustering for using linkage
        cm = sns.clustermap(
            data=w_df, row_linkage=linkage, col_cluster=False, yticklabels=1,
            cbar_kws={"ticks":[0.0, 0.5, 1.0]}, linewidths=0.1, linecolor='grey'
            )

        # Assign cluster ranks
        model_index_order = [model_space_report_df.ix[int(t.get_text())]['model_idx'] for t in cm.ax_heatmap.get_yticklabels()]

        sorterIndex = dict(zip(model_index_order, range(len(model_index_order))))
        model_space_report_df['cluster_rank'] = model_space_report_df['model_idx'].map(sorterIndex)
        model_space_report_df['cluster_rank'] = model_space_report_df['cluster_rank'].values + 1
        model_space_report_df.sort_values(['cluster_rank'], inplace=True)

        hm = cm.ax_heatmap.get_position()

        cm.ax_heatmap.set_position([hm.x0, hm.y0, hm.width*0.25, hm.height])
        cm.ax_heatmap.set(yticklabels=[])
        cm.ax_heatmap.set(xticklabels=[])

        cm.ax_heatmap.set(xlabel='')
        cm.ax_heatmap.set(ylabel='')
        cm.ax_heatmap.tick_params(left=False, bottom=False, right=False)
        # cm.cax.set_visible(False)

        col = cm.ax_col_dendrogram.get_position()
        cm.ax_col_dendrogram.set_position([col.x0, col.y0, col.width*0.25, col.height])


        col = cm.ax_col_dendrogram.get_position()
        hm = cm.ax_heatmap.get_position()

        output_path = combined_analysis_output_dir + "NMF_dendrogram_" + Path(file).stem + "_.pdf"
        plt.savefig(output_path, dpi=500, bbox_inches='tight', transparent=True)
        plt.close()


        # Sort model space report df
        sorterIndex = dict(zip(model_index_order, range(len(model_index_order))))
        model_space_report_df['cluster_rank'] = model_space_report_df['model_idx'].map(sorterIndex)
        model_space_report_df['cluster_rank'] = model_space_report_df['cluster_rank'].values + 1
        model_space_report_df.sort_values(['cluster_rank'], inplace=True)

        # print(model_space_report_df.cluster_rank)

        bar_color = sns.xkcd_palette(["amber"])[0]
        width_inches = 77.5*4 / 25.4
        height_inches = 19*4 / 25.4

        fig, ax = plt.subplots(figsize=(width_inches, height_inches))
        sns.barplot(x='cluster_rank', y='model_marginal_mean', data=model_space_report_df, ax=ax, dodge=False, color=bar_color, linewidth=0)
        ax.tick_params(left=True, bottom=False, right=False)

        ax.errorbar(model_space_report_df.cluster_rank, 
                    model_space_report_df['model_marginal_mean'], 
                    yerr=model_space_report_df['model_marginal_std'], fmt=',', color='black', alpha=1,
                    label=None, elinewidth=0.1)

        ax.unicode_minus = True

        if hide_x_ticks:
            ax.set(xticklabels=[])
            ax.set(xlabel='')
            ax.legend().remove()    
        
        else:
            ax.set(xticklabels=model_space_report_df['model_idx'])
            ax.set(xlabel='Model')
            ax.legend()

        ax.set_ylabel('')
        # ax.set(xlim=(-0.5,None))
        ax.set(ylim=(-0))
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["left"].set_alpha(0.5)
        ax.spines["bottom"].set_alpha(0.5)
        ax.tick_params(labelsize=15)
        ax.margins(x=0)
        ax.margins(y=0)

        output_path = combined_analysis_output_dir + "dend_order_posterior_" + Path(file).stem + "_.pdf"
        plt.savefig(output_path, dpi=500, bbox_inches='tight', transparent=True)



def main():
    adj_mat_dir = "/Users/behzakarkaria/Documents/UCL/barnes_lab/PhD_project/research_code/model_sel_doe/input_files/input_files_two_species_5/adj_matricies/"
    combined_analysis_output_dir = "/Volumes/Samsung_T5/BK_manu_data_backup/output/two_species_5_SMC_0b/experiment_analysis/"
    NMF_dir = combined_analysis_output_dir + "nmf_analysis/"

    # adj_mat_dir = "/Users/behzakarkaria/Documents/UCL/barnes_lab/PhD_project/research_code/model_sel_doe/input_files/input_files_three_species_7/adj_matricies/"
    # combined_analysis_output_dir = "/Volumes/Samsung_T5/BK_manu_data_backup/output/three_species_7_SMC_5/experiment_analysis/"
    # make_hierarchical_cluster_with_NMF_weights(combined_analysis_output_dir, NMF_dir, drop_eqless=-1)

    for levels in [5]:
        make_hierarchical_clustering(combined_analysis_output_dir, adj_mat_dir, drop_eqless=-1, hide_x_ticks=True, max_level=levels, use_bar=True, plot_error=False, average_clusters=True, log_scale=True)

    exit()
    make_model_network(combined_analysis_output_dir, adj_mat_dir, use_adjacent_neighbour=True)
    exit()
    exit()
    exit()
    make_model_network(combined_analysis_output_dir, adj_mat_dir, use_nearest_neighbour=True)

    make_model_network(combined_analysis_output_dir, adj_mat_dir, use_adjacent_neighbour=True, colour_by_motif=True)
    make_model_network(combined_analysis_output_dir, adj_mat_dir, use_nearest_neighbour=True, colour_by_motif=True)

if __name__ == "__main__":
    main()




