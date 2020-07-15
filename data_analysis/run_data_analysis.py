import argparse
import yaml
from shutil import copy
from glob import glob
from pathlib import Path

from . import data_utils
from . import data_analysis_ABCSMC as ABC_DA
from . import NMF_analysis
from . import nearest_neighbours
from . import network_vis

def main(data_analysis_config):
    experiment_name = data_analysis_config['experiment_name']
    inputs_dir = data_analysis_config['inputs_folder']
    output_dir = data_analysis_config['output_folder']

    combined_analysis_output_dir = output_dir + experiment_name + "/experiment_analysis/"
    data_utils.make_folder(combined_analysis_output_dir)

    adj_mat_dir = inputs_dir + "adj_matricies/"

    exp_dir = output_dir + experiment_name + "/"
    clean_exp_repeat_dirs = []
    final_pop_dirs = []
    first_pop_dirs = []
    
    # Find experiments that have finished (they have model space report written in final population)
    finished_exp_final_population_dirs = ABC_DA.find_finished_experiments(exp_dir)

    # Get final population directories of each repeat experiment
    final_pop_dirs = ABC_DA.find_latest_pop_dirs(exp_dir)


    model_space_report_df = ABC_DA.write_combined_model_space_report(finished_exp_final_population_dirs, combined_analysis_output_dir) 
    model_space_report_df = ABC_DA.write_combined_model_space_with_motif_counts(finished_exp_final_population_dirs, combined_analysis_output_dir, adj_mat_dir, window=0, normalise=False, plot=False)

    network_vis.make_hierarchical_clustering(combined_analysis_output_dir, adj_mat_dir, drop_eqless=-1, hide_x_ticks=True, max_level=1000, 
        use_bar=True, plot_error=True, average_clusters=False, log_scale=False)

    ABC_DA.generate_marginal_probability_distribution(finished_exp_final_population_dirs, 
        combined_analysis_output_dir, hide_x_ticks=True, show_median=False, show_BF=False, drop_eqless=-1)

    ABC_DA.compare_top_models_by_parts(combined_analysis_output_dir, adj_mat_dir, 
    	combined_analysis_output_dir, drop_eqless=-1, 
    	show_median=False, hide_x_ticks=False)

    nearest_neighbours.get_motif_neighbours(combined_analysis_output_dir, load_pickle=False, remove_zero_change=True)
    NMF_analysis.nmf_decomposition(combined_analysis_output_dir, adj_mat_dir)
    ABC_DA.combine_model_params(finished_exp_final_population_dirs, combined_analysis_output_dir)
    ABC_DA.plot_all_model_param_distributions(combined_analysis_output_dir, inputs_dir, combined_analysis_output_dir + "dist_plots/")




