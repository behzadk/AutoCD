from .model_space import ModelSpace
from . import algorithm_utils as alg_utils

from . import population_modules

import numpy as np
import os
import csv
from timeit import default_timer as timer
import sys
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pylab
import matplotlib.pyplot as plt
from . import plotting
from scipy.optimize import fsolve
import copy
import pickle
import math
import time

from pyunicorn.timeseries import RecurrencePlot, RecurrenceNetwork
import seaborn as sns


## Routines to run optimisation/model selection algorithms.
# All routines are initialised with the same information, these are accessed by member functions tp
# run different algorithms.
class ABC:
    ## Constructor
    # Arguements for ABC object initialisation containing all information needed to perform ABC rejection or ABC SMC
    # algorithms. The input model list should be a list of Model objects.
    #
    # @param t_0 initial time point
    # @param t_end final time point
    # @param dt time step size
    # @param exp_num experiment number used as suffix for recording data
    # @param model_list list of Model objects in the model space
    # @param population_size Number of particles that must be accepted in a population
    # @param n_sims_batch Number of simulations to conduct in batch
    # @param fit_species List of species indexes that will be fitted.
    # @param initial_epsilon Final distance thresholds
    # @param initial_epsilon Initial distance thresholds
    # @param distance_function_mode Defines which distance objective to run
    # @param n_distances Number of distances per fit species
    # @param abs_tol Absolute tolerance of error stepper
    # @param rel_tol Relative tolerance of error stepper
    # @param out_dir Output directory
    def __init__(self, t_0, t_end, dt, exp_num,
                 model_list, population_size, n_sims_batch,
                 fit_species, final_epsilon, initial_epsilon, distance_function_mode, n_distances, abs_tol, rel_tol,
                 out_dir):
        self.t_0 = t_0
        self.t_end = t_end
        self.dt = dt
        self.model_list = model_list

        self.exp_num = exp_num

        self.population_size = population_size
        self.population_accepted_count = 0
        self.population_total_simulations = 0
        self.n_sims_batch = n_sims_batch
        self.fit_species = fit_species
        self.distance_function_mode = distance_function_mode
        self.n_distances = n_distances

        self.population_accepted_particle_distances = []
        self.population_accepted_particles_count = 0

        self.population_model_refs = []
        self.population_judgements = []

        self.population_accepted_particles = []

        self.abs_tol = abs_tol
        self.rel_tol = rel_tol

        # self.epsilon = [100, 10, 1e4]
        # self.epsilon = [0.01, 0.001, 1e-5]

        self.final_epsilon = final_epsilon

        self.current_epsilon = initial_epsilon
        print("Initial epsilon: ", self.current_epsilon)

        # Init model space
        self.model_space = ModelSpace(model_list)

        self.out_dir = out_dir
        self.population_number = 0
        self.finished = False

    ## Saves ABC object as pickle.
    # named checkpoint.pickle
    # @param output_dir Output directory
    def save_object_pickle(self, output_dir):
        pickle_name = output_dir + "checkpoint.pickle"
        with open(pickle_name, 'wb') as handle:
            pickle.dump(self, handle, protocol=-1)

    ## Plots accepted particles.
    # For a given batch of simulated particles, the distance judgements are iterated. If a particle has been accepted
    # it will be plotted.
    #
    # @param out_dir Output directory
    # @param pop_num population number used to generate plot name
    # @param batch_num batch number used to generate plot name
    # @param part_judgements list of bools dictating whether a particle was accepted or rejected
    # @param init_states list of initial states for each particle
    # @param model_refs list containing the model refs for each particle.
    def plot_accepted_particles(self, out_dir, pop_num, batch_num, part_judgements, init_states, model_refs):
        out_path = out_dir + "Population_" + str(pop_num) + "_batch_" + str(batch_num) + "_accepted_plots.pdf"

        if sum(part_judgements) == 0:
            return 0

        # Make new pdf
        with PdfPages(out_path) as pdf:
            # Iterate all particles in batch
            for sim_idx, is_accepted in enumerate(part_judgements):
                if is_accepted:
                    state_list = self.pop_obj.get_particle_state_list(sim_idx)
                    time_points = self.pop_obj.get_timepoints_list()

                else:
                    continue

                try:
                    state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))

                except(ValueError):
                    # print(len(state_list)/len(init_states[sim_idx]))
                    time_points = range(int(len(state_list) / len(init_states[sim_idx])))
                    state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))

                model_ref = model_refs[sim_idx]

                plot_species = self.fit_species
                error_msg = self.pop_obj.get_particle_integration_error(sim_idx)

                plotting.plot_simulation(pdf, sim_idx, model_ref, state_list, time_points, plot_species, error_msg)

        # pdf.close()

    ## Plots all particles.
    # For a given batch of simulated particles, all particles are plotted
    # @param out_dir Output directory
    # @param pop_num population number used to generate plot name
    # @param batch_num batch number used to generate plot name
    # @param init_states list of initial states for each particle
    # @param model_refs list containing the model refs for each particle.
    def plot_all_particles(self, out_dir, pop_num, batch_num, init_states, model_refs):
        out_path = out_dir + "/simulation_plots/Population_" + str(pop_num) + "_batch_" + str(
            batch_num) + "all_plots.pdf"

        # Make new pdf
        pdf = PdfPages(out_path)

        negative_count = 0
        # Iterate all particles in batch
        for sim_idx, m_ref in enumerate(model_refs):
            state_list = self.pop_obj.get_particle_state_list(sim_idx)
            time_points = self.pop_obj.get_timepoints_list()

            try:
                state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))

            except(ValueError):
                # print(len(state_list)/len(init_states[sim_idx]))
                time_points = range(int(len(state_list) / len(init_states[sim_idx])))
                state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))

            if np.min(state_list) < 0 or np.isnan(state_list).any():
                negative_count += 1

            model_ref = model_refs[sim_idx]
            error_msg = self.pop_obj.get_particle_integration_error(sim_idx)

            plot_species = [i for i in self.fit_species]
            # plot_species = [i for i in range(np.shape(state_list)[1])]

            plotting.plot_simulation(pdf, sim_idx, model_ref, state_list, time_points, plot_species, error_msg)

        print("negative count: ", negative_count)
        pdf.close()

    ## Writes particle state list to csv.
    # For a given batch of simulated particles, the complete states of all simulations are written to a .csv file
    # @param out_dir Output directory
    # @param pop_num population number used to generate file name
    # @param batch_num batch number used to generate file name
    # @param init_states list of initial states for each particle
    # @param model_refs list containing the model refs for each particle.
    def write_all_particle_state_lists(self, out_dir, pop_num, batch_num, init_states, model_refs):
        for sim_idx, m_ref in enumerate(model_refs):
            out_path = out_dir + "/simulation_states/Population_" + str(pop_num) + "_batch_" + \
                       str(batch_num) + "_idx_" + str(sim_idx) + "_state.csv"

            state_list = self.pop_obj.get_particle_state_list(sim_idx)
            time_points = self.pop_obj.get_timepoints_list()

            tp = int(len(state_list) / len(init_states[sim_idx]))
            state_list = np.reshape(state_list, (tp, len(init_states[sim_idx])))
            np.savetxt(out_path, state_list, delimiter=',')

    ## Writes accepted particle distances to csv.
    # For a given batch of simulated particles, the distances are written to .csv. One .csv contains all accepted
    # particle distances for a population.
    # @param out_dir Output directory
    # @param model_refs list containing the model refs for each particle.
    # @param part_judgements list of bools dictating whether a particle was accepted or rejected
    # @param distances list of distances for each particle and each species being fit.
    def write_accepted_particle_distances(self, out_dir, model_refs, part_judgments, distances):
        out_path = out_dir + "distances.csv"
        with open(out_path, 'a') as out_csv:
            wr = csv.writer(out_csv)
            for idx, is_accepted in enumerate(part_judgments):
                record_vals = [idx, model_refs[idx]]
                if is_accepted:
                    for n in self.fit_species:
                        for d in distances[idx][n]:
                            record_vals.append(d)

                    wr.writerow(record_vals)

    ## Writes all particle distances to csv.
    # For a given batch of simulated particles, the distances are written to .csv. One .csv contains all accepted
    # particle distances for a population.
    # @param out_dir Output directory
    # @param model_refs list containing the model refs for each particle.
    # @param batch_num batch number used to generate file name
    # @param pop_num population number used to generate file name
    # @param part_judgements list of bools dictating whether a particle was accepted or rejected
    # @param distances list of distances for each particle and each species being fit.
    def write_particle_distances(self, out_dir, model_refs, batch_num, pop_num, part_judgments, distances,
                                 only_accepted=False):
        out_path = out_dir + "distances.csv"

        # If file doesn't exist, write header
        if not os.path.isfile(out_path):
            col_header = ['sim_idx', 'batch_idx', 'population_num', 'exp_num', 'model_ref', 'Accepted', 'integ_error']
            idx = 1
            for n in self.fit_species:
                for d in self.current_epsilon:
                    col_header.append('d' + str(idx))
                    idx += 1

            with open(out_path, 'a') as out_csv:
                wr = csv.writer(out_csv, quoting=csv.QUOTE_NONNUMERIC)
                wr.writerow(col_header)

        # Write distances
        with open(out_path, 'a') as out_csv:
            wr = csv.writer(out_csv, quoting=csv.QUOTE_NONNUMERIC)

            for idx, m_ref in enumerate(model_refs):
                if only_accepted and part_judgments[idx] == False:
                    continue

                error_msg = self.pop_obj.get_particle_integration_error(idx)
                row_vals = [idx, batch_num, pop_num, self.exp_num, m_ref, part_judgments[idx], error_msg]

                for n_idx, n in enumerate(self.fit_species):
                    for d in distances[idx][n_idx]:
                        row_vals.append(d)

                wr.writerow(row_vals)

    ## Writes parameters and initial species for a given list of particles.
    # @param out_dir Output directory
    # @param batch_num batch number used to generate file name
    # @param simulated_particles list of simulated particle objects
    # @param input_params list of input parameters for each particle
    # @param input_init_species list of initial species values for each particle
    # @param judgement_array list of bools dictating whether a particle was accepted or rejected
    # @param particle_weights list of particle weights.
    def write_particle_params(self, out_dir, batch_num, simulated_particles,
                              input_params, input_init_species, judgement_array, particle_weights):
        for m in self.model_space._model_list:
            out_path = out_dir + "model_" + str(m.get_model_ref()) + "_all_params"

            # Add header if file does not yet exist
            if not os.path.isfile(out_path):
                col_header = ['sim_idx', 'batch_idx', 'Accepted', 'particle_weight'] + list(
                    sorted(m._params_prior, key=str.lower)) + list(
                    m._init_species_prior)
                with open(out_path, 'a') as out_csv:
                    wr = csv.writer(out_csv)
                    wr.writerow(col_header)

            # Write data
            with open(out_path, 'a') as out_csv:
                wr = csv.writer(out_csv)
                for idx, particle in enumerate(simulated_particles):
                    if m is particle:
                        wr.writerow(
                            [idx] + [batch_num] + [judgement_array[idx]] + [particle_weights[idx]] + input_params[idx] +
                            input_init_species[idx])

    def write_population_particle_params(self, out_dir):
        for m in self.model_space._model_list:
            out_path = out_dir + "model_" + str(m.get_model_ref()) + "_population_all_params"

            # Add header if file does not yet exist
            if not os.path.isfile(out_path):
                col_header = ['sim_idx', 'batch_idx', 'population_num', 'exp_num', 'particle_weight'] + list(
                    sorted(m._params_prior, key=str.lower)) + list(
                    m._init_species_prior)
                with open(out_path, 'a') as out_csv:
                    wr = csv.writer(out_csv)
                    wr.writerow(col_header)

            # Write data
            with open(out_path, 'a') as out_csv:
                wr = csv.writer(out_csv)
                for idx, particle in enumerate(self.population_accepted_particles):
                    if m._model_ref is particle.curr_model._model_ref:
                        wr.writerow(
                            [particle.sim_idx] + [particle.batch_idx] + [self.population_number] + [self.exp_num] + [
                                particle.curr_weight] + particle.curr_params + particle.curr_init_state)

    ## Writes the population distance thresholds.
    # @param out_dir Output directory
    # @param epsilon current distance thresholds
    def write_epsilon(self, out_dir, epsilon):
        out_path = out_dir + "epsilon.txt"

        if not os.path.isfile(out_path):
            col_header = ['e_' + str(idx) for idx, _ in enumerate(epsilon)]

            with open(out_path, 'a') as out_csv:
                wr = csv.writer(out_csv)
                wr.writerow(col_header)
                wr.writerow(epsilon)

    def calculate_separtion_coeffs(self, model_refs, input_params, init_states):
        init_theta_states = [copy.deepcopy(state) for state in init_states]

        # Peturb by init species theta_0
        for theta_state in init_theta_states:
            perturb_idx = np.random.choice(self.fit_species)
            perturb_idx = 2
            theta_state[perturb_idx] += self.theta_0

        separation_coefficients = [self.theta_0 for c in range(len(model_refs))]

        curr_states = [copy.deepcopy(state) for state in init_states]
        curr_theta_states = [copy.deepcopy(state) for state in init_theta_states]
        
        # Simulate step
        for x in range(self.sep_length):
            if len(curr_states) < self.n_sims_batch:
                adjusted_sims_batch = len(curr_states)

            else:
                adjusted_sims_batch = self.n_sims_batch

            pop_obj = population_modules.Population(adjusted_sims_batch, self.t_0, self.dt*2,
                                                         self.dt, curr_states, input_params, model_refs,
                                                         self.fit_species, self.abs_tol, self.rel_tol)

            pop_obj_theta = population_modules.Population(adjusted_sims_batch, self.t_0, self.dt*2,
                                                         self.dt, curr_theta_states, input_params, model_refs,
                                                         self.fit_species, self.abs_tol, self.rel_tol)
            pop_obj.generate_particles()
            pop_obj_theta.generate_particles()

            pop_obj.simulate_particles()
            pop_obj_theta.simulate_particles()

            # Get step separations
            for sim_idx in range(len(model_refs)):
                final_vals = pop_obj.get_particle_final_species_values(sim_idx)
                final_vals_theta = pop_obj_theta.get_particle_final_species_values(sim_idx)

                # Calculate separation
                theta_1 = 0
                for y_a, y_b in zip(final_vals, final_vals_theta):
                    theta_1 += (y_a - y_b)**2 

                theta_1 = theta_1**0.5
                theta_1 = theta_1

                # Set next step init states
                species_idx = 0
                for y_a, y_b in zip(final_vals, final_vals_theta):
                    # Set state for next sim
                    curr_states[sim_idx][species_idx] = y_a

                    # Adjusted state of y_theta for next sim
                    # curr_theta_states[sim_idx][species_idx] = y_a + (self.theta_0 / theta_1) * (y_b - y_a)
                    curr_theta_states[sim_idx][species_idx] = y_a + (self.theta_0 * (y_b - y_a)) / theta_1

                    species_idx += 1

                if x > self.n_skip:
                    # Add separation coeff to total
                    separation_coefficients[sim_idx] += np.log2(abs(theta_1/self.theta_0))

        # Make into mean
        separation_coefficients = [sep / (self.sep_length - self.n_skip) for sep in separation_coefficients]

        return separation_coefficients, init_theta_states

    def get_RQA_distances(self, model_refs, init_states, pop_num, batch_num, plot_RQA=False):
        distances = []

        # Iterate all particles in batch
        for sim_idx, m_ref in enumerate(model_refs):
            state_list = self.pop_obj.get_particle_state_list(sim_idx)
            time_points = self.pop_obj.get_timepoints_list()

            state_list = state_list[50000:]
            time_points = time_points[50000:]

            try:
                state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))


            except(ValueError):
                # print(len(state_list)/len(init_states[sim_idx]))
                time_points = range(int(len(state_list) / len(init_states[sim_idx])))
                state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))

            # state_list = state_list[40000:]
            # time_points = time_points[40000:]

            #  Settings for the recurrence plot
            RR = 0.05
            # Distance metric in phase space ->
            # Possible choices ("manhattan","euclidean","supremum")
            METRIC = "supremum"

            #  Generate a recurrence plot object with fixed recurrence rate RR
            rp = RecurrencePlot(state_list, metric=METRIC, 
                normalize=True, recurrence_rate=RR, silence_level=3)
            det = rp.determinism(l_min=2)
            entropy = rp.diag_entropy(l_min=2)
            laminarity = rp.laminarity(v_min=2)

            distances.append([det, entropy, laminarity])

            if plot_RQA:
                out_path_template = self.out_dir + "Population_" + str(self.population_number) + "/simulation_plots/RQA" + "_batch_" + str(self.batch_num) + "_idx_#SIM_IDX#_sep.pdf"
                out_path = out_path_template.replace('#SIM_IDX#', str(sim_idx))
                fig = plt.figure()
                plt.matshow(rp.recurrence_matrix())
                plt.savefig(out_path, dpi=500)
                plt.close()
                plt.clf()

        return distances


    def plot_RQA(self, model_refs, init_states, pop_num, batch_num):
        # Make new pdf
        negative_count = 0

        # Iterate all particles in batch
        for sim_idx, m_ref in enumerate(model_refs):
            state_list = self.pop_obj.get_particle_state_list(sim_idx)
            time_points = self.pop_obj.get_timepoints_list()

            try:
                state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))


            except(ValueError):
                # print(len(state_list)/len(init_states[sim_idx]))
                time_points = range(int(len(state_list) / len(init_states[sim_idx])))
                state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))

            # state_list = state_list[40000:]
            # time_points = time_points[40000:]

            #  Settings for the recurrence plot
            RR = 0.05
            # Distance metric in phase space ->
            # Possible choices ("manhattan","euclidean","supremum")
            METRIC = "supremum"

            #  Generate a recurrence plot object with fixed recurrence rate RR
            rp = RecurrencePlot(state_list,
                metric=METRIC, normalize=False, recurrence_rate=RR)
            print(rp.recurrence_rate())
            det = rp.determinism(l_min=2)
            entropy = rp.diag_entropy(l_min=2)
            rp_mat = rp.recurrence_matrix()
            print(rp_mat)
            print(det)
            print(entropy)
            out_path_template = self.out_dir + "Population_" + str(self.population_number) + "/simulation_plots/RQA" + "_batch_" + str(self.batch_num) + "_idx_#SIM_IDX#_sep.pdf"
            out_path = out_path_template.replace('#SIM_IDX#', str(sim_idx))
            fig = plt.figure()
            plt.matshow()
            plt.savefig(out_path, dpi=500)
            plt.close()
            plt.clf()


    def plot_separation(self, model_refs, states, theta_states, params, pop_num, batch_num):
        out_path_template = self.out_dir + "Population_" + str(self.population_number) + "/simulation_plots/Population_" + str(pop_num) + "_batch_" + str(batch_num) + "_idx_#SIM_IDX#_sep.pdf"

        if len(model_refs) < self.n_sims_batch:
            adjusted_sims_batch = len(model_refs)

        else:
            adjusted_sims_batch = self.n_sims_batch


        # Simulate to attractor
        pop_obj = population_modules.Population(adjusted_sims_batch, self.plot_t_0, self.plot_t_end,
                                                     self.dt, states, params, model_refs,
                                                     self.fit_species, self.abs_tol, self.rel_tol)

        # Simulate to attractor
        theta_pop_obj = population_modules.Population(adjusted_sims_batch, self.plot_t_0, self.plot_t_end,
                                                     self.dt, theta_states, params, model_refs,
                                                     self.fit_species, self.abs_tol, self.rel_tol)

        pop_obj.generate_particles()
        theta_pop_obj.generate_particles()

        pop_obj.simulate_particles()
        theta_pop_obj.simulate_particles()


        for sim_idx in range(len(model_refs)):
            state_list = pop_obj.get_particle_state_list(sim_idx)
            theta_state_list = theta_pop_obj.get_particle_state_list(sim_idx)

            n_species = len(states[sim_idx])

            time_points = pop_obj.get_timepoints_list()

            try:
                state_list = np.reshape(state_list, (len(time_points), n_species))

            except(ValueError):
                # print(len(state_list)/len(init_states[sim_idx]))
                idxs = range(int(len(state_list) / n_species))
                state_list = np.reshape(state_list, (len(idxs), n_species))
                time_points = time_points[:len(idxs)]


            try:
                theta_state_list = np.reshape(theta_state_list, (len(time_points), n_species))

            except(ValueError):
                # print(len(state_list)/len(init_states[sim_idx]))
                idxs = range(int(len(theta_state_list) / n_species))
                theta_state_list = np.reshape(theta_state_list, (len(idxs), n_species))
                time_points = time_points[:len(idxs)]



            model_ref = model_refs[sim_idx]
            plot_species = self.fit_species

            print("stat_list: ", np.shape(state_list))
            print("theta_state_list: ", np.shape(theta_state_list))

            # print(np.shape(time_points))
            time_points = time_points[self.n_skip_plot:]
            state_list = state_list[self.n_skip_plot:]
            theta_state_list = theta_state_list[self.n_skip_plot:]

            if len(state_list) < 1 or len(theta_state_list) < 1:
                print("stat_list: ", np.shape(state_list))
                print("theta_state_list: ", np.shape(theta_state_list))
                return 0

            elif len(state_list) != len(theta_state_list):
                print("stat_list: ", np.shape(state_list))
                print("theta_state_list: ", np.shape(theta_state_list))
                return 0

            else:
                out_path = out_path_template.replace('#SIM_IDX#', str(sim_idx))
                with PdfPages(out_path) as pdf:
                    plotting.plot_separation(pdf, sim_idx, model_ref, state_list, theta_state_list, time_points)

    def run_chaos_separation_model_selection_ABC_SMC(self, alpha=0.5):
        while not self.finished:
            folder_name = self.out_dir + "Population_" + str(self.population_number) + "/"

            try:
                os.mkdir(folder_name)
            except FileExistsError:
                pass

            sim_params_folder = folder_name + 'model_sim_params/'
            try:
                os.mkdir(sim_params_folder)
            except FileExistsError:
                pass

            try:
                os.mkdir(folder_name + 'simulation_plots')
            except FileExistsError:
                pass

            try:
                os.mkdir(folder_name + 'simulation_states')

            except FileExistsError:
                pass

            if self.current_epsilon == self.final_epsilon:
                print("Running final epsilon")
                self.finished = True

            # Reset for new population
            self.batch_num = 0

            while self.population_accepted_count < self.population_size:
                print("")
                start_time = time.time()
                start_time_sampling = time.time()

                # 1. Sample from model space
                if self.population_number == 0:
                    particles = self.model_space.sample_particles_from_prior(self.n_sims_batch)
                    init_weights = 1 / len(particles)
                    for sim_idx, p in enumerate(particles):
                        p.curr_weight = init_weights
                        p.prev_weight = init_weights
                        p.batch_idx = self.batch_num
                        p.sim_idx = sim_idx

                else:
                    print("Sampling new population")
                    # Sample and perturb particles from previous population
                    particles = self.model_space.sample_particles_from_previous_population(self.n_sims_batch)
                    for sim_idx, p in enumerate(particles):
                        p.batch_idx = self.batch_num
                        p.sim_idx = sim_idx

                end_time_sampling = time.time()
                print("Particle sampling time elapsed: ", end_time_sampling - start_time_sampling)
                

                init_states = [copy.deepcopy(p.curr_init_state) for p in particles]
                input_params = [copy.deepcopy(p.curr_params) for p in particles]
                model_refs = [copy.deepcopy(p.curr_model._model_ref) for p in particles]
                particle_weights = [p.prev_weight for p in particles]

                particle_models = [p.curr_model for p in particles]

                alg_utils.rescale_parameters(input_params, init_states, particle_models)

                self.pop_obj = population_modules.Population(self.n_sims_batch, self.t_0, self.t_end,
                                                             self.dt, init_states, input_params, model_refs,
                                                             self.fit_species, self.abs_tol, self.rel_tol)

                self.pop_obj.generate_particles()
                start_time_sim = time.time()
                self.pop_obj.simulate_particles()
                end_time_sim = time.time()

                # print("Particle simulation time elapsed: ", end_time_sim - start_time_sim)

                # 3. Calculate survival distances for population
                self.pop_obj.calculate_particle_distances(2)
                self.pop_obj.accumulate_distances()

                batch_distances = self.pop_obj.get_flattened_distances_list()
                batch_distances = np.reshape(batch_distances,
                                 (self.n_sims_batch, len(self.fit_species), 1))

                RQA_distances = self.get_RQA_distances(model_refs, init_states, self.population_number, self.batch_num, plot_RQA=True)
                # RQA_distances = self.get_RQA_distances(sep_coeff_model_refs, sep_coeff_init_state, self.population_number, self.batch_num, plot_RQA=False)
                

                # Extract state and params of surviving simulations
                sep_coeff_init_state = []
                sep_coeff_init_params = []
                sep_coeff_model_refs = []
                for sim_idx in range(len(model_refs)):
                    final_vals = self.pop_obj.get_particle_final_species_values(sim_idx)
                    sep_coeff_init_state.append(final_vals)
                    sep_coeff_init_params.append(input_params[sim_idx])
                    sep_coeff_model_refs.append(model_refs[sim_idx])

                separation_coefficients, init_theta_states = self.calculate_separtion_coeffs(sep_coeff_model_refs, sep_coeff_init_params, sep_coeff_init_state)
                self.plot_separation(sep_coeff_model_refs, sep_coeff_init_state, init_theta_states, 
                    sep_coeff_init_params, self.population_number, self.batch_num)

                n_spaces = self.population_size - len(self.population_accepted_particles)


                batch_part_judgements = [True for x in range(self.n_sims_batch)]

                # If we have more accepted than spaces, set some to false
                if n_spaces < sum(batch_part_judgements):
                    num_trim = sum(batch_part_judgements) - n_spaces
                    trimmed_judgements = []

                    for judgement in batch_part_judgements:
                        if judgement:
                            if sum(trimmed_judgements) < n_spaces:
                                trimmed_judgements.append(True)

                            else:
                                trimmed_judgements.append(False)
                        else: 
                            trimmed_judgements.append(False)

                    batch_part_judgements = trimmed_judgements

                accepted_particles = [p for p, judgement in zip(particles, batch_part_judgements) if judgement]

                combined_distances = []
                comp_exl_count = 0
                for particle_idx, p in enumerate(particles):
                    p_dist = list(np.concatenate( [batch_distances[particle_idx].reshape(-1), [separation_coefficients[particle_idx]], RQA_distances[particle_idx] ] )) * len(self.fit_species)
                    p_dist = np.reshape(p_dist, (len(self.fit_species), -1))
                    combined_distances.append(p_dist)

                self.population_accepted_particle_distances += [part_d for part_d, judgement in
                                                zip(combined_distances, batch_part_judgements)
                                                if judgement]

                self.population_accepted_particles = self.population_accepted_particles + accepted_particles
                self.population_total_simulations += len(particles)

                # print("Writing data")
                self.write_epsilon(folder_name, self.current_epsilon)


                start_time_write_distance = time.time()

                if self.final_epsilon == self.current_epsilon:
                    # self.plot_accepted_particles(folder_name, self.population_number, self.batch_num, batch_part_judgements, init_states, model_refs)
                    self.write_particle_distances(folder_name, model_refs, self.batch_num, self.population_number,
                                                  batch_part_judgements, combined_distances, only_accepted=True)

                end_time_write_distance = time.time()

                # print("Write distance time elapsed: ", end_time_write_distance - start_time_write_distance)

                self.population_accepted_count += sum(batch_part_judgements)
                # self.population_total_simulations += len(model_refs)

                self.model_space.update_population_sample_data_v2(model_refs, batch_part_judgements)

                print("Population: ", self.population_number, "Accepted particles: ", self.population_accepted_count,
                      "Total simulations: ", self.population_total_simulations)

                del self.pop_obj

                sys.stdout.flush()
                self.batch_num += 1
                end_time = time.time()

                print("Batch time elapsed: ", end_time - start_time)

            # Trim to population size
            self.population_accepted_particles = self.population_accepted_particles[: self.population_size]
            self.model_space.accepted_particles = self.population_accepted_particles

            if self.population_number == 0:
                for p in self.model_space.accepted_particles:
                    p.curr_weight = 1

            else:
                self.model_space.compute_particle_weights()

            self.model_space.normalize_particle_weights()
            self.model_space.update_model_marginals()

            if self.final_epsilon == self.current_epsilon:
                self.save_object_pickle(folder_name)
                self.write_population_particle_params(sim_params_folder)

            self.model_space.prepare_next_population()

            print("Generating model kernels")
            self.model_space.generate_model_kernels(self.population_accepted_particles, self.population_number)
            print("Generating aux info")

            self.model_space.generate_kernel_aux_info()
            self.model_space.count_dead_models()

            print("Generating model space report")
            print(folder_name)

            self.model_space.model_space_report(folder_name, use_sum=False)

            if self.current_epsilon != self.final_epsilon:
                self.current_epsilon = alg_utils.update_epsilon(self.current_epsilon, self.final_epsilon,
                                                                self.population_accepted_particle_distances, alpha)

                print("Current epsilon: ", self.current_epsilon)
                print("Starting new population... ")
                print("")

                self.population_number += 1
                self.population_accepted_count = 0
                self.population_total_simulations = 0
                self.population_accepted_particle_distances = []
                self.population_model_refs = []
                self.population_judgements = []
                self.population_accepted_particles = []
