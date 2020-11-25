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
import matplotlib
# matplotlib.use('Agg')
from . import plotting
from scipy.optimize import fsolve
import scipy as sp
import copy
import pickle
import math
import time

from pyunicorn.timeseries import RecurrencePlot, RecurrenceNetwork
import seaborn as sns

import multiprocessing

from itertools import product
import threading

plt.rcParams['figure.figsize'] = [15, 10]

font = {'size'   : 15, }
axes = {'labelsize': 'medium', 'titlesize': 'medium'}

sns.set_context("talk")
sns.set_style("white")
matplotlib.rc('font', **font)
matplotlib.rc('axes', **axes)


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
    def plot_all_particles(self, out_dir, pop_num, batch_num, init_states, model_refs, batch_judgements, inverse_max_exponents):
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

            if error_msg != '':
                continue

            plot_species = [i for i in self.fit_species]
            # plot_species = [i for i in range(np.shape(state_list)[1])]

            plotting.plot_simulation(pdf, sim_idx, model_ref, state_list, time_points, 
                plot_species, batch_judgements[sim_idx], 1/(inverse_max_exponents[sim_idx][0]) -1, error_msg)

            plotting.plot_LV_four_species(pdf, sim_idx, model_ref, state_list, time_points, plot_species, error_msg)

        print("negative count: ", negative_count)
        pdf.close()

    def plot_lorenz(self, out_dir, pop_num, batch_num, init_states, model_refs):
        out_path = out_dir + "/simulation_plots/Population_" + str(pop_num) + "_batch_" + str(batch_num) + "lorenz_3D.pdf"

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

            plotting.plot_lorenz(pdf, sim_idx, model_ref, state_list, time_points, plot_species, error_msg)

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
                if error_msg == '':
                    row_vals = [idx, batch_num, pop_num, self.exp_num, m_ref, part_judgments[idx], error_msg]

                    for n_idx, n in enumerate(self.fit_species):
                        for d in distances[idx][n_idx]:
                            row_vals.append(d)

                    wr.writerow(row_vals)

    def write_particle_data(self, out_dir, model_refs, batch_num, pop_num, data, batch_part_judgements, particle_weights):
        out_path = out_dir + "distances.csv"

        # If file doesn't exist, write header
        if not os.path.isfile(out_path):
            col_header = ['sim_idx', 'batch_idx', 'population_num', 'exp_num', 'model_ref', 'integ_error', 'accepted', 'particle_weight']
            idx = 1
            for d in data[0]:
                col_header.append('d' + str(idx))
                idx += 1

            with open(out_path, 'a') as out_csv:
                wr = csv.writer(out_csv, quoting=csv.QUOTE_NONNUMERIC)
                wr.writerow(col_header)

        # Write distances
        with open(out_path, 'a') as out_csv:
            wr = csv.writer(out_csv, quoting=csv.QUOTE_NONNUMERIC)

            for idx, m_ref in enumerate(model_refs):
                error_msg = self.pop_obj.get_particle_integration_error(idx)
                if error_msg == '':
                    row_vals = [idx, batch_num, pop_num, self.exp_num, m_ref, error_msg, batch_part_judgements[idx], particle_weights[idx]]

                    sim_data = data[idx]
                    for d in sim_data:
                        row_vals.append(d)

                    wr.writerow(row_vals)


    def write_particle_chaos_params(self, out_dir, batch_num, pop_num, simulated_particles,
                              input_params, input_init_species, batch_part_judgements):
        for m in self.model_space._model_list:
            out_path = out_dir + "model_" + str(m.get_model_ref()) + "_all_params.csv"

            # Add header if file does not yet exist
            if not os.path.isfile(out_path):
                col_header = ['sim_idx', 'batch_idx', 'population_num', 'exp_num', 'model_ref', 'integ_error', 'accepted'] + list(
                    sorted(m._params_prior, key=str.lower)) + list(
                    m._init_species_prior)
                with open(out_path, 'a') as out_csv:
                    wr = csv.writer(out_csv)
                    wr.writerow(col_header)

            # Write data
            with open(out_path, 'a') as out_csv:
                wr = csv.writer(out_csv)
                for idx, particle in enumerate(simulated_particles):
                    error_msg = self.pop_obj.get_particle_integration_error(idx)
                    if error_msg == '':
                        if m._model_ref is particle.curr_model._model_ref:

                            wr.writerow(
                                [idx, batch_num, pop_num, self.exp_num, m._model_ref, error_msg, batch_part_judgements[idx]] + input_params[idx] +
                                input_init_species[idx])


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
            out_path = out_dir + "model_" + str(m.get_model_ref()) + "_all_params.csv"

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

    def wolf_LE_est(self, model_refs, input_params, init_states, all_exponents_dict, all_timepoint_exponents_dict, sim_idx):
        err = self.pop_obj.get_particle_integration_error(sim_idx)
        if err:
            all_exponents_dict[sim_idx] = [np.nan for x in range(len(init_states[sim_idx]))]
            all_timepoint_exponents_dict[sim_idx] = [np.nan]
            return 0

        w = np.eye(len(init_states[sim_idx]), dtype=np.float64)

        exponents = [0 for x in range(len(init_states[sim_idx]))]
        time_point_exponents = [[] for x in range(len(init_states[sim_idx]))]
        num_species = len(init_states[sim_idx])

        # Get states
        state_list = self.pop_obj.get_particle_state_list(sim_idx)
        time_points = self.pop_obj.get_timepoints_list()
        state_list = np.reshape(state_list, (len(time_points), num_species))

        skip_percent = int(len(time_points) * 0.1)

        I_num_species = np.eye(num_species) 

        # Skip transient
        # state_list = state_list[self.n_skip:]
        # time_points = time_points[self.n_skip:]

        for t_idx, t in enumerate(time_points):
            state = list(state_list[t_idx])

            J = self.pop_obj.get_particle_jacobian(state, sim_idx)
            J = np.reshape(J, [num_species, num_species])

            # print(J)
            # exit()
            J = I_num_species + J * self.dt

            # w += np.inner(w, J)

            w = np.matmul(J, w)
            # print(w)

            # w, orth = self.gram_schmidt_2(w)
            w, orth  = np.linalg.qr(w)

            # exit()

            # w, orth  = np.linalg.qr(J)

            # print("")
            # print(J)
            # print(orth)

            # w, orth = np.linalg.qr(w)
            # print("w", w)
            # exit()

            if t_idx < skip_percent:
                continue

            orth = orth.diagonal()

            for e_idx, e in enumerate(exponents):
                log_orth = np.log2(abs(orth[e_idx]))
                exponents[e_idx] += log_orth
                # time_point_exponents[e_idx].append(exponents[e_idx] / ((t_idx + 1) *self.dt))


        for e_idx, e in enumerate(exponents):
            exponents[e_idx] = exponents[e_idx] / ( (len(time_points) - skip_percent) * self.dt)

        all_exponents_dict[sim_idx] = exponents
        all_timepoint_exponents_dict[sim_idx] = time_point_exponents


    def sprott_max_LE_est(self, model_refs, input_params, init_states):
        init_theta_states = [copy.deepcopy(state) for state in init_states]
        perturb_idx = 2

        # Peturb by init species theta_0
        for theta_state in init_theta_states:
            perturb_idx = np.random.choice(self.fit_species)
            # perturb_idx = 2
            theta_state[perturb_idx] += self.theta_0

        separation_coefficients = [self.theta_0 for c in range(len(model_refs))]

        curr_states = [copy.deepcopy(state) for state in init_states]
        curr_theta_states = [copy.deepcopy(state) for state in init_theta_states]

        dt = self.dt
        time_points = np.arange(self.t_0, self.t_end, dt)
        skip_percent = int(len(np.arange(self.t_0, self.t_end, dt)) * 0.1)

        # Simulate step
        for x in range(len(time_points)):
            if len(curr_states) < self.n_sims_batch:
                adjusted_sims_batch = len(curr_states)

            else:
                adjusted_sims_batch = self.n_sims_batch

            pop_obj = population_modules.Population(adjusted_sims_batch, self.t_0, dt*2,
                                                         self.dt, curr_states, input_params, model_refs,
                                                         self.fit_species, self.abs_tol, self.rel_tol)

            pop_obj_theta = population_modules.Population(adjusted_sims_batch, self.t_0, dt*2,
                                                         self.dt, curr_theta_states, input_params, model_refs,
                                                         self.fit_species, self.abs_tol, self.rel_tol)
            pop_obj.generate_particles()
            pop_obj_theta.generate_particles()

            pop_obj.simulate_particles()
            pop_obj_theta.simulate_particles()

            # Get step separations
            for sim_idx in range(len(model_refs)):
                err = self.pop_obj.get_particle_integration_error(sim_idx)
                if err:
                    separation_coefficients[sim_idx] = np.nan

                final_vals = pop_obj.get_particle_final_species_values(sim_idx)
                final_vals_theta = pop_obj_theta.get_particle_final_species_values(sim_idx)

                # Calculate separation
                theta_1 = 0
                species_idx = 0
                for y_a, y_b in zip(final_vals, final_vals_theta):
                    theta_1 += (y_a - y_b)**2

                    species_idx += 1

                    if species_idx > self.fit_species[-1]:
                        break

                theta_1 = theta_1**0.5


                species_idx = 0
                for y_a, y_b in zip(final_vals, final_vals_theta):
                    # Set state for next sim
                    curr_states[sim_idx][species_idx] = y_a

                    # Adjusted state of y_theta for next sim
                    # curr_theta_states[sim_idx][species_idx] = y_a + (self.theta_0 / theta_1) * (y_b - y_a)
                    # print(x, sim_idx, theta_1)
                    curr_theta_states[sim_idx][species_idx] = y_a + (self.theta_0 * (y_b - y_a)) / theta_1

                    species_idx += 1


                if x > skip_percent:
                    # Add separation coeff to total
                    separation_coefficients[sim_idx] += np.log2(abs(theta_1/self.theta_0))

        # Make into mean
        separation_coefficients = [sep / ((len(time_points) - skip_percent) * dt) for sep in separation_coefficients]

        return separation_coefficients, init_theta_states

    def get_RQA_distances(self, model_refs, init_states, pop_num, batch_num, plot_RQA=False):
        distances = []

        # Iterate all particles in batch
        for sim_idx, m_ref in enumerate(model_refs):
            err = self.pop_obj.get_particle_integration_error(sim_idx)
            if err:
                distances.append([0.0, 0.0, 0.0])
                continue

            state_list = self.pop_obj.get_particle_state_list(sim_idx)
            time_points = self.pop_obj.get_timepoints_list()

            try:
                state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))


            except(ValueError):
                # print(len(state_list)/len(init_states[sim_idx]))
                time_points = range(int(len(state_list) / len(init_states[sim_idx])))
                state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))
            

            skip_percent = int(len(time_points) * 0.5)

            state_list = state_list[-skip_percent:]

            time_points = time_points[-skip_percent:]

            # state_list = state_list[5000:]
            # time_points = time_points[5000:]


            # state_list = state_list[40000:]
            # time_points = time_points[40000:]

            #  Settings for the recurrence plot
            RR = 0.05
            # Distance metric in phase space ->
            # Possible choices ("manhattan","euclidean","supremum")
            METRIC = "euclidean"

            #  Generate a recurrence plot object with fixed recurrence rate RR
            rp = RecurrencePlot(state_list, metric=METRIC, 
                normalize=True, recurrence_rate=RR, silence_level=3, sparse_rqa=False)

            recurrence_rate = rp.recurrence_rate()
            recurrence_prob =  rp.recurrence_probability()

            np.fill_diagonal(rp.R, 0.0)

            det = rp.determinism(l_min=2)
            entropy = rp.diag_entropy(l_min=2)
            laminarity = rp.laminarity(v_min=2)
            max_diag = rp.max_diaglength()
            div = 1/ max_diag
            distances.append([div, entropy])

            if plot_RQA:
                out_path_template = self.out_dir + "Population_" + str(self.population_number) + "/simulation_plots/RQA" + "_batch_" + str(self.batch_num) + "_idx_#SIM_IDX#_sep.pdf"
                out_path = out_path_template.replace('#SIM_IDX#', str(sim_idx))
                fig = plt.figure()
                plt.matshow(rp.recurrence_matrix())
                plt.savefig(out_path, dpi=150)
                plt.close()
                plt.clf()

        return distances

    def gram_schmidt_2(self, A):
        # https://gist.github.com/iizukak/1287876
        """Orthogonalize a set of vectors stored as the columns of matrix A."""
        # Get the number of vectors.
        n = A.shape[1]
        orth = [0 for x in range(n)]
        for j in range(n):
            # To orthogonalize the vector in column j with respect to the
            # previous vectors, subtract from it its projection onto
            # each of the previous vectors.
            for k in range(j):
                A[:, j] -= np.dot(A[:, k], A[:, j]) * A[:, j]

            orth[j] = np.linalg.norm(A[:, j])
            A[:, j] = A[:, j] / np.linalg.norm(A[:, j])

        return A, orth

    def gram_schmidt(self, mat):
        basis = []
        n_columns = np.shape(mat)[1]
        orth = [0 for x in range(n_columns)]
        for idx in range(n_columns):
            v = mat[idx]
            w = v - np.sum( np.dot(v,b)*b  for b in basis )
            
            orth[idx] = np.linalg.norm(w)

            basis.append(w/np.linalg.norm(w))

        return np.array(basis), orth

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
            det = rp.determinism(l_min=2)
            entropy = rp.diag_entropy(l_min=2)
            rp_mat = rp.recurrence_matrix()
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
            err = self.pop_obj.get_particle_integration_error(sim_idx)
            print(sim_idx, err)
            if err:
                continue

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

    def plot_time_LE_data(self, timepoint_exponents_data, model_refs, output_dir):
        sns.set_context("talk")
        sns.set_style("white")


        colours = ['#e30b17', '#1d71b8', '#73ba65', '#ac7cb5', '#706f6f']
        for sim_idx in range(len(model_refs)):
            print(np.shape(timepoint_exponents_data))
            sim_tp_exponents = timepoint_exponents_data[sim_idx]
            t_idx = 0
            per_time_exponents = []

            width_inches = 95*4 / 25.4
            height_inches = (51*4 / 25.4) / 2

            fig, ax = plt.subplots(figsize=(width_inches, height_inches))

            for idx, exponent_list in enumerate(sim_tp_exponents):
                exponent_list = exponent_list

                sns.lineplot(range(len(exponent_list)), exponent_list, color=colours[idx])
                plt.plot(range(len(exponent_list)), exponent_list)
                print(exponent_list[-1])


            ax.set(ylabel='Average LE')
            ax.set(xlabel='Time index')
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.spines["left"].set_alpha(0.5)
            ax.spines["bottom"].set_alpha(0.5)
            fig.tight_layout()

            output_path = output_dir + "average_LE_idx_" + str(sim_idx) + ".pdf"

            plt.savefig(output_path, dpi=500)

    def calculate_LE_spectra(self, model_refs, init_states, pop_num, batch_num):
        all_particle_exponents = []

        #Iterate all simulated particles
        for sim_idx in range(len(model_refs)):
            num_species = len(init_states[sim_idx])

            exponents = [0 for i in range(num_species)]

            w = np.eye(num_species)

            # Get trajectory
            state_list = self.pop_obj.get_particle_state_list(sim_idx)
            time_points = self.pop_obj.get_timepoints_list()

            print(np.shape(state_list))

            try:
                state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))


            except(ValueError):
                # print(len(state_list)/len(init_states[sim_idx]))
                time_points = range(int(len(state_list) / len(init_states[sim_idx])))
                state_list = np.reshape(state_list, (len(time_points), len(init_states[sim_idx])))
            
            state_list = state_list[self.n_skip:]
            time_points = time_points[self.n_skip:]

            dt = time_points[1] - time_points[0]

            for idx, t in enumerate(time_points):
                jac_inputs = list(state_list[idx])
                particle_jacobian = np.reshape(self.pop_obj.get_particle_jacobian(jac_inputs, sim_idx), [num_species, num_species])
                # J = np.eye(num_species) + np.dot(particle_jacobian, dt)
                J = particle_jacobian
                # J = np.eye(num_species) + particle_jacobian

                w = np.dot(J, w)
                
                for e_idx, e in enumerate(exponents):
                    exponents[e_idx] += (np.linalg.norm(w[:, e_idx]))

                # Orthonormal matrix w 
                w = self.gram_schmidt(w)

                # # # # # # Re normalise into orthogonal vectors
                # for norm_idx in range(num_species):
                #     w[:, norm_idx] = w[:, norm_idx] / np.linalg.norm(w[:, norm_idx])


            # Exponents as average of time
            for e_idx, e in enumerate(exponents):
                exponents[e_idx] = np.log(exponents[e_idx]) / (len(time_points) * self.dt)

            # Eigenvalues of longtime product matrix
            w_eig, _ = np.linalg.eig(w)

            print(exponents)
            all_particle_exponents.append(exponents)

        return all_particle_exponents

    def make_four_strain_interaction_masks(self):
        n, m = 4, 4

        all_mats = product([1, 0], repeat=n*m)
        all_mats = np.reshape(list(all_mats), (-1, n, m))

        # Keep only matricies with valid trace
        valid_trace_mats = []
        for x in all_mats:
            if np.trace(x) != 4:
                continue

            else:
                valid_trace_mats.append(x)

        non_identical_mats = []
        adj_mat_sums = []

        for idx, mat in enumerate(valid_trace_mats):
            # print(idx)
            exists = False
            candidate_idxs = np.argwhere(np.array(adj_mat_sums) == np.sum(mat))

            for i, _ in enumerate(range(n)):
                for j, _ in enumerate(range(n)):
                    if i == j :
                        continue
                    permuted_matrix = np.copy(mat)

                    # Swap rows
                    permuted_matrix[i] = np.copy(mat[j])
                    permuted_matrix[j] = np.copy(mat[i])

                    # Swap columns
                    col_i = np.copy(permuted_matrix[:, i])
                    col_j = np.copy(permuted_matrix[:, j])
                    permuted_matrix[:, i] = col_j
                    permuted_matrix[:, j] = col_i
                    # permuted_matrix = np.reshape(permuted_matrix, (-1))


                    for k in candidate_idxs:
                        if np.array_equal(non_identical_mats[k[0]], permuted_matrix):
                            exists = True
                            break

                    if exists:
                        break

                if exists:
                    break

            if not exists:
                non_identical_mats.append(mat)
                adj_mat_sums.append(np.sum(mat))

        return non_identical_mats


    def apply_four_species_param_mask(self, mask_list, particle):
        interaction_params = []

        mask_example = np.eye(4)

        for idx, id in enumerate(sorted(particle.curr_model._params_prior, key=str.lower)):
            if "m_" in id:
                interaction_params.append(particle.curr_params[idx])

        random_mask_idx = np.random.choice(range(len(mask_list)))
        interaction_params = np.array(interaction_params).reshape(4, 4)
        masked_interaction_params  = interaction_params * mask_list[random_mask_idx]
        

        masked_interaction_params = np.reshape(masked_interaction_params, (-1))

        int_param_idx = 0
        for idx, id in enumerate(sorted(particle.curr_model._params_prior, key=str.lower)):
            if "m_" in id:
                particle.curr_params[idx] = masked_interaction_params[int_param_idx]
                int_param_idx += 1


    def run_chaos_separation_model_selection_ABC_SMC(self, alpha=0.5):
        # interactions_mask_list = self.make_four_strain_interaction_masks()
        # print("Numbe of masks: ", interactions_mask_list)

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
            LHC_idx = 0

            while self.population_accepted_count < self.population_size:
                print("")
                start_time = time.time()
                start_time_sampling = time.time()

                # 1. Sample from model space
                if self.population_number == 0:
                    # particles = self.model_space.sample_particles_from_prior(self.n_sims_batch)
                    # csv_path = '/home/behzad/Documents/AutoCD/output/sprott_small_SMC/experiment_analysis/all_strange_params.csv'
                    # particles = self.model_space.sample_particles_from_csv(self.n_sims_batch, 0, csv_path)
                    csv_path = '/home/behzad/Documents/AutoCD/output/gen_LV_four_chem_3_rej_1/experiment_analysis/all_strange_params.csv'
                    csv_path = '/home/behzad/Documents/AutoCD/output/gen_LV_four_chem_3_rej/sprott_reevaluation/experiment_analysis/all_strange_params.csv'
                    csv_path = '/home/behzad/Documents/AutoCD/output/three_species_7_chaos_rej/experiment_analysis/all_strange_params.csv'

                    particles = self.model_space.sample_particles_from_ordered_csv(self.n_sims_batch, self.population_total_simulations, csv_path)

                    self.theta_0 =1e-10

                    # for p in particles:
                    #     self.apply_four_species_param_mask(interactions_mask_list, p)

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
                
                # exit()

                init_states = [copy.deepcopy(p.curr_init_state) for p in particles]
                input_params = [copy.deepcopy(p.curr_params) for p in particles]

                model_refs = [copy.deepcopy(p.curr_model._model_ref) for p in particles]
                particle_weights = [p.prev_weight for p in particles]

                particle_models = [p.curr_model for p in particles]

                # alg_utils.rescale_parameters(input_params, init_states, particle_models)
                self.pop_obj = population_modules.Population(self.n_sims_batch, self.t_0, self.t_end,
                                                             self.dt, init_states, input_params, model_refs,
                                                             self.fit_species, self.abs_tol, self.rel_tol)

                self.pop_obj.generate_particles()
                start_time_sim = time.time()
                self.pop_obj.simulate_particles()


                if 0:
                    max_jobs_running = 4

                    start_time_sim = time.time()
                    manager = multiprocessing.Manager()                
                    all_exponents_dict = manager.dict()
                    all_timepoint_exponents_dict = manager.dict()

                    jobs = []
                    jobs_running = 0
                    for sim_idx in range(len(model_refs)):
                        p = multiprocessing.Process(target=self.wolf_LE_est, args=(model_refs, input_params, init_states, all_exponents_dict, all_timepoint_exponents_dict, sim_idx))
                        jobs.append(p)
                        p.start()

                        jobs_running += 1

                        if jobs_running >= max_jobs_running:
                            while jobs_running >= max_jobs_running:
                                jobs_running = 0
                                # print(len(jobs))

                                finished_jobs = []
                                for p in jobs:
                                    if p.is_alive():
                                        jobs_running += 1
                                    else:
                                        # pass
                                        p.join()
                                        jobs.remove(p)

                    for proc in jobs:
                        proc.join()
                        # jobs.remove(proc)

                    end_time_sim = time.time()
                    print("LE time elapsed: ", end_time_sim - start_time_sim)
                    all_exponents_data = [all_exponents_dict[sim_idx] for sim_idx in range(len(model_refs))]
                    all_exponents_timepoint_data = [all_timepoint_exponents_dict[sim_idx] for sim_idx in range(len(model_refs))]

                    del jobs
                # self.plot_time_LE_data(all_exponents_timepoint_data, model_refs, folder_name)

                self.pop_obj.calculate_particle_distances(2)
                self.pop_obj.accumulate_distances()

                batch_distances = self.pop_obj.get_flattened_distances_list()
                batch_distances = np.reshape(batch_distances,
                                 (self.n_sims_batch, len(self.fit_species)))
                
                # RQA_dist = self.get_RQA_distances(model_refs, init_states, self.population_number, self.batch_num, plot_RQA=False)
                sprott_exponents = self.sprott_max_LE_est(model_refs, input_params, init_states)[0]

                data = []
                inverse_max_exponents = []
                for sim_idx, e in enumerate(sprott_exponents):
                    inverse_max_exponents.append([1/(e+1)])
                    data.append([e] + list(batch_distances[sim_idx]))
                    
                batch_part_judgements = alg_utils.check_distances_LE(inverse_max_exponents, epsilon_array=self.current_epsilon)

                error_messages = [self.pop_obj.get_particle_integration_error(sim_idx) for sim_idx, _ in enumerate(model_refs)]
                # print(error_messages)
                # batch_part_judgements = [True if i == '' else False for i in error_messages]

                all_errors = True

                non_error_sims = 0
                for msg in error_messages:
                    if msg == '':
                        non_error_sims += 1
                        all_errors = False
                
                if not all_errors:
                    # pass
                    # self.plot_lorenz(folder_name, self.population_number, self.batch_num, init_states, model_refs)
                    # if sum(batch_part_judgements) != 0:
                    self.plot_all_particles(folder_name, self.population_number, self.batch_num, init_states, 
                        model_refs, batch_part_judgements, inverse_max_exponents)
                    # self.plot_four_species(folder_name, self.population_number, self.batch_num, init_states, 
                    #     model_refs, batch_part_judgements, inverse_max_exponents)
                # exit()
                self.batch_num += 1


                n_spaces = self.population_size - self.population_accepted_count
                # # If we have more accepted than spaces, set some to false
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

                self.write_particle_data(folder_name, model_refs, self.batch_num, self.population_number, data, batch_part_judgements, particle_weights)
                self.write_particle_chaos_params(sim_params_folder, self.batch_num, self.population_number, particles,
                              input_params, init_states, batch_part_judgements)

                end_time_sim = time.time()
                print("Batch time elapsed: ", end_time_sim - start_time_sim)

                del self.pop_obj

                self.population_accepted_particle_distances += [part_d for part_d, judgement in
                                                zip(inverse_max_exponents, batch_part_judgements)
                                                if judgement]

                self.population_accepted_particles = self.population_accepted_particles + accepted_particles
                self.population_total_simulations += len(particles)
                self.population_accepted_count += sum(batch_part_judgements)

                # # print("Writing data")
                # self.write_epsilon(folder_name, self.current_epsilon)


                # start_time_write_distance = time.time()

                # if self.final_epsilon == self.current_epsilon:
                #     # self.plot_accepted_particles(folder_name, self.population_number, self.batch_num, batch_part_judgements, init_states, model_refs)
                #     self.write_particle_distances(folder_name, model_refs, self.batch_num, self.population_number,
                #                                   batch_part_judgements, combined_distances, only_accepted=False)

                # end_time_write_distance = time.time()

                # # self.population_total_simulations += len(model_refs)

                # self.model_space.update_population_sample_data_v2(model_refs, batch_part_judgements)

                print("Population: ", self.population_number, "Accepted particles: ", self.population_accepted_count,
                      "Total simulations: ", self.population_total_simulations)
            
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
                self.current_epsilon = alg_utils.update_epsilon_LE(self.current_epsilon, self.final_epsilon,
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
