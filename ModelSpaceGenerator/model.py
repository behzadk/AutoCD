from . import equation_builder
import sympy
from sympy.printing.cxxcode import cxxcode
import numpy as np
import pandas as pd
import csv
from . import species
import os
from . import utils
from collections import OrderedDict


class Model:
    def __init__(self, model_idx, strain_list):
        self.idx = model_idx
        self.strains = strain_list

        self.strain_ids = self.get_strain_species()
        self.AHL_ids = self.get_AHL_species()
        self.microcin_ids = self.get_microcin_species()
        self.substrate_ids = self.get_substrate_species()
        self.antitoxin_ids = self.get_antitoxin_species()
        self.immunity_ids = self.get_immunity_species()
        self.toxin_ids = self.get_toxin_species()

        self.all_ids = self.microcin_ids + self.AHL_ids + self.substrate_ids + self.strain_ids + self.immunity_ids + self.toxin_ids

        self.diff_eqs = OrderedDict()
        self.symbolic_equations = None
        self.jac = None

        self.params_list = []
        self.species_list = []

        # Model at this stage can be split into terms
        # Substrate rhs
        # Cell number rhs
        # AHL rhs

    # edges show interactions between species from j (column) to i (row)
    def generate_adjacency_matrix(self, max_sub, max_AHL, max_mic, max_strains, max_antitoxin, max_immunity, max_toxin):
        # Cell1, Cell2, M1, M2, AHL1, AHL2, V1, V2
        total_species = max_AHL + max_mic + max_strains + max_sub + max_antitoxin + max_immunity + max_toxin
        # total_species = len(self.AHL_ids) + len(self.microcin_ids) + len(self.strain_ids)
        adjacency_matrix = np.zeros([total_species, total_species])

        strain_init_idx = 0
        substrate_init_idx = strain_init_idx + max_strains
        microcin_init_idx = substrate_init_idx + max_sub
        AHL_init_idx = microcin_init_idx + max_mic
        antitoxin_init_idx = AHL_init_idx + max_AHL
        immunity_init_idx = antitoxin_init_idx + max_antitoxin
        toxin_init_idx = immunity_init_idx + max_immunity

        all_microcin_objects = []
        all_microcin_ids = []
        all_AHLs = []
        all_substrates = []
        all_antitoxin_objects = []
        all_antitoxin_ids = []

        all_toxin_ids = []
        all_toxin_objects = []

        all_immunity_ids = []
        all_immunity_objects = []

        # Collect species objects
        for strain in self.strains:
            for microcin in strain.microcins:
                if microcin not in all_microcin_objects:
                    all_microcin_objects.append(microcin)

                if microcin.id not in all_microcin_ids:
                    all_microcin_ids.append(microcin.id)

            for AHL in strain.AHLs:
                if AHL not in all_AHLs:
                    all_AHLs.append(AHL)

            for s_dependence in strain.substrate_dependences:
                if s_dependence not in all_substrates:
                    all_substrates.append(s_dependence)

            for s_production in strain.substrate_production:
                if s_production not in all_substrates:
                    all_substrates.append(s_production)

            for v in strain.antitoxins:
                if v not in all_antitoxin_objects:
                    all_antitoxin_objects.append(v)

                if v.id not in all_antitoxin_ids:
                    all_antitoxin_ids.append(v.id)

            for t in strain.toxins:
                if t not in all_toxin_objects:
                    all_toxin_objects.append(t)

                if t.id not in all_toxin_ids:
                    all_toxin_ids.append(t.id)

            for i in strain.immunity:
                if i not in all_immunity_objects:
                    all_immunity_objects.append(i)

                if i.id not in all_immunity_ids:
                    all_immunity_ids.append(i.id)

        # Fill strain sensitivities to microcin sensitivity is a negative interaction from
        # microcin to strain. (i strain row, j mic col )
        for idx_strain, strain in enumerate(self.strains):
            for sens in strain.sensitivities:
                for idx_mic, mic_id in enumerate(all_microcin_ids):
                    if sens == mic_id:
                        from_node = idx_mic + microcin_init_idx
                        to_node = idx_strain + strain_init_idx

                        adjacency_matrix[to_node, from_node] = -1

        # Fill strain substrate dependency and production
        for idx_strain, strain in enumerate(self.strains):
            # Dependency
            for sub_dependence in strain.substrate_dependences:
                sub_indexs = [all_substrates.index(x) for i, x in enumerate(all_substrates) if x == sub_dependence]
                for s_idx in sub_indexs:
                    from_node = s_idx + substrate_init_idx
                    to_node = idx_strain + strain_init_idx

                    adjacency_matrix[to_node, from_node] = 1

        # Fill strain toxin sensitivity
        for idx_strain, strain in enumerate(self.strains):

            # Production
            for sub_production in strain.substrate_production:
                sub_indexs = [all_substrates.index(x) for i, x in enumerate(all_substrates) if x == sub_production]
                for s_idx in sub_indexs:
                    from_node = idx_strain + strain_init_idx
                    to_node = s_idx + substrate_init_idx

                    adjacency_matrix[to_node, from_node] = 1

        # Fill strain production of microcin, AHL, substrate and antitoxin
        for idx_strain, strain in enumerate(self.strains):
            # Microcin production
            for idx_strain_mic, strain_mic in enumerate(strain.microcins):
                mic_produced_idx = [all_microcin_ids.index(x.id) for i, x in enumerate(all_microcin_objects) if
                                    x == strain_mic]

                for idx_mic in mic_produced_idx:
                    from_node = idx_strain + strain_init_idx
                    to_node = idx_mic + microcin_init_idx

                    adjacency_matrix[to_node, from_node] = 1

            # Antitoxin production
            for strain_antitoxin in strain.antitoxins:
                v_produced_idx = [all_antitoxin_ids.index(x.id) for i, x in enumerate(all_antitoxin_objects) if
                                  x == strain_antitoxin]

                for v_idx in v_produced_idx:
                    from_node = idx_strain + strain_init_idx
                    to_node = v_idx + antitoxin_init_idx
                    adjacency_matrix[to_node, from_node] = 1

            # Immuity production
            for strain_immunity in strain.immunity:
                i_produced_idx = [all_immunity_ids.index(x.id) for i, x in enumerate(all_immunity_objects) if
                                  x == strain_immunity]

                for i_idx in i_produced_idx:
                    from_node = idx_strain + strain_init_idx
                    to_node = i_idx + immunity_init_idx
                    adjacency_matrix[to_node, from_node] = 1

            # Toxin production and sensitivity
            for strain_toxin in strain.toxins:
                t_produced_idx = [all_toxin_ids.index(x.id) for i, x in enumerate(all_toxin_objects) if
                                  x == strain_toxin]

                for t_idx in t_produced_idx:
                    from_node = idx_strain + strain_init_idx
                    to_node = t_idx + toxin_init_idx
                    adjacency_matrix[to_node, from_node] = 1

                for t_idx in t_produced_idx:
                    from_node = t_idx + toxin_init_idx
                    to_node = idx_strain + strain_init_idx
                    adjacency_matrix[to_node, from_node] = -1

            # Antitoxin inhibition of toxin
            for v_idx, v in enumerate(all_antitoxin_ids):
                from_node = v_idx + antitoxin_init_idx

                for toxin_idx, toxin in enumerate(all_toxin_ids):
                    if v.split('_')[-1] == toxin:
                        to_node = toxin_idx + toxin_init_idx
                        adjacency_matrix[to_node, from_node] = -1

            # Immunity inhibition of mic
            for i_idx, i in enumerate(all_immunity_ids):
                from_node = i_idx + immunity_init_idx

                for mic_idx, mic in enumerate(all_microcin_ids):
                    if i.split('_')[-1] == mic:
                        to_node = mic_idx + microcin_init_idx
                        adjacency_matrix[to_node, from_node] = -1

            # AHL production
            for idx_strain_AHL, strain_AHL in enumerate(strain.AHLs):
                AHL_produced_idx = [i for i, x in enumerate(all_AHLs) if x == strain_AHL]
                for idx_AHL in AHL_produced_idx:
                    from_node = idx_strain + strain_init_idx
                    to_node = idx_AHL + AHL_init_idx
                    adjacency_matrix[to_node, from_node] = 1

            # AHL mic interactions. from AHL to mic
            for mic in all_microcin_objects:
                idx_mic = all_microcin_ids.index(mic.id)

                # Repressors
                try:  #
                    if mic.AHL_repressors is np.nan:
                        continue

                    repressor = mic.AHL_repressors[0]
                    AHL_repressor_idx = [i for i, x in enumerate(all_AHLs) if x == repressor]

                    for idx_AHL in AHL_repressor_idx:
                        from_node = idx_AHL + AHL_init_idx
                        to_node = idx_mic + microcin_init_idx

                        adjacency_matrix[to_node, from_node] = -1

                except(IndexError):
                    pass

                # Inducers
                try:
                    if mic.AHL_repressors is np.nan:
                        continue

                    inducer = mic.AHL_inducers[0]
                    AHL_inducer_idx = [i for i, x in enumerate(all_AHLs) if x == inducer]
                    for idx_AHL in AHL_inducer_idx:
                        from_node = idx_AHL + AHL_init_idx
                        to_node = idx_mic + microcin_init_idx

                        adjacency_matrix[to_node, from_node] = 1

                except(IndexError):
                    pass

            # AHL antitoxin interactions from AHL to V
            for v in all_antitoxin_objects:
                v_idx = all_antitoxin_ids.index(v.id)

                # Repressors
                try:
                    if v.AHL_repressors is np.nan:
                        continue

                    repressor = v.AHL_repressors[0]
                    AHL_repressor_idx = [i for i, x in enumerate(all_AHLs) if x == repressor]

                    for idx_AHL in AHL_repressor_idx:
                        from_node = idx_AHL + AHL_init_idx
                        to_node = v_idx + antitoxin_init_idx

                        adjacency_matrix[to_node, from_node] = -1

                except(IndexError):
                    pass

                # Inducers
                try:
                    if v.AHL_repressors is np.nan:
                        continue

                    inducer = v.AHL_inducers[0]
                    AHL_inducer_idx = [i for i, x in enumerate(all_AHLs) if x == inducer]
                    for idx_AHL in AHL_inducer_idx:
                        from_node = idx_AHL + AHL_init_idx
                        to_node = v_idx + antitoxin_init_idx

                        adjacency_matrix[to_node, from_node] = 1

                except(IndexError):
                    pass

            # AHL toxin interactions from AHL to V
            for t in all_toxin_objects:
                t_idx = all_toxin_ids.index(t.id)

                # Repressors
                try:
                    if t.AHL_repressors is np.nan:
                        continue

                    repressor = t.AHL_repressors[0]
                    AHL_repressor_idx = [i for i, x in enumerate(all_AHLs) if x == repressor]

                    for idx_AHL in AHL_repressor_idx:
                        from_node = idx_AHL + AHL_init_idx
                        to_node = t_idx + toxin_init_idx

                        adjacency_matrix[to_node, from_node] = -1

                except(IndexError):
                    pass

                # Inducers
                try:
                    if t.AHL_repressors is np.nan:
                        continue

                    inducer = t.AHL_inducers[0]
                    AHL_inducer_idx = [i for i, x in enumerate(all_AHLs) if x == inducer]
                    for idx_AHL in AHL_inducer_idx:
                        from_node = idx_AHL + AHL_init_idx
                        to_node = t_idx + toxin_init_idx

                        adjacency_matrix[to_node, from_node] = 1

                except(IndexError):
                    pass

            # AHL immunity interactions from AHL to I
            for i in all_immunity_objects:
                i_idx = all_immunity_ids.index(i.id)

                # Repressors
                try:
                    if i.AHL_repressors is np.nan:
                        continue

                    repressor = i.AHL_repressors[0]
                    AHL_repressor_idx = [i for i, x in enumerate(all_AHLs) if x == repressor]

                    for idx_AHL in AHL_repressor_idx:
                        from_node = idx_AHL + AHL_init_idx
                        to_node = i_idx + immunity_init_idx

                        adjacency_matrix[to_node, from_node] = -1

                except(IndexError):
                    pass

                # Inducers
                try:
                    if i.AHL_repressors is np.nan:
                        continue

                    inducer = i.AHL_inducers[0]
                    AHL_inducer_idx = [i for i, x in enumerate(all_AHLs) if x == inducer]
                    for idx_AHL in AHL_inducer_idx:
                        from_node = idx_AHL + AHL_init_idx
                        to_node = i_idx + immunity_init_idx

                        adjacency_matrix[to_node, from_node] = 1

                except(IndexError):
                    pass

        self.adjacency_matrix = adjacency_matrix

    def get_strain_species(self):
        strain_id_list = []
        for strain in self.strains:
            strain_id_list.append(strain.id)

        return list(set(strain_id_list))

    def get_microcin_species(self):
        microcin_id_list = []
        for strain in self.strains:
            strain_microcins = strain.microcins
            for m in strain_microcins:
                microcin_id_list.append(m.id)

        return list(set(microcin_id_list))

    def get_AHL_species(self):
        AHL_id_list = []
        for strain in self.strains:
            strain_AHLs = strain.AHLs
            for a in strain_AHLs:
                AHL_id_list.append(a.id)

        AHL_id_list = list(set(AHL_id_list))

        sorting_func = lambda x: int(x)

        return list(set(AHL_id_list))

    def get_substrate_species(self):
        substrate_id_list = []
        for strain in self.strains:
            strain_susbtrates = strain.substrate_dependences
            for s in strain_susbtrates:
                substrate_id_list.append(s.id)

        for strain in self.strains:
            strain_susbtrates = strain.substrate_production
            for s in strain_susbtrates:
                substrate_id_list.append(s.id)

        return list(set(substrate_id_list))

    def get_antitoxin_species(self):
        antitoxin_list = []
        for strain in self.strains:
            strain_antitoxins = strain.antitoxins
            for v in strain_antitoxins:
                antitoxin_list.append(v.id)

        return list(set(antitoxin_list))

    def get_immunity_species(self):
        immunity_list = []
        for strain in self.strains:
            strain_immunities = strain.immunity
            for i in strain_immunities:
                immunity_list.append(i.id)

        return list(set(immunity_list))

    def get_toxin_species(self):
        toxin_list = []
        for strain in self.strains:
            strain_toxins = strain.toxins
            for t in strain_toxins:
                toxin_list.append(t.id)

        return list(set(toxin_list))

    def is_legal(self):
        required_microcin = []
        required_AHL = []
        required_sub = []

        # Legal requirements
        for s in self.strains:
            if len(s.substrate_dependences) == 0:
                return False

            required_sub += s.substrate_dependences

            # Load AHLs which have an activity on microcin or antitoxin expression
            for m in s.microcins:
                if m.AHL_inducers is not np.nan:
                    for a in m.AHL_inducers:
                        required_AHL.append(a)
                if m.AHL_repressors is not np.nan:
                    for a in m.AHL_repressors:
                        required_AHL.append(a)

            for v in s.antitoxins:
                if v.AHL_inducers is not np.nan:
                    for a in v.AHL_inducers:
                        required_AHL.append(a)
                if v.AHL_repressors is not np.nan:
                    for a in v.AHL_repressors:
                        required_AHL.append(a)

            for i in s.immunity:
                if i.AHL_inducers is not np.nan:
                    for a in i.AHL_inducers:
                        required_AHL.append(a)
                if i.AHL_repressors is not np.nan:
                    for a in i.AHL_repressors:
                        required_AHL.append(a)

            for t in s.toxins:
                if t.AHL_inducers is not np.nan:
                    for a in t.AHL_inducers:
                        required_AHL.append(a)
                if t.AHL_repressors is not np.nan:
                    for a in t.AHL_repressors:
                        required_AHL.append(a)

            for m_sens in s.sensitivities:
                required_microcin += [m_sens]

        for a in required_AHL:
            if a.id not in self.AHL_ids:
                return False

        for m in required_microcin:
            if m not in self.microcin_ids:
                return False

        # Remove redundancies
        # Remove models where AHL has no action
        required_AHL_ids = [a.id for a in required_AHL]
        for a_expressed in self.AHL_ids:
            if a_expressed not in required_AHL_ids:
                return False

        # Expression of microcin which no strain is sensitive to
        system_sensitivities = []
        for s in self.strains:
            system_sensitivities = system_sensitivities + s.sensitivities

        for m in self.microcin_ids:
            if m not in system_sensitivities:
                return False

        # Remove models where a produced substrate is not consumed
        for strain in self.strains:
            for sub in strain.substrate_production:
                if sub not in required_sub:
                    return False

        # Remove models where a produced substrate is also a dependency
        for strain in self.strains:
            for sub in strain.substrate_production:
                if sub in strain.substrate_dependences:
                    return False

        # Remove models where a substrate dependency does not exist
        # for strain in self.strains:
        #     for sub in strain.substrate_dependences:
        #         if sub not in self.substrate_ids:
        #             return False

        for strain in self.strains:
            if len(strain.substrate_dependences) == 0:
                return False

        # Remove models where antitoxin has no cognate toxin
        for strain in self.strains:
            for v in strain.antitoxins:
                if v.id.split('_')[-1] not in [t.id for t in strain.toxins]:
                    return False

        # Remove models where immunity has no cognate microcin
        for strain in self.strains:
            for i in strain.immunity:
                if i.id.split('_')[-1] not in [i.id for i in strain.immunity]:
                    return False

        # Remove models 

        return True

    def build_equations(self):

        # For each strain
        for n in self.strains:
            dN_dt = equation_builder.gen_strain_growth_diff(n.id, self.strains)
            self.diff_eqs.update(dN_dt)

        # For each substrate
        for s in self.substrate_ids:
            dS_dt = equation_builder.gen_diff_eq_substrate(s, self.strains)
            self.diff_eqs.update(dS_dt)

        # For each microcin
        for b in self.microcin_ids:
            dB_dt = equation_builder.gen_microcin_diff_eq(b, self.strains)
            self.diff_eqs.update(dB_dt)

        # For each AHL
        for a in self.AHL_ids:
            dA_dt = equation_builder.gen_AHL_diff_eq(a, self.strains)
            self.diff_eqs.update(dA_dt)

        # For each antitoxin
        for v in self.antitoxin_ids:
            dV_dt = equation_builder.gen_diff_eq_antitoxin(v, self.strains)
            self.diff_eqs.update(dV_dt)

        # For each immunity
        for i in self.immunity_ids:
            dI_dt = equation_builder.gen_diff_eq_immunity(i, self.strains)
            self.diff_eqs.update(dI_dt)

        # For each immunity
        for t in self.toxin_ids:
            dT_dt = equation_builder.gen_toxin_diff_eq(t, self.strains)
            self.diff_eqs.update(dT_dt)

    def build_jacobian(self):
        species_names = list(self.diff_eqs.keys())
        order = sympy.symbols(species_names)
        J = self.symbolic_equations.jacobian(order)

        self.jac = J

    def build_symbolic_equations(self):
        species_names = list(self.diff_eqs.keys())
        # print(list(self.diff_eqs.keys()))
        order = sympy.symbols(species_names)
        zeros_list = [0 for i in range(len(order))]
        symbolic_equations = sympy.Matrix(zeros_list)

        for idx, eq_key in enumerate(species_names):
            symbolic_equations[idx] = sympy.sympify(self.diff_eqs[eq_key], locals=locals())
            # print(eq_key)

        self.symbolic_equations = symbolic_equations

    def extract_species(self):
        self.species_list = list(self.diff_eqs.keys())

    def extract_params(self):
        all_params = []

        for eq in self.symbolic_equations:
            free_symbols = eq.free_symbols
            for symbol in free_symbols:
                if str(symbol) not in self.species_list:
                    all_params.append(str(symbol))

        all_params = sorted(list(set(all_params)), key=str.lower)  # Alphabetical order!

        self.params_list = all_params

    def config_data(self):
        # N, S, B, A
        for N in self.strains:
            print(N.id)
            for S in N.substrate_dependences:
                print(S.id)

            for B in N.microcins:
                print(B.config_idx)

            for A in N.AHLs:
                print(A.id)

    def write_adj_matrix(self, output_dir, mic_ids, AHL_ids, strain_ids, substrate_ids, antitoxin_ids, immunity_ids,
                         toxin_ids):
        adj_mat_dir = output_dir + "adj_matricies/"
        utils.make_folder(adj_mat_dir)

        adj_mat_path = adj_mat_dir + 'model_' + str(self.idx) + '_adj_mat.csv'

        new_mic_ids = []
        for idx, i in enumerate(mic_ids):
            new_mic_ids.append("B_" + i)

        new_AHL_ids = []
        for idx, i in enumerate(AHL_ids):
            new_AHL_ids.append("A_" + i)

        new_strain_ids = []
        for idx, i in enumerate(strain_ids):
            new_strain_ids.append("N_" + i)

        new_substrate_ids = []
        for idx, i in enumerate(substrate_ids):
            new_substrate_ids.append("S_" + i)

        new_antitoxin_ids = []
        for idx, i in enumerate(antitoxin_ids):
            new_antitoxin_ids.append("V_" + i)

        new_immunity_ids = []
        for idx, i in enumerate(immunity_ids):
            new_immunity_ids.append("I_" + i)

        new_toxin_ids = []
        for idx, i in enumerate(toxin_ids):
            new_toxin_ids.append("T_" + i)

        adj_matrix = self.adjacency_matrix

        with open(adj_mat_path, 'w') as f:
            w = csv.writer(f)
            adj_mat_species = new_strain_ids + new_substrate_ids + new_mic_ids + new_AHL_ids + new_antitoxin_ids + new_immunity_ids + new_toxin_ids
            header = [None] + adj_mat_species
            w.writerow(header)

            for row_idx in range(len(adj_matrix)):
                w.writerow([adj_mat_species[row_idx]] + adj_matrix[row_idx].tolist())

    ##
    # Writes upper and lower boundaries for uniform priors to a csv, separately for
    # parameters and species using a csv containing default parameters.
    ##
    def write_prior_parameter_dict(self, default_params_path, output_dir):
        sim_inputs_dir = output_dir + "input_files/"
        utils.make_folder(sim_inputs_dir)

        sim_params_path = sim_inputs_dir + 'params_' + str(self.idx) + ".csv"

        model_parameters = self.params_list
        default_params = pd.read_csv(default_params_path)
        default_params.sort_values('parameter', inplace=True)

        prior_dict = OrderedDict()
        for idx, row in default_params.iterrows():
            p = row['parameter']

            # Adds param to dict if not linked to a particular species
            if p in model_parameters:
                prior_dict[p] = [row['lower_bound'], row['upper_bound']]
                continue

            # Adds param to dict if it is linked to a particular species, identified by the id tag
            p = row['parameter'] + '_#ID#'

            for id in self.all_ids:
                param_species = p.replace('#ID#', id)
                if param_species in model_parameters:
                    prior_dict[param_species] = [row['lower_bound'], row['upper_bound']]
                    continue

        if len(model_parameters) != len(prior_dict):
            params_missing = [i for i in model_parameters if i not in list(prior_dict.keys())]
            raise RuntimeError('Mismatch in length of prior dict and model parameters.', 'Prior: ',
                               len(prior_dict), 'Params needed: ', len(model_parameters), params_missing)

        with open(sim_params_path, 'w') as f:  # Just use 'w' mode in 3.x
            w = csv.writer(f)
            for key, value in prior_dict.items():
                w.writerow([key, value[0], value[1]])

    def write_init_species_dict(self, default_init_species_path, output_dir):
        sim_inputs_dir = output_dir + "input_files/"
        utils.make_folder(sim_inputs_dir)

        sim_species_path = sim_inputs_dir + 'species_' + str(self.idx) + ".csv"

        model_species = self.species_list
        default_species = pd.read_csv(default_init_species_path)

        prior_dict = OrderedDict()
        for idx, row in default_species.iterrows():
            p = row['species']

            # Adds param to dict if not linked to a particular species
            if p in model_species:
                prior_dict[p] = [row['lower_bound'], row['upper_bound']]
                continue

            # Adds param to dict if it is linked to a particulr species, identified by the id tag
            p = row['species'] + '_#ID#'
            for id in self.all_ids:
                param_species = p.replace('#ID#', id)

                if param_species in model_species:
                    prior_dict[param_species] = [row['lower_bound'], row['upper_bound']]
                    continue

        if len(model_species) != len(prior_dict):
            species_missing = [i for i in model_species if i not in list(model_species.keys())]
            raise RuntimeError('Mismatch in length of prior dict and model species.', 'Prior: ',
                               len(prior_dict), 'Params needed: ', len(model_parameters), species_missing)

        with open(sim_species_path, 'w') as f:  # Just use 'w' mode in 3.x
            w = csv.writer(f)
            for key, value in prior_dict.items():
                w.writerow([key, value[0], value[1]])

    def write_python_equations(self, output_path):
        py_eqs_dir = output_path + "py_eqs_txt_files/"
        utils.make_folder(py_eqs_dir)

        model_py_eqs_path = py_eqs_dir + "model_" + str(self.idx) + "_eqs.py"

        with open(model_py_eqs_path, 'w') as txt_file:  # Just use 'w' mode in 3.x
            for eq in self.diff_eqs:
                res = "\t" + "d" + eq + " = " + str(self.diff_eqs[eq]) + "\n\n"
                res = res.replace("^", "**")
                txt_file.write(res)
