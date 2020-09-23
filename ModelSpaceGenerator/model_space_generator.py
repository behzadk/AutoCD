from . import species
import pandas as pd
import itertools
from . import model
import numpy as np
from . import utils
from .cpp_output import Cpp_source_output
from .cpp_output import Cpp_header_output


def generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids, antitoxin_ids,
                                 immunity_ids, toxin_ids, output_dir):
    utils.make_folder(output_dir)

    for idx, m in enumerate(model_list):
        m.write_adj_matrix(output_dir, microcin_ids, AHL_ids, strain_ids, substrate_ids, antitoxin_ids, immunity_ids,
                           toxin_ids)


def generate_simulation_files(model_list, params_path, init_species_path, output_dir):
    utils.make_folder(output_dir)

    for idx, m in enumerate(model_list):
        m.build_equations()
        m.build_symbolic_equations()
        m.build_jacobian()
        m.extract_species()
        m.extract_params()
        m.write_python_equations(output_dir)

        m.write_prior_parameter_dict(params_path, output_dir)
        m.write_init_species_dict(init_species_path, output_dir)

    cpp_out = Cpp_source_output(model_list)
    cpp_out.write_source_file(output_dir)
    header = Cpp_header_output(model_list)
    header.write_header_file(output_dir)


##
# Generates microcin objects for all combinations, given a list of AHL objects and information on
# if microcin can be induced, repressed, or constitutively expressed. Currently assumes only one AHl can mediate
# expression at any one time.
#
# @param microcin_ids - List of possible ids for microcin, a proxy for the number of independent microcin
# @param AHL_objects - A list containing AHL objects
# @param microcin_induced - Determines whether microcin can be induced by AHL
# @param microcin_repressed - Determines whether microcin can be repressed by AHL
# @param microcin_constitutive - Determines whether microcin can be constitutively expressed
##
def generate_microcin_combinations(microcin_ids, AHL_objects, microcin_induced=False,
                                   microcin_repressed=False, microcin_constitutive=False):
    microcin_config_idx = 0
    microcin_objects = []
    microcin_config_data = []  # Microcin reference contains - config index, microcin_id, inducer_id, repressor_id, constitutive

    if microcin_induced is True:
        for AHL in AHL_objects:
            for m_id in microcin_ids:
                microcin_objects.append(species.Microcin(microcin_config_idx, m_id, [AHL], []))
                config_data = [microcin_config_idx, m_id, AHL.id, np.nan]
                microcin_config_data.append(config_data)
                microcin_config_idx += 1

    if microcin_repressed is True:
        for AHL in AHL_objects:
            for m_id in microcin_ids:
                microcin_objects.append(species.Microcin(microcin_config_idx, m_id, [], [AHL]))

                config_data = [microcin_config_idx, m_id, np.nan, AHL.id]
                microcin_config_data.append(config_data)

                microcin_config_idx += 1

    if microcin_constitutive is True:
        for m_id in microcin_ids:
            microcin_objects.append(
                species.Microcin(microcin_config_idx, m_id, np.nan, np.nan, constitutive_expression=True))
            microcin_config_idx += 1

    microcin_config_df = pd.DataFrame(columns=["microcin_idx", "microcin_id", "inducer_id", "repressor_id"],
                                      data=microcin_config_data)

    return microcin_objects, microcin_config_df


##
# Generates antitoxin objects for all combinations, given a list of AHL objects and information on
# if antitoxin can be induced, repressed, or constitutively expressed. Currently assumes only one AHl can mediate
# expression at any one time.
#
# @param antitoxin - List of possible ids for antitoxin, a proxy for the number of independent antitoxin
# @param AHL_objects - A list containing AHL objects
# @param antitoxin - Determines whether antitoxin can be induced by AHL
# @param antitoxin_repressed - Determines whether antitoxin can be repressed by AHL
# @param antitoxin_constitutive - Determines whether antitoxin can be constitutively expressed
##
def generate_antitoxin_combinations(antitoxin_ids, AHL_objects, antitoxin_induced=False,
                                    antitoxin_repressed=False, antitoxin_constitutive=False):
    antitoxin_config_idx = 0
    antitoxin_objects = []
    antitoxin_config_data = []  # Microcin reference contains - config index, microcin_id, inducer_id, repressor_id, constitutive

    if antitoxin_induced is True:
        for AHL in AHL_objects:
            for v_id in antitoxin_ids:
                antitoxin_objects.append(species.Antitoxin(antitoxin_config_idx, v_id, [AHL], []))
                config_data = [antitoxin_config_idx, v_id, AHL.id, np.nan]
                antitoxin_config_data.append(config_data)
                antitoxin_config_idx += 1

    if antitoxin_repressed is True:
        for AHL in AHL_objects:
            for v_id in antitoxin_ids:
                antitoxin_objects.append(species.Antitoxin(antitoxin_config_idx, v_id, [], [AHL]))

                config_data = [antitoxin_config_idx, v_id, np.nan, AHL.id]
                antitoxin_config_data.append(config_data)

                antitoxin_config_idx += 1

    if antitoxin_constitutive is True:
        for v_id in antitoxin_ids:
            antitoxin_objects.append(
                species.Antitoxin(antitoxin_config_idx, v_id, np.nan, np.nan, constitutive_expression=True))
            antitoxin_config_idx += 1

    antitoxin_config_df = pd.DataFrame(columns=["antitoxin_idx", "antitoxin_id", "inducer_id", "repressor_id"],
                                       data=antitoxin_config_data)

    return antitoxin_objects, antitoxin_config_df


##
# Generates immunity objects for all combinations, given a list of AHL objects and information on
# if immunity can be induced, repressed, or constitutively expressed. Currently assumes only one AHl can mediate
# expression at any one time.
#
# @param immunity_ids - List of possible ids for immunity, a proxy for the number of independent immunity
# @param AHL_objects - A list containing AHL objects
# @param immunity_induced - Determines whether immunity can be induced by AHL
# @param immunity_repressed - Determines whether immunity can be repressed by AHL
# @param immunity_constitutive - Determines whether immunity can be constitutively expressed
##
def generate_immunity_combinations(immunity_ids, AHL_objects, immunity_induced=False,
                                   immunity_repressed=False, immunity_constitutive=False):
    immunity_config_idx = 0
    immunity_objects = []
    immunity_config_data = []  # Microcin reference contains - config index, immunity_id, inducer_id, repressor_id, constitutive

    if immunity_induced is True:
        for AHL in AHL_objects:
            for i_id in immunity_ids:
                immunity_objects.append(species.Immunity(immunity_config_idx, i_id, [AHL], []))
                config_data = [immunity_config_idx, i_id, AHL.id, np.nan]
                immunity_config_data.append(config_data)
                immunity_config_idx += 1

    if immunity_repressed is True:
        for AHL in AHL_objects:
            for i_id in immunity_ids:
                immunity_objects.append(species.Immunity(immunity_config_idx, i_id, [], [AHL]))

                config_data = [immunity_config_idx, i_id, np.nan, AHL.id]
                immunity_config_data.append(config_data)

                immunity_config_idx += 1

    if immunity_constitutive is True:
        for i_id in immunity_ids:
            immunity_objects.append(
                species.Immunity(immunity_config_idx, i_id, np.nan, np.nan, constitutive_expression=True))
            immunity_config_idx += 1

    immunity_config_df = pd.DataFrame(columns=["immunity_idx", "immunity_id", "inducer_id", "repressor_id"],
                                      data=immunity_config_data)

    return immunity_objects, immunity_config_df


##
# Generates toxin objects for all combinations, given a list of AHL objects and information on
# if toxin can be induced, repressed, or constitutively expressed. Currently assumes only one AHl can mediate
# expression at any one time.
#
# @param toxin_ids - List of possible ids for toxin, a proxy for the number of independent toxin
# @param AHL_objects - A list containing AHL objects
# @param toxin_induced - Determines whether toxin can be induced by AHL
# @param toxin_repressed - Determines whether toxin can be repressed by AHL
# @param toxin_constitutive - Determines whether toxin can be constitutively expressed
##
def generate_toxin_combinations(toxin_ids, AHL_objects, toxin_induced=False,
                                toxin_repressed=False, toxin_constitutive=False):
    toxin_config_idx = 0
    toxin_objects = []
    toxin_config_data = []  # Microcin reference contains - config index, toxin_id, inducer_id, repressor_id, constitutive

    if toxin_induced is True:
        for AHL in AHL_objects:
            for t_id in toxin_ids:
                toxin_objects.append(species.Toxin(toxin_config_idx, t_id, [AHL], []))
                config_data = [toxin_config_idx, t_id, AHL.id, np.nan]
                toxin_config_data.append(config_data)
                toxin_config_idx += 1

    if toxin_repressed is True:
        for AHL in AHL_objects:
            for t_id in toxin_ids:
                toxin_objects.append(species.Toxin(toxin_config_idx, t_id, [], [AHL]))

                config_data = [toxin_config_idx, t_id, np.nan, AHL.id]
                toxin_config_data.append(config_data)

                toxin_config_idx += 1

    if toxin_constitutive is True:
        for t_id in toxin_ids:
            toxin_objects.append(species.Toxin(toxin_config_idx, t_id, np.nan, np.nan, constitutive_expression=True))
            toxin_config_idx += 1

    toxin_config_df = pd.DataFrame(columns=["toxin_idx", "toxin_id", "inducer_id", "repressor_id"],
                                   data=toxin_config_data)

    return toxin_objects, toxin_config_df


def remove_identical_part_lists(parts_list):
    unique_part_combos = []

    for candidate in parts_list:
        unique = True
        for existing_combo in unique_part_combos:
            if len(existing_combo) != len(candidate):
                continue

            match = True
            for p in candidate:
                if p not in existing_combo:
                    match = False

            if match:
                unique = False
                break

        if unique:
            unique_part_combos.append(candidate)

    return unique_part_combos


##
# Generates strain objects for all combinations, given a list of microcin objects and substrate objects
##
class model_space():
    def __init__(self, strain_ids, microcin_objects, AHL_objects, substrate_objects, antitoxin_objects,
                 immunity_objects, toxin_objects,
                 max_microcin_parts, max_AHL_parts, max_substrate_dependencies, max_antitoxins, max_immunity,
                 max_toxins, max_microcin_sensitivities=1):
        self.strain_ids = strain_ids
        self.strain_objects = []
        self.microcin_objects = microcin_objects
        self.AHL_objects = AHL_objects
        self.substrate_objects = substrate_objects
        self.antitoxin_objects = antitoxin_objects
        self.immunity_objects = immunity_objects
        self.toxin_objects = toxin_objects

        self.microcin_ids = list(set([m.id for m in microcin_objects]))

        self.part_combinations = []
        self.models_list = []

        # Maximum parts for each strain
        self.max_microcin_parts = max_microcin_parts
        self.max_AHL_parts = max_AHL_parts
        self.max_substrate_parts = max_substrate_dependencies
        self.max_microcin_sensitivities = max_microcin_sensitivities
        self.max_antitoxins = max_antitoxins
        self.max_immunity = max_immunity
        self.max_toxins = max_toxins

    def generate_part_combinations(
            self, strain_max_microcin, strain_max_AHL, strain_max_sub_dependencies,
            strain_max_microcin_sens, strain_max_sub_production, strain_max_antitoxin,
            strain_max_immunity, strain_max_toxin):
        # Construct possible combinations for each part. [None] added to AHL production, microcin production and
        # microcin sensitivity represent empty part. The purpose of this is so we generate combinations with one or more
        # of each part.

        microcin_production_lists = [list(i for i in m if i != None) for m in
                                     itertools.combinations(
                                         self.microcin_objects + [None for n in range(strain_max_microcin - 1)],
                                         strain_max_microcin)]

        AHL_production_lists = [list(i for i in a if i != None) for a in
                                itertools.combinations(self.AHL_objects + [None for n in range(strain_max_AHL - 1)],
                                                       strain_max_AHL)]

        substrate_dependencies_list = [list(i for i in s if i != None) for s in
                                       itertools.combinations(self.substrate_objects + [None for n in range(
                                           strain_max_sub_dependencies - 1)], strain_max_sub_dependencies)]

        microcin_sensitivities_list = [list(i for i in m_id if i != None) for m_id in
                                       itertools.combinations(
                                           self.microcin_ids + [None for n in range(strain_max_microcin_sens - 1)],
                                           strain_max_microcin_sens)]

        substrate_production_list = [list(i for i in s if i != None) for s in
                                     itertools.combinations(
                                         self.substrate_objects + [None for n in range(strain_max_sub_production - 1)],
                                         strain_max_sub_production)]

        antitoxin_list = [list(i for i in v if i != None) for v in
                          itertools.combinations(
                              self.antitoxin_objects + [None for n in range(strain_max_antitoxin - 1)],
                              strain_max_antitoxin)]

        immunity_list = [list(i for i in v if i != None) for v in
                         itertools.combinations(self.immunity_objects + [None for n in range(strain_max_immunity - 1)],
                                                strain_max_immunity)]

        toxin_list = [list(i for i in v if i != None) for v in
                      itertools.combinations(self.toxin_objects + [None for n in range(strain_max_toxin - 1)],
                                             strain_max_toxin)]

        microcin_production_lists = remove_identical_part_lists(microcin_production_lists)
        AHL_production_lists = remove_identical_part_lists(AHL_production_lists)
        substrate_dependencies_list = remove_identical_part_lists(substrate_dependencies_list)
        substrate_production_list = remove_identical_part_lists(substrate_production_list)
        microcin_sensitivities_list = remove_identical_part_lists(microcin_sensitivities_list)
        antitoxin_list = remove_identical_part_lists(antitoxin_list)
        immunity_list = remove_identical_part_lists(immunity_list)
        toxin_list = remove_identical_part_lists(toxin_list)

        # Append empty list representing no production or sensitivity, only necessary if more than two max parts
        if strain_max_AHL >= 1:
            AHL_production_lists.append([])

        if strain_max_microcin >= 1:
            microcin_production_lists.append([])

        if strain_max_microcin_sens >= 1:
            microcin_sensitivities_list.append([])

        if strain_max_sub_dependencies >= 1:
            substrate_dependencies_list.append([])

        if strain_max_sub_production >= 1:
            substrate_production_list.append([])

        if strain_max_antitoxin >= 1:
            antitoxin_list.append([])

        if strain_max_immunity >= 1:
            immunity_list.append([])

        if strain_max_toxin >= 1:
            toxin_list.append([])

        # for sub in microcin_production_lists:
        #     for s in sub:
        #         print(s.id)
        #     print("")

        # exit()

        # Generate all different combinations of parts
        for m in microcin_production_lists:
            for a in AHL_production_lists:
                for s in substrate_dependencies_list:
                    for sensi in microcin_sensitivities_list:
                        for s_prod in substrate_production_list:
                            for v in antitoxin_list:
                                for i in immunity_list:
                                    for t in toxin_list:
                                        self.part_combinations.append([m, a, s, sensi, s_prod, v, i, t])

        keep_combinations = []
        for x in self.part_combinations:
            if len(x[2]) == 0:
                continue

            if x in keep_combinations:
                print(x)
                print("yes!")
                continue

            else:
                keep_combinations.append(x)

        self.part_combinations = keep_combinations
        # exit()

        return self.part_combinations

    def remove_direct_symmetries(self):
        all_adj_mats = []
        keep_idx_0 = []
        adj_mat_sums = []

        # Check for direct matches where two adjacency matrices match
        for idx, model in enumerate(self.models_list):

            keep_model = True

            mat_sum = model.adjacency_matrix.sum()
            candidate_idxs = np.argwhere(np.array(adj_mat_sums) == mat_sum)

            # Check if candidate is equal to any matrix in any already stored
            for i in candidate_idxs:
                if np.array_equal(all_adj_mats[i[0]], model.adjacency_matrix):
                    keep_model = False
                    break

            if keep_model:
                keep_idx_0.append(idx)
                all_adj_mats.append(model.adjacency_matrix)
                adj_mat_sums.append(mat_sum)

        self.models_list = [self.models_list[i] for i in keep_idx_0]

    def remove_indirect_symmetries(self):
        strain_init_idx = 0
        microcin_init_idx = strain_init_idx + len(self.strain_ids) + len(self.substrate_objects)
        AHL_init_idx = microcin_init_idx + self.max_microcin_parts

        strains_index_range = range(strain_init_idx, len(self.strain_ids))
        mic_index_range = range(microcin_init_idx, microcin_init_idx + self.max_microcin_parts)
        AHL_index_range = range(AHL_init_idx, AHL_init_idx + self.max_AHL_parts)

        # Flip strain columns and remove symmetrical strains
        keep_idx = []
        clean_adj_mats = []
        adj_mat_sums = []

        for idx, model in enumerate(self.models_list):
            permute_strains = list(itertools.permutations(strains_index_range))
            original_config = permute_strains[0]
            candidate_adj_mat = model.adjacency_matrix

            mat_sum = candidate_adj_mat.sum()
            candidate_idxs = np.argwhere(np.array(adj_mat_sums) == mat_sum)

            match = False
            for perm in permute_strains[1:]:
                new_adj = candidate_adj_mat.copy()

                # Shuffle first configuration to new configuration
                new_adj.T[[original_config]] = new_adj.T[[perm]]
                new_adj[[original_config]] = new_adj[[perm]]

                for adj_idx in candidate_idxs:
                    if np.array_equal(clean_adj_mats[adj_idx[0]], new_adj):
                        match = True
                        break
                if match:
                    break

            if match is False:
                keep_idx.append(idx)
                clean_adj_mats.append(candidate_adj_mat)
                adj_mat_sums.append(mat_sum)

        self.models_list = [self.models_list[i] for i in keep_idx]

    def remove_symmetries(self):
        all_adj_mats = []
        keep_idx_0 = []

        # Check for direct matches where two adjacency matrices match
        for idx, model in enumerate(self.models_list):

            keep_model = True

            # Check if candidate is equal to any matrix in any already stored
            for x in all_adj_mats:
                if np.array_equal(x, model.adjacency_matrix):
                    keep_model = False
                    break

            if keep_model:
                keep_idx_0.append(idx)
                all_adj_mats.append(model.adjacency_matrix)

        clean_stage_1_adj_mats = []

        for idx, model in enumerate(self.models_list):
            if idx in keep_idx_0:
                clean_stage_1_adj_mats.append(model.adjacency_matrix)

        strain_init_idx = 0
        microcin_init_idx = strain_init_idx + len(self.strain_ids) + len(self.substrate_objects)
        AHL_init_idx = microcin_init_idx + self.max_microcin_parts

        strains_index_range = range(strain_init_idx, len(self.strain_ids))
        mic_index_range = range(microcin_init_idx, microcin_init_idx + self.max_microcin_parts)
        AHL_index_range = range(AHL_init_idx, AHL_init_idx + self.max_AHL_parts)

        # Flip strain columns and remove symmetrical strains
        keep_idx_1 = []

        # Shuffles the strain columns and compares to
        for idx in keep_idx_0:
            permute_strains = list(itertools.permutations(strains_index_range))
            original_config = permute_strains[0]
            model_adj = self.models_list[idx].adjacency_matrix

            match = False
            for perm in permute_strains[1:]:
                new_adj = model_adj.copy()
                # Shuffle first configuration to new configuration
                new_adj.T[[original_config]] = new_adj.T[[perm]]
                new_adj[[original_config]] = new_adj[[perm]]

                if any(np.array_equal(x, new_adj) for x in clean_stage_1_adj_mats):
                    match = True
                    break

            if match is False:
                keep_idx_1.append(idx)

        self.models_list = [self.models_list[i] for i in keep_idx_1]

    def generate_models(self):
        model_idx = 0

        system_combinations = itertools.combinations(self.part_combinations, len(self.strain_ids))
        total_sys = 0
        for sys in system_combinations:
            model_strains = []
            for idx, N_id in enumerate(self.strain_ids):
                new_strain = species.Strain(N_id, *sys[idx])
                model_strains.append(new_strain)

            new_model = model.Model(model_idx, model_strains)

            if new_model.is_legal():
                self.models_list.append(new_model)
                model_idx += 1
                new_model.generate_adjacency_matrix(self.max_substrate_parts, self.max_AHL_parts,
                                                    self.max_microcin_parts, len(self.strain_ids), self.max_antitoxins,
                                                    self.max_immunity, self.max_toxins)
                total_sys += 1

    def reset_model_indexes(self):
        for idx, model in enumerate(self.models_list):
            model.idx = idx

    def spock_manu_model_filter(self):
        keep_list = []

        # keep models with only one strain engineered
        for model in self.models_list:
            if not any(model.substrate_ids):
                continue

            keep = True

            # Only one strain is producing stuff
            if sum(model.adjacency_matrix[:, 0]) == 0 or sum(model.adjacency_matrix[:, 1]) == 0:
                pass
            else:
                keep = False

            # neither strain produces glucose
            if sum(model.adjacency_matrix[2]) != 0:
                keep = False

            # Both strains use glucose
            if sum(model.adjacency_matrix[:, 2]) != 2:
                keep = False

            if keep:
                keep_list.append(model)
            # if sum(model.adjacency_matrix[:, 1]) == 0 and sum(model.adjacency_matrix[2]) == 0 and sum(model.adjacency_matrix[:, 2]) == 2 and model.adjacency_matrix[0][3] == 0:
            #     keep_list.append(model)

        self.models_list = keep_list

    def one_species_filter(self):
        keep_list = []

        for model in self.models_list:
            if sum(model.adjacency_matrix[:, 0]) == 3 and sum(model.adjacency_matrix[:, 2]) == -1:
                keep_list.append(model)

        self.models_list = keep_list

    def aux_filter(self):
        keep_list = []

        for model in self.models_list:
            substrate_production = []
            substrate_dependencies = []

            for strain in model.strains:
                for s in strain.substrate_dependences:
                    if s.id != 'glu':
                        substrate_dependencies.append(s)
                for s in strain.substrate_production:
                    if s.id != 'glu':
                        substrate_production.append(s)

            keep = True
            for s in substrate_dependencies:
                if s not in substrate_production:
                    keep = False

            if keep:
                keep_list.append(model)

        self.models_list = keep_list

    def max_immunity_filter(self, max_immunity):
        keep_list = []
        for model in self.models_list:
            keep_model = True

            # Get all microcins of a model
            microcins = []
            for strain in model.strains:
                microcins = microcins + [m.id for m in strain.microcins]

            microcins = list(set(microcins))

            for strain in model.strains:
                n_immunity = len(microcins) - len(strain.sensitivities)

                if n_immunity > max_immunity:
                    keep_model = False

            if keep_model:
                keep_list.append(model)

        self.models_list = keep_list

    # Remove models containing strains that are 
    # sensitive to their own constitutive microcin
    def self_sensitivity_filter(self, remove_constitutive_only=False):
        keep_list = []
        for model in self.models_list:
            keep_model = True

            for strain in model.strains:
                # Get constitutive microcins of a model

                if remove_constitutive_only:
                    expressed_microcins = [m.id for m in strain.microcins if m.constitutive_expression]

                else:
                    expressed_microcins = [m.id for m in strain.microcins]

                for m in expressed_microcins:
                    if m in strain.sensitivities:
                        keep_model = False

            if keep_model:
                keep_list.append(model)

        self.models_list = keep_list

    def one_predator_two_prey_filter(self):
        keep_list = []

        for model in self.models_list:
            keep_model = True
            n_microcin_producers = 0

            # Only one strain produces microcins
            for strain in model.strains:
                if len(strain.microcins) > 0:
                    n_microcin_producers += 1

            if n_microcin_producers > 1:
                keep_model = False
                continue

            # All species sensitive to at least one microcin
            for strain in model.strains:
                if len(strain.sensitivities) < 0:
                    keep_model = False

            if keep_model == True:
                keep_list.append(model)
        
        self.models_list = keep_list

    def generate_model_reference_table(self, max_microcin_parts, max_AHL_parts,
                                       max_substrate_dependencies, max_microcin_sensitivities):
        # Make column headers
        cell_prefix = 'cell_IDX_'
        microcin_prefix = 'M_IDX'
        AHL_prefix = 'AHL_IDX'
        substrate_prefix = 'S_IDX'
        sensitivity_prefix = 'Sens_IDX'

        models_datasheet_cols = []

        for cell_idx, i in enumerate(self.strain_ids):
            for m_idx, m in enumerate(range(max_microcin_parts)):
                cell = cell_prefix.replace('IDX', str(cell_idx))
                j = microcin_prefix.replace('IDX', str(m_idx))
                models_datasheet_cols.append(cell + j)

            for a_idx, a in enumerate(range(max_AHL_parts)):
                cell = cell_prefix.replace('IDX', str(cell_idx))
                j = AHL_prefix.replace('IDX', str(a_idx))
                models_datasheet_cols.append(cell + j)

            for s_idx, s in enumerate(range(max_substrate_dependencies)):
                cell = cell_prefix.replace('IDX', str(cell_idx))
                j = substrate_prefix.replace('IDX', str(s_idx))
                models_datasheet_cols.append(cell + j)

            for sens_idx, sens in enumerate(range(max_microcin_sensitivities)):
                cell = cell_prefix.replace('IDX', str(cell_idx))
                j = sensitivity_prefix.replace('IDX', str(sens_idx))
                models_datasheet_cols.append(cell + j)

        models_datasheet = pd.DataFrame(columns=models_datasheet_cols)

        # NaN's used to fill an empty cell, where no action takes place
        for model_idx, model in enumerate(self.models_list):
            model_data = []
            for strain_idx, strain in enumerate(model.strains):

                for idx_m, m in enumerate(range(max_microcin_parts)):
                    try:
                        model_data.append(strain.microcins[idx_m].config_idx)
                    except(IndexError):
                        model_data.append(np.nan)

                for idx_a, a in enumerate(range(max_AHL_parts)):
                    try:
                        model_data.append(strain.AHLs[idx_a].id)
                    except(IndexError):
                        model_data.append(np.nan)

                for idx_s, s in enumerate(range(max_substrate_dependencies)):
                    try:
                        model_data.append(strain.substrate_dependences[idx_s].id)
                    except(IndexError):
                        model_data.append(np.nan)

                for idx_sens, sens in enumerate(range(max_microcin_sensitivities)):
                    try:
                        model_data.append(strain.sensitivities[idx_sens])

                    except(IndexError):
                        model_data.append(np.nan)

            models_datasheet.loc[model.idx] = model_data

        return models_datasheet
