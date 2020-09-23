from .species import Microcin
from .species import AHL
from .species import Strain
from .species import Substrate

from .model import Model
from .cpp_output import Cpp_source_output
from .cpp_output import Cpp_header_output

from .model import Model
from . import model_space_generator
import argparse
import yaml


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--MSG_config', type=str)
    parser.add_argument('--ABC_config', type=str)
    parser.add_argument('--DA_config', type=str)
    parser.add_argument("--exp_num", required=False, default=0)



    args = parser.parse_args()

    config_yaml_path = args.MSG_config

    with open(config_yaml_path, 'r') as yaml_file:
        experiment_config = yaml.load(yaml_file, Loader=yaml.FullLoader)

    # Set strain IDs
    strain_ids = experiment_config['strain_ids']
    use_one_pred_two_prey_filter = experiment_config['use_one_pred_two_prey_filter']

    # Set substrate parts
    substrate_ids = experiment_config['substrate_ids']
    S_glu = Substrate(substrate_ids[0])
    substrate_objects = [Substrate(S) for S in substrate_ids]

    # Set AHL/QS parts
    AHL_ids = experiment_config['AHL_ids']
    AHL_objects = [AHL(QS) for QS in AHL_ids]

    # Set bacteriocin parts
    microcin_ids = experiment_config['bacteriocin_ids']
    # Unpack expression options
    microcin_induced = experiment_config['bacteriocin_induced']
    microcin_repressed = experiment_config['bacteriocin_repressed']
    microcin_constitutive = experiment_config['bacteriocin_constitutive']

    # Generate microcin expression objects from AHLs and microcins
    microcin_objects, microcin_configs_df = model_space_generator.generate_microcin_combinations(microcin_ids,
                                                                                                 AHL_objects,
                                                                                                 microcin_induced=microcin_induced,
                                                                                                 microcin_repressed=microcin_repressed,
                                                                                                 microcin_constitutive=microcin_constitutive)
    # Features currently turned off
    antitoxin_ids = immunity_ids = toxin_ids = []
    antitoxin_objects = immunity_objects = toxin_objects = []

    # Numerical maximum number of parts in a model
    max_substrate_parts = len(substrate_objects)
    max_microcin_parts = len(microcin_ids)
    max_AHL_parts = len(AHL_objects)
    max_strains_parts = len(strain_ids)
    max_toxin_parts = len(toxin_ids)
    max_antitoxins = len(antitoxin_ids)
    max_immunity_parts = len(immunity_ids)

    max_microcin_sensitivities = experiment_config['strain_max_bacteriocin_sens']

    model_space = model_space_generator.model_space(strain_ids=strain_ids, microcin_objects=microcin_objects,
                                                    AHL_objects=AHL_objects,
                                                    substrate_objects=substrate_objects,
                                                    antitoxin_objects=antitoxin_objects,
                                                    immunity_objects=immunity_objects,
                                                    toxin_objects=toxin_objects, max_microcin_parts=max_microcin_parts,
                                                    max_AHL_parts=max_AHL_parts,
                                                    max_substrate_dependencies=max_substrate_parts,
                                                    max_antitoxins=max_antitoxins,
                                                    max_immunity=max_immunity_parts, max_toxins=max_toxin_parts,
                                                    max_microcin_sensitivities=max_microcin_sensitivities)

    # Unpack part combination options
    strain_max_microcin = experiment_config['strain_max_bacteriocin']
    strain_max_AHL = experiment_config['strain_max_AHL']
    strain_max_microcin_sens = experiment_config['strain_max_bacteriocin_sens']

    # These options are currently turned off
    strain_max_sub_dependencies = 1
    strain_max_sub_production = 0
    strain_max_antitoxin = 0
    strain_max_immunity = 0
    strain_max_toxin = 0

    print("Generating part combinations")
    part_combos = model_space.generate_part_combinations(
        strain_max_microcin=strain_max_microcin, strain_max_AHL=strain_max_AHL,
        strain_max_sub_dependencies=strain_max_sub_dependencies,
        strain_max_microcin_sens=strain_max_microcin_sens, strain_max_sub_production=strain_max_sub_production,
        strain_max_antitoxin=strain_max_antitoxin,
        strain_max_immunity=strain_max_immunity, strain_max_toxin=strain_max_toxin
    )

    default_params_path = experiment_config['default_parameters_path']
    default_init_species_path = experiment_config['default_initial_species_path']
    output_dir = experiment_config['output_dir'] + experiment_config['model_space_name'] + '/'

    print("Generating models... \n")
    model_list = model_space.generate_models()

    print("Removing direct symmeteries... \n")
    model_space.remove_direct_symmetries()

    print("Removing indirect symmeteries... \n")
    model_space.remove_indirect_symmetries()

    if use_one_pred_two_prey_filter:
        model_space.one_predator_two_prey_filter()

    # model_space.remove_symmetries()
    print("Total number of models: ", len(model_space.models_list))

    model_space.reset_model_indexes()

    model_list = model_space.models_list

    print("Generating adjacency matricies... \n")
    model_space_generator.generate_adjacency_matricies(model_list, substrate_ids, microcin_ids, AHL_ids, strain_ids,
                                                       antitoxin_ids, immunity_ids, toxin_ids, output_dir)

    print("Building source and header files... \n")
    model_space_generator.generate_simulation_files(model_list, default_params_path, default_init_species_path,
                                                    output_dir)


if __name__ == "__main__":
    main()
