"""
Following flags for a species type
#N# - strain
#S# - substrate
#B# - bacteriocin
#A# - AHL

The flags are replaced with a index
"""


class Strain:
    def __init__(self, strain_id, microcin_expression, AHL_expression, substrate_dependences, microcin_sensitivities,
                 substrate_production, antitoxins, immunity_expression, toxin_expression):
        self.id = strain_id
        self.substrate_dependences = substrate_dependences
        self.microcins = microcin_expression
        self.AHLs = AHL_expression
        self.sensitivities = microcin_sensitivities
        self.substrate_production = substrate_production
        self.antitoxins = antitoxins
        self.immunity = immunity_expression
        self.toxins = toxin_expression

        # self.diff_eqs = {}


class Microcin:
    def __init__(self, config_idx, microcin_id, AHL_inducer_list, AHL_repressor_list, constitutive_expression=False):
        self.config_idx = config_idx
        self.id = microcin_id
        self.AHL_inducers = AHL_inducer_list
        self.AHL_repressors = AHL_repressor_list
        self.constitutive_expression = constitutive_expression


class AHL:
    def __init__(self, AHL_id):
        self.id = AHL_id


class Substrate:
    def __init__(self, substrate_id):
        self.id = substrate_id


class Antitoxin:
    def __init__(self, config_idx, antitoxin_id, AHL_inducer_list, AHL_repressor_list, constitutive_expression=False):
        self.config_idx = config_idx
        self.id = antitoxin_id
        self.AHL_inducers = AHL_inducer_list
        self.AHL_repressors = AHL_repressor_list
        self.constitutive_expression = constitutive_expression


class Immunity:
    def __init__(self, config_idx, immunity_id, AHL_inducer_list, AHL_repressor_list, constitutive_expression=False):
        self.config_idx = config_idx
        self.id = immunity_id
        self.AHL_inducers = AHL_inducer_list
        self.AHL_repressors = AHL_repressor_list
        self.constitutive_expression = constitutive_expression


class Toxin:
    def __init__(self, config_idx, toxin_id, AHL_inducer_list, AHL_repressor_list, constitutive_expression=False):
        self.config_idx = config_idx
        self.id = toxin_id
        self.AHL_inducers = AHL_inducer_list
        self.AHL_repressors = AHL_repressor_list
        self.constitutive_expression = constitutive_expression
