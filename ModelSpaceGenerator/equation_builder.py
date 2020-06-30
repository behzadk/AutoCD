import sympy
import numpy as np

funcs = {
    # Strain growth rate
    'mu_#N#': '( mu_max_#N# * S_#S# / ( K_mu_#S# + S_#S# ) )',

    # Induction of bacteriocin expression by AHL
    'k_b_ind_#B#': '( (C_extra * A_#A#)^n_A_B_#B# / ( K_A_B_#B#^n_A_B_#B# + (C_extra * A_#A#)^n_A_B_#B# ) )',

    # Repression of bacteriocin expression by AHL
    'k_b_repr_#B#': '( K_A_B_#B#^n_A_B_#B# / ( K_A_B_#B#^n_A_B_#B# + (C_extra * A_#A#)^n_A_B_#B# ) )',
    # 'k_b_repr_#B#': '( 1 / 1 + ( KB_#B# + A_#A# )^nB_#B# )',

    # Production of an AHL species
    'A_production': 'kA_#A# * N_#N# * C_OD / C_extra',

    # Function defining sensitivity to microcin
    'omega_B': ' omega_max * (C_extra * B_#B#)^n_omega / ( k_omega_B_#B#^n_omega + (C_extra *  B_#B#)^n_omega ) ',

    # Function defining sensitivity to microcin
    'omega_T': '( k_omega_T_#T# * T_#T# )',

    # Induction of antitoxin expression by AHL
    'k_v_ind_#V#': '( (C_extra *  A_#A#)^n_A_V_#V# / ( K_A_V_#V#^n_A_V_#V#  + (C_extra * A_#A#)^n_A_V_#V#) )',

    # Repression of antitoxin expression by AHL
    'k_v_repr_#V#': '( K_A_V_#V#^n_A_V_#V#  / ( K_A_V_#V#^n_A_V_#V#  + (C_extra * A_#A#)^n_A_V_#V# ) )',
    # 'k_v_repr_#V#': '( 1  / 1 + ( A_#A# / kV_#V#)^nV_#V# )'

    # Impact of immunity protein
    'I_immunity': ' ( k_I_#I#^nI_#I# / ( k_I_#I#^nI_#I# + C_intra * I_#I#^nI_#I# ) ) ',

    # Induction of immunity expression by AHL
    'k_i_ind_#I#': '( (C_extra * A_#A#)^n_A_I_#I# / ( K_A_I_#I#^n_A_I_#I#  + (C_extra * A_#A#)^n_A_I_#I#) )',

    # Repression of immunity expression by AHL
    'k_i_repr_#I#': '( K_A_I_#I#^n_A_I_#I#  / ( K_A_I_#I#^n_A_I_#I#  + (C_extra *  A_#A#)^n_A_I_#I# ) )',

    # Induction of toxin expression by AHL
    'k_t_ind_#T#': '( (C_extra * A_#A#)^n_A_T_#T# / ( K_A_T_#T#^n_A_T_#T#  + (C_extra * A_#A#)^n_A_T_#T#) )',

    # Repression of toxin expression by AHL
    'k_t_repr_#T#': '( K_A_T_#T#^n_A_T_#T#  / ( K_A_T_#T#^n_A_T_#T#  + (C_extra * A_#A#)^n_A_T_#T# ) )',

    # Annihilation of toxin by antitoxin
    'k_ann_T_#T#_V_#V#': '( k_TV_ann * T_#T# * V_#V# )'
}

# Base eqs contain the description of species before interactions with other species.
# This usually consists of the dilution term which we can apply production terms to.
base_eqs = {
    'N_#N#': '( - D * N_#N# )',
    'S_#S#': '( D * ( S0_#S# - S_#S# ) )',
    'B_#B#': '( - D * B_#B# )',
    'A_#A#': '( - D * A_#A# )',
    'V_#V#': ' ',
    'I_#I#': ' ',
    'T_#T#': ' '
}


def gen_strain_growth_diff(strain_id, strain_list):
    dN_dt = base_eqs['N_#N#']

    for strain in strain_list:
        if strain.id is strain_id:

            # Substrate dependent growth and dilution
            dN_dt = dN_dt + ' + N_#N# '
            for s in strain.substrate_dependences:
                dN_dt = dN_dt + ' * ' + funcs['mu_#N#'].replace('#S#', s.id)

            # Strain sensitive to microcin
            for m in strain.sensitivities:
                dN_dt = dN_dt + ' - ( ' + funcs['omega_B']
                dN_dt = dN_dt.replace('#B#', m)

                # Apply immunity function if cognate immunity is present.
                for imm in strain.immunity:
                    if imm.id.split('_')[-1] == m:
                        dN_dt = dN_dt + ' * ' + funcs['I_immunity']
                        dN_dt = dN_dt.replace('#I#', imm.id)
                        break

                dN_dt = dN_dt + ' ) '
                dN_dt = dN_dt + ' * N_#N# '

            for t in strain.toxins:
                dN_dt = dN_dt + ' - ( ' + funcs['omega_T']
                dN_dt = dN_dt.replace('#T#', t.id)

                dN_dt = dN_dt + ' ) '
                dN_dt = dN_dt + ' * N_#N# '

    dN_dt = dN_dt.replace('#N#', strain_id)
    N_key = 'N_#N#'.replace('#N#', strain_id)

    return {N_key: dN_dt}


def gen_diff_eq_antitoxin(antitoxin_id, strain_list):
    dV_dt = base_eqs['V_#V#']

    for strain in strain_list:

        for v in strain.antitoxins:
            if v.id is antitoxin_id:

                # Production term
                dV_dt = dV_dt + ' + ' + ' kV_max_#V# '
                if v.AHL_inducers is not np.nan:
                    # Induction terms
                    for a in v.AHL_inducers:
                        dV_dt = dV_dt + ' * ' + funcs['k_v_ind_#V#'].replace('#A#', a.id)

                if v.AHL_repressors is not np.nan:
                    for a in v.AHL_repressors:
                        dV_dt = dV_dt + ' * ' + funcs['k_v_repr_#V#'].replace('#A#', a.id)

                # Get growth rate of strain producing to create degradation term
                dV_dt = dV_dt + ' - ( V_#V# '
                for s in strain.substrate_dependences:
                    dV_dt = dV_dt + ' * ( ' + funcs['mu_#N#'].replace('#S#', s.id) + ' ) '

                dV_dt = dV_dt + ' / 2 ) '
                dV_dt = dV_dt.replace('#N#', strain.id)

                # Add annihilation term for cognate toxin
                for toxin in strain.toxins:
                    if toxin.id.split('_')[-1] == v.id:
                        dV_dt = dV_dt + ' - ' + funcs['k_ann_T_#T#_V_#V#']
                        dV_dt = dV_dt.replace('#T#', toxin.id)
                        dV_dt = dV_dt.replace('#V#', v.id)
                        break

    dV_dt = dV_dt.replace('#V#', antitoxin_id)
    V_key = 'V_#V#'.replace('#V#', antitoxin_id)

    return {V_key: dV_dt}


def gen_diff_eq_immunity(immunity_id, strain_list):
    dI_dt = base_eqs['I_#I#']

    for strain in strain_list:

        for i in strain.immunity:
            if i.id is immunity_id:
                dI_dt = dI_dt + ' + ' + ' kI_max_#I# / C_intra '

                if i.AHL_inducers is not np.nan:

                    # Induction terms
                    for a in i.AHL_inducers:
                        dI_dt = dI_dt + ' * ' + funcs['k_i_ind_#I#'].replace('#A#', a.id)

                if i.AHL_repressors is not np.nan:
                    for a in i.AHL_repressors:
                        dI_dt = dI_dt + ' * ' + funcs['k_i_repr_#I#'].replace('#A#', a.id)

                # Get growth rate of strain producing
                dI_dt = dI_dt + ' - ( I_#I# * C_intra '

                grwth_dilution = ''
                for s in strain.substrate_dependences:
                    grwth_dilution = grwth_dilution + ' * ( ' + funcs['mu_#N#'].replace('#S#', s.id)

                dI_dt = dI_dt + grwth_dilution
                dI_dt = dI_dt + ' / 2 ) ) '
                dI_dt = dI_dt.replace('#N#', strain.id)

                dI_dt = dI_dt + ' * 1  / C_intra '
    dI_dt = dI_dt.replace('#I#', immunity_id)
    I_key = 'I_#I#'.replace('#I#', immunity_id)

    return {I_key: dI_dt}


def gen_diff_eq_substrate(substrate_id, strain_list):
    dS_dt = base_eqs['S_#S#'].replace('#S#', substrate_id)

    # Term defining consumption of substrate by a strain
    strain_growth_rate = funcs['mu_#N#']
    strain_consumption = strain_growth_rate + ' * N_#N# * C_OD / g_#N# '
    strain_production = ' N_#N# * p_#S# * C_OD'

    # Sum of all consumption by strains
    for strain in strain_list:
        substrates = strain.substrate_dependences
        for consume_s in substrates:
            if consume_s.id is substrate_id:
                dS_dt = dS_dt + ' - ' + strain_consumption
                dS_dt = dS_dt.replace('#N#', strain.id)

    for strain in strain_list:
        for produce_s in strain.substrate_production:
            if produce_s.id is substrate_id:
                dS_dt = dS_dt + ' + ' + strain_production
                dS_dt = dS_dt.replace('#N#', strain.id)

    dS_dt = dS_dt.replace('#S#', substrate_id)
    S_key = 'S_#S#'.replace('#S#', substrate_id)
    return {S_key: dS_dt}


def gen_AHL_diff_eq(AHL_id, strain_list):
    dA_dt = base_eqs['A_#A#']
    production_term = funcs['A_production']

    for strain in strain_list:
        for a in strain.AHLs:
            if a.id is AHL_id:
                dA_dt = dA_dt + ' + ' + production_term.replace('#N#', strain.id)

    dA_dt = dA_dt.replace('#A#', AHL_id)
    AHL_key = 'A_#A#'.replace('#A#', AHL_id)

    return {AHL_key: dA_dt}


def gen_microcin_diff_eq(microcin_id, strain_list):
    dB_dt = base_eqs['B_#B#']

    for strain in strain_list:

        for b in strain.microcins:
            if b.id is microcin_id:
                dB_dt = dB_dt + ' + ' + ' kB_max_#B# '

                if b.AHL_inducers is not np.nan:
                    # Induction terms
                    for a in b.AHL_inducers:
                        dB_dt = dB_dt + ' * ' + funcs['k_b_ind_#B#'].replace('#A#', a.id)

                if b.AHL_repressors is not np.nan:
                    for a in b.AHL_repressors:
                        dB_dt = dB_dt + ' * ' + funcs['k_b_repr_#B#'].replace('#A#', a.id)

                dB_dt = dB_dt + ' * N_#N# * C_OD / C_extra '
                dB_dt = dB_dt.replace('#N#', strain.id)

    dB_dt = dB_dt.replace('#B#', microcin_id)
    B_key = 'B_#B#'.replace('#B#', microcin_id)

    return {B_key: dB_dt}


def gen_toxin_diff_eq(toxin_id, strain_list):
    dT_dt = base_eqs['T_#T#']

    for strain in strain_list:
        for t in strain.toxins:
            if t.id is toxin_id:

                # Add production terms
                dT_dt = dT_dt + ' + ' + ' kT_max_#T# '

                if t.AHL_inducers is not np.nan:
                    # Induction terms
                    for a in t.AHL_inducers:
                        dT_dt = dT_dt + ' * ' + funcs['k_t_ind_#T#'].replace('#A#', a.id)

                if t.AHL_repressors is not np.nan:
                    for a in t.AHL_repressors:
                        dT_dt = dT_dt + ' * ' + funcs['k_t_repr_#T#'].replace('#A#', a.id)

                # Get growth rate of strain producing
                dT_dt = dT_dt + ' - ( T_#T# '
                for s in strain.substrate_dependences:
                    dT_dt = dT_dt + ' * ( ' + funcs['mu_#N#'].replace('#S#', s.id)

                dT_dt = dT_dt + ' / 2 ) )'
                dT_dt = dT_dt.replace('#N#', strain.id)

                # Add annihilation term for cognate toxin
                for antitoxin in strain.antitoxins:
                    if antitoxin.id.split('_')[-1] == t.id:
                        dT_dt = dT_dt + ' - ' + funcs['k_ann_T_#T#_V_#V#']
                        dT_dt = dT_dt.replace('#T#', t.id)
                        dT_dt = dT_dt.replace('#V#', antitoxin.id)
                        break

    dT_dt = dT_dt.replace('#T#', toxin_id)
    T_key = 'T_#T#'.replace('#T#', toxin_id)

    # print(dT_dt)
    # print(strain_list[0].substrate_dependences)
    return {T_key: dT_dt}
