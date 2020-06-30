import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib as mpl

plt.rcParams['figure.figsize'] = [15, 10]

font = {'size': 30, }
axes = {'labelsize': 'medium', 'titlesize': 'medium'}


def plot_simulation(out_pdf, sim_idx, model_ref, state, time_points, plot_species_idx, error_msg=None):
    plt_title = "Model idx: " + str(model_ref) + " Sim_idx: " + str(sim_idx)

    if error_msg is not None and error_msg is not np.nan and error_msg is not "":
        plt_title = plt_title + " error: " + error_msg

        # for col in range(np.shape(state)[1]):
        #     print(col)
        #     print(np.min(state[:, col]))
        # print("")

    sns.set_context("talk")
    sns.set_style("white")
    mpl.rc('font', **font)
    mpl.rc('axes', **axes)
    ## for Palatino and other serif fonts use:
    mpl.rc('font', **{'family': 'serif', 'serif': ['Palatino']})
    # plt.rcParams['text.usetex'] = True
    plt.rcParams['axes.unicode_minus'] = False

    fig = plt.figure()
    for i in plot_species_idx:
        sns.lineplot(x=time_points, y=state[:, i], label=str(i))
        # print("max: ", np.max(state[:, i]))
        # print("min: ", np.min(state[:, i]))

    # plt.axhline(0.001, ls='--')
    # plt.ylim(10**-4, 10**1.5)
    plt.yscale('log')
    plt.legend()
    plt.title(plt_title)
    out_pdf.savefig()
    plt.close()
