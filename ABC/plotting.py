import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams['figure.figsize'] = [15, 10]

font = {'size': 30, }
axes = {'labelsize': 'medium', 'titlesize': 'medium'}


def plot_simulation(out_pdf, sim_idx, model_ref, state, time_points, plot_species_idx, judgement, max_LE, error_msg=None):
    plt_title = "Model idx: " + str(model_ref) + " Sim_idx: " + str(sim_idx) + "Accepted: " + str(judgement) + "MLE: " + str(np.around(max_LE, 4))

    if error_msg is not None and error_msg is not np.nan and error_msg is not "":
        plt_title = plt_title + " error: " + error_msg
        return 0

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

    plot_skip = int(len(time_points) * 0.1)
    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4

    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    for i in plot_species_idx:
        sns.lineplot(x=time_points[plot_skip:], y=state[:, i][plot_skip:], label=str(i))
        # print("max: ", np.max(state[:, i]))
        # print("min: ", np.min(state[:, i]))


    # plt.axhline(0.001, ls='--')
    # plt.ylim(10**-4, 10**1.5)
    # plt.yscale('symlog')

    ax.get_legend().remove()

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.tick_params(labelsize=15)
    # ax.margins(x=0)
    # ax.margins(y=0)
    fig.tight_layout()

    plt.title(plt_title)
    out_pdf.savefig()
    plt.close()
    plt.clf()

def plot_separation(out_pdf, sim_idx, model_ref, sol, sol_theta, t, log_timeseries=False):
    width_inches = 95*4 / 25.4
    height_inches = 95*4 / 25.4
    fig, axes = plt.subplots(figsize=(width_inches,height_inches), nrows=4, ncols=2)
    plt_title = "Model idx: " + str(model_ref) + " Sim_idx: " + str(sim_idx)

    species_0_min = np.min([np.min(sol[:, 0]), np.min(sol_theta[:, 0])])
    species_1_min = np.min([np.min(sol[:, 1]), np.min(sol_theta[:, 1])])
    species_2_min = np.min([np.min(sol[:, 2]), np.min(sol_theta[:, 2])])
    species_3_min = np.min([np.min(sol[:, 3]), np.min(sol_theta[:, 3])])
    
    species_0_max = np.max([np.max(sol[:, 0]), np.max(sol_theta[:, 0])])
    species_1_max = np.max([np.max(sol[:, 1]), np.max(sol_theta[:, 1])])
    species_2_max = np.max([np.max(sol[:, 2]), np.max(sol_theta[:, 2])])
    species_3_max = np.max([np.max(sol[:, 3]), np.max(sol_theta[:, 3])])

    ax_row = axes[0]
    # Plot time series
    for idx, data in enumerate([sol, sol_theta]):
        ax = ax_row[idx]
        sns.lineplot(x=t, y=data[:, 0], ax=ax, estimator=None)
        sns.lineplot(x=t, y=data[:, 1], ax=ax, estimator=None)
        sns.lineplot(x=t, y=data[:, 2], ax=ax, estimator=None)
        sns.lineplot(x=t, y=data[:, 3], ax=ax, estimator=None)

        if log_timeseries:
            ax.set_yscale('log')

    ax_row = axes[1]

    # Plot prey_1 predator phase
    for idx, data in enumerate([sol, sol_theta]):
        ax = ax_row[idx]
        # sns.lineplot(x=data[:, 0], y=data[:, 2], ax=ax, estimator=None)
        ax.plot(data[:, 0], data[:, 2])
        x_min = species_0_min - 0.0001
        x_max = species_0_max + 0.0001

        y_min = species_2_min - 0.0001
        y_max = species_2_max + 0.0001

        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_min, y_max])

    ax_row = axes[2]
    # Plot prey_2 predator phase

    for idx, data in enumerate([sol, sol_theta]):
        ax = ax_row[idx]
        # sns.lineplot(x=data[:, 1], y=data[:, 2], ax=ax, estimator=None, ci=None)
        ax.plot(data[:, 1], data[:, 2])
        # ax.fill_between(data[:, 1], data[:, 2], color="red", alpha=0.3)
        x_min = species_1_min - 0.0001
        x_max = species_1_max + 0.0001

        y_min = species_2_min - 0.0001
        y_max = species_2_max + 0.0001
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_min, y_max])

    ax_row = axes[3]
    # Plot prey_2 predator phase
    for idx, data in enumerate([sol, sol_theta]):
        ax = ax_row[idx]
        # sns.lineplot(x=data[:, 1], y=data[:, 2], ax=ax, estimator=None, ci=None)
        ax.plot(data[:, 2], data[:, 3])
        # ax.fill_between(data[:, 1], data[:, 2], color="red", alpha=0.3)
        x_min = species_2_min - 0.0001
        x_max = species_2_max + 0.0001

        y_min = species_3_min - 0.0001
        y_max = species_3_max + 0.0001

        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_min, y_max])



    for ax_row in axes:
        for ax in ax_row:
            ax.unicode_minus = True

            ax.set_ylabel('')
            # ax.set(xlim=(-0.5,None))
            # ax.set(ylim=(-0))
            ax.spines["right"].set_visible(False)
            ax.spines["top"].set_visible(False)
            ax.spines["left"].set_alpha(0.5)
            ax.spines["bottom"].set_alpha(0.5)
            # ax.tick_params(labelsize=15)
            # ax.margins(x=0)
            # ax.margins(y=0)

    out_pdf.savefig()

    # plt.title(str(sep_coeff))

    for idx, data in enumerate([sol, sol_theta]):
        fig, axes = plt.subplots(figsize=(width_inches,height_inches))
        ax = Axes3D(fig)
        ax.plot(data[:, 0], data[:, 1], data[:, 2])
        out_pdf.savefig()
        plt.close('all')
        plt.clf()


def plot_lorenz(out_pdf, sim_idx, model_ref, state, time_points, plot_species_idx, error_msg=None):
    colours = ['#e30b17', '#1d71b8', '#73ba65', '#ac7cb5', '#706f6f']
    sns.set_context("talk")
    sns.set_style("white")


    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4

    plot_skip = int(len(time_points) * 0.1)


    fig = plt.figure(figsize=(width_inches, height_inches))
    ax = fig.gca(projection='3d')
    ax.plot(state[:, 0][plot_skip:], state[:, 1][plot_skip:], state[:, 2][plot_skip:], color=colours[1])
    ax.tick_params(labelsize=0)
    ax.set(xticklabels=[])
    ax.set(xlabel='')

    ax.set(yticklabels=[])
    ax.set(ylabel='')
    ax.set(zticklabels=[])
    ax.set(zlabel='')

    fig.tight_layout()

    out_pdf.savefig()
    plt.close()
    plt.clf()

def plot_LV_four_species(out_pdf, sim_idx, model_ref, state, time_points, plot_species_idx, error_msg=None):
    colours = ['#e30b17', '#1d71b8', '#73ba65', '#ac7cb5', '#706f6f']
    sns.set_context("talk")
    sns.set_style("white")


    width_inches = 95*4 / 25.4
    height_inches = 51*4 / 25.4

    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    plot_skip = int(len(time_points) * 0.1)
    ax.plot(state[:, 0][plot_skip:], state[:, 1][plot_skip:], label=str(0))
   
    out_pdf.savefig()
    plt.close()

    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    plot_skip = int(len(time_points) * 0.1)
    ax.plot(state[:, 1][plot_skip:], state[:, 2][plot_skip:], label=str(0))
   
    out_pdf.savefig()
    plt.close()


    fig, ax = plt.subplots(figsize=(width_inches, height_inches))
    plot_skip = int(len(time_points) * 0.1)
    ax.plot(state[:, 2][plot_skip:], state[:, 3][plot_skip:], label=str(0))
   
    out_pdf.savefig()
    plt.close()

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_alpha(0.5)
    ax.spines["bottom"].set_alpha(0.5)
    ax.tick_params(labelsize=15)
    # ax.margins(x=0)
    # ax.margins(y=0)
    fig.tight_layout()

    out_pdf.savefig()
    plt.close()
    plt.clf()
