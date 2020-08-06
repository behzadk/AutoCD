import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
def main():
    n = 50
    
    m_vals = np.linspace(0.1, 1, n)
    k1_vals = np.linspace(1, 20, n)

    mat = np.zeros([n, n])
    print(mat)

    i = 0
    for m in m_vals:
        j = 0
        for k1 in k1_vals:
            l = (k1 * m)**2 - (4 * k1)

            if l > 0:
                mat[i][j] = 1

            j += 1

        i += 1

    fig, ax = plt.subplots(figsize=(13,10)) 

    sns.heatmap(mat, yticklabels=m_vals, xticklabels=k1_vals, ax=ax)
    # sns.heatmap(binned_grid.statistic, xticklabels=K_1_bins, yticklabels=m_bins, ax=ax)

    ax.set_xticklabels(k1_vals)
    ax.set_yticklabels(m_vals)
    ax.set_facecolor("yellow")
    # ax.set_xticks(K_1_bins)
    # ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    # ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    # ax.invert_xaxis()
    ax.invert_yaxis()

    plt.show()
    plt.matshow(mat)
    plt.show() 


if __name__ == "__main__":
    main()