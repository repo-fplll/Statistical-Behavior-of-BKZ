"""Plot the Head and Tail of covariances."""

import matplotlib.pyplot as plt


def marker(i):
    """Different markers."""
    if i == 0:
        return "^"
    if i == 1:
        return "s"
    if i == 2:
        return "v"
    if i == 3:
        return "o"
    if i == 4:
        return ">"
    if i == 5:
        return "D"
    if i == 6:
        return "<"
    if i == 7:
        return "p"


def classify(l_beta, l_n):
    """Plot the Head and Tail from the covariance matrix."""
    for beta in l_beta:
        plt.figure(figsize=(8, 5))
        plt.title("Covariances " + "in BKZ_" + str(beta) + " basis", fontsize=20)
        index = 0
        for n in l_n:
            M = [[] for j in range(10)]
            f = open("cov_mat_" + str(n) + "_" + str(beta))
            cnt = 0
            for line in f:
                if line == "":
                    break
                Line = line.split()
                for j in range(cnt + 1, min(cnt + 10, n - 1)):
                    M[j - cnt] += [float(Line[j])]
                cnt += 1
            L = M[1]
            m = len(L)
            plt.plot(range(1, m / 2), L[0:m / 2 - 1], marker(index) + "b-", markersize=5, linewidth=1)
            plt.plot(range(138 - m + m / 2 - 1, 138), L[m / 2 - 1:m], marker(index) + "b-", markersize=5, linewidth=1, label='dim=' + str(n))
            L = M[2]
            m = len(L)
            plt.plot(range(1, m / 2), L[0:m / 2 - 1], marker(index) + "r:", markersize=5, linewidth=1)
            plt.plot(range(138 - m + m / 2 - 1, 138), L[m / 2 - 1:m], marker(index) + "r:", markersize=5, linewidth=1)
            index += 1

        plt.xlabel("$i$", fontsize=18)
        if beta == 2:
            plt.plot([beta, beta], [-0.0015, 0.0005], "k--")
            plt.annotate('Cov($r_i$,$r_{i+2}$)', fontsize=18, xy=(90, -0.00005), xytext=(80, 0.00025), arrowprops=dict(color='k', arrowstyle="->"))
            plt.annotate('Cov($r_i$,$r_{i+1}$)', fontsize=18, xy=(95, -0.0006), xytext=(80, -0.0004), arrowprops=dict(color='k', arrowstyle="->"))
            plt.plot([140 - 2 - beta, 140 - 2 - beta], [-0.0015, 0.0005], "k--")
        else:
            plt.ylim([-0.002, 0.0008])
            plt.plot([beta, beta], [-0.001, 0.0007], "k--")
            plt.annotate('Cov($r_i$,$r_{i+2}$)', fontsize=18, xy=(90, 0.00007), xytext=(80, 0.00037), arrowprops=dict(color='k', arrowstyle="->"))
            plt.annotate('Cov($r_i$,$r_{i+1}$)', fontsize=18, xy=(90, -0.0006), xytext=(80, -0.0003), arrowprops=dict(color='k', arrowstyle="->"))
            plt.plot([140 - 2 - beta, 140 - 2 - beta], [-0.001, 0.0007], "k--")
        plt.legend(loc='lower left', ncol=2, prop={'size': 18})
        plt.savefig("ht_cov_" + str(beta) + ".eps")

N = [140, 120, 100, 80]
Beta = [2]
classify(Beta, N)
