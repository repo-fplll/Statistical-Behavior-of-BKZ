"""
Plot the Head and Tail.

Input the file 'n_b' consisting of the mean and std of r_1,...r_{n-1} of n-dim BKZ_b bases
"""

import matplotlib.pyplot as plt
from math import sqrt


def stat(l):
    """Return mean and standard deviation of elements in list l."""
    ans = 0.0
    ave = sum(l) / len(l)
    for x in l:
        ans += (x - ave)**2
    ans /= len(l)
    return (ave, sqrt(ans))


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


def plot_ht(l_n, l_beta):
    """Plot head and tail."""
    MAX = max(l_n)
    for beta in l_beta:
        fig = plt.figure(figsize=(8, 5))
        plt.title("Head and Tail in " + "BKZ_" + str(beta) + " basis", fontsize=20)
        i = 0
        for n in l_n:
            h_Av = []
            h_Std = []
            t_Av = []
            t_Std = []
            f = open(str(n) + "_" + str(beta))
            cnt = 0
            for line in f:
                if line == "":
                    break
                data = line.split()
                if cnt < n / 2 - 1:
                    h_Av += [float(data[0])]
                    h_Std += [float(data[1])]
                elif cnt >= n / 2 - 1:
                    t_Av += [float(data[0])]
                    t_Std += [float(data[1])]
                cnt += 1
            plt.plot(range(1, n / 2), h_Av, "r" + marker(i) + "-", markersize=5, linewidth=1)
            plt.plot(range(MAX - 1 - n + n / 2 + 1, MAX), t_Av, "r" + marker(i) + "-", markersize=5, linewidth=1, label='Average for n=' + str(n))
            i += 1
        i = 0
        for n in l_n:
            h_Av = []
            h_Std = []
            t_Av = []
            t_Std = []
            f = open(str(n) + "_" + str(beta))
            cnt = 0
            for line in f:
                if line == "":
                    break
                data = line.split()
                if cnt < n / 2 - 1:
                    h_Av += [float(data[0])]
                    h_Std += [float(data[1])]
                elif cnt >= n / 2 - 1:
                    t_Av += [float(data[0])]
                    t_Std += [float(data[1])]
                cnt += 1
            plt.plot(range(1, n / 2), h_Std, "b" + marker(i + 4) + "-", markersize=5, linewidth=1)
            plt.plot(range(MAX - 1 - n + n / 2 + 1, MAX), t_Std, "b" + marker(i + 4) + ":", markersize=5, linewidth=1, label='Std for n=' + str(n))
            i += 1
        plt.legend(loc=1, ncol=2, prop={'size': 16})
        plt.ylim(-0.02, 0.14)
        plt.plot([beta, beta], [0, 0.125], "k--")
        plt.plot([MAX - beta, MAX - beta], [0, 0.125], "k--")
        plt.xticks([0, 20, 40, 60, 80, 100, 120, 140])
        ax = fig.add_subplot(1, 1, 1)
        ax.set_xticklabels([0, 20, 40, 60, "n-60", "n-40", "n-20", "n"])
        plt.savefig("ht_in_bkz_" + str(beta) + '.eps')
        plt.close()

N = [140, 120, 100, 80]
Beta = [2]
plot_ht(N, Beta)
