"""
Calculate and plot covariance matrix.

Input files locate in "dist_r", which is calculated by "/Data Process/Gen_Statistics/gen_R.py". The files "n_b_i" are calculated from "calculate_r/gen_R.py".
"""

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

import matplotlib


def stat(l):
    """Return mean and standard deviation of elements in list l."""
    ans = 0.0
    ave = sum(l) / len(l)
    for x in l:
        ans += (x - ave)**2
    ans /= len(l)
    return (ave, sqrt(ans))


def gen_cov(n, beta, i, j):
    """Return covariance between r_i and r_j."""
    data1 = open("dist_r/" + str(n) + "_" + str(beta) + "_" + str(i))
    data2 = open("dist_r/" + str(n) + "_" + str(beta) + "_" + str(j))
    L1 = []
    L2 = []
    for line in data1:
        if line == "":
            break
        L1 += [float(line)]
    for line in data2:
        if line == "":
            break
        L2 += [float(line)]
    e1 = sum(L1) / len(L1)
    e2 = sum(L2) / len(L2)
    s = 0.0
    for i in range(len(L1)):
        s += (L1[i] - e1) * (L2[i] - e2)
    s /= len(L1)
    return s

# Plot Matrix
Beta = [2]
n = 120
for beta in Beta:
    Matrix = [[] for i in range(n - 1)]
    f = open('cov_mat_' + str(n) + '_' + str(beta), 'w')
    for row in range(1, n):
        Row = []
        for col in range(1, n):
            if col != row:
                Row += [str(gen_cov(n, beta, row, col)) + str(' ')]
                Matrix[row - 1] += [gen_cov(n, beta, row, col)]
            else:
                Row += [str(0) + str(' ')]
                Matrix[row - 1] += [0]
        f.writelines(Row)
        f.write('\n')
    plt.figure(figsize=(8, 5))
    plt.title("Covariance Matrix for BKZ_" + str(beta) + " Basis", fontsize=20)
    X = np.array(Matrix)
    cmap = matplotlib.cm.gray_r
    plt.pcolor(X, cmap=cmap)
    plt.colorbar()
    plt.savefig("no_diag_cov_mat_" + str(n) + "_" + str(beta) + ".eps")
