"""Plot the c(beta)."""

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


def plot_c(l_n):
    """Plot c(beta)."""
    plt.figure(figsize=(8, 5))
    plt.xlabel("$\\beta$", fontsize=18)
    plt.xlim([0, 32])
    plt.xticks([2, 3, 4, 5] + range(6, 32, 2))
    plt.title("Experimental measure of $c(\\beta)$", fontsize=20)
    index = 0
    for n in l_n:
        N = []
        C = []
        f = open("fulldata_" + str(n))
        for line in f:
            if line == "":
                break
            data = line.split()
            N += [int(data[0])]
            C += [float(data[1])]
        plt.plot(N, C, marker(index) + "-", label='Dim=' + str(n))
        index += 1
    plt.legend(loc="lower right", prop={'size': 18})
    plt.savefig("c.eps")


def plot_var_rhf(n):
    """Plot (v(beta)+2c(beta))/3."""
    n = 100
    plt.figure(figsize=(8, 5))
    f = open("fulldata_" + str(n))
    N = []
    C = []
    for line in f:
        if line == "":
            break
        data = line.split()
        N += [int(data[0])]
        C += [float(data[1])]

    f = open("v_" + str(n))
    V = []
    for line in f:
        if line == "":
            break
        data = line.split()
        V += [float(data[1])**2]

    L = [(2 * C[i] + V[i]) / 3 for i in range(len(C))]
    print L[0]
    plt.plot(N[1:], L[1:], "^-", markersize=6, label="$\\frac{v(\\beta)+2c(\\beta)}{3}$")
    plt.xlim([0, 32])
    plt.xticks(N)
    plt.xlabel("$\\beta$", fontsize=18)
    plt.legend(loc='lower right', prop={'size': 35})
    plt.savefig("varrhf.eps")

N = [140, 120, 100, 80]
plot_c(N)
plot_var_rhf(100)
