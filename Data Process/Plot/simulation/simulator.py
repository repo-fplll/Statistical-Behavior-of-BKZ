"""
Chen and Nguyen's BKZ Simulator.

The python code is adapted from Michael Walter's implementation.
"""

import argparse
# import numpy as np
import matplotlib.pyplot as plt
from math import log, gamma, pi, sqrt
from copy import copy


def parse_commandline():
    """
    Parse command line arguments and check for consistency.

    returns: parameters
    rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description='BKZ 2.0')
    parser.add_argument('-b', type=int, help="Blocksize")
    parser.add_argument('-c', type=float, help="C")
    args = parser.parse_args()
    return args


def gauss_est(n):
    """Gaussian Heuristic Prediction."""
    return (gamma(n / 2. + 1) / (pi**(n / 2.)))**(1. / n)


def norm(l1, l2):
    """The norm of the difference between two vectors."""
    l = len(l1)
    L = [l1[i] - l2[i] for i in range(l)]
    s = 0.0
    for i in range(l):
        s += L[i]**2
    return s


def simulator(l, beta, tour, rk, delta):
    """
    BKZ Simulator.

    l: Input lattice profile, i.e. log(||b_i*||), here log is base-2 logarithm.
    beta: Blocksize.
    tour: Maximum Tour Number.
    rk: Average Profile of a BKZ_n basis of an n-dimensional unitary volume lattice. Here we set n = 45.
    delta: Reduced parameter in BKZ algorithm. Here we set delta = 0.99.

    Return a simulated profile and tour number.
    """
    n = len(l)
    l1 = copy(l)
    l2 = copy(l)
    c = [rk[-i] - sum(rk[-i:]) / i for i in range(1, 46)]
    c += [log(gauss_est(d)) / log(2.0) for d in xrange(46, beta + 1)]

    for j in xrange(tour):
        phi = True
        for k in xrange(n - min(45, beta)):
            d = min(beta, n - k)
            f = k + d
            logV = sum(l1[:f]) - sum(l2[:k])
            lma = logV / d + c[d - 1]
            if phi:
                if lma < l1[k] + log(sqrt(delta)) / log(2.0):
                    l2[k] = lma
                    phi = False
            else:
                l2[k] = lma

        # early termination
        if phi or norm(l1, l2) == 0:
            return l1, j
        else:
            d = min(45, beta)
            logV = sum(l1) - sum(l2[:-d])

            if beta < 45:
                tmp = sum(rk[-beta:]) / beta
                rk1 = [r - tmp for r in rk[-beta:]]
            else:
                rk1 = rk

            for k, r in zip(range(n - d, n), rk1):
                l2[k] = logV / d + r
            l1 = copy(l2)

    return l1, tour


"""
Average profile of 45-dim BKZ_45 bases with delta = 0.99, sample = 500.
"""

rk = [
    0.7822502258076371450812303010025061666965,
    0.7714030343026966884956152625818504020572,
    0.7458414716124702192878714868129463866353,
    0.7155255753322982692488096745364600792527,
    0.686347773029179819559431052766740322113,
    0.655241542364430551437237681966507807374,
    0.6199560206576186227600544498272938653827,
    0.5888465732297837873510104600427439436316,
    0.5556778156310949179363944949727738276124,
    0.5213747134029462215742967146070441231132,
    0.4854167760379433094186651942436583340168,
    0.449886369660688262284153893233451526612,
    0.4165263287630752708778913984133396297693,
    0.3791376431804212021603461835184134542942,
    0.3419367546283794967454383595395484007895,
    0.3070517154118544977947635743475984781981,
    0.2684530232983071418306764144290355034173,
    0.2304625959537567186918849415633303578943,
    0.1910913902906814648696709468822518829256,
    0.1514901699615026076095603002613643184304,
    0.113846290691147024320972747091218479909,
    0.07825131750415731390901896702416706830263,
    0.03471801835688430662441028573539369972423,
    -0.005001018248190601295406765292028694602777,
    -0.04726216122386221740736544916217098943889,
    -0.08748098995283861561977900578312983270735,
    -0.1249237780153876001265622619484929600731,
    -0.1704775632285835118562289380861329846084,
    -0.2122450046865790029837484098607092164457,
    -0.2589717676567012405453027668045251630247,
    -0.3004558489374823115891643965369439683855,
    -0.3452633452730204111080780648990185000002,
    -0.3876698289547565456913957859796937555075,
    -0.4322739546369419594284977392817381769419,
    -0.4816234027277476181794213516695890575647,
    -0.5251211578216516406847347298025852069259,
    -0.5673073869699161719815094784280518069863,
    -0.6138492136256009645833842114370781928301,
    -0.6661554773089147785114505495585035532713,
    -0.7093281446019751301790279285341966897249,
    -0.7496492099454789559409562116343295201659,
    -0.795611619983795148591809720528544858098,
    -0.8369569117373906741441658141411608085036,
    -0.875182693595705829547881648977636359632,
    -0.8979226599764336680387089018040569499135,
]


def run(l, beta, c, delta):
    """
    Run BKZ simulator.

    l, beta and delta follow the previous meaning.
    c: a parameter to control the tour number.
    Tour number = c * (n/b)^2 * max(log(n) + loglog(||b_i*|| / vol^{1/n})), as suggested by  Hanrot, Pujol and Stehle.
    """
    n = len(l)
    max_log = log(l[0] * n - sum(l))
    max_tour = int(n**2 * 1.0 / beta**2 * max_log * c) + 1
    ll, tour = simulator(l, beta, max_tour, rk, delta)
    return ll, tour, max_tour


def marker(i):
    """Different markers."""
    if i == 0:
        return "^-"
    if i == 1:
        return "s-"
    if i == 2:
        return "v-"
    if i == 3:
        return "o-"
    if i == 4:
        return ">-"
    if i == 5:
        return "D-"
    if i == 6:
        return "<-"
    if i == 7:
        return "p-"


def color(i):
    """Different colors."""
    if i == 0:
        return "b"
    if i == 1:
        return "g"
    if i == 2:
        return "r"
    if i == 3:
        return "c"


def gen_stat(c=10, l_n=[600, 400], l_beta=range(6, 60, 4)):
    """Calculate the statistics of e, d, d'."""
    for n in l_n:
        f = open("sim_stat_" + str(n), "w")
        l = [10 * n - 10] + [-10 for i in range(1, n)]
        l1 = [i for i in l]
        for beta in l_beta:
            l1, n1, m1 = run(l1, beta, c, 0.99)
            print "Maximum and Predicted Tour Number : ", m1, n1

            # Profile.
            profile = open("sim_profile_" + str(n) + "_" + str(beta), "w")
            L = [str(item - sum(l1) / len(l1)) + "\n" for item in l1]
            profile.writelines(L)

            # Profile.
            profile = open("sim_" + str(n) + "_" + str(beta), "w")
            L = [str(l1[i - 1] - l1[i]) + "\n" for i in range(1, n)]
            profile.writelines(L)

            # Head and tail length.
            h = max(15, beta)
            t = max(15, beta)

            # Statitics.
            r = [(l1[i - 1] - l1[i]) * log(2.0) for i in range(1, n)]
            e = sum(r[h: n - t - 1]) / len(r[h: n - t - 1])
            d = sum(r[0: h]) - (h + 0.5) * e
            _d = 0
            for i in range(1, h + 1):
                _d += (i * 1.0 / 2) * (r[i - 1] - e)
            for i in range(1, t + 1):
                _d += (i * 1.0 / 2) * (r[n - 1 - i] - e)
            f.write(str(beta) + " " + str(e) + " " + str(d) + " " + str(_d) + "\n")


def gen_sum(c=10, l_n=[600, 400], l_beta=range(6, 60, 4)):
    """Calculate the statistics of sum_head, weightsum_head, weightsum_tail."""
    for n in l_n:
        f = open("sim_sum_" + str(n), "w")
        l = [10 * n - 10] + [-10 for i in range(1, n)]
        l1 = [i for i in l]
        for beta in l_beta:
            l1, n1, m1 = run(l1, beta, c, 0.99)
            print "Maximum and Predicted Tour Number : ", m1, n1

            # Head and tail length.
            h = max(15, beta)
            t = max(15, beta)

            # Statitics.
            r = [(l1[i - 1] - l1[i]) * log(2.0) for i in range(1, n)]
            sum_head = sum(r[0: h])
            weightsum_head = 0
            weightsum_tail = 0
            for i in range(1, h + 1):
                weightsum_head += (i * 0.5) * r[i - 1]
            for i in range(1, t + 1):
                weightsum_tail += (i * 0.5) * r[n - 1 - i]
            f.write(str(beta) + " " + str(sum_head) + " " + str(weightsum_head) + " " + str(weightsum_tail) + "\n")


def plot_t(l_n=[600, 400], l_beta=range(6, 60, 4)):
    """Plot tail."""
    MAX = max(l_n)
    for beta in l_beta:
        plt.figure(figsize=(8, 5))
        plt.title("Tail in " + "Simualted BKZ_" + str(beta) + " basis", fontsize=20)
        i = 0
        s = 0
        for n in l_n:
            h_Av = []
            t_Av = []
            f = open("sim_" + str(n) + "_" + str(beta))
            cnt = 0
            for line in f:
                if line == "":
                    break
                data = line.split()
                if cnt < n / 2 - 1:
                    h_Av += [float(data[0])]
                elif cnt >= n / 2:
                    t_Av += [float(data[0])]
                cnt += 1
            if i == 0:
                s = (sum(h_Av[2 * max(15, beta):]) + sum(t_Av[0: len(t_Av) - 2 * max(15, beta)])) / (len(h_Av[2 * max(15, beta):]) + len(t_Av[0: len(t_Av) - 2 * max(15, beta)]))
            plt.plot(range(100 + MAX - 1 - n + n / 2 + 1, MAX - 1), t_Av[100:], color(i) + marker(i) + "-", markeredgecolor=color(i), markersize=3, linewidth=1, label='Average for dim ' + str(n))
            i += 1
        plt.legend(loc=1, prop={'size': 16})
        plt.ylim(-0.02, 0.14)
        plt.plot([100 + MAX - 1 - MAX + MAX / 2 + 1 - 10, MAX - 1], [s, s], "k--")
        plt.plot([MAX - beta, MAX - beta], [0, 0.125], "k--")
        plt.savefig("sim_tail_in_bkz_" + str(beta) + '.eps')
        plt.close()


def plot_h(l_n=[600, 400], l_beta=range(6, 60, 4)):
    """Plot head."""
    MAX = max(l_n)
    for beta in l_beta:
        plt.figure(figsize=(8, 5))
        plt.title("Head in " + "Simualted BKZ_" + str(beta) + " basis", fontsize=20)
        i = 0
        s = 0
        for n in l_n:
            h_Av = []
            t_Av = []
            f = open("sim_" + str(n) + "_" + str(beta))
            cnt = 0
            for line in f:
                if line == "":
                    break
                data = line.split()
                if cnt < n / 2 - 1:
                    h_Av += [float(data[0])]
                elif cnt >= n / 2:
                    t_Av += [float(data[0])]
                cnt += 1
            if i == 0:
                s = (sum(h_Av[2 * max(15, beta):]) + sum(t_Av[0: len(t_Av) - 2 * max(15, beta)])) / (len(h_Av[2 * max(15, beta):]) + len(t_Av[0: len(t_Av) - 2 * max(15, beta)]))
            plt.plot(range(1, n / 2 - 100), h_Av[0: len(h_Av) - 100], color(i) + marker(i), markeredgecolor=color(i), markersize=3, linewidth=1, label='Average for dim ' + str(n))
            i += 1
        plt.legend(loc=1, prop={'size': 16})
        plt.ylim(-0.02, 0.14)
        plt.plot([1, MAX / 2 - 100 + 10], [s, s], "k--")
        plt.savefig("sim_head_in_bkz_" + str(beta) + '.eps')
        plt.close()


def plot_ht(l_n=[600, 400], l_beta=range(6, 60, 4)):
    """Plot head and tail."""
    MAX = max(l_n)
    for beta in l_beta:
        plt.figure(figsize=(8, 5))
        plt.title("Head and tail in " + "Simualted BKZ_" + str(beta) + " basis", fontsize=20)
        i = 0
        s = 0
        for n in l_n:
            h_Av = []
            t_Av = []
            f = open("sim_" + str(n) + "_" + str(beta))
            cnt = 0
            for line in f:
                if line == "":
                    break
                data = line.split()
                if cnt < n / 2 - 1:
                    h_Av += [float(data[0])]
                elif cnt >= n / 2:
                    t_Av += [float(data[0])]
                cnt += 1
            if i == 0:
                s = (sum(h_Av[2 * max(15, beta):]) + sum(t_Av[0: len(t_Av) - 2 * max(15, beta)])) / (len(h_Av[2 * max(15, beta):]) + len(t_Av[0: len(t_Av) - 2 * max(15, beta)]))
            plt.plot(range(1, n / 2), h_Av[0: len(h_Av)], color(i) + marker(i), markeredgecolor=color(i), markersize=3, linewidth=1)
            plt.plot(range(MAX - 1 - n + n / 2 + 1, MAX - 1), t_Av, color(i) + marker(i) + "-", markeredgecolor=color(i), markersize=3, linewidth=1, label='Dim=' + str(n))
            i += 1
        plt.legend(loc=1, prop={'size': 16})
        plt.ylim(-0.02, 0.14)
        plt.plot([1, MAX - 1 + 10], [s, s], "k--")
        plt.plot([MAX - beta, MAX - beta], [0, 0.125], "k--")
        plt.xlim([1, MAX - 1 + 10])
        plt.savefig("sim_ht_in_bkz_" + str(beta) + '.eps')
        plt.close()


def plot_comparison(l_n=[400]):
    """Plot the comparison."""
    index = 0
    N0 = []
    E0 = []
    f = open("r_100")
    for line in f:
        if line == "":
            break
        data = line.split()
        N0 += [int(data[0])]
        E0 += [float(data[1])]
    plt.figure('fig')
    plt.plot(N0, E0, marker(index), label="BKZ")
    index += 1

    N1 = []
    E1 = []
    f = open("0.25_r_180")
    for line in f:
        if line == "":
            break
        data = line.split()
        N1 += [int(data[0])]
        E1 += [float(data[1])]
    plt.figure('fig')
    plt.plot(N1, E1, marker(index), label="BKZ 2.0")
    index += 1

    for n in l_n:
        N = []
        E = []
        D = []
        _D = []
        f = open("sim_stat_" + str(n))
        for line in f:
            if line == "":
                break
            data = line.split()
            N += [int(data[0])]
            E += [float(data[1])]
            D += [float(data[2])]
            _D += [float(data[3])]

        plt.figure('fig')
        plt.plot(N, E, marker(index), label="Simulation")
        index += 1

    plt.figure('fig')
    Gauss = [2 * log(gauss_est(i)) / (i - 1) for i in range(14, 64, 2)]
    plt.plot(range(14, 64, 2), Gauss, marker(index), label="$ln(GH(\\beta)^{\\frac{2}{\\beta-1}})$")
    plt.ylabel("$e(\\beta)$", fontsize=18)
    plt.legend(loc="lower right", prop={'size': 18})
    plt.xticks(range(6, 66, 4))
    plt.xlabel("$\\beta$", fontsize=18)
    plt.savefig("E_sim.eps")


def plot_sum_comparison(l_n=[400]):
    """Plot the comparison on sum."""
    index = 0
    N0 = []
    SH0 = []
    WH0 = []
    WT0 = []
    f = open("sum_head_100")
    for line in f:
        if line == "":
            break
        data = line.split()
        N0 += [int(data[0])]
        SH0 += [float(data[1])]
    f = open("weightsum_head_100")
    for line in f:
        if line == "":
            break
        data = line.split()
        WH0 += [float(data[1])]
    f = open("weightsum_tail_100")
    for line in f:
        if line == "":
            break
        data = line.split()
        WT0 += [float(data[1])]
    plt.figure('fig1')
    plt.plot(N0, SH0, marker(index), label="BKZ")
    plt.figure('fig2')
    plt.plot(N0, WH0, marker(index), label="BKZ")
    plt.figure('fig3')
    plt.plot(N0, WT0, marker(index), label="BKZ")
    index += 1

    N1 = []
    SH1 = []
    WH1 = []
    WT1 = []
    f = open("0.25_sum_head_180")
    for line in f:
        if line == "":
            break
        data = line.split()
        N1 += [int(data[0])]
        SH1 += [float(data[1])]
    f = open("0.25_weightsum_head_180")
    for line in f:
        if line == "":
            break
        data = line.split()
        WH1 += [float(data[1])]
    f = open("0.25_weightsum_tail_180")
    for line in f:
        if line == "":
            break
        data = line.split()
        WT1 += [float(data[1])]
    plt.figure('fig1')
    plt.plot(N1, SH1, marker(index), label="BKZ 2.0")
    plt.figure('fig2')
    plt.plot(N1, WH1, marker(index), label="BKZ 2.0")
    plt.figure('fig3')
    plt.plot(N1, WT1, marker(index), label="BKZ 2.0")
    index += 1

    for n in l_n:
        N = []
        SH = []
        WH = []
        WT = []
        f = open("sim_sum_" + str(n))
        for line in f:
            if line == "":
                break
            data = line.split()
            N += [int(data[0])]
            SH += [float(data[1])]
            WH += [float(data[2])]
            WT += [float(data[3])]

        plt.figure('fig1')
        plt.plot(N, SH, marker(index), label="Simulation")
        plt.figure('fig2')
        plt.plot(N, WH, marker(index), label="Simulation")
        plt.figure('fig3')
        plt.plot(N, WT, marker(index), label="Simulation")
        index += 1

    plt.figure('fig1')
    plt.ylabel("$s^{(h)}(\\beta)$", fontsize=18)
    plt.legend(loc="lower right", prop={'size': 18})
    plt.xticks(range(6, 66, 4))
    plt.xlabel("$\\beta$", fontsize=18)
    plt.savefig("SH_sim.eps")
    plt.figure('fig2')
    plt.ylabel("$w^{(h)}(\\beta)$", fontsize=18)
    plt.xticks(range(6, 66, 4))
    plt.xlabel("$\\beta$", fontsize=18)
    plt.legend(loc="lower right", prop={'size': 18})
    plt.savefig("WH_sim.eps")
    plt.figure('fig3')
    plt.ylabel("$w^{(t)}(\\beta)$", fontsize=18)
    plt.xticks(range(6, 66, 4))
    plt.xlabel("$\\beta$", fontsize=18)
    plt.legend(loc="lower right", prop={'size': 18})
    plt.savefig("WT_sim.eps")


if __name__ == '__main__':
    gen_stat(0.25, [400], range(6, 66, 4))
    plot_ht([400, 300, 200], [30, 62])
    gen_sum(0.25, [400], range(6, 66, 4))
    plot_comparison()
    plot_sum_comparison()
