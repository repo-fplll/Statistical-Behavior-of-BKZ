"""Calculate statistics of r_i's."""


from math import sqrt


def stat(l):
    """Return mean and standard deviation of elements in list l."""
    ans = 0.0
    ave = sum(l) / len(l)
    for x in l:
        ans += (x - ave)**2
    ans /= len(l)
    return (ave, sqrt(ans))


def middle_r(l_n, h, l_beta):
    """Calculate e and v from middle r_i's, i.e. max(h, beta) < i < n - max(h, beta)."""
    for n in l_n:
        f = open("r_" + str(n), "w")
        Av = []
        Std = []
        for beta in l_beta:
            L = []
            for i in range(max(h, beta) + 1, n - max(h, beta)):
                data = open("dist_r/" + str(n) + "_" + str(beta) + "_" + str(i))
                for line in data:
                    if line == "":
                        break
                    L += [float(line)]
            av, std = stat(L)
            Av += [av]
            Std += [std * std]
            f.write(str(beta) + " " + str(av) + " " + str(std) + "\n")


def sum_head(l_n, h, l_beta):
    """Calculate sh. h(beta) = max(h, beta)."""
    for n in l_n:
        f = open("sum_head_" + str(n), "w")
        Av = []
        for beta in l_beta:
            L = []
            for i in range(1, max(h, beta) + 1):
                data = open("dist_r/" + str(n) + "_" + str(beta) + "_" + str(i))
                for line in data:
                    if line == "":
                        break
                    L += [float(line)]
            av, std = stat(L)
            Av += [max(h, beta) * av]
            f.write(str(beta) + " " + str(max(h, beta) * av) + "\n")


def weightsum_head(l_n, h, l_beta):
    """Calculate wh. h(beta) = max(h, beta)."""
    for n in l_n:
        f = open("weightsum_head_" + str(n), "w")
        Av = []
        for beta in l_beta:
            L = []
            for i in range(1, max(h, beta) + 1):
                data = open("dist_r/" + str(n) + "_" + str(beta) + "_" + str(i))
                for line in data:
                    if line == "":
                        break
                    L += [float(line) * i * 0.5]
            av, std = stat(L)
            Av += [max(h, beta) * av]
            f.write(str(beta) + " " + str(max(h, beta) * av) + "\n")


def weightsum_tail(l_n, h, l_beta):
    """Calculate wt. t(beta) = max(h, beta)."""
    for n in l_n:
        f = open("weightsum_tail_" + str(n), "w")
        Av = []
        for beta in l_beta:
            L = []
            for i in range(n - max(h, beta), n):
                data = open("dist_r/" + str(n) + "_" + str(beta) + "_" + str(i))
                for line in data:
                    if line == "":
                        break
                    L += [float(line) * (n - i) * 0.5]
            av, std = stat(L)
            Av += [max(h, beta) * av]
            f.write(str(beta) + " " + str(max(h, beta) * av) + "\n")


def d(l_n, h, l_beta):
    """Calculate d. h(beta) = max(h, beta)."""
    for n in l_n:
        f = open("d_" + str(n), "w")
        for beta in l_beta:
            E = []
            for i in range(max(h, beta) + 1, n - max(h, beta)):
                data = open("dist_r/" + str(n) + "_" + str(beta) + "_" + str(i))
                for line in data:
                    if line == "":
                        break
                    E += [float(line)]
            e, v = stat(E)
            L = []
            for i in range(1, max(h, beta) + 1):
                data = open("dist_r/" + str(n) + "_" + str(beta) + "_" + str(i))
                for line in data:
                    if line == "":
                        break
                    L += [float(line)]
            av, std = stat(L)
            f.write(str(beta) + " " + str(max(h, beta) * av - (max(h, beta) + 0.5) * e) + "\n")


def _d(l_n, h, l_beta):
    """Calculate d'. h(beta) t(beta) = max(h, beta)."""
    for n in N:
        f = open("_d_" + str(n), "w")
        for beta in Beta:
            fin = open("stat_r/" + str(n) + "_" + str(beta))
            R = []
            for line in fin:
                if line == "":
                    break
                data = line.split()
                R += [float(data[0])]
            e = sum(R[max(h, beta):n - max(h, beta) - 1]) / len(R[max(h, beta):n - max(h, beta) - 1])
            term = 0
            for i in range(1, max(h, beta) + 1):
                term += i * (R[i - 1]) / 2
                term += i * (R[n - 1 - i]) / 2
            term -= ((max(h, beta) * (max(h, beta) + 1)) / 2) * e
            f.write(str(beta) + " " + str(term) + "\n")

N = [100]
Beta = [2]
middle_r(N, 15, Beta)
sum_head(N, 15, Beta)
weightsum_head(N, 15, Beta)
weightsum_tail(N, 15, Beta)
d(N, 15, Beta)
_d(N, 15, Beta)
