"""
Calculate (r_1,...r_{n-1}) and their statistics from (|b_1*|^2,...,|b_n*|^2).

Input file locates in "profile"
Output file locates in "dist_r" and "stat_r"
The data of r_j stored in the file named 'n_b_j' in "dist_r", where n is the dimension of lattices and b denotes the blocksize
The file 'n_b' in "stat_r" stored the mean and std of r_1,...,r_{n-1}
"""


from math import sqrt, log


def stat(l):
    """Return mean and standard deviation of elements in list l."""
    ans = 0.0
    ave = sum(l) / len(l)
    for x in l:
        ans += (x - ave)**2
    ans /= len(l)
    return (ave, sqrt(ans))


def gen_r(l_n, l_beta):
    """Return r_i from profile file."""
    for n in l_n:
        for beta in l_beta:
            f = open("profile/" + str(n) + "_" + str(beta))
            L = [open("dist_r/" + str(n) + "_" + str(beta) + "_" + str(i), "a") for i in range(1, n)]
            for line in f:
                if line == "":
                    break
                gs = line.split()
                for i in range(1, n):
                    L[i - 1].write(str((log(float(gs[i - 1])) - log(float(gs[i]))) / 2.))
                    L[i - 1].write("\n")


def gen_stat(l_n, l_beta):
    """Return the statistics of r_i's."""
    for n in l_n:
        for beta in l_beta:
            f = open("stat_r/" + str(n) + "_" + str(beta), "w")
            for i in range(1, n):
                data = open("dist_r/" + str(n) + "_" + str(beta) + "_" + str(i))
                L = []
                for line in data:
                    if line == "":
                        break
                    L += [float(line)]
                av, std = stat(L)
                f.write(str(av) + " " + str(std) + "\n")

n = [120]
block = range(2, 4, 2)
gen_r(n, block)
gen_stat(n, block)
