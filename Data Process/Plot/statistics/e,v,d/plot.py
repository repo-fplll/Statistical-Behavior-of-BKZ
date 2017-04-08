"""Plot e(beta) and d(beta)."""

import matplotlib.pyplot as plt
from math import sqrt, gamma, pi


def gauss_est(n):
    """Gaussian Heuristic Prediction."""
    return (gamma(n / 2. + 1) / (pi**(n / 2.)))**(1. / n)


def stat(l):
    """Mean and standard deviation of all element in a list."""
    ans = 0.0
    ave = sum(l) / len(l)
    for x in l:
        ans += (x - ave)**2
    ans /= len(l)
    return (ave, sqrt(ans))


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

"""BKZ: dim = 100."""
f = open("r_100")
N = []
E = []
Var = []
for line in f:
    if line == "":
        break
    L = line.split()
    N += [int(L[0])]
    E += [float(L[1])]
    Var += [float(L[2])**2]

"""BKZ 2.0: dim = 180."""
C = [0.25, 2.0, 8.0]
N180 = [[] for c in C]
E180 = [[] for c in C]
Var180 = [[] for c in C]
for i in range(len(C)):
    f = open(str(C[i]) + "_r_180")
    for line in f:
        if line == "":
            break
        L = line.split()
        N180[i] += [int(L[0])]
        E180[i] += [float(L[1])]
        Var180[i] += [float(L[2])**2]

index = 0
plt.figure(figsize=(8, 5))
plt.xlabel("$\\beta$", fontsize=18)
plt.ylabel("$e(\\beta)$", fontsize=18)
plt.plot(N, E, marker(index), label='BKZ: n=100')

for i in range(len(C)):
    plt.plot(N180[i], E180[i], marker(index), label='BKZ 2.0: n=180, C=' + str(C[i]))
    index += 1

# Gauss = [2 * log(gauss_est(i)) / (i - 1) for i in range(22, 64, 2)]
# plt.plot(range(22, 64, 2), Gauss, "k-", label="log(GH($\\beta$))*2/($\\beta$-1)")

plt.legend(loc='upper right', prop={'size': 14})
plt.xticks(range(2, 66, 4))
plt.savefig('E.eps')
plt.close()

index = 0
plt.figure(figsize=(8, 5))
plt.xlabel("$\\beta$", fontsize=18)
plt.ylabel("$v(\\beta$)", fontsize=18)
plt.plot(N, Var, marker(index), label='BKZ: n=100')
index += 1

for i in range(len(C)):
    plt.plot(N180[i], Var180[i], marker(index), label='BKZ 2.0: n=180, C=' + str(C[i]))
    index += 1

plt.legend(loc='upper right', prop={'size': 15})
plt.xticks(range(2, 66, 4))
plt.savefig('var.eps')
plt.close()

f = open("d_100")
N = []
D = []
for line in f:
    if line == "":
        break
    L = line.split()
    N += [int(L[0])]
    D += [float(L[1])]

N180 = [[] for c in C]
D180 = [[] for c in C]
for i in range(len(C)):
    f = open(str(C[i]) + "_d_180")
    for line in f:
        if line == "":
            break
        L = line.split()
        N180[i] += [int(L[0])]
        D180[i] += [float(L[1])]

index = 0
plt.figure(figsize=(8, 5))
plt.xlabel("$\\beta$", fontsize=18)
plt.ylabel("$d(\\beta)$", fontsize=18)
plt.plot(N, D, marker(index), label='BKZ: n=100')
index += 1

for i in range(len(C)):
    plt.plot(N180[i], D180[i], marker(index), label='BKZ 2.0: n=180, C=' + str(C[i]))
    index += 1

plt.legend(loc='upper right', prop={'size': 10})
plt.xticks(range(2, 66, 4))
plt.savefig('D.eps')
plt.close()
