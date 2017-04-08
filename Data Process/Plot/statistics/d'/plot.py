"""Plot d'(beta)."""

import matplotlib.pyplot as plt


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
f = open("_d_100")
N = []
E = []
for line in f:
    if line == "":
        break
    L = line.split()
    N += [int(L[0])]
    E += [float(L[1])]

"""BKZ 2.0: dim = 180."""
C = [0.25, 2.0, 8.0]
N180 = [[] for c in C]
E180 = [[] for c in C]
for i in range(len(C)):
    f = open("_d_180_" + str(C[i]))
    for line in f:
        if line == "":
            break
        L = line.split()
        N180[i] += [int(L[0])]
        E180[i] += [float(L[1])]

index = 0
plt.figure(figsize=(8, 5))
plt.xlabel("$\\beta$", fontsize=18)
plt.ylabel("$d'(\\beta)$", fontsize=18)
plt.plot(N, E, marker(index), label='BKZ: n=100')

for i in range(len(C)):
    plt.plot(N180[i], E180[i], marker(index), label='BKZ 2.0: n=180, C=' + str(C[i]))
    index += 1

plt.legend(loc='upper right', prop={'size': 14})
plt.ylim([-3.0, 2.0])
plt.xticks(range(2, 66, 4))
plt.savefig('_D.eps')
plt.close()
