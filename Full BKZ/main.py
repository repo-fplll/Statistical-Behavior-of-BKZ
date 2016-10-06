"""
Run Schnorr & Euchner' BKZ (without abortion).

Requirement: All lattices to be reduced should be stored in the folder 'lattice' and named 'n_i' respectively, where n is the dimension and i is the index.
Also we should set a folder named 'data' to record the results.
"""

from fpylll import *
from simple_bkz import BKZReduction as BKZ
import time


def run(n, a, b, c):
    """
    Run BKZ.

    n: dimension,
    a, b, c represent that the index set of lattices to reduce is range(a, b, c).
    """
    Block = range(2, 32, 2)
    file_list = [open('data/' + str(n) + '_' + str(i), 'a') for i in Block]
    for s in range(a, b, c):
        B = IntegerMatrix.from_file('lattice/' + str(n) + "_" + str(s))
        bkz_B = BKZ(B)
        for i in range(len(Block)):
            beta = Block[i]
            t = time.clock()
            bkz_B.__call__(beta, False)
            print "Time for blocksize = " + str(beta) + ": ", time.clock() - t, "s"
            gso_B = GSO.Mat(B)
            gso_B.update_gso()
            L = [str(gso_B.get_r(j, j)) + "  " for j in range(n)]
            file_list[i].writelines(L)
            file_list[i].write("\n")
run(100, 0, 4, 1)
