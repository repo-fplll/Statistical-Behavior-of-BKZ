"""

BKZ 2.0.

strategy "95.json".
fixed tour number=c*(n/b)^2*(log n + loglogmax(||b_i*||/vol^{1/n})).
enumeration radius = 0.99 ||b_i*|| for each local block L_[i,j].

Requirement: All lattices to be reduced should be stored in the folder 'lattice' and named 'n_i' respectively, where n is the dimension and i is the index.
Also we should set folder named 'basis' and 'data' to record the results.
"""

import argparse
import time
from math import log
# from math import ceil, sqrt, log, exp, gamma, pi

# import argparse
# from random import uniform, randint
# import sys
# import numpy as np
# from copy import copy

from fpylll import IntegerMatrix, GSO
from fpylll.algorithms.bkz2 import BKZReduction as BKZ2
from fpylll import BKZ as BKZ
# from fpylll.algorithms.simple_bkz import BKZReduction as SimpleBKZ
# from fpylll.algorithms.simple_dbkz import DBKZReduction as SimpleDualBKZ
# from fpylll.util import set_random_seed


def parse_commandline():
    """
    Parse command line arguments and check for consistency.

    returns: parameters
    rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description='BKZ 2.0')
    parser.add_argument('-n', type=int, help="Dimension")
    parser.add_argument('-c', type=float, help="C")
    parser.add_argument('-l', type=int, help="Range_Lower")
    parser.add_argument('-u', type=int, help="Range_Upper")
    parser.add_argument('-s', type=int, help="Range_Step")
    args = parser.parse_args()
    return args


def run(n, c, a, b, d):
    """
    Run BKZ2.0.

    n: dimension,
    c: tour number parameter,
    a,b,d represents that the index set of lattices to reduce is range(a,b,d).
    """
    Block = range(2, 66, 4)
    for beta in Block:
        for s in range(a, b, d):
            if beta == 2:
                B = IntegerMatrix.from_file('lattice/' + str(n) + '_' + str(s))
            else:
                B = IntegerMatrix.from_file('basis/C=' + str(c) + '/' + str(n) + '_bkz_' + str(beta - 4) + '_' + str(s))
            """Compute tour number."""
            Y = 0
            gso_B = GSO.Mat(B, float_type="ld")
            gso_B.update_gso()
            if beta == 2:
                Y = log(10 * n**2) / log(2)
            else:
                L = [0.5 * log(gso_B.get_r(j, j)) / log(2.) for j in range(n)]
                vol = sum(L)
                max_gso = max(L) * n
                Y = max_gso - vol
                Y = log(Y) / log(2.)

            t = time.clock()
            X = c * ((n**2) / beta**2) * Y + 1
            """BKZ 2.0."""
            param = BKZ.Param(block_size=beta, min_success_probability=0.95,
                              strategies="95.json", flags=BKZ.MAX_LOOPS, max_loops=X)
            if beta <= 6:
                BKZ2(gso_B)(params=param)
            else:
                BKZ2(B)(params=param)
            gso_B = GSO.Mat(B)
            gso_B.update_gso()
            L = [str(gso_B.get_r(j, j)) + "  " for j in range(n)]
            data_file = open('data/C=' + str(c) + '/' + str(a) + 'to' + str(b) + '_' + str(n) + '_' + str(beta), 'a')
            data_file.writelines(L)
            data_file.write("\n")
            data_file.close()
            print "C = " + str(c) + "  Time for blocksize = " + str(beta) + ": ", time.clock() - t, "s"
            file = open('basis/C=' + str(c) + '/' + str(n) + '_bkz_' + str(beta) + '_' + str(s), 'w')
            file.write('[')
            file.writelines(str(B))
            file.write(']\n')
            file.close()


def main():
    """Main Function."""
    args = parse_commandline()
    n, c, a, b, d = args.n, args.c, args.l, args.u, args.s
    run(n, c, a, b, d)

if __name__ == '__main__':
    main()
