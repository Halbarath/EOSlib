#!/usr/bin/env python3


from __future__ import print_function
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize
import numpy as np

from subprocess import PIPE, Popen
import sys
import os


def P_of_rhot(rho, T, input="aneos-dunite-collins-2014.input"):
    cmd = "./calcPressureRhoTMANEOS {:} {:} {:}".format(rho, T, input)

    process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    out = process.communicate()[0]

    P = float(out.split()[5])
    return P

def rho_of_pt(P, T, input="aneos-dunite-collins-2014.input"):
        func = lambda rho: P - P_of_rhot(rho, T, input)

        rho_min = 1e-25
        rho_max = 1e3

        sol = optimize.root_scalar(func, method='bisect', bracket=(rho_min, rho_max))

        return sol.root

rho = 19.0
T= 1.0
P = P_of_rhot(rho, T, input="aneos-iron-stewart-2020.input")

print("rho= {:15.7E} T= {:15.7E} P= {:15.7E}".format(rho, T, P))

rho_inv = rho_of_pt(P, T, input="aneos-iron-stewart-2020.input")


print("rho= {:15.7E} T= {:15.7E} P= {:15.7E} rho_inv= {:15.7E}".format(rho, T, P, rho_inv))


