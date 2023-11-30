#!/usr/bin/env python3

from scipy import optimize
import numpy as np

from subprocess import PIPE, Popen

def EOSPofRhoT(iMat, rho, T):
    cmd = "./callEOSlibforMixing {:} {:} {:}".format(iMat, rho, T)
    process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    out = process.communicate()[0]
    P = float(out.split()[5])
    return P

def EOSUofRhoT(iMat, rho, T):
    cmd = "./callEOSlibforMixing {:} {:} {:}".format(iMat, rho, T)
    process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    out = process.communicate()[0]
    U = float(out.split()[7])
    return U

def EOSPUofRhoT(iMat, rho, T):
    cmd = "./callEOSlibforMixing {:} {:} {:}".format(iMat, rho, T)
    process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    out = process.communicate()[0]
    P = float(out.split()[5])
    U = float(out.split()[7])
    return P, U

def calcRhoOfPT(iMat, P, T, rho_min = 1e-4, rho_max = 1e3):
    func = lambda rho: P - EOSPofRhoT(iMat, rho, T)
    sol = optimize.root_scalar(func, method='bisect', bracket=(rho_min, rho_max))
    return sol.root

if __name__ == '__main__':
    iMat = 52
    rho = 5.0
    T = 300.0

    rho_min = 1e-4
    rho_max = 999.0
    P = EOSPofRhoT(iMat, rho, T)
    print(P)
    rho_new = calcRhoOfPT(iMat, P, T, rho_min, rho_max)
    print(rho_new)
    u = EOSUofRhoT(iMat, rho_new, T)
    print(u)
