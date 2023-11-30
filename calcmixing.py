#!/usr/bin/env python3

from scipy import optimize
import numpy as np
import sys

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

def calcRhoOfPT(iMat, P, T):
    cmd = "./callEOSlibRhoofPT {:} {:} {:}".format(iMat, P, T)
    process = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    out = process.communicate()[0]
    return float(out.split()[1])

def calcRhoMixOfPT(iMatX, iMatY, X, Y, P, T):
    rhoMixInv = X / calcRhoOfPT(iMatX,P,T) + Y / calcRhoOfPT(iMatY,P,T)
    print(1.0 / rhoMixInv)
    return 1.0 / rhoMixInv

if __name__ == '__main__':
    iMatX = int(sys.argv[1])
    iMatY = int(sys.argv[2])
    filename = sys.argv[3]
    X = 0.5
    Y = 1.0 - X
    P_min = 1e6
    P_max = 1e16

    rho_axis = np.logspace(0,2,101)
    T_axis = np.logspace(3,5,101)

    f = open(filename, "w")
    for T in T_axis:
        for rho in rho_axis:
            print("Calculating at rho = {:}, T = {:}".format(rho, T))
            func = lambda P: rho - calcRhoMixOfPT(iMatX, iMatY, X, Y, P, T)
            sol = optimize.root_scalar(func, method='brenth', bracket=(P_min, P_max))
            P = sol.root
            rhoX = calcRhoOfPT(iMatX, P, T)
            rhoY = calcRhoOfPT(iMatY, P, T)
            umix = X * EOSUofRhoT(iMatX, rhoX, T) + Y * EOSUofRhoT(iMatY, rhoY, T)
            f.write("{:} {:} {:} {:}\n".format(rho, T, P/1e10, umix/1e10))
            f.flush()
    f.close()
