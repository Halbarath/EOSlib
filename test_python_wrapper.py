#!/usr/bin/env python3
import eoslib

dKpcUnit = 2.06701e-13;
dMsolUnit = 4.80438e-08;

fe = eoslib.Material(63, dKpcUnit, dMsolUnit)
P, c, T = fe.PCTofRhoU(20.0, 0.1)
print(P, c, T)

gran = eoslib.Material(1, 0.0, 0.0)
gran.initIsentropicLookup()
u = gran.UofRhoT(2.7,1500)
rho2 = 100
u2 = gran.Isentropic(2.7,u,rho2)
T2 = gran.TofRhoU(rho2,u2)
print(T2)
