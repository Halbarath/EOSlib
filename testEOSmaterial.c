/*  This file is part of EOSlib.
 *  Copyright (c) 2020-2021 Thomas Meier & Christian Reinhardt
 *
 *  EOSlib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  EOSlib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with EOSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "EOSlib.h"

int main(int argc, char *argv[])
{
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	//double dKpcUnit = -1.0;
	//double dMsolUnit = -1.0;
	
	// Version check
	if (EOS_VERSION_MAJOR != 1) {
        fprintf(stderr, "EOS library has the wrong version (%s)\n", EOS_VERSION_TEXT);
        exit(1);
    }
	
	// Initialization
	EOSMATERIAL *material1;
	EOSMATERIAL *material2;

	material1 = EOSinitMaterial(54, dKpcUnit, dMsolUnit, 0);
	material2 = EOSinitMaterial(6, dKpcUnit, dMsolUnit, 0);

	EOSinitIsentropicLookup(material1, 0);
	EOSinitIsentropicLookup(material2, 0);
	
	printf("minSoundSpeed mat1 %.15e\n", material1->minSoundSpeed);
	printf("minSoundSpeed mat2 %.15e\n", material2->minSoundSpeed);
	
	// Test p(rho,u)
	double rho = 8/material1->ANEOSmaterial->CodeUnitstoCGSforRho;
	double u = 1e12/material1->ANEOSmaterial->CodeUnitstoCGSforU;
	
	double P1 = EOSPofRhoU(material1, rho, u);
	printf("Pressure %.15e\n", P1);
	double P2 = EOSPofRhoU(material2, rho, u);
	printf("Pressure %.15e\n", P2);
	
	// Test c(rho,u) and p_and_c(rho,u)
	double c1 = 0;
	double c2 = 0;
	c1 = EOSCofRhoU(material1, rho, u);
	printf("Sound speed %.15e\n", c1);
	c2 = EOSCofRhoU(material2, rho, u);
	printf("Sound speed %.15e\n", c2);
	
	P1 = EOSPCofRhoU(material1, rho, u, &c1);
	printf("Pressure %.15e\n", P1);
	printf("Sound speed %.15e\n", c1);
	P2 = EOSPCofRhoU(material2, rho, u, &c2);
	printf("Pressure %.15e\n", P2);
	printf("Sound speed %.15e\n", c2);
	
	// Test rho0
	printf("rho0 %.15e\n", material1->rho0);
	printf("rho0 %.15e\n", material2->rho0);
	
	// Test isentropic evolution
	double rhonew = 8.1/material1->ANEOSmaterial->CodeUnitstoCGSforRho;
	double unew1 = EOSIsentropic(material1, rho, u, rhonew);
	double unew2 = EOSIsentropic(material2, rho, u, rhonew);
	printf("old u %.15e, new u %.15e\n", u, unew1);
	printf("old u %.15e, new u %.15e\n", u, unew2);

	// Test T(rho,u)
	double T1 = EOSTofRhoU(material1, rho, u);
	double T2 = EOSTofRhoU(material2, rho, u);
	printf("Temperature %.15e\n", T1);
	printf("Temperature %.15e\n", T2);

	// Test u(rho,T)
	double u1 = EOSUofRhoT(material1, rho, T1);
	double u2 = EOSUofRhoT(material2, rho, T2);
	printf("u %.15e\n", u1);
	printf("u %.15e\n", u2);

	// Test rho(p,T)
	double rho1 = EOSRhoofPT(material1, P1, T1);
	double rho2 = EOSRhoofPT(material2, P2, T2);
	printf("rho %.15e\n", rho1);
	printf("rho %.15e\n", rho2);

	// Test derivatives
	double dPdRho1 = EOSdPdRho(material1, rho, u);
	double dPdRho2 = EOSdPdRho(material2, rho, u);
	double dPdU1 = EOSdPdU(material1, rho, u);
	double dPdU2 = EOSdPdU(material2, rho, u);
	double dUdRho1 = EOSdUdRho(material1, rho, u);
	double dUdRho2 = EOSdUdRho(material2, rho, u);
	printf("derivatives: dPdRho %.15e, dPdU %.15e, dUdRho %.15e\n", dPdRho1, dPdU1, dUdRho1);
	printf("derivatives: dPdRho %.15e, dPdU %.15e, dUdRho %.15e\n", dPdRho2, dPdU2, dUdRho2);

	// Test solveBC
	double rhosolved;
	double usolved;
	int ret = EOSSolveBC(material2, material1, rho, u, &rhosolved, &usolved);
	printf("returnvalue %d\n", ret);
	printf("rhosolved %.15e\n", rhosolved);
	printf("usolved %.15e\n", usolved);
	P1 = EOSPofRhoU(material2, rho, u);
	P2 = EOSPofRhoU(material1, rhosolved, usolved);
	T1 = EOSTofRhoU(material2, rho, u);
	T2 = EOSTofRhoU(material1, rhosolved, usolved);
	printf("T: original %.15e, solved %.15e\n", T1, T2);
	printf("P: original %.15e, solved %.15e\n", P1, P2);

	// Test Woolfson correction
	double woolfsoncoeff = EOSWoolfsonCoeff(material1, material2, P1, T1);
	printf("woolfsoncoeff %.15e\n", woolfsoncoeff);

	// Test EOSisbelowColdCurve
	int ret1 = EOSisbelowColdCurve(material1, 25, 0.1);
	int ret2 = EOSisbelowColdCurve(material1, 25, 60);
	printf("isbelowColdCurve: true = %d, false = %d\n", ret1, ret2);

	ret1 = EOSisbelowColdCurve(material2, 25, 0.1);
	ret2 = EOSisbelowColdCurve(material2, 25, 60);
	printf("isbelowColdCurve: true = %d, false = %d\n", ret1, ret2);

	// Test coldU
	printf("cold u = %g\n",EOSUCold(material1, rho));
	printf("cold u = %g\n",EOSUCold(material2, rho));

	// Test printing
	EOSPrintMat(material1,stderr);
	EOSPrintMat(material2,stderr);

    // Test isInTable
    printf("Test neg density mat1: %d\n",EOSIsInTable(material1,2.7e-4,u1));
    printf("Test neg density mat2: %d\n",EOSIsInTable(material2,0.9e-4,u1));

    printf("Test large density mat1: %d\n",EOSIsInTable(material1,275,u1));
    printf("Test large density mat2: %d\n",EOSIsInTable(material2,200,u1));

    printf("Test small u mat1: %d\n",EOSIsInTable(material1,rho,EOSUCold(material1,rho)*0.9999));
    printf("Test small u mat2: %d\n",EOSIsInTable(material2,rho,EOSUCold(material2,rho)*0.9999));

    printf("Test large u mat1: %d\n",EOSIsInTable(material1,10,100000));
    printf("Test large u mat2: %d\n",EOSIsInTable(material2,10,100000));

    printf("Test ok parameters mat1: %d\n",EOSIsInTable(material1,rho,u));
    printf("Test ok parameters mat2: %d\n",EOSIsInTable(material2,rho,u));

	// Finalize
	EOSfinalizeMaterial(material1);
	EOSfinalizeMaterial(material2);
	return 0;
}
