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
	
	
	EOSMATERIAL *material1;
	EOSMATERIAL *material2;
	
	material1 = EOSinitMaterial(54, dKpcUnit, dMsolUnit, 0);
	material2 = EOSinitMaterial(6, dKpcUnit, dMsolUnit, 0);
	
	double rho = 8/material1->ANEOSmaterial->CodeUnitstoCGSforRho;
	double u = 1e12/material1->ANEOSmaterial->CodeUnitstoCGSforU;
	
	double P1 = EOSPofRhoU(material1, rho, u);
	printf("Pressure %.15e\n", P1);
	double P2 = EOSPofRhoU(material2, rho, u);
	printf("Pressure %.15e\n", P2);
	
	
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
	
	
	printf("rho0 %.15e\n", material1->rho0);
	printf("rho0 %.15e\n", material2->rho0);
	
	printf("conversion %.15e\n", material1->ANEOSmaterial->CodeUnitstoCGSforC);
	EOSfinalizeMaterial(material1);
	EOSfinalizeMaterial(material2);
	return 0;
}
