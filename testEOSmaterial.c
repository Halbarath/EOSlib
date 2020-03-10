#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "EOSlib.h"

int main(int argc, char *argv[])
{
	double dKpcUnit = -1.0;
	double dMsolUnit = -1.0;
	
	EOSMATERIAL *material;
	
	material = EOSinitMaterial(54, dKpcUnit, dMsolUnit, 0);
	
	double P = EOSPofRhoU(material, 8, 1e12);
	printf("Pressure %.15e\n", P);
	
	EOSfinalizeMaterial(material);
	return 0;
}
