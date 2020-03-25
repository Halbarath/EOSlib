#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "EOSlib.h"

int main(int argc, char *argv[])
{
	if (argc != 4) {
        fprintf(stderr,"Usage: calcTemperature <iMat> <rho> <u>\n");
        exit(1);
    }
	
	int iMat = atoi(argv[1]);
	double rho = atof(argv[2]);
	double u = atof(argv[3]);
	
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	
	// Version check
	if (EOS_VERSION_MAJOR != 1) {
        fprintf(stderr, "EOS library has the wrong version (%s)\n", EOS_VERSION_TEXT);
        exit(1);
    }
	
	// Initialization
	EOSMATERIAL *material;

	int flag = 1;

	material = EOSinitMaterial(iMat, dKpcUnit, dMsolUnit);
	EOSinitIsentropicLookup(material, &flag);
	
	double T = EOSTofRhoU(material, rho, u);
	
	fprintf(stdout, "Temperature: %.15e\n", T);
	
	// Finalize
	EOSfinalizeMaterial(material);
	return 0;
}
