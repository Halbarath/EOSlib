#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "EOSlib.h"

int main(int argc, char *argv[])
{
	if (argc != 3) {
        fprintf(stderr,"Usage: calcColdcurveEnergy <iMat> <rho>\n");
        exit(1);
    }
	
	int iMat = atoi(argv[1]);
	double rho = atof(argv[2]);
	
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
	
	double u = EOSUofRhoT(material, rho, 1e-4);
	
	fprintf(stdout, "Energy on cold curve: %.15e\n", u);
	
	// Finalize
	EOSfinalizeMaterial(material);
	return 0;
}
