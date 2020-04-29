#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "EOSlib.h"

int main(int argc, char *argv[])
{
	if (argc != 2) {
        fprintf(stderr,"Usage: calcColdcurveEnergy <iMat>\n");
        exit(1);
    }
	
	int iMat = atoi(argv[1]);
	
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
	
	// Version check
	if (EOS_VERSION_MAJOR != 1) {
        fprintf(stderr, "EOS library has the wrong version (%s)\n", EOS_VERSION_TEXT);
        exit(1);
    }
	
	// Initialization
	EOSMATERIAL *material;

	material = EOSinitMaterial(iMat, dKpcUnit, dMsolUnit, 0);
	EOSinitIsentropicLookup(material, 0);
	
    int nRho = 500;
    double rhoMin = 1e-3;
    double rhoMax = 100;

	double deltaLogRho = (log(rhoMax)-log(rhoMin))/(nRho-1);
	for (int i=0; i<nRho; i++)
	{
		double rho = exp(log(rhoMin)+i*deltaLogRho);
        fprintf(stdout, "%.15e %.15e\n", rho, EOSUCold(material,rho));
	}
	
	// Finalize
	EOSfinalizeMaterial(material);
	return 0;
}
