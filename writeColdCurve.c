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
    if (iMat == EOSIDEALGAS) {
        struct igeosParam param;
        param.dConstGamma = 5.0/3.0;
        param.dMeanMolMass = 1.0;
        material = EOSinitMaterial(iMat, dKpcUnit, dMsolUnit, &param);
    } else {
        material = EOSinitMaterial(iMat, dKpcUnit, dMsolUnit, NULL);
    }

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
