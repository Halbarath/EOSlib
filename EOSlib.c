/*
 * Equation Of State library EOSlib
 * EOS material
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "EOSlib.h"

EOSMATERIAL *EOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit, int featureset)
{
	EOSMATERIAL *material;
	material = (EOSMATERIAL *) calloc(1, sizeof(EOSMATERIAL));
	material->iMat = iMat;

	if (iMat>=0 && iMat<50)
	{
		//Tillotson material number
		material->matType = 0;
		material->tillmaterial = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
		if (featureset > 0)
		{
			// do lookup initialization here
		}
	} else if (iMat>=50 && iMat <100)
	{
		//ANEOS material number
		material->matType = 1;
		material->ANEOSmaterial = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);
	}
	
	return material;
}

void EOSfinalizeMaterial(EOSMATERIAL *material)
{
	switch(material->matType)
	{
		case 0:
			tillFinalizeMaterial(material->tillmaterial);
			break;
		case 1:
			ANEOSfinalizeMaterial(material->ANEOSmaterial);
			break;
		default:
		assert(0);
	}

	free(material);
}


double EOSPofRhoU(EOSMATERIAL *material, double rho, double u)
{
	double P = 0;
	
	switch(material->matType)
	{
		case 0:
			P = eosPressure(material->tillmaterial, rho, u);
			break;
		case 1:
			P = ANEOSPofRhoU(material->ANEOSmaterial, rho, u);
			break;
		default:
		assert(0);
	}
	
	return P;
}
