/*
 * Equation Of State library EOSlib
 * EOS material
 *
 */

#include <stdlib.h>
#include <stdio.h>
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
	if (material->tillmaterial != NULL)
	{
		tillFinalizeMaterial(material->tillmaterial);
	}
	if (material->ANEOSmaterial != NULL)
	{
		ANEOSfinalizeMaterial(material->ANEOSmaterial);
	}
	free(material);
}


double EOSPofRhoU(EOSMATERIAL *material, double p, double u)
{
	double P = 0;
	
	if (material->tillmaterial != NULL)
	{
	} else if (material->ANEOSmaterial !
	
	return 1.0;
}