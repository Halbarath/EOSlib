/*
 * Equation Of State library EOSlib
 *
 * This is a wrapper to be used in Gasoline and ballic
 * It wraps around the Tillotson material library and the ANEOSmaterial library
 * but can be extended with additional material libraries
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "EOSlib.h"

EOSMATERIAL *EOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit, const void * additional_data)
{
	EOSMATERIAL *material;
	material = (EOSMATERIAL *) calloc(1, sizeof(EOSMATERIAL));
	material->iMat = iMat;

	if (iMat == iMatIdealGas)
	{
		// not implemented
	} else if (iMat>=iMatTillotsonmin && iMat<=iMatTillotsonmax)
	{
		//Tillotson material number
		material->matType = EOSTillotson;
		material->tillmaterial = tillInitMaterial(iMat, dKpcUnit, dMsolUnit);
		if (additional_data != NULL)
		{
			printf("called with non-NULL additional_data\n");
			// do lookup initialization here
		}
		material->rho0 = material->tillmaterial->rho0;
	} else if (iMat>=iMatANEOSmin && iMat <=iMatANEOSmax)
	{
		//ANEOS material number
		material->matType = EOSANEOS;
		material->ANEOSmaterial = ANEOSinitMaterial(iMat, dKpcUnit, dMsolUnit);
		material->rho0 = ANEOSgetRho0(material->ANEOSmaterial);
	}
	
	return material;
}

void EOSfinalizeMaterial(EOSMATERIAL *material)
{
	switch(material->matType)
	{
		case EOSIdealGas:
		// not implemented
			break;
		case EOSTillotson:
			tillFinalizeMaterial(material->tillmaterial);
			break;
		case EOSANEOS:
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
		case EOSIdealGas:
		// not implemented
			break;
		case EOSTillotson:
			P = tillPressure(material->tillmaterial, rho, u);
			break;
		case EOSANEOS:
			P = ANEOSPofRhoU(material->ANEOSmaterial, rho, u);
			break;
		default:
		assert(0);
	}
	
	return P;
}

double EOSCofRhoU(EOSMATERIAL *material, double rho, double u)
{
	double c = 0;
	double P = 0;
	
	switch(material->matType)
	{
		case EOSIdealGas:
		// not implemented
			break;
		case EOSTillotson:
			P = tillPressureSound(material->tillmaterial, rho, u, &c);
			break;
		case EOSANEOS:
			c = ANEOSCofRhoU(material->ANEOSmaterial, rho, u);
			break;
		default:
		assert(0);
	}
	
	return c;
}

double EOSPCofRhoU(EOSMATERIAL *material, double rho, double u, double *c)
{
	double P = 0;
	double T = 0;
	
	switch(material->matType)
	{
		case EOSIdealGas:
		// not implemented
			break;
		case EOSTillotson:
			P = tillPressureSound(material->tillmaterial, rho, u, c);
			break;
		case EOSANEOS:
			T = ANEOSTofRhoU(material->ANEOSmaterial, rho, u);
			P = ANEOSPofRhoT(material->ANEOSmaterial, rho, T);
			*c = ANEOSCofRhoT(material->ANEOSmaterial, rho, T);
			break;
		default:
		assert(0);
	}
	
	return P;
}
