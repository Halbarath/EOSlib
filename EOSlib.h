/*
 * Equation Of State library EOSlib
 * Header of EOSlib.c
 *
 */
 
#include "../tillotson/tillotson.h"
#include "../ANEOSmaterial/ANEOSmaterial.h"
 
typedef struct EOSmaterial
{
	int iMat; // Material number
	int matType; // Material type, 0: Tillotson, 1: ANEOS
	TILLMATERIAL *tillmaterial; // Pointer to tillotson material
	ANEOSMATERIAL *ANEOSmaterial; // Pointer to ANEOS material
	
} EOSMATERIAL;

// Initialization and finalization

EOSMATERIAL *EOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit, int featureset);
void EOSfinalizeMaterial(EOSMATERIAL *material);

// Access functions

// Main functions
double EOSPofRhoU(EOSMATERIAL *material, double p, double u);
double EOSCofRhoU(EOSMATERIAL *material, double p, double u);
double EOSIsentropic(EOSMATERIAL *material, double rho1, double u1, double rho2);

// T<-->U conversion functions
double EOSTofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSUofRhoT(EOSMATERIAL *material, double rho, double T);

// Derivatives
double EOSdPdRho(EOSMATERIAL *material, double rho, double u);
double EOSdPdRho(EOSMATERIAL *material, double rho, double u);
double EOSdUdRho(EOSMATERIAL *material, double rho, double u);

// Boundary condition solver
int EOSSolveBC(EOSMATERIAL *material1, EOSMATERIAL *material2, double rho1, double u1, double *prho2, double *pu2);
