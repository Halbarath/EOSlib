/*
 * Equation Of State library EOSlib
 *
 * This is a wrapper to be used in Gasoline and ballic
 * It wraps around the Tillotson material library and the ANEOSmaterial library
 * but can be extended with additional material libraries
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
double EOSPofRhoU(EOSMATERIAL *material, double rho, double u); // used in standard Gasoline
double EOSCofRhoU(EOSMATERIAL *material, double rho, double u); // used in standard Gasoline
double EOSIsentropic(EOSMATERIAL *material, double rho1, double u1, double rho2); // used in Gasoline with ISPH

// Inverse functions
double EOSTofRhoU(EOSMATERIAL *material, double rho, double u); // used in ballic
double EOSUofRhoT(EOSMATERIAL *material, double rho, double T); // used in ballic
double EOSRhoofPT(EOSMATERIAL *material, double p, double T); // used in Gasoline using Woolfson correction


// Derivatives
double EOSdPdRho(EOSMATERIAL *material, double rho, double u); // used in ballic
double EOSdPdRho(EOSMATERIAL *material, double rho, double u); // used in ballic
double EOSdUdRho(EOSMATERIAL *material, double rho, double u); // used in ballic

// Boundary condition solver
int EOSSolveBC(EOSMATERIAL *material1, EOSMATERIAL *material2, double rho1, double u1,
 double *prho2, double *pu2); // used in ballic
