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

#define EOS_VERSION_TEXT    "1.0.0"
#define EOS_VERSION_MAJOR   1
#define EOS_VERSION_MINOR   0
#define EOS_VERSION_PATCH   0

#define EOSIdealGas  0
#define EOSTillotson 1
#define EOSANEOS     2
#define iMatIdealGas 0
#define iMatTillotsonmin 1
#define iMatTillotsonmax 50
#define iMatANEOSmin 51
#define iMatANEOSmax 100

 
typedef struct EOSmaterial
{
	int iMat; // Material number
	int matType; // Material type, 0: Tillotson, 1: ANEOS
	double rho0; // reference density
	TILLMATERIAL *tillmaterial; // Pointer to tillotson material
	ANEOSMATERIAL *ANEOSmaterial; // Pointer to ANEOS material
	
} EOSMATERIAL;

// Initialization and finalization

EOSMATERIAL *EOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit, const void * additional_data);
void EOSfinalizeMaterial(EOSMATERIAL *material);

// Access functions

// Main functions
double EOSPofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSCofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSPCofRhoU(EOSMATERIAL *material, double rho, double u, double *c);
double EOSIsentropic(EOSMATERIAL *material, double rho1, double u1, double rho2);

// Inverse functions
double EOSTofRhoU(EOSMATERIAL *material, double rho, double u);
double EOSUofRhoT(EOSMATERIAL *material, double rho, double T);
double EOSRhoofPT(EOSMATERIAL *material, double p, double T);


// Derivatives
double EOSdPdRho(EOSMATERIAL *material, double rho, double u);
double EOSdPdU(EOSMATERIAL *material, double rho, double u);
double EOSdUdRho(EOSMATERIAL *material, double rho, double u);

// Boundary condition solver
int EOSSolveBC(EOSMATERIAL *material1, EOSMATERIAL *material2, double rho1, double u1,
 double *prho2, double *pu2);
