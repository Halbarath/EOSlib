/*
 * Equation Of State library EOSlib
 *
 * This is a wrapper to be used in Gasoline and ballic
 * It wraps around the Tillotson material library and the ANEOSmaterial library
 * but can be extended with additional material libraries
 *
 */

#ifndef EOSLIB_HINCLUDED
#define EOSLIB_HINCLUDED
 
#include "../tillotson/tillotson.h"
#include "../ANEOSmaterial/ANEOSmaterial.h"

#define EOS_VERSION_TEXT    "1.0.0"
#define EOS_VERSION_MAJOR   1
#define EOS_VERSION_MINOR   0
#define EOS_VERSION_PATCH   0

#define WOOLFSON_MIN_PRESSURE 1e-3

#define EOS_N_MATERIAL_MAX 100

#define EOSIDEALGAS  0
#define EOSTILLOTSON 1
#define EOSANEOS     2
#define iMATIDEALGAS 0
#define iMATTILLOTSONMIN 1
#define iMATTILLOTSONMAX 50
#define iMATANEOSMIN 51
#define iMATANEOSMAX 100

#define EOS_TRUE 1
#define EOS_FALSE 0

#define EOS_SUCCESS 0
#define EOS_FAIL 1
 
typedef struct EOSmaterial
{
	int iMat; // Material number
	int matType; // Material type, 0: Tillotson, 1: ANEOS
	double rho0; // reference density
	int bEntropyTableInit; // flag to signal if the entropy table is initialized
	double minSoundSpeed; // sound speed at reference values
    char MatString[256];
	TILLMATERIAL *tillmaterial; // Pointer to tillotson material
	ANEOSMATERIAL *ANEOSmaterial; // Pointer to ANEOS material
	
} EOSMATERIAL;

// Initialization and finalization

EOSMATERIAL *EOSinitMaterial(int iMat, double dKpcUnit, double dMsolUnit, const void * additional_data);
void EOSinitIsentropicLookup(EOSMATERIAL *material, const void * additional_data);
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

double EOSRhoofUT(EOSMATERIAL *material, double u, double T);

int EOSisbelowColdCurve(EOSMATERIAL *material, double rho, double u);
int EOSIsInTable(EOSMATERIAL *material, double rho, double u);

// Derivatives
double EOSdPdRho(EOSMATERIAL *material, double rho, double u);
double EOSdPdU(EOSMATERIAL *material, double rho, double u);
double EOSdUdRho(EOSMATERIAL *material, double rho, double u);

double EOSUCold(EOSMATERIAL *material, double rho);

// Boundary condition solver
int EOSSolveBC(EOSMATERIAL *material1, EOSMATERIAL *material2, double rho1, double u1,
 double *prho2, double *pu2);
 
// Woolfson correction
double EOSWoolfsonCoeff(EOSMATERIAL *material1, EOSMATERIAL *material2, double P, double T);
#endif
