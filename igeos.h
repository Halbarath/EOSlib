/*
 * Header file for igeos.c.
 */
#ifndef IGEOS_HINCLUDED
#define IGEOS_HINCLUDED

/* Version. */
#define IGEOS_VERSION_TEXT    "1.0.1"
#define IGEOS_VERSION_MAJOR   1
#define IGEOS_VERSION_MINOR   0
#define IGEOS_VERSION_PATCH   1

/*
 * Error codes.
 */
#define IGEOS_SUCCESS 0
#define IGEOS_FAIL   -1

typedef struct igeosMat {
    /* Units. */
    double dKpcUnit;
    double dMsolUnit;
    double dGasConst;
    double dErgPerGmUnit;
    double dGmPerCcUnit;
    double dSecUnit;

    int iMat;
    
    /* Constant for ideal gas EOS. */
    double dConstGamma;
    double dMeanMolMass;

    /* The reference density required by ballic. */
    double rho0;

    /* Specific heat capacity. */
    double cv;
} IGEOSMAT; 

IGEOSMAT *igeosInitMat(int iMat, double dConstGamma, double dMeanMolMass, double dKpcUnit, double dMsolUnit);

int igeosFinalizeMat(IGEOSMAT *Mat);

void igeosPrintMat(IGEOSMAT *Mat, FILE *fp);

double igeosPofRhoU(IGEOSMAT *Mat, double rho, double u);
double igeosPofRhoT(IGEOSMAT *Mat, double rho, double T);
double igeosCofRhoU(IGEOSMAT *Mat, double rho, double u);
double igeosPCofRhoU(IGEOSMAT *Mat, double rho, double u, double *c);
double igeosIsentropicU(IGEOSMAT *Mat, double rho1, double u1, double rho2);
double igeosTofRhoU(IGEOSMAT *Mat, double rho, double u);
double igeosUofRhoT(IGEOSMAT *Mat, double rho, double T);
double igeosRhoofPT(IGEOSMAT *Mat, double P, double T);
double igeosRhoofUT(IGEOSMAT *Mat, double U, double T);
double igeosdPdRho(IGEOSMAT *Mat, double rho, double u);
double igeosdPdU(IGEOSMAT *Mat, double rho, double u);

#endif
