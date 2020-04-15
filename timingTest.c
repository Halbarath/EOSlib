#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "EOSlib.h"

int main(int argc, char *argv[])
{
	double dKpcUnit = 2.06701e-13;
	double dMsolUnit = 4.80438e-08;
    
    double MSOLG = 1.99e33;        /* solar mass in grams */
    double GCGS = 6.67e-8;         /* G in cgs */
    double KPCCM = 3.085678e21;    /* kiloparsec in centimeters */
	
	// unit of u is erg/g = cm^2/s^2
	double CodeUnitstoCGSforU = GCGS*dMsolUnit*MSOLG/(dKpcUnit*KPCCM);
	// unit of rho is g/cm^3
	double CodeUnitstoCGSforRho = (dMsolUnit*MSOLG)/pow(dKpcUnit*KPCCM,3.0);
	// unit of p is erg/cm^3=1/(cm*s^2)
	// combine the upper two
	double CodeUnitstoCGSforP = CodeUnitstoCGSforU*CodeUnitstoCGSforRho;
	// unit of c is cm/s
	double CodeUnitstoCGSforC = dKpcUnit*KPCCM*sqrt((CodeUnitstoCGSforRho*GCGS));

	// Version check
	if (EOS_VERSION_MAJOR != 1) {
        fprintf(stderr, "EOS library has the wrong version (%s)\n", EOS_VERSION_TEXT);
        exit(1);
    }
	
	// Initialization
	EOSMATERIAL *material;

	material = EOSinitMaterial(55, dKpcUnit, dMsolUnit, 0);
    EOSinitIsentropicLookup(material, 0);
    
    for (int i = 0; i < 100000; i++)
    {
        //double T = pow(10.0,(double)rand()/RAND_MAX*3.0 + 3.0);
		//double rho = 1/CodeUnitstoCGSforRho*pow(10.0,(double)rand()/RAND_MAX*0.5 + 1.0);
        double T = pow(10.0,(double)rand()/RAND_MAX*9.0 - 3.5);
		double rho = 1/CodeUnitstoCGSforRho*pow(10.0,(double)rand()/RAND_MAX*4.0 - 2.5);
        double u = EOSUofRhoT(material, rho, T);
        double T1 = EOSTofRhoU(material, rho, u);
        double P = EOSPofRhoU(material, rho, u);
        double c = EOSCofRhoU(material, rho, u);
        P = EOSPCofRhoU(material, rho, u, &c);
        double u2 = EOSIsentropic(material, rho, u, rho*0.999999);
        //double rho2 = EOSRhoofPT(material, P, T);
        //double rho3 = EOSRhoofUT(material, u, T);
        //double dpdrho = EOSdPdRho(material, rho, u);
        double dpdu = EOSdPdU(material, rho, u);
        double dudrho = EOSdUdRho(material, rho, u);
        double ucold = EOSUCold(material, rho);
    }
    EOSfinalizeMaterial(material);    
	return 0;
}
