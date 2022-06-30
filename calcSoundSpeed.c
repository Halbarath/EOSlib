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
        fprintf(stderr,"Usage: calcSoundSpeed <iMat>\n");
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
    
    char axesFilename[256] = "axes.in";
	
	FILE *fp;
	char str[1000];
	int nRho;
	int nT;

	fp = fopen(axesFilename, "r");
	if (fp == NULL){
        fprintf(stderr,"Could not open file %s",axesFilename);
    }
	if (fgets(str, 1000, fp) != NULL)
	{
	nRho = (int) strtol(str, (char **)NULL, 10);
	}
	if (fgets(str, 1000, fp) != NULL)
	{
	nT = (int) strtol(str, (char **)NULL, 10);
	}
	
	double *rhoAxis = (double *)malloc(nRho * sizeof(double));
	double *TAxis = (double *)malloc(nT * sizeof(double));
	
	for (int i=0; i<nRho; i++)
	{
		if (fgets(str, 1000, fp) != NULL)
		{
			rhoAxis[i] = (double) strtod(str, (char **)NULL);
		}
	}

	for (int i=0; i<nT; i++)
	{
		if (fgets(str, 1000, fp) != NULL)
		{
			TAxis[i] = (double) strtod(str, (char **)NULL);
		}
	}
	
	fclose(fp);
    
    double **csArray = (double **)malloc(sizeof(double*)*nT);
    for (int i=0; i<nT; i++)
	{
		csArray[i] = (double *)malloc(nRho * sizeof(double)); 
	}
	
	// Initialization
	EOSMATERIAL *material;
    if (iMat == MAT_IDEALGAS) {
        struct igeosParam param;
        param.dConstGamma = 5.0/3.0;
        param.dMeanMolMass = 1.0;
        material = EOSinitMaterial(iMat, dKpcUnit, dMsolUnit, &param);
    } else {
        material = EOSinitMaterial(iMat, dKpcUnit, dMsolUnit, NULL);
    }

	EOSinitIsentropicLookup(material, 0);
		
    double cmin = material->minSoundSpeed;
    
    for (int i = 0; i<nT; i++)
	{
		if ( (i%50)==0){
			fprintf(stderr, "Iterating over T for %d of %d\n",i,nT);
		}
		for (int j = 0; j<nRho; j++)
		{
            double T = TAxis[i];
			double rho = rhoAxis[j];
            double u = EOSUofRhoT(material, rho, T);
            double cs = EOSCofRhoU(material, rho, u);
            if (cs < cmin)
            {
                cs = cmin;
            }
            csArray[i][j] = cs;
        }
    }
	
	fp = fopen("soundspeed.txt", "w");

	for (int i = 0; i<nT; i++)
	{
		for (int j = 0; j<nRho; j++)
		{
			fprintf(fp," %.15e", csArray[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	// Finalize
	EOSfinalizeMaterial(material);
	return 0;
}
