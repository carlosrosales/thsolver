/******************************************************************************
* File      : thermalPostProcessColCubat_tria6.c                              *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the temperature T at the required point Xin for the case using   *
* the cubature solver. Avoid points in the original surfaces where the field  *
* was calculated, as no consideration for the singular case is made. Excludes *
* points within the specifiec nColumns from the calculation.                  *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
*                                                                             *
* domainSize    : number of nodes in the domain integral evaluation           *
* nColumns      : number of columns to avoid in the calculations              *
* xColumns[][3] : positions of the columns. xColumns[0][2] = y for column 1   *
* integralNodes : array of structures with the info on the integral eval nodes*
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* This file is part of thSolver.                                              *
*                                                                             *
* thSolver is free software; you can redistribute it and/or modify            *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* thSolver is distributed in the hope that it will be useful,                 *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with thSolver; if not, write to the Free Software                     *
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
******************************************************************************/

#include "constants.h"
#include "thermalPostProcessColCubat_tria6.h"

int thermalPostProcessColCubat_tria6(double **trapSize, double **Xinner, 
                                     double **mNodes, unsigned int **mElems, 
                                     double *vB, double **cellNodes, 
                                     unsigned int **cells, double **vMatParam, 
                                     unsigned int nColumns, double **xColumns, 
                                     double *b)
{
	FILE   *fT;
	unsigned int count, i, j, k, m, Tcount;
	double dx, dy, dz, size, sx, sy, sz, T, Tdev, Tmax, Tmean, Tmin;
	double Xin[3];
	double *Tstats;

	/* Initialize */
	T = 0.0;
	sx = trapSize[0][1] - trapSize[0][0];
	sy = trapSize[1][1] - trapSize[1][0];
	sz = trapSize[2][1] - trapSize[2][0];
	fprintf(file_log,"\n\tthermalPostProcessColCubat_tria6(): ");
    fprintf(file_log,"doing temperature calculation...");
	fT = fopen("temperature.dat","w");
	fprintf(fT,"x y z T\n");

	/* Use NTxNTxNT points to get T statistics */
	dx = sx/((double)NT - 1.0);
	dy = sy/((double)NT - 1.0);
	dz = sz/((double)NT - 1.0);
	size = NT*NT*NT;
	Tstats = doubleVector(size,1);

	/* Temperature at required points */
	for(i = 0; i < nInternalPoints; i++){
	
	   	/* Evaluation point's coordinates */
	   	Xin[0] = Xinner[i][0];
	   	Xin[1] = Xinner[i][1];
	   	Xin[2] = Xinner[i][2];
	   	
   		count = 0;
   		for(j = 0; j < nColumns; j++){
			if(nColumnType == 1){
				if(sqrt(pow(Xin[0] - xColumns[j][0],2.0) + 
				pow(Xin[1] - xColumns[j][1],2.0)) > xColumns[j][2]) count++;
			}
			else{
				if((fabs(Xin[0] - xColumns[j][0])> xColumns[j][2]) && 
				(fabs(Xin[1] - xColumns[j][1]) > xColumns[j][2])) count++;
			}
   		}
   		if(count == nColumns){
			temperatureCubat_tria6(Xin,mNodes,mElems,vB,cellNodes,cells,b,&T);
			fprintf(fT,"%le\t%le\t%le\t%le\n",Xin[0],Xin[1],Xin[2],T);
		}
		else fprintf(fT,"%le\t%le\t%le\t%le\n",Xin[0],Xin[1],Xin[2],T0);
	}
	fclose(fT);

	/* Temperature statistics */
	Tmax = -T0;
	Tmin = 1.0E+04;
	Tcount = 0;
	Tmean = 0.0;
	for(i = 0; i < NT; i++){
		Xin[0] = trapSize[0][0] + i*dx;
		for(j = 0; j < NT; j++){
			Xin[1] = trapSize[1][0] + j*dy;
			for(k = 0; k < NT; k++){
				Xin[2] = trapSize[2][0] + k*dz;
				
				count = 0;
				for(m = 0; m < nColumns; m++){
	   				if(nColumnType == 1){
				        if(sqrt(pow(Xin[0] - xColumns[m][0],2.0) + 
				        pow(Xin[1] - xColumns[m][1],2.0)) > xColumns[m][2]) count++;
			        }
			        else{
				        if((fabs(Xin[0] - xColumns[m][0])> xColumns[m][2]) && 
				        (fabs(Xin[1] - xColumns[m][1]) > xColumns[m][2])) count++;
			        }
   				}
   				if(count == nColumns){
					temperatureCubat_tria6(Xin,mNodes,mElems,vB,cellNodes,
					cells,b,&T);
					Tstats[Tcount] = T;
					Tmean += T;
					Tcount++;
					if(T > Tmax) Tmax = T;
					else if(T < Tmin) Tmin = T;
				}
			}
		}
	}
	Tmean = Tmean/Tcount;
	for(i = 0; i < Tcount; i++) Tdev += pow(Tstats[i] - Tmean,2.0);
	Tdev = sqrt(Tdev/(Tcount - 1.0));
	fT = fopen("tempStats.dat","w");	
	fprintf(fT,"Tmax = %le\nTmin = %le\n",Tmax,Tmin);
	fprintf(fT,"Tmean = %le\nTdev = %le\n",Tmean,Tdev);
	fclose(fT);
	
	free(Tstats);

	return 0;
}
