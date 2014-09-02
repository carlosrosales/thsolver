/******************************************************************************
* File      : thermalPostProcess_tria3.c                                      *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the temperature T at the required point Xin for the case using   *
* the direct solver. Avoid points in the original surfaces where the field    *
* was calculated, as no consideration for the singular case is made.          *
* Works for linear interpolation in triangular elements (3-noded triangles).  *
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
#include "thermalPostProcessCubat_tria3.h"

int thermalPostProcess_tria3(int *nInterp, double **TrapSize, double **Xinner, 
                             double **mNodes, unsigned int **mElems, 
                             double *vB, double **vMatParam, double ***bi)
{
	FILE   *fT;
	unsigned int i, j, k, m, nxi, nyi, nzi, Tcount;
	double b0, dxi, dyi, dzi, dxStats, dyStats, dzStats, size, sx, sy, sz, T, 
	       Tdev, Tmax, Tmean, Tmin;
	double Xin[3];
	double *xi, *yi, *zi, *Tstats;

	/* Initialize */
	T = 0.0;
	nxi = nInterp[0];
	nyi = nInterp[1];
	nzi = nInterp[2];
	sx = TrapSize[0][1] - TrapSize[0][0];
	sy = TrapSize[1][1] - TrapSize[1][0];
	sz = TrapSize[2][1] - TrapSize[2][0];
	b0 = sx*sy*sz*vMatParam[1][0]/(vMatParam[1][1]*nxi*nyi*nzi*pi4);
	dxi = sx/((double)nxi - 1.0);
	dyi = sy/((double)nyi - 1.0);
	dzi = sz/((double)nzi - 1.0);
	xi = doubleVector(nxi,0);
	yi = doubleVector(nyi,0);
	zi = doubleVector(nzi,0);
	for(i = 0; i < nxi; i++) xi[i] = TrapSize[0][0] + i*dxi;
	for(i = 0; i < nyi; i++) yi[i] = TrapSize[1][0] + i*dyi;
	for(i = 0; i < nzi; i++) zi[i] = TrapSize[2][0] + i*dzi;
	fprintf(file_log,"\n\tthermalPostProcess_tria3(): ");
	fprintf(file_log,"doing temperature calculation...");
	fT = fopen("temperature.dat","w");
	fprintf(fT,"x y z T\n");
	
	/* Use NTxNTxNT points to get temperature statistics */
	dxStats = sx/((double)NT - 1.0);
	dyStats = sy/((double)NT - 1.0);
	dzStats = sz/((double)NT - 1.0);
	size = NT*NT*NT;
	Tstats = doubleVector(size,1);
	

	/* Calculate temperature at the required points */
	for(i = 0; i < nInternalPoints; i++){
	   	Xin[0] = Xinner[i][0];
	   	Xin[1] = Xinner[i][1];
	   	Xin[2] = Xinner[i][2];
		temperature_tria3(nInterp,b0,Xin,mNodes,mElems,vB,xi,yi,zi,bi,&T);
		fprintf(fT,"%le	%le	%le	%le\n",Xin[0],Xin[1],Xin[2],T);
	}
	fclose(fT);
	
	/* Get temperature statistics */
	Tmax = -T0;
	Tmin = 1.0E+04;
	Tcount = 0;
	Tmean = 0.0;
	for(i = 0; i < NT; i++){
		Xin[0] = TrapSize[0][0] + i*dxStats;
		for(j = 0; j < NT; j++){
			Xin[1] = TrapSize[1][0] + j*dyStats;
			for(k = 0; k < NT; k++){
				Xin[2] = TrapSize[2][0] + k*dzStats;	
				temperature_tria3(nInterp,b0,Xin,mNodes,mElems,vB,xi,yi,zi,bi,&T);
				Tstats[Tcount] = T;
				Tmean += T;
				Tcount++;
				if(T > Tmax) Tmax = T;
				else if(T < Tmin) Tmin = T;
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
	
	free(xi);
	free(yi);
	free(zi);
	free(Tstats);

	return 0;
}
