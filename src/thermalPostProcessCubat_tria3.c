/******************************************************************************
* File      : thermalPostProcessCubat_tria3.c                                 *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the temperature T at the required point Xin for the case using   *
* the cubature solver. Avoid points in the original surfaces where the field  *
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

int thermalPostProcessCubat_tria3(double **trapSize, double **Xinner, 
                                  double **mNodes, unsigned int **mElems, 
                                  double *vB, double **cellNodes,
                                  unsigned int **cells, double **vMatParam, 
                                  double *b)
{
    FILE    *fT;
    unsigned int i, j, k, m, Tcount;
    double  dx, dy, dz, size, sx, sy, sz, T, Tdev, Tmax, Tmean, Tmin;
    double  Xin[3];
    double  *Tstats;

    /* Initialize */
    T = 0.0;
    sx = trapSize[0][1] - trapSize[0][0];
    sy = trapSize[1][1] - trapSize[1][0];
    sz = trapSize[2][1] - trapSize[2][0];
    fprintf(file_log,"\n\tthermalPostProcessCubat_tria3(): ");
    fprintf(file_log,"doing temperature calculation...");
    fT = fopen("temperature.dat","w");
    fprintf(fT,"x\ty\tz\tT\n");
        
    /* Use NTxNTxNT points to get T statistics */
    dx = sx/((double)NT - 1.0);
    dy = sy/((double)NT - 1.0);
    dz = sz/((double)NT - 1.0);
    size = NT*NT*NT;
    Tstats = doubleVector(size,1);

    /* Temperature at required points */
    for(i = 0; i < nInternalPoints; i++){     
        Xin[0] = Xinner[i][0];
        Xin[1] = Xinner[i][1];
        Xin[2] = Xinner[i][2];        
        temperatureCubat_tria3(Xin,mNodes,mElems,vB,cellNodes,cells,b,&T);
        fprintf(fT,"%le\t%le\t%le\t%le\n",Xin[0],Xin[1],Xin[2],T);
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
                temperatureCubat_tria3(Xin,mNodes,mElems,vB,cellNodes,cells,
                b,&T);
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
    
    free(Tstats);

    return 0;
}
