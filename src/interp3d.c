/******************************************************************************
* File      : bodyForceAdapt_tria3.c                                          *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* This function takes 3 columns x[nx], y[ny], and z[nz] and produces a new    *
* file	with the nInterp[] interpolated values of the remaining columns       *
* (trilinear interpolation). Vectors x, y, z, must be monotously increasing.  *
* Interpolation is based on a unit cube with the origin at node 1:            *
*	      7--------8                                                          *
*	     /|       /|                                                          *
*	    5------- 6 |                                                          *
*	    | 3------|-4                                                          *
*	    | /      | /                                                          *
*	    1--------2                                                            *
* Subindex i indicates interpolated quantity (xi, nxi, fi, ...) and no        *
* subindex indicates original input (x, nx, f, ...).                          *
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double *doubleVector(unsigned int LENGTH, unsigned int INIT);

int interp3d(int *n, int *nInterp, double *x, double *y, double *z, 
             double ***f, double ***fi)
{
	int i, j, k, m, mx, my, mz;
	unsigned int nx, ny, nz, nxi, nyi, nzi;
	double dx, dy, dz, dxi, dyi, dzi, factor, L1, L2, L3, L1m, L2m, L3m;
	double M[8], ftmp[8];
	double *xi, *yi, *zi;
	
	/* INITIALIZE */
	nx = n[0];
	ny = n[1];
	nz = n[2];
	nxi = nInterp[0];
	nyi = nInterp[1];
	nzi = nInterp[2];
	dx = 1.0/(x[1] - x[0]);
	dy = 1.0/(y[1] - y[0]);
	dz = 1.0/(z[1] - z[0]);
	dxi = (x[nx-1] - x[0])/((double)nxi - 1.0);
	dyi = (y[ny-1] - y[0])/((double)nyi - 1.0);
	dzi = (z[nz-1] - z[0])/((double)nzi - 1.0);
	xi = doubleVector(nxi,0);
	yi = doubleVector(nyi,0);
	zi = doubleVector(nzi,0);
		
	/* GET INTERPOLATION POSITIONS */
	for(i = 0; i < nxi; i++) xi[i] = x[0] + i*dxi;
	for(i = 0; i < nyi; i++) yi[i] = y[0] + i*dyi;
	for(i = 0; i < nzi; i++) zi[i] = z[0] + i*dzi;
	
	/* IDENTIFY CUBE WHERE EACH INTERPOLATION POINT LIES */
	for(i = 0; i < nxi; i++)
	{
		for(m = 0; m < nx-1; m++) if(xi[i] >= x[m] && xi[i] <= x[m+1]) mx = m;
		L1 = (xi[i] - x[mx])*dx;
		L1m = 1.0 - L1;
		for(j = 0; j < nyi; j++)
		{
			for(m = 0; m < ny-1; m++) if(yi[j] >= y[m] && yi[j] <= y[m+1]) my = m;
			L2 = (yi[j] - y[my])*dy;
			L2m = 1.0 - L2;
			for(k = 0; k < nzi; k++)
			{
				for(m = 0; m < nz-1; m++) if(zi[k] >= z[m] && zi[k] <= z[m+1]) mz = m;
				
				/* OBTAIN f AT CORNERS OF CUBE CONTAINING POINT (xi,yi,zi) */
				ftmp[0] = f[mx][my][mz];
				ftmp[1] = f[mx+1][my][mz];
				ftmp[2] = f[mx][my+1][mz];
				ftmp[3] = f[mx+1][my+1][mz];
				ftmp[4] = f[mx][my][mz+1];
				ftmp[5] = f[mx+1][my][mz+1];
				ftmp[6] = f[mx][my+1][mz+1];
				ftmp[7] = f[mx+1][my+1][mz+1];
							
				/* NORMALIZE (xi,yi,zi) WITH RESPECT TO CUBE */
				L3 = (zi[k] - z[mz])*dz;
				L3m = 1.0 - L3;
	
				/* GET SHAPE FUNCTIONS @ (xtmp,ytmp,ztmp) */
				M[0] = L1m*L2m*L3m;
				M[1] = L1*L2m*L3m;
				M[2] = L1m*L2*L3m;
				M[3] = L1*L2*L3m;
				M[4] = L1m*L2m*L3;
				M[5] = L1*L2m*L3;
				M[6] = L1m*L2*L3;
				M[7] = L1*L2*L3;
											
				fi[i][j][k] = M[0]*ftmp[0] + M[1]*ftmp[1] + M[2]*ftmp[2] + M[3]*ftmp[3] + M[4]*ftmp[4] + M[5]*ftmp[5] + M[6]*ftmp[6] +M[7]*ftmp[7];
			}
		}
	}
		
free(xi);
free(yi);
free(zi);

return 0;
}
