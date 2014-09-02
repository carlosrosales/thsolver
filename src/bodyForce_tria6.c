/******************************************************************************
* File      : bodyForceCubat_tria6.c                                          *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the body force term at the required points and returns it in     *
* array vDomain[]. The evaluation is direct in a N[0]xN[1]xN[2] mesh that is  *
* interpolated to a nInterp[0]xnInterp[1]xnInterp[2] mesh. 100 points are     *
* randomly sampled to estimate the error.                                     *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
*                                                                             *
* N[]   : number of points in each direction for electric field calculation.  *
* nInterp[] : number of interpolated values in each direction.                *
* trapSize[0][0] == xmin  |  trapSize[0][1] == xmax                           *
* trapSize[1][0] == ymin  |  trapSize[1][1] == ymax                           *
* trapSize[2][0] == zmin  |  trapSize[2][1] == zmax                           *
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
#include "bodyForce_tria6.h"

int bodyForce_tria6(int *N, int *nInterp, double **trapSize,
                    double **mDomainNodes, unsigned int **mDomainElems, 
                    double *vDomainSolution, double **mNodes, 
                    unsigned int **mElems, double **vMatParam, double *vDomain,
                    unsigned int *vBCType, double ***bi)
{
	long seed;
	int  index[100][3];
	unsigned int bc, i, j, k, m, nx, ny, nz, nxi, nyi, nzi;
	double max, nr, b0, dr, dx, dy, dz, dxi, dyi, dzi, sx, sy, sz;
	double Xin[3], E[6], normal[3];
	double *domainTmp, *x, *y, *z, *xi, *yi, *zi;
	double ***b;

	/* Initialize */
	nx = N[0];
	ny = N[1];
	nz = N[2];
	nxi = nInterp[0];
	nyi = nInterp[1];
	nzi = nInterp[2];
	sx = trapSize[0][1] - trapSize[0][0];
	sy = trapSize[1][1] - trapSize[1][0];
	sz = trapSize[2][1] - trapSize[2][0];
	dx = sx/((double)nx - 1.0);
	dy = sy/((double)ny - 1.0);
	dz = sz/((double)nz - 1.0);
	dxi = sx/((double)nxi - 1.0);
	dyi = sy/((double)nyi - 1.0);
	dzi = sz/((double)nzi - 1.0);
	b0 = sx*sy*sz*vMatParam[1][0]/(vMatParam[1][1]*nxi*nyi*nzi*pi4);
	x = doubleVector(nx,0);
	y = doubleVector(ny,0);
	z = doubleVector(nz,0);
	xi = doubleVector(nxi,0);
	yi = doubleVector(nyi,0);
	zi = doubleVector(nzi,0);
	b = doubleTensor(nx,ny,nz,0);
	domainTmp = doubleVector(nNodes,1);
	
    /* Calculate |E|^2 in the evaluation grid nodes */
	for(k = 0; k < nz; k++){
		Xin[2] = trapSize[2][0]+dz*k;
		for(j = 0; j < ny; j++){
			Xin[1] = trapSize[1][0]+dy*j;
			for(i = 0; i < nx; i++){
				Xin[0] = trapSize[0][0]+dx*i;
				field_tria6(Xin,mDomainNodes,mDomainElems,vDomainSolution,E);
				b[i][j][k] = E[0]*E[0] + E[1]*E[1] + E[2]*E[2] + E[3]*E[3] + 
				E[4]*E[4] + E[5]*E[5];
			}
		}
	}

	/* Interpolate b values to improve integral accuracy */
	for(i = 0; i < nx; i++) x[i] = trapSize[0][0] + i*dx;
	for(i = 0; i < ny; i++) y[i] = trapSize[1][0] + i*dy;
	for(i = 0; i < nz; i++) z[i] = trapSize[2][0] + i*dz;
	for(i = 0; i < nxi; i++) xi[i] = trapSize[0][0] + i*dxi;
	for(i = 0; i < nyi; i++) yi[i] = trapSize[1][0] + i*dyi;
	for(i = 0; i < nzi; i++) zi[i] = trapSize[2][0] + i*dzi;
	interp3d(N,nInterp,x,y,z,b,bi);

	/* Domain integral for each node */
	for(m = 0 ; m < nNodes; m++){
	
		Xin[0] = mNodes[m][0];
		Xin[1] = mNodes[m][1];
		Xin[2] = mNodes[m][2];
		bc = vBCType[m];
		if(bc == 6) bc = 0;
		
		if(bc == 0) getNormal_tria6(m+1,mNodes,mElems,normal);		
		
		for(i = 0; i < nxi; i++){
			dx = Xin[0] - xi[i];
			for(j = 0; j < nyi; j++){
				dy = Xin[1] - yi[j];
				for(k = 0; k < nzi; k++){
					dz = Xin[2] - zi[k];
					dr = sqrt(dx*dx + dy*dy + dz*dz);
					
					if(bc == 1) vDomain[m] -= bi[i][j][k]/dr;
					else{
						nr = dx*normal[0] + dy*normal[1] + dz*normal[2];
						vDomain[m] += bi[i][j][k]*nr/(dr*dr*dr);
					}
				}					
			}
		}
		vDomain[m] *= b0;
		domainTmp[m] *= b0;
	}
	
	free(x);
	free(y);
	free(z);
	free(xi);
	free(yi);
	free(zi);
	free(domainTmp);
	freeDoubleTensor(b,nx,ny);

	return 0;
}
