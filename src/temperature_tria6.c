/******************************************************************************
* File      : temperatureAdaptive_tria6.c                                     *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the temperature T at the required point Xin. Avoid points in the *
* original surfaces where the field was calculated, as no consideration for   *
* the singular case is made.                                                  *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
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
#include "temperature_tria6.h"

int temperature_tria6(int *nInterp, double b0, double *Xin, double **mNodes, 
    unsigned int **mElems, double *vB, double *xi, double *yi, double *zi, 
    double ***bi, double *T)
{
	const unsigned int ELEMS = nElems, NODES_IN_ELEM = 6;
	unsigned int currentNode, i, j, k, SinNode; 
	double A, dr, dx, dy, dz, Tdomain, Ttmp;
	double normal[3], SubT[6];
	double X[6][3];
	int nxi, nyi, nzi;

	/* Initialize */
	A = 1.0/pi4;
	Ttmp = 0.0;
	Tdomain = 0.0;
	nxi = nInterp[0];
	nyi = nInterp[1];
	nzi = nInterp[2];

	for(i = 0; i < ELEMS; i++){
	
	   	/* Coordinates of this element's nodes */
	   	for(j = 0; j < NODES_IN_ELEM; j++){
		   	currentNode = mElems[i][j] - 1;
			X[j][0] = mNodes[currentNode][0];
			X[j][1] = mNodes[currentNode][1];
			X[j][2] = mNodes[currentNode][2];
   		}

   		/* Test for singular case */
   		for(j = 0; j < NODES_IN_ELEM; j++){
			dx = X[j][0] - Xin[0];
			dy = X[j][1] - Xin[1];
			dz = X[j][2] - Xin[2];
			if(dx == 0.0 && dy == 0.0 && dz == 0.0){
    			SinNode = j + 1;
    			getNormal_tria6(mElems[i][j],mNodes,mElems,normal);
    			break;
			}
   		}

   		if(SinNode == 0) intG_tria6(X,Xin,SubT);
   		else intSingularG_tria6(SinNode,X,Xin,SubT);

    	/* Add contributions to temperature from all j nodes in element i */
    	for(j = 0; j < NODES_IN_ELEM; j++) Ttmp += SubT[j]*vB[mElems[i][j]-1];
	}
	
	/* Contribution from domain integral */	
	for(i = 0; i < nxi; i++){
		dx = Xin[0] - xi[i];
		for(j = 0; j < nyi; j++){
			dy = Xin[1] - yi[j];
			for(k = 0; k < nzi; k++){
				dz = Xin[2] - zi[k];
				dr = sqrt(dx*dx + dy*dy + dz*dz);
				if(dr != 0.0) Tdomain -= bi[i][j][k]/dr;
			}
		}
	}

	*T = A*Ttmp - b0*Tdomain + T0;

	return 0;
}
