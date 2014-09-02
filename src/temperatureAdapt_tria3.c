/******************************************************************************
* File      : temperatureAdapt_tria3.c                                        *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the temperature T at the required point Xin. Avoid points in the *
* original surfaces where the field was calculated, as no consideration for   *
* the singular case is made.                                                  *
* Works for linear interpolation in triangular elements (3-noded triangles)   *
*                                                                             *
* domainSize    : number of nodes in the domain integral evaluation           *
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
#include "temperatureAdapt_tria3.h"

int temperatureAdapt_tria3(double *Xin, double **mNodes, unsigned int **mElems,
                           double *vB, evalNode *integralNodes, double *T)
{
	const unsigned int ELEMS = nElems, NODES_IN_ELEM = 3;
	unsigned int currentNode, i, j, SinNode; 
	double A, dr, dx, dy, dz, Ttmp, Tdomain;
	double SubT[3];
	double X[3][3];

	/* Initialize */
	A = 1.0/pi4;
	Ttmp = 0.0;
	Tdomain = 0.0;
 
	for(i = 0; i < ELEMS; i++){
	
	   	/* Element node's coordinates */
	   	for(j = 0; j < NODES_IN_ELEM; j++){
		   	currentNode = mElems[i][j] - 1;
			X[j][0] = mNodes[currentNode][0];
			X[j][1] = mNodes[currentNode][1];
			X[j][2] = mNodes[currentNode][2];
   		}

   		/* Check for singular case */
   		SinNode = 0;
   		for(j = 0; j < NODES_IN_ELEM; j++){
			dx = X[j][0] - Xin[0];
			dy = X[j][1] - Xin[1];
			dz = X[j][2] - Xin[2];
			if(dx == 0 && dy == 0.0 && dz == 0.0){
    			SinNode = j + 1;
    			break;
			}
   		}

   		/* Contribution from all j nodes in element i */
   		if(SinNode == 0) intG_tria3(X,Xin,SubT);
   		else intSingularG_tria3(SinNode,X,Xin,SubT);

	    /* Add total contribution to temperature from element i */
	    for(j = 0; j < NODES_IN_ELEM; j++) Ttmp += A*SubT[j]*vB[mElems[i][j]-1];
	}
	
	/* ADD CONTRIBUTION FROM DOMAIN INTEGRAL */		
	for(i = 0; i < domainSize; i++){
    	if(integralNodes[i].level > -1){
        	dx = Xin[0] - integralNodes[i].x;
        	dy = Xin[1] - integralNodes[i].y;
        	dz = Xin[2] - integralNodes[i].z;
            dr = sqrt( dx*dx + dy*dy + dz*dz );
            if(dr != 0.0) Tdomain -= integralNodes[i].Esq/dr;
		}
	}
    *T = A*Ttmp - Tdomain + T0;

	return 0;
}
