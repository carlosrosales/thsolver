/******************************************************************************
* File      : field_tria6.c                                                   *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the Electric Field at the required point Xin[] and returns its   *  
* real and imaginary parts in mE[]:                                           *
*                                                                             *
* mE[0] -> Re(Ex)	|	mE[1] -> Im(Ex)                                       *
* mE[2] -> Re(Ey)	|	mE[3] -> Im(Ey)                                       *
* mE[4] -> Re(Ez)	|	mE[5] -> Im(Ez)                                       *
*                                                                             *
* This version uses nDomainElems so that the correct mesh and number of       *
* elements are used in the calculation of the electric field. This is for     *
* compatibility with the temperature calculation, which actually uses nElems  *
* and nNodes for its own mesh.                                                *
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
******************************************************************************/

#include "constants.h"
#include "field_tria6.h"

int field_tria6(double *Xin, double **mNodes, unsigned int **mElems, double *vB,
                double *mE)
{
	const unsigned int ELEMS = nDomainElems, DIM = 3, FLAG = 1, NODES_IN_ELEM = 6;
	unsigned int i, j, k, currentNode, SinNode, test;
	double A, B, dx, dy, dz, normal[3], SubE[18], X[6][3];

	/* INITIALIZATION */
	B = 0.5/eps0;
	A = B/pi2;
	for(i = 0; i < 6; i++) mE[i] = 0.0;

	for(i = 0; i < ELEMS; i++)
	{
	   	/* GET COORDINATES OF ELEMENT NODES */
	   	for(j = 0; j < NODES_IN_ELEM; j++)
	   	{
		   	currentNode = mElems[i][j]-1;
			X[j][0] = mNodes[currentNode][0];
			X[j][1] = mNodes[currentNode][1];
			X[j][2] = mNodes[currentNode][2];
   		}

   		/* IF Xin IN SURFACE USE SINGULAR INTEGRATION */
   		test = 0;
   		for(j = 0; j < NODES_IN_ELEM; j++)
   		{
			dx = X[j][0] - Xin[0];
			dy = X[j][1] - Xin[1];
			dz = X[j][2] - Xin[2];
			if(dx == 0.0 && dy == 0.0 && dz == 0.0 && test == 0)
			{
	    		test = 1;
    			SinNode = j+1;
    			getNormal_tria6(mElems[i][j],mNodes,mElems,normal);
			}
   		}

   		/* DO REGULAR INTEGRAL FOR ALL j NODES IN ELEMENT i */
   		if(test == 0) intH_tria6(FLAG,X,Xin,SubE,normal);
   		else intSingularH_tria6(FLAG,SinNode,X,Xin,SubE,normal);

   		/* ADD CONTRIBUTION TO FIELD FROM ELEMENT i */
   		for(j = 0; j < NODES_IN_ELEM; j++)          
   		{
		   	currentNode = (mElems[i][j] - 1)*2;
			for(k = 0; k < DIM; k++)
			{
	    		mE[k*2] += A*SubE[j*DIM+k]*vB[currentNode];		/* Re(E) */
    			mE[k*2+1] += A*SubE[j*DIM+k]*vB[currentNode+1];	/* Im(E) */
			}
			if(test == 1)
			{
				for(k = 0; k < DIM; k++)
				{
					mE[k*2] += B*vB[currentNode]*normal[k];
					mE[k*2+1] += B*vB[currentNode+1]*normal[k];
				}
				test = 2;
			}
   		}
	}

	/* SUCCESSFUL EXIT */
	return 0;
}
