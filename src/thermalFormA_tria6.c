/******************************************************************************
* File      : bodyForceAdapt_tria6.c                                          *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Forms the coefficient matrix mA by assembling submatrices of the element,   *
* SubA, for the IBEM thermal problem.                                         *
*                                                                             *
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
#include "thermalFormA_tria6.h"

int thermalFormA_tria6(double **mNodes, unsigned int **mElems, 
                       unsigned int *vBCType, unsigned int **vInterfaces, 
                       double **vMatParam, double **mA, double **mBC, 
                       double *vB, double *vDomain)
{
	const unsigned int FLAG = 0, NODES_IN_ELEM = 6;
	unsigned int currentNode, i, interfaceID, j, k, mat1, mat2, SinNode, bc;  
	double B, Kdif, Ktot;
	double Xeval[3], SubA[6], normal[3];
	double X[6][3];
	double *A, *C;
	
	/* DEFINE AUXILIARY QUANTITIES FOR THE EQS AT THE CONDUCTORS */
	B = 1.0/pi4;

	/* DEFINE AUXILIARY COEFFICIENT FOR DIFFERENT INTERFACES */
	if(nInterfaces){
		A = doubleVector(nInterfaces,0);
		C = doubleVector(nInterfaces,0);
		for(i = 0; i < nInterfaces; i++){
			mat1 = vInterfaces[i][0] - 1;
			mat2 = vInterfaces[i][1] - 1;
			Kdif = vMatParam[mat1][1] - vMatParam[mat2][1];
			Ktot = vMatParam[mat1][1] + vMatParam[mat2][1];
			C[i] = Kdif/Ktot;
			A[i] = B*C[i];
		}
	}

	for(i = 0; i < nNodes; i++)
	{
    	/* GET EVALUATION NODE */
    	Xeval[0] = mNodes[i][0];
    	Xeval[1] = mNodes[i][1];
    	Xeval[2] = mNodes[i][2];
    	bc = vBCType[i];
    	if(bc == 6) bc = 0;

	    /* OPERATIONS NEEDED FOR INTERFACES ONLY */
    	if(bc == 0)
    	{
			/* GET INTERFACE AND NORMAL FOR EVALUATION NODE */
			interfaceID = (unsigned int)mBC[i][1] - 1;
			getNormal_tria6(i+1,mNodes,mElems,normal);

			/* ADD DIAGONAL TERMS TO mA */
			mA[i][i] -= 0.5;
    	}  

    	/* DO INTEGRALS FOR THE MATRIX COEFFICIENTS */
    	for(j = 0; j < nElems; j++)
    	{	
			/* FIND OUT IF THERE IS SINGULARITY AND WHERE IN THE ELEMENT */
        	SinNode = 0;
        	for(k = 0; k < NODES_IN_ELEM; k++)
	    	{
	      		currentNode = mElems[j][k] - 1;
				if(currentNode == i) SinNode = k + 1;
	      		X[k][0] = mNodes[currentNode][0];
	      		X[k][1] = mNodes[currentNode][1];
	      		X[k][2] = mNodes[currentNode][2];
        	}

        	/* TEMPERATURE GIVEN */
    		if(bc == 1)
    		{
	    		/* CALL INTEGRATION ROUTINES */
	    		if(SinNode == 0) intG_tria6(X,Xeval,SubA);
	    		else intSingularG_tria6(SinNode,X,Xeval,SubA);

	    		/* ASSEMBLE COEFFICIENTS */
	    		for(k = 0; k < NODES_IN_ELEM; k++)
	    		{		
					currentNode = mElems[j][k] - 1;
        			mA[i][currentNode] += B*SubA[k];
	    		}
    		}
    	
    		/* THERMAL INTERFACE */
	    	else
	    	{
		    	/* CALL INTEGRATION ROUTINES */
		    	if(SinNode == 0) intH_tria6(FLAG,X,Xeval,SubA,normal);
	    		else intSingularH_tria6(FLAG,SinNode,X,Xeval,SubA,normal);

	    		/* ASSEMBLE COEFFICIENTS */
		    	for(k = 0; k < NODES_IN_ELEM; k++)
		    	{
			    	currentNode = mElems[j][k] - 1;
	 		    	mA[i][currentNode] += A[interfaceID]*SubA[k];
	    		}
	     	}
    	}

    	/* ASSEMBLE RIGHT HAND SIDE VECTOR */
    	switch(bc)
    	{
    		case 0: vB[i] = C[interfaceID]*vDomain[i]; break;
    		case 1: vB[i] = mBC[i][0] + vDomain[i] - T0; break;
    	}
	}

	/* FREE DYNAMICALLY ALLOCATED MEMORY IF NECESSARY */
	if(nInterfaces)
	{
		free(A);
		free(C);
	}

	/* SUCCESSFUL EXIT */
	return 0;
}
