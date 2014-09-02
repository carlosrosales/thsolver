/******************************************************************************
* File      : bodyForceCubat_tria6.c                                          *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the body force term at the required points and returns it in     *
* array vDomain[]. The evaluation is done by qaussian cubature in hexaedral   *
* cells (cubes) with 1 to 216 points depending on the value of NG (the number *
* of integration points is NG*NG*NG). Standard version, no adaptive meshing,  *
* no adaptive integration.                                                    *
* Assumes cells are given in the form         8-----7                         *
* shown on the right. The electric field     /|    /|                         *
* is calculated using quadratic elements    5-----6 |                         *
* in the boundary (6-noded triangles) and   | 4---|-3                         *
* then interpolated using trilinear         |/    |/                          *
* functions at the cubature points.         1-----2                           *
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
#include "bodyForceCubat_tria6.h"

int bodyForceCubat_tria6(double **cellNodes, unsigned int **cells, 
                         double **mDomainNodes, unsigned int **mDomainElems, 
                         double *vDomainSolution, double **mNodes, 
                         unsigned int **mElems, double **vMatParam, 
                         double *vDomain, unsigned int *vBCType, double *b)
{
    unsigned int bc, currentCell, i, j, k, l, m, p;
    double  A, b0, bLocal, nr, r, dx, dx1, dx2, dx3, dy, dy1, dy2, dy3, dz, dz1,
            dz2, dz3, J, jx, jy, jz;
    double  E[6], L[3], N[8], normal[3], W[3], Xin[3], XL[3];
    double  X[8][3];

    /* Initialize */
    b0 = vMatParam[1][0]/(vMatParam[1][1]*pi4);
       
    /* Calculate |E|^2 in the evalation grid nodes */
    for(i = 0; i < nCellNodes; i++){
        Xin[0] = cellNodes[i][0];
        Xin[1] = cellNodes[i][1];
        Xin[2] = cellNodes[i][2];    
        field_tria6(Xin,mDomainNodes,mDomainElems,vDomainSolution,E);
        b[i] = b0*(E[0]*E[0] + E[1]*E[1] + E[2]*E[2] + E[3]*E[3] + E[4]*E[4] + 
        E[5]*E[5]);
    }

    /* Calculate domain integral contribution for each boundary node */
    for(p = 0 ; p < nNodes; p++){
        Xin[0] = mNodes[p][0];
        Xin[1] = mNodes[p][1];
        Xin[2] = mNodes[p][2];
        bc = vBCType[p];
        if(bc == 0 || bc == 6) getNormal_tria6(p+1,mNodes,mElems,normal);

        for(i = 0; i < nCells; i++){
            
            /* Vertexes of this integration cell */
            for(j = 0; j < 8; j++){
                currentCell = cells[i][j] - 1;
                 X[j][0] = cellNodes[currentCell][0];
                X[j][1] = cellNodes[currentCell][1];
                X[j][2] = cellNodes[currentCell][2];
            }
            
            for(k = 0; k < NG; k++){
                for(l = 0; l < NG; l++){
                    for(m = 0; m < NG; m++){
                        switch(NG){
                            case 1: 
                                    L[0] = G1[k][0];
                                    L[1] = G1[l][0];
                                    L[2] = G1[m][0];
                                    W[0] = G1[k][1];
                                    W[1] = G1[l][1];
                                    W[2] = G1[m][1];
                                    break;
                            case 2: 
                                    L[0] = G2[k][0];
                                    L[1] = G2[l][0];
                                    L[2] = G2[m][0];
                                    W[0] = G2[k][1];
                                    W[1] = G2[l][1];
                                    W[2] = G2[m][1];
                                    break;
                            case 3: 
                                    L[0] = G3[k][0];
                                    L[1] = G3[l][0];
                                    L[2] = G3[m][0];
                                    W[0] = G3[k][1];
                                    W[1] = G3[l][1];
                                    W[2] = G3[m][1];
                                    break;
                            case 4: 
                                    L[0] = G4[k][0];
                                    L[1] = G4[l][0];
                                    L[2] = G4[m][0];
                                    W[0] = G4[k][1];
                                    W[1] = G4[l][1];
                                    W[2] = G4[m][1];
                                    break;
                            case 5: 
                                    L[0] = G5[k][0];
                                    L[1] = G5[l][0];
                                    L[2] = G5[m][0];
                                    W[0] = G5[k][1];
                                    W[1] = G5[l][1];
                                    W[2] = G5[m][1];
                                    break;
                            default: 
                                    L[0] = G6[k][0];
                                    L[1] = G6[l][0];
                                    L[2] = G6[m][0];
                                    W[0] = G6[k][1];
                                    W[1] = G6[l][1];
                                    W[2] = G6[m][1];
                        }
                                    	
                        /* Shape functions for trilinear interpolation */
                    	N[0] = 0.125*(1.0 - L[0])*(1.0 - L[1])*(1.0 - L[2]);
                    	N[1] = 0.125*(1.0 + L[0])*(1.0 - L[1])*(1.0 - L[2]);
                    	N[2] = 0.125*(1.0 + L[0])*(1.0 + L[1])*(1.0 - L[2]);
                    	N[3] = 0.125*(1.0 - L[0])*(1.0 + L[1])*(1.0 - L[2]);
                    	N[4] = 0.125*(1.0 - L[0])*(1.0 - L[1])*(1.0 + L[2]);
                    	N[5] = 0.125*(1.0 + L[0])*(1.0 - L[1])*(1.0 + L[2]);
                    	N[6] = 0.125*(1.0 + L[0])*(1.0 + L[1])*(1.0 + L[2]);
                    	N[7] = 0.125*(1.0 - L[0])*(1.0 + L[1])*(1.0 + L[2]);
                    	
                    	/* Get integral evaluation point and interpolated */
                    	/* field intensity at this position  */
                    	XL[0] = 0.0;
                        XL[1] = 0.0;
                        XL[2] = 0.0;
                        bLocal = 0.0;
                    	for(j = 0; j < 8; j++){
                        	XL[0] += N[j]*X[j][0];
                        	XL[1] += N[j]*X[j][1];
                        	XL[2] += N[j]*X[j][2];
                            bLocal += b[cells[i][j] - 1]*N[j];
                    	}

                    	/* Derivatives of (x,y,z) with respect to (L1,L2,L3) */
                    	/* except for a 0.125 factor in each (added in J) */
                    	dx1 = (1.0 - L[1])*(L[2]*(X[5][0] - X[4][0] + X[0][0] - 
                    	X[1][0]) + X[1][0] - X[0][0] + X[5][0] - X[4][0]) + 
                    	(1.0 + L[1])*(L[2]*(X[6][0] - X[7][0] + X[3][0] - 
                    	X[2][0]) + X[6][0] - X[7][0] + X[2][0] - X[3][0]);
                    	
                    	dx2 = (1.0 - L[0])*(L[2]*(X[7][0] - X[4][0] + X[0][0] - 
                    	X[3][0]) + X[7][0] - X[4][0] + X[3][0] - X[0][0]) + 
                    	(1.0 + L[0])*(L[2]*(X[6][0] - X[5][0] + X[1][0] - 
                    	X[2][0]) + X[6][0] - X[5][0] + X[2][0] - X[1][0]);
                    	
                    	dx3 = (1.0 - L[0])*(L[1]*(X[7][0] - X[3][0] + X[0][0] - 
                    	X[4][0]) + X[7][0] - X[3][0] - X[0][0] + X[4][0]) + 
                    	(1.0 + L[0])*(L[1]*(X[6][0] - X[2][0] + X[1][0] - 
                    	X[5][0]) + X[6][0] - X[2][0] + X[5][0] - X[1][0]);
                    	
                    	dy1 = (1.0 - L[1])*(L[2]*(X[5][1] - X[4][1] + X[0][1] - 
                    	X[1][1]) + X[1][1] - X[0][1] + X[5][1] - X[4][1]) + 
                    	(1.0 + L[1])*(L[2]*(X[6][1] - X[7][1] + X[3][1] - 
                    	X[2][1]) + X[6][1] - X[7][1] + X[2][1] - X[3][1]);
                    	
                    	dy2 = (1.0 - L[0])*(L[2]*(X[7][1] - X[4][1] + X[0][1] - 
                    	X[3][1]) + X[7][1] - X[4][1] + X[3][1] - X[0][1]) + 
                    	(1.0 + L[0])*(L[2]*(X[6][1] - X[5][1] + X[1][1] - 
                    	X[2][1]) + X[6][1] - X[5][1] + X[2][1] - X[1][1]);
                    	
                    	dy3 = (1.0 - L[0])*(L[1]*(X[7][1] - X[3][1] + X[0][1] - 
                    	X[4][1]) + X[7][1] - X[3][1] - X[0][1] + X[4][1]) + 
                    	(1.0 + L[0])*(L[1]*(X[6][1] - X[2][1] + X[1][1] - 
                    	X[5][1]) + X[6][1] - X[2][1] + X[5][1] - X[1][1]);
                    	
                    	dz1 = (1.0 - L[1])*(L[2]*(X[5][2] - X[4][2] + X[0][2] - 
                    	X[1][2]) + X[1][2] - X[0][2] + X[5][2] - X[4][2]) + 
                    	(1.0 + L[1])*(L[2]*(X[6][2] - X[7][2] + X[3][2] - 
                    	X[2][2]) + X[6][2] - X[7][2] + X[2][2] - X[3][2]);
                    	
                    	dz2 = (1.0 - L[0])*(L[2]*(X[7][2] - X[4][2] + X[0][2] - 
                    	X[3][2]) + X[7][2] - X[4][2] + X[3][2] - X[0][2]) + 
                    	(1.0 + L[0])*(L[2]*(X[6][2] - X[5][2] + X[1][2] - 
                    	X[2][2]) + X[6][2] - X[5][2] + X[2][2] - X[1][2]);
                    	
                    	dz3 = (1.0 - L[0])*(L[1]*(X[7][2] - X[3][2] + X[0][2] - 
                    	X[4][2]) + X[7][2] - X[3][2] - X[0][2] + X[4][2]) + 
                    	(1.0 + L[0])*(L[1]*(X[6][2] - X[2][2] + X[1][2] - 
                    	X[5][2]) + X[6][2] - X[2][2] + X[5][2] - X[1][2]);
 
                        jx = dx1*(dy2*dz3 - dy3*dz2);
                        jy = dy1*(dx3*dz2 - dx2*dz3);
                        jz = dz1*(dx2*dy3 - dx3*dy2);                   	
                    	J = 0.125*0.125*0.125*(jx + jy + jz);
                        
                    	dx = Xin[0] - XL[0];
                        dy = Xin[1] - XL[1];
                        dz = Xin[2] - XL[2];
                        A = W[0]*W[1]*W[2]*J;
            
                        r = sqrt(dx*dx + dy*dy + dz*dz);
                        if(bc == 1) vDomain[p] -= A*bLocal/r;
                        else{
                            nr = dx*normal[0] + dy*normal[1] + dz*normal[2];
                            vDomain[p] += A*bLocal*nr/(r*r*r);
                        }
                    }
                }
            }                 
        }
    }

    return 0;
}
