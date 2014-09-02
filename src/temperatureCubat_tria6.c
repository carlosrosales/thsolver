/******************************************************************************
* File      : temperatureCubat_tria6.c                                        *
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
#include "temperatureCubat_tria6.h"

int temperatureCubat_tria6(double *Xin, double **mNodes, unsigned int **mElems, 
    double *vB,  double **cellNodes, unsigned int **cells, double *b, double *T)
{
    const unsigned int ELEMS = nElems, NODES_IN_ELEM = 6;
    unsigned int currentCell, currentNode, i, j, k, l, m, SinNode; 
    double  A, B, bLocal, dx, dx1, dx2, dx3, dy, dy1, dy2, dy3, dz, dz1, dz2,
            dz3, J, jx, jy, jz, Ttmp, Tdomain;
    double  L[3], N[8], normal[3], SubT[6], W[3], XL[3];
    double  X[6][3], XC[8][3];

    /* Initialize */
    A = 1.0/pi4;
    Ttmp = 0.0;
    Tdomain = 0.0;

    for(i = 0; i < ELEMS; i++){
        
        /* Coordinates of this element's nodes */
        for(j = 0; j < NODES_IN_ELEM; j++){
            currentNode = mElems[i][j] - 1;
            X[j][0] = mNodes[currentNode][0];
            X[j][1] = mNodes[currentNode][1];
            X[j][2] = mNodes[currentNode][2];
        }

        /* Test for singular case */
        SinNode = 0;
        for(j = 0; j < NODES_IN_ELEM; j++){
            dx = X[j][0] - Xin[0];
            dy = X[j][1] - Xin[1];
            dz = X[j][2] - Xin[2];
            if(dx == 0.0 && dy == 0.0 && dz == 0.0){                SinNode = j + 1;
                getNormal_tria6(mElems[i][j],mNodes,mElems,normal);
                break;
            }
        }

        if(SinNode == 0) intG_tria6(X,Xin,SubT);
        else intSingularG_tria6(SinNode,X,Xin,SubT);

        /* Add contributions to temperature from all j nodes in element i */
        for(j = 0; j < NODES_IN_ELEM; j++) Ttmp += SubT[j]*vB[mElems[i][j] - 1];
    }
    
    /* Contribution from domain integral */
    for(i = 0; i < nCells; i++){
            
        /* Vertexes of this integration cell */
        for(j = 0; j < 8; j++){
            currentCell = cells[i][j] - 1;
            XC[j][0] = cellNodes[currentCell][0];
            XC[j][1] = cellNodes[currentCell][1];
            XC[j][2] = cellNodes[currentCell][2];
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
                        XL[0] += N[j]*XC[j][0];
                        XL[1] += N[j]*XC[j][1];
                        XL[2] += N[j]*XC[j][2];
                        bLocal += b[cells[i][j] - 1]*N[j];
                    }
                    
                    /* Derivatives of (x,y,z) with respect to (L1,L2,L3) */
                    /* except for a 0.125 factor in each (added in J) */
                    dx1 = (1.0 - L[1])*(L[2]*(XC[5][0] - XC[4][0] + XC[0][0] - 
                    XC[1][0]) + XC[1][0] - XC[0][0] + XC[5][0] - XC[4][0]) + 
                    (1.0 + L[1])*(L[2]*(XC[6][0] - XC[7][0] + XC[3][0] - 
                    XC[2][0]) + XC[6][0] - XC[7][0] + XC[2][0] - XC[3][0]);
                    
                    dx2 = (1.0 - L[0])*(L[2]*(XC[7][0] - XC[4][0] + XC[0][0] - 
                    XC[3][0]) + XC[7][0] - XC[4][0] + XC[3][0] - XC[0][0]) + 
                    (1.0 + L[0])*(L[2]*(XC[6][0] - XC[5][0] + XC[1][0] - 
                    XC[2][0]) + XC[6][0] - XC[5][0] + XC[2][0] - XC[1][0]);
                    
                    dx3 = (1.0 - L[0])*(L[1]*(XC[7][0] - XC[3][0] + XC[0][0] - 
                    XC[4][0]) + XC[7][0] - XC[3][0] - XC[0][0] + XC[4][0]) + 
                    (1.0 + L[0])*(L[1]*(XC[6][0] - XC[2][0] + XC[1][0] - 
                    XC[5][0]) + XC[6][0] - XC[2][0] + XC[5][0] - XC[1][0]);
                    
                    dy1 = (1.0 - L[1])*(L[2]*(XC[5][1] - XC[4][1] + XC[0][1] - 
                    XC[1][1]) + XC[1][1] - XC[0][1] + XC[5][1] - XC[4][1]) + 
                    (1.0 + L[1])*(L[2]*(XC[6][1] - XC[7][1] + XC[3][1] - 
                    XC[2][1]) + XC[6][1] - XC[7][1] + XC[2][1] - XC[3][1]);
                    
                    dy2 = (1.0 - L[0])*(L[2]*(XC[7][1] - XC[4][1] + XC[0][1] - 
                    XC[3][1]) + XC[7][1] - XC[4][1] + XC[3][1] - XC[0][1]) + 
                    (1.0 + L[0])*(L[2]*(XC[6][1] - XC[5][1] + XC[1][1] - 
                    XC[2][1]) + XC[6][1] - XC[5][1] + XC[2][1] - XC[1][1]);
                    
                    dy3 = (1.0 - L[0])*(L[1]*(XC[7][1] - XC[3][1] + XC[0][1] - 
                    XC[4][1]) + XC[7][1] - XC[3][1] - XC[0][1] + XC[4][1]) + 
                    (1.0 + L[0])*(L[1]*(XC[6][1] - XC[2][1] + XC[1][1] - 
                    XC[5][1]) + XC[6][1] - XC[2][1] + XC[5][1] - XC[1][1]);
                    
                    dz1 = (1.0 - L[1])*(L[2]*(XC[5][2] - XC[4][2] + XC[0][2] - 
                    XC[1][2]) + XC[1][2] - XC[0][2] + XC[5][2] - XC[4][2]) + 
                    (1.0 + L[1])*(L[2]*(XC[6][2] - XC[7][2] + XC[3][2] - 
                    XC[2][2]) + XC[6][2] - XC[7][2] + XC[2][2] - XC[3][2]);
                    
                    dz2 = (1.0 - L[0])*(L[2]*(XC[7][2] - XC[4][2] + XC[0][2] - 
                    XC[3][2]) + XC[7][2] - XC[4][2] + XC[3][2] - XC[0][2]) + 
                    (1.0 + L[0])*(L[2]*(XC[6][2] - XC[5][2] + XC[1][2] - 
                    XC[2][2]) + XC[6][2] - XC[5][2] + XC[2][2] - XC[1][2]);
                    
                    dz3 = (1.0 - L[0])*(L[1]*(XC[7][2] - XC[3][2] + XC[0][2] - 
                    XC[4][2]) + XC[7][2] - XC[3][2] - XC[0][2] + XC[4][2]) + 
                    (1.0 + L[0])*(L[1]*(XC[6][2] - XC[2][2] + XC[1][2] - 
                    XC[5][2]) + XC[6][2] - XC[2][2] + XC[5][2] - XC[1][2]);
                        
                    jx = dx1*(dy2*dz3 - dy3*dz2);
                    jy = dy1*(dx3*dz2 - dx2*dz3);
                    jz = dz1*(dx2*dy3 - dx3*dy2);                   	
                    J = 0.125*0.125*0.125*(jx + jy + jz);

                    dx = Xin[0] - XL[0];
                    dy = Xin[1] - XL[1];
                    dz = Xin[2] - XL[2];
                    B = W[0]*W[1]*W[2]*J;
            
                    Tdomain += bLocal*B/sqrt(dx*dx + dy*dy + dz*dz);
                }
            }
        }                 
    }
       
    *T = A*Ttmp + Tdomain + T0;

    return 0;
}
