/******************************************************************************
* File      : temperatureCubat_tria6.h                                        *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
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

/* Funcion prototypes */
int getNormal_tria6(unsigned int nNodeID, double **mNodes, 
                    unsigned int **mElems, double *normal);

int intG_tria6(double X[6][3], double *Xeval, double *Int);

int intSingularG_tria6(unsigned int SinNode, double X[6][3], double *Xeval, 
                       double *Int);

/* Global variables */
extern unsigned int nCellNodes, nCells, nElems, NG;
extern const double G1[1][2], G2[2][2], G3[3][2], G4[4][2], G5[5][2], G6[6][2];
