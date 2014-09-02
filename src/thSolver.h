/******************************************************************************
* File      : thSolver.h                                                      *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

/* Common function prototypes */
void errorHandler(char errorText[]);

void comFilter(char *FileName);

void freeDoubleMatrix(double **M, unsigned int ROWS);

void freeDoubleTensor(double ***T, unsigned int ROWS, unsigned int COLS);

void freeUintMatrix(unsigned int **M, unsigned int ROWS);

void solverGMRES(unsigned int preCond, unsigned int nInit, unsigned int nNodes, 
                 double **Amatrix, double *Rhs);

unsigned int *uintVector(unsigned int LENGTH, unsigned int INIT);

unsigned int **uintMatrix(unsigned int ROWS, unsigned int COLS, 
                            unsigned int INIT);
                            int elemType();

int thermalFormA_tria3(double **mNodes, unsigned int **mElems, 
                       unsigned int *vBCType, unsigned int **vInterfaces, 
                       double **vMatParam, double **mA, double **mBC, 
                       double *vB, double *vDomain);
    
int thermalFormA_tria6(double **mNodes, unsigned int **mElems, 
                       unsigned int *vBCType, unsigned int **vInterfaces, 
                       double **vMatParam, double **mA, double **mBC, 
                       double *vB, double *vDomain);
                       
int gaussBksb(unsigned int N, double **A, double *B);

double *doubleVector(unsigned int LENGTH, unsigned int INIT);

double **doubleMatrix(unsigned int ROWS, unsigned int COLS, 
                       unsigned int INIT);

double ***doubleTensor(unsigned int ROWS, unsigned int COLS, 
                        unsigned int DEPTH, unsigned int INIT);

/* Function Prototypes for Direct evaluation version */                      
int bodyForce_tria3(int *N, int *nInterp, double **trapSize, 
                    double **mDomainNodes, unsigned int **mDomainElems, 
                    double *vDomainSolution, double **mNodes, 
                    unsigned int **mElems, double **vMatParam, double *vDomain,
                    unsigned int *vBCType, double ***bi);
    
int bodyForce_tria6(int *N, int *nInterp, double **trapSize, 
                    double **mDomainNodes, unsigned int **mDomainElems, 
                    double *vDomainSolution, double **mNodes, 
                    unsigned int **mElems, double **vMatParam, double *vDomain,
                    unsigned int *vBCType, double ***bi);
                    
int thermalPostProcess_tria3(int *nInterp, double **trapSize, double **Xinner,
                             double **mNodes, unsigned int **mElems, 
                             double *vB, double **vMatParam, double ***bi);
    
int thermalPostProcess_tria6(int *nInterp, double **trapSize, double **Xinner,
                             double **mNodes, unsigned int **mElems, 
                             double *vB, double **vMatParam, double ***bi);
                             
int thermalPostProcessCol_tria3(int *nInterp, double **trapSize,
                                double **Xinner, double **mNodes, 
                                unsigned int **mElems, double *vB,
                                double **vMatParam, unsigned int nColumns, 
                                double **xColumns, double ***bi);
                                
int thermalPostProcessCol_tria6(int *nInterp, double **trapSize, 
                                double **Xinner, double **mNodes, 
                                unsigned int **mElems, double *vB,
                                double **vMatParam, unsigned int nColumns, 
                                double **xColumns, double ***bi);

/* Function Prototypes for Cubat version */                     
int bodyForceCubat_tria3(double **cellNodes, unsigned int **cells, 
                         double **mDomainNodes, unsigned int **mDomainElems, 
                         double *vDomainSolution, double **mNodes, 
                         unsigned int **mElems, double **vMatParam, 
                         double *vDomain, unsigned int *vBCType, double *b);
                         
int bodyForceCubat_tria6(double **cellNodes, unsigned int **cells,
                         double **mDomainNodes, unsigned int **mDomainElems, 
                         double *vDomainSolution, double **mNodes, 
                         unsigned int **mElems, double **vMatParam, 
                         double *vDomain, unsigned int *vBCType, double *b);
    
int thermalPostProcessCubat_tria3(double **trapSize, double **Xinner,
                                  double **mNodes, unsigned int **mElems, 
                                  double *vB, double **cellNodes, 
                                  unsigned int **cells, double **vMatParam, 
                                  double *b);
    
int thermalPostProcessCubat_tria6(double **trapSize, double **Xinner,
                                  double **mNodes, unsigned int **mElems, 
                                  double *vB, double **cellNodes,
                                  unsigned int **cells, double **vMatParam, 
                                  double *b);
    
int thermalPostProcessColCubat_tria3(double **trapSize, double **Xinner,
                                     double **mNodes, unsigned int **mElems, 
                                     double *vB, double **cellNodes,
                                     unsigned int **cells, double **vMatParam, 
                                     unsigned int nColumns, double **xColumns, 
                                     double *b);
    
int thermalPostProcessColCubat_tria6(double **trapSize, double **Xinner,
                                     double **mNodes, unsigned int **mElems, 
                                     double *vB, double **cellNodes,
                                     unsigned int **cells, double **vMatParam, 
                                     unsigned int nColumns, double **xColumns, 
                                     double *b);
                                     
/* Function Prototypes for Adapt version */
int bodyForceAdapt_tria3(int *N, double **trapSize, double maxError, 
                         double **mDomainNodes, unsigned int **mDomainElems, 
                         double *vDomainSolution, double **mNodes, 
                         unsigned int **mElems, double **vMatParam, 
                         double *vDomain, unsigned int *vBCType, 
                         evalNode *integralNodes);
                         
int bodyForceAdapt_tria6(int *N, double **trapSize, double maxError,
                         double **mDomainNodes, unsigned int **mDomainElems, 
                         double *vDomainSolution, double **mNodes, 
                         unsigned int **mElems, double **vMatParam, 
                         double *vDomain, unsigned int *vBCType, 
                         evalNode *integralNodes);
                         
int thermalPostProcessAdapt_tria3(double **trapSize, double **Xinner,
                                  double **mNodes, unsigned int **mElems, 
                                  double *vB, evalNode *integralNodes);
                                  
int thermalPostProcessAdapt_tria6(double **trapSize, double **Xinner,
                                  double **mNodes, unsigned int **mElems, 
                                  double *vB, evalNode *integralNodes);
                                  
int thermalPostProcessColAdapt_tria3(double **trapSize, double **Xinner,
                                     double **mNodes, unsigned int **mElems,
                                     double *vB, unsigned int nColumns,
                                     double **xColumns, evalNode *integralNodes);
                                     
int thermalPostProcessColAdapt_tria6(double **trapSize, double **Xinner,
                                     double **mNodes, unsigned int **mElems, 
                                     double *vB, unsigned int nColumns,
                                     double **xColumns, evalNode *integralNodes); 


/* Define global variables */
FILE    *file_log;
char    cElemType[6];
double  localError;
unsigned int  domainSize, nCellNodes, nCells, nDim, nDomain, nDomainElems, 
                nElems, nElemType, NG, nInternalPoints, nInterfaces, nMats, 
                nNodes, nNodesInElem, nColumnType;
