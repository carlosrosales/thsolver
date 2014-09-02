/******************************************************************************
* File		: thSolver.c                                                      *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Main function for Poisson 3D thermal problems using triangular elements.    *
* Gaussian cubature, direct evaluation and adaptive direct evaluation are     *
* available as options for the volume integral of the source term.            *
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
#include "gaussData.h"
#include "thSolver.h"

int main()
{
	FILE   *fIn, *fAux;
	time_t currentTime;
	struct tm *localTime;
	char  cFilename[]="input.bem", cBuffer[32], cBCFile[32],
	       cCellNodesFile[32], cCellsFile[32], cDomain[32], cDomainElems[32], 
	       cElemFile[32], cNodeFile[32], cInternalPointsFile[32], cSection[20], 
	       cTmp[32], domainSolver[10], solver[10];
	int nxi, nyi, nzi, t1m, t2m, t3m, t4m, t_totm;
	int n[3], nInterp[3];
	unsigned int ANALYSIS, i, j, nBuffer, nColumns, nInit, preCond;
	unsigned int *vBCType, *indx;
    unsigned int **cells, **vInterfaces, **mElems, **mDomainElems;
    double exchanges, maxError, t1, t1domain, t2, t3, t4, t_tot, t1s, t2s, t3s,
            t4s, t_tots;
    double *b, *vB, *vDomain, *vDomainSolution;
    double **cellNodes, **vMatParam, **mBC, **mNodes, **Xinner, **mA, 
            **mDomainNodes, **trapSize, ***bi, **xColumns;
    evalNode integralNodes[MAXNODES];

	/* Open and initialize log file */
	file_log = fopen("bem.log","a");
        fprintf(file_log,"\nthSolver(): Started @ ");
	currentTime = time(NULL);
	localTime = localtime(&currentTime);
	fputs(asctime(localTime),file_log);

	/* Clean & read input file */
	comFilter(cFilename);
	fIn=fopen(cFilename,"r");

	/* Node and element sections */
	fscanf(fIn,"%s %d %s",cSection,&nNodes,cNodeFile);
	fscanf(fIn,"%s %d %s %s",cSection,&nElems,cElemType,cElemFile);
	for(i = 0; i < 6; i++) cElemType[i] = tolower(cElemType[i]);

	/* Choose number of nodes in element from cElementType */
	nElemType = elemType();
	if(nElemType != 5 && nElemType != 6)
	   errorHandler("Error: Unknown element type specified.\n");

	/* Materials and interfaces sections */
	fscanf(fIn,"%s %d",cSection,&nMats);
	vMatParam = doubleMatrix(nMats,2,0);
	for(i = 0; i < nMats; i++)
        fscanf(fIn,"%d %le %le", &nBuffer, &vMatParam[i][0], &vMatParam[i][1]);
	fscanf(fIn,"%s %d",cSection,&nInterfaces);
	if(nInterfaces > 0){
		vInterfaces = uintMatrix(nInterfaces,2,0);
    	for(i = 0; i < nInterfaces; i++)
            fscanf(fIn,"%d %d %d",&nBuffer,&vInterfaces[i][0],&vInterfaces[i][1]);
	}

	/* Problems section */
	trapSize = doubleMatrix(3,2,0);
	fscanf(fIn,"%s",cSection);
	fscanf(fIn,"%d %d",&nDomain,&nDomainElems);
	fscanf(fIn,"%le %le %le",&trapSize[0][0],&trapSize[1][0],&trapSize[2][0]); /* Min Size */
	fscanf(fIn,"%le %le %le",&trapSize[0][1],&trapSize[1][1],&trapSize[2][1]); /* Max Size */
    fscanf(fIn,"%s",domainSolver);
    for(i = 0; i < strlen(domainSolver); i++)
        domainSolver[i] = toupper(domainSolver[i]);
    if(strcmp(domainSolver,"DIRECT") == 0){
    	fscanf(fIn,"%d %d %d",&n[0],&n[1],&n[2]);
    	fscanf(fIn,"%d %d %d",&nInterp[0],&nInterp[1],&nInterp[2]);
        nxi = nInterp[0];
        nyi = nInterp[1];
        nzi = nInterp[2];
        bi = doubleTensor(nxi,nyi,nzi,0);
        domainSize = (n[0] - 1)*(n[1] - 1)*(n[2] - 1);
    }
    else if(strcmp(domainSolver,"CUBATURE") == 0){
    	fscanf(fIn,"%d %d %d",&NG,&nCellNodes,&nCells);
        fscanf(fIn,"%s %s",cCellNodesFile,cCellsFile);
        cellNodes = doubleMatrix(nCellNodes,3,0);
	    cells = uintMatrix(nCells,8,0);
	    b = doubleVector(nCellNodes,1);
	    comFilter(cCellNodesFile);
        fAux = fopen(cCellNodesFile,"r");
        for(i = 0; i < nCellNodes; i++){
            fscanf(fAux,"%d ",&nBuffer);
            fscanf(fAux,"%le %le %le ",&cellNodes[i][0],&cellNodes[i][1],&cellNodes[i][2]);
        }
        fclose(fAux);   
        comFilter(cCellsFile);
        fAux = fopen(cCellsFile,"r");
        for(i = 0; i < nCells; i++){
            fscanf(fAux,"%d ",&nBuffer);
            for(j = 0; j < 8; j++) fscanf(fAux,"%d ",&cells[i][j]);
        }
        fclose(fAux);
    }
    else if(strcmp(domainSolver,"ADAPTIVE") == 0){
        fscanf(fIn,"%d %d %d",&n[0],&n[1],&n[2]);
        fscanf(fIn,"%le %le",&maxError,&localError);
        domainSize = (n[0] - 1)*(n[1] - 1)*(n[2] - 1);
    }
    else errorHandler("Error: Unknown domain integral solver specified.\n");
	fscanf(fIn,"%s %s",cDomain,cDomainElems);
	fscanf(fIn,"%s",cBCFile);
	vBCType = uintVector(nNodes,1);
	mBC = doubleMatrix(nNodes,3,1);
	vDomainSolution = doubleVector(nDomain*2,0);
	mDomainNodes = doubleMatrix(nDomain,3,0);
	mDomainElems = uintMatrix(nDomainElems,nNodesInElem,0);
	
    /* Clean and read domain files */
    /*************************************************************************** 
    * NOTE: mDomain contains the data for the nDomain points from the electric * 
    * calculation. The first three elements per row are the (x,y,z) coordinates* 
    * of the nodes, and the next two are the real and imaginary parts of the   *
    * source density s(x,y,z). mDomainElems contains the connectivity for the  *
    * mesh in the electrical problem. It is assumed that both the electrical   *
    * and the thermal problems are solved with the same type of element.       *
    ***************************************************************************/
	comFilter(cDomain);
	fAux = fopen(cDomain,"r");
	for(i = 0; i < nDomain; i++){
		fscanf(fAux,"%le %le",&mDomainNodes[i][0],&mDomainNodes[i][1]);
        fscanf(fAux,"%le",&mDomainNodes[i][2]);
        fscanf(fAux,"%le %le",&vDomainSolution[2*i],&vDomainSolution[2*i+1]);
    }
	fclose(fAux);

	comFilter(cDomainElems);
	fAux = fopen(cDomainElems,"r");
	for(i = 0; i < nDomainElems; i++){
	    fscanf(fAux,"%d",&nBuffer); 		/* Discard index */
	    for(j = 0; j < nNodesInElem; j++) fscanf(fAux,"%d",&mDomainElems[i][j]);
	}
	fclose(fAux);

   /* Identify BC type and store as integer */
	comFilter(cBCFile);
	fAux = fopen(cBCFile,"r");
	for(j = 0; j < nNodes; j++){
	    fscanf(fAux,"%d",&nBuffer);     /* Discard index */
	    fscanf(fAux,"%s",cBuffer);
	    if(strncmp(cBuffer,"C",1) == 0 || strncmp(cBuffer,"c",1) == 0){
        	strncpy(cTmp,&cBuffer[1],strlen(cBuffer)-1);
        	nBuffer = atoi(cTmp) + 5;
    	}
    	else nBuffer = atoi(cBuffer);
	        
    	vBCType[j] = nBuffer;

    	/* Decide how many values to read */
		if((nBuffer > 2) && (nBuffer < 6))
        	fscanf(fAux,"%le %le %le",&mBC[j][0],&mBC[j][1],&mBC[j][2]);
    	else if(nBuffer == 0 || nBuffer == 6)
	    	fscanf(fAux,"%le %le",&mBC[j][0],&mBC[j][1]);
		else
	    	fscanf(fAux,"%le",&mBC[j][0]);
    }
	fclose(fAux);

    /* Read analysis section - default solver is Gauss-Jordan */
	fscanf(fIn,"%s %s ",cBuffer,solver);
	for(i = 0; i < strlen(solver); i++) solver[i] = toupper(solver[i]);
	if(strcmp(solver,"GMRES") == 0) fscanf(fIn,"%d %d",&preCond,&nInit);
	
	/* Read internal points section */
	fscanf(fIn,"%s %d %s",cBuffer,&nInternalPoints,cInternalPointsFile); 
	Xinner = doubleMatrix(nInternalPoints,nDim,0);
	
	/* Read columns positions */
	fscanf(fIn,"%s %d",cBuffer,&nColumns);
	if(nColumns){
		fscanf(fIn,"%d",&nColumnType);
		xColumns = doubleMatrix(nColumns,3,0);
		for(i = 0; i < nColumns; i++){
			fscanf(fIn,"%le %le %le",&xColumns[i][0],&xColumns[i][1],&xColumns[i][2]);
		}
	}
	else xColumns = doubleMatrix(1,1,1);
	fclose(fIn);

	/* Allocate memory for nodes and elements */
	mNodes = doubleMatrix(nNodes,nDim,0);
	mElems = uintMatrix(nElems,nNodesInElem,0);

	/* Clean and read nodes file */
	comFilter(cNodeFile);
	fAux = fopen(cNodeFile,"r");
	for(i = 0; i < nNodes; i++){
	    fscanf(fAux,"%d",&nBuffer);    /* Discard index */
	    for(j = 0; j < nDim; j++) fscanf(fAux,"%le",&mNodes[i][j]);
	}
	fclose(fAux);

	/* Clean and read elements file */
    comFilter(cElemFile);
    fAux = fopen(cElemFile,"r");
    for(i = 0; i < nElems; i++){
        fscanf(fAux,"%d",&nBuffer);     /* Discard index */
        for(j = 0; j < nNodesInElem; j++) fscanf(fAux,"%d",&mElems[i][j]);
	}
	fclose(fAux);

	/* Clean and read internal points file */
	comFilter(cInternalPointsFile);
	fAux = fopen(cInternalPointsFile,"r");
	for(i = 0; i < nInternalPoints; i++){
    		fscanf(fAux,"%d",&nBuffer);   /* Discard index */
    		for(j = 0; j < nDim; j++) fscanf(fAux,"%le",&Xinner[i][j]);
	}
	fclose(fAux);

	/* Record program status */
	fprintf(file_log,"\nUsing %d Nodes and %d ",nNodes,nElems);
	fprintf(file_log,"%s Elems to solve Thermal Problem",cElemType);
	fprintf(file_log,"\nNumber of interfaces: %d",nInterfaces);
	fprintf(file_log,"\nNumber of different materials: %d",nMats);
    fprintf(file_log,"\nSystem solver is %s.",solver);
    fprintf(file_log,"\nDomain integral solved using %s method",domainSolver);
	t1 = (double)clock()/CLOCKS_PER_SEC;

	/****************************************************************/
	/***     INPUT DATA OBTAINED, SOLUTION OF PROBLEM FOLLOWS     ***/
	/****************************************************************/

    /* Allocate memory for the coefficient matrix and the Rhs vector */
	mA = doubleMatrix(nNodes,nNodes,1);
	vB = doubleVector(nNodes,1);
	vDomain = doubleVector(nNodes,1);
	
    /* Form coefficient matrix */
	if(nElemType == 5){
		fprintf(file_log,"\n\nthSolver(): calling bodyForce_tria3() ...");
        if(strcmp(domainSolver,"DIRECT") == 0){
            bodyForce_tria3(n,nInterp,trapSize,mDomainNodes,mDomainElems,
            vDomainSolution,mNodes,mElems,vMatParam,vDomain,vBCType,bi);
        }
        else if(strcmp(domainSolver,"CUBATURE") == 0){
            bodyForceCubat_tria3(cellNodes,cells,mDomainNodes,mDomainElems,
            vDomainSolution,mNodes,mElems,vMatParam,vDomain,vBCType,b);
        }
        else{
            bodyForceAdapt_tria3(n,trapSize,maxError,mDomainNodes,mDomainElems,
            vDomainSolution,mNodes,mElems,vMatParam,vDomain,vBCType,
            integralNodes);
        }
		
		t1domain = (double)clock()/CLOCKS_PER_SEC - t1;
		if(t1domain > 60){
    		t2m = (int)(t1domain/60.0);
	    	t2s = t1domain-t2m*60;
    		fprintf(file_log,"\nTime in domain integral:");
    		fprintf(file_log,"\t\t%d m %2.1f s",t2m,t2s);
		}
		else{
            fprintf(file_log,"\nTime in domain integral:");
            fprintf(file_log,"\t\t%2.1f seconds",t1domain);
        }
		
    	fprintf(file_log,"\nthSolver(): calling thermalFormA_tria3() ...");
		thermalFormA_tria3(mNodes,mElems,vBCType,vInterfaces,vMatParam,mA,mBC,
		vB,vDomain);
	}
	else{
		fprintf(file_log,"\nthSolver(): calling bodyForce_tria6() ...");
        if(strcmp(domainSolver,"DIRECT") == 0){
		  bodyForce_tria6(n,nInterp,trapSize,mDomainNodes,mDomainElems,
		  vDomainSolution,mNodes,mElems,vMatParam,vDomain,vBCType,bi);
        }
        else if(strcmp(domainSolver,"CUBATURE") == 0){
            bodyForceCubat_tria6(cellNodes,cells,mDomainNodes,mDomainElems,
            vDomainSolution,mNodes,mElems,vMatParam,vDomain,vBCType,b);
        }
        else{
            bodyForceAdapt_tria6(n,trapSize,maxError,mDomainNodes,mDomainElems,
            vDomainSolution,mNodes,mElems,vMatParam,vDomain,vBCType,
            integralNodes);
        }
		
		t1domain = (double)clock()/CLOCKS_PER_SEC - t1;
		if(t1domain > 60){
    		t2m = (int)(t1domain/60.0);
	    	t2s = t1domain-t2m*60;
    		fprintf(file_log,"\nTime in domain integral:");
    		fprintf(file_log,"\t\t%d m %2.1f s",t2m,t2s);
		}
		else{
            fprintf(file_log,"\nTime in domain integral:");
            fprintf(file_log,"\t\t%2.1f seconds",t1domain);
        }
		
		fprintf(file_log,"\nthSolver(): calling formMatrix_tria6() ...");
    	thermalFormA_tria6(mNodes,mElems,vBCType,vInterfaces,vMatParam,mA,mBC,
    	vB,vDomain);
	}
	t2 = (double)clock()/CLOCKS_PER_SEC-t1;

	/* Solve linear system */
    if(strcmp(solver,"GMRES") == 0){
        fprintf(file_log,"\nthSolver(): calling solverGMRES() ...");
        solverGMRES(preCond,nInit,nNodes,mA,vB);
    }
    else{
        fprintf(file_log,"\nthSolver(): calling gaussBksb() ...");
        gaussBksb(nNodes,mA,vB);
    }
	t3 = (double)clock()/CLOCKS_PER_SEC-t2-t1;

	/**************************************************************/
	/***      SOLUTION OBTAINED, POST-PROCESSING FOLLOWS        ***/
	/**************************************************************/

    /* Temperature at the required internal points */
    if(nElemType == 5){
        fprintf(file_log,"\nthSolver(): calling postProcess_tria3() ...");
        if(strcmp(domainSolver,"DIRECT") == 0){
            if(nColumns != 0){
                thermalPostProcessCol_tria3(nInterp,trapSize,Xinner,mNodes,
                mElems,vB,vMatParam,nColumns,xColumns,bi);
            }
            else{
                thermalPostProcess_tria3(nInterp,trapSize,Xinner,mNodes,
                mElems,vB,vMatParam,bi);
            }
        }
        else if(strcmp(domainSolver,"CUBATURE") == 0){
            if(nColumns != 0){
                thermalPostProcessColCubat_tria3(trapSize,Xinner,mNodes,
                mElems,vB,cellNodes,cells,vMatParam,nColumns,xColumns,b);
            }
            else{
                thermalPostProcessCubat_tria3(trapSize,Xinner,mNodes,
                mElems,vB,cellNodes,cells,vMatParam,b);
            }
        }
        else{
            if(nColumns != 0){
                thermalPostProcessColAdapt_tria3(trapSize,Xinner,mNodes,mElems,
                vB,nColumns,xColumns,integralNodes);
            }
            else{
                thermalPostProcessAdapt_tria3(trapSize,Xinner,mNodes,mElems,vB,
                integralNodes);
            }
	}
    }
    else{
        fprintf(file_log,"\nthSolver(): calling postProcess_tria6() ...");
        if(strcmp(domainSolver,"DIRECT") == 0){
            if(nColumns != 0){
                thermalPostProcessCol_tria6(nInterp,trapSize,Xinner,mNodes,
                mElems,vB,vMatParam,nColumns,xColumns,bi);
            }
            else{
                thermalPostProcess_tria6(nInterp,trapSize,Xinner,mNodes,
                mElems,vB,vMatParam,bi);
            }
        }
        else if(strcmp(domainSolver,"CUBATURE") == 0){
            if(nColumns != 0){
                thermalPostProcessColCubat_tria6(trapSize,Xinner,mNodes,
                mElems,vB,cellNodes,cells,vMatParam,nColumns,xColumns,b);
            }
            else{
                thermalPostProcessCubat_tria6(trapSize,Xinner,mNodes,
                mElems,vB,cellNodes,cells,vMatParam,b);
            }
        }
        else{
            if(nColumns != 0){
                thermalPostProcessColAdapt_tria6(trapSize,Xinner,mNodes,mElems,
                vB,nColumns,xColumns,integralNodes);
            }
            else{
                thermalPostProcessAdapt_tria6(trapSize,Xinner,mNodes,mElems,vB,
                integralNodes);
            }
        }
    }

	/* Save solution vector to file */
	fAux = fopen("solution.dat","w");
	for(i = 0; i < nNodes; i++){
	   fprintf(fAux,"%le %le %le",mNodes[i][0],mNodes[i][1],mNodes[i][2]);
	   fprintf(fAux,"%le\n",vB[i]);
    }
	fclose(fAux);
	t4 = (double)clock()/CLOCKS_PER_SEC-t3-t2-t1;

	/*******************************************************************/
	/***  POST-PROCESSING FINISHED, ANALYSIS OF PERFORMANCE FOLLOWS  ***/
	/*******************************************************************/

	fprintf(file_log,"\n\n*** thSolver(): performance analysis ***");
	if(t1 > 60){
	    t1m = (int)(t1/60.0);
	    t1s = t1-t1m*60.0;
	    fprintf(file_log,"\nTIME READING INPUT: \t\t%d m %2.1f s",t1m,t1s);
	}
	else fprintf(file_log,"\nTIME READING INPUT: \t\t%2.1f seconds",t1);
	
	if(t2 > 60){
    	t2m = (int)(t2/60.0);
	    t2s = t2-t2m*60;
    	fprintf(file_log,"\nTIME IN ASSEMBLY: \t\t%d m %2.1f s",t2m,t2s);
	}
	else fprintf(file_log,"\nTIME IN ASSEMBLY: \t\t%2.1f seconds",t2);
	
	if(t3 > 60){
	    t3m = (int)(t3/60.0);
	    t3s = t3-t3m*60.0;
	    fprintf(file_log,"\nTIME IN SOLVER: \t\t%d m %2.1f s",t3m,t3s);
	}
	else  fprintf(file_log,"\nTIME IN SOLVER: \t\t%2.1f seconds",t3);
	
	if(t4 > 60){
	    t4m = (int)(t4/60.0);
	    t4s = t4-t4m*60.0;
	    fprintf(file_log,"\nTIME IN POST-PROCESSING: \t%d m %2.1f s",t4m,t4s);
	}
	else fprintf(file_log,"\nTIME IN POST-PROCESSING: \t%2.1f seconds",t4);
	
	t_tot = t1+t2+t3+t4;
	if(t_tot > 60){
	    t_totm = (int)(t_tot/60.0);
	    t_tots = t_tot-t_totm*60.0;
	    fprintf(file_log,"\nTOTAL EXECUTION TIME: \t%d m %2.1f s",t_totm,t_tots);
	}
	else fprintf(file_log,"\nTOTAL EXECUTION TIME: \t\t%2.1f seconds",t_tot);

	/************************************************************************/
	/***  PERFORMANCE ANALYSIS FINISHED, CLOSE ALL FILES AND FREE MEMORY  ***/
	/************************************************************************/

	/* Free dynamically allocated memory */
	free(vB);
	free(vDomain);
	free(vDomainSolution);
	free(vBCType);
	freeUintMatrix(mElems,nElems);
	freeUintMatrix(mDomainElems,nDomainElems);
	freeDoubleMatrix(mA,nNodes);
	freeDoubleMatrix(mBC,nNodes);
	freeDoubleMatrix(trapSize,3);
	freeDoubleMatrix(mNodes,nNodes);
	freeDoubleMatrix(vMatParam,nMats);
	freeDoubleMatrix(mDomainNodes,nDomain);
	freeDoubleTensor(bi,nxi,nyi);
	if(nInterfaces != 0) freeUintMatrix(vInterfaces,nInterfaces);
    
    fprintf(file_log,"\n\nthSolver(): Finished @ ");
	currentTime = time(NULL);
	localTime = localtime(&currentTime);
	fputs(asctime(localTime),file_log);
	fprintf(file_log,"\n\n**************************");
	fprintf(file_log,"********************************\n\n");
	fclose(file_log);

	return 0;
}
