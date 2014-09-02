/******************************************************************************
* File      : bodyForceAdapt_tria6.c                                          *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the body force term at the required points and returns it in     *    
* array vDomain[nNodes], together with the info on the chosen nodes for the   *
* integral evaluation. Evaluation points for electric field are chosen using  *
* an adaptive meshing based on cell subdivision. Integral quality is          *
* controlled indirectly by the quality of the evaluation of the electric      *
* electric field.                                                             *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
*                                                                             *   
* N[3]          : initial number of points in each direction for E-field      *
*                 calculation.                                                *
* maxError      : maximum allowed error in E-field evaluation.                *
* vDomain[]     : array that returns the domain integral for each node.       *
* domainSize    : number of nodes evaluated in the domain.                    *
* integralNodes : array of structures containing the info on the eval nodes.  *
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
#include "bodyForceAdapt_tria6.h"

int bodyForceAdapt_tria6(int *N, double **trapSize, double maxError, 
                         double **mDomainNodes, unsigned int **mDomainElems, 
                         double *vDomainSolution, double **mNodes, 
                         unsigned int **mElems, double **vMatParam, 
                         double *vDomain, unsigned int *vBCType, 
                         evalNode *integralNodes)         
{
    FILE    *fp;
    int level;
    unsigned int bc, count, currentSize, i, j, k, last, m, maxReached, 
                   newNodes, nx, ny, nz, size, sublevel;
    double average, b0, cellVolume, dr, dx, dy, dz, DX, DY, DZ, Lx, Ly, Lz, sx,
            sy, sz, nr, oldSum, passError, sum, trapVolume, xc, yc, zc;
    double  E[6], normal[3], Xin[3];

    /* Initialize */
    nx = N[0] - 1;
    ny = N[1] - 1;
    nz = N[2] - 1;
    size = nx*ny*nz;
    sx = trapSize[0][1] - trapSize[0][0];
    sy = trapSize[1][1] - trapSize[1][0];
    sz = trapSize[2][1] - trapSize[2][0];
    trapVolume = sx*sy*sz;
    cellVolume = trapVolume/size;
    DX = sx/nx;
    DY = sy/ny;
    DZ = sz/nz;
    Lx = 0.5*(sx - DX);
    Ly = 0.5*(sy - DY);
    Lz = 0.5*(sz - DZ);

    /* Create initial mesh */
    fp = fopen("original-domain-mesh.dat","w");
    for(k = 0; k < nz; k++){
        for(j = 0; j < ny; j++){
            for(i = 0; i < nx; i++){
                m = k*ny*nx + j*nx + i;
                integralNodes[m].level = 0;
                Xin[0] = xc = integralNodes[m].x = - Lx + i*DX;
                Xin[1] = yc = integralNodes[m].y = - Ly + j*DY;
                Xin[2] = zc = integralNodes[m].z = - Lz + k*DZ;
                integralNodes[m].vol = cellVolume;
                field_tria6(Xin,mDomainNodes,mDomainElems,vDomainSolution,E);
                integralNodes[m].Esq = E[0]*E[0] + E[1]*E[1] + E[2]*E[2] + 
                E[3]*E[3] + E[4]*E[4] + E[5]*E[5];
                fprintf(fp,"%le %le %le %le\n",xc,yc,zc,integralNodes[m].Esq);
            }
        }
    }
    fclose(fp);

    /* Adaptive remeshing */
    fp = fopen("domain-info.dat","w");
    sublevel = 0;
    maxReached = 0;
    do{
        /* Calculate refinement criteria */
        sum = 0.0;
        count = 0;
        for(i = 0; i < size; i++){
            if(integralNodes[i].level > -1){
                sum += integralNodes[i].Esq*integralNodes[i].vol;
                count++;
            }
        }
        average = sum/count; 

        /* Save info about this subdivision pass */
        if( sublevel > 0 ){
            passError = fabs(fabs(sum - oldSum)/sum);
            fprintf(fp,"PASS #%d\n----------\nError\t\t:  %le\n",sublevel,passError);
            fprintf(fp,"New Nodes\t: %d\nTotal Nodes\t: %d\n\n",newNodes,count);
            if( passError < maxError ) break;
        }
        
        /* Refine mesh as necessary */
        newNodes = 0;
        for(i = 0; i < size; i++){
            currentSize = size + newNodes;
            if(integralNodes[i].level > -1){//exclude previously subdivided cells
                if( (integralNodes[i].Esq*integralNodes[i].vol/average) > localError ){
                    if((size + newNodes + 8) < MAXNODES){ // subdivide this cell
                        xc = integralNodes[i].x;
                        yc = integralNodes[i].y;
                        zc = integralNodes[i].z;
                        level = integralNodes[i].level + 1;
                        dx = 0.5*DX/pow(2.0,level);
                        dy = 0.5*DY/pow(2.0,level);
                        dz = 0.5*DZ/pow(2.0,level);
                        last = currentSize;

                        integralNodes[last].x = xc + dx;
                        integralNodes[last].y = yc + dy;
                        integralNodes[last].z = zc + dz;
            
                        last++;
                        integralNodes[last].x = xc + dx;
                        integralNodes[last].y = yc - dy;
                        integralNodes[last].z = zc + dz;
            
                        last++;
                        integralNodes[last].x = xc - dx;
                        integralNodes[last].y = yc + dy;
                        integralNodes[last].z = zc + dz;
            
                        last++;
                        integralNodes[last].x = xc - dx;
                        integralNodes[last].y = yc - dy;
                        integralNodes[last].z = zc + dz;

                        last++;
                        integralNodes[last].x = xc + dx;
                        integralNodes[last].y = yc + dy;
                        integralNodes[last].z = zc - dz;
            
                        last++;
                        integralNodes[last].x = xc + dx;
                        integralNodes[last].y = yc - dy;
                        integralNodes[last].z = zc - dz;
 
                        last++;
                        integralNodes[last].x = xc - dx;
                        integralNodes[last].y = yc + dy;
                        integralNodes[last].z = zc - dz;

                        last++;
                        integralNodes[last].x = xc - dx;
                        integralNodes[last].y = yc - dy;
                        integralNodes[last].z = zc - dz;

                        for(j = currentSize; j <= last; j++){
                            integralNodes[j].level = level;
                            integralNodes[j].vol = 0.125*integralNodes[i].vol;
                            Xin[0] = xc = integralNodes[j].x;
                            Xin[1] = yc = integralNodes[j].y;
                            Xin[2] = zc = integralNodes[j].z;
                            field_tria6(Xin,mDomainNodes,mDomainElems,
                            vDomainSolution,E);
                            integralNodes[j].Esq = E[0]*E[0] + E[1]*E[1] + 
                            E[2]*E[2] + E[3]*E[3] + E[4]*E[4] + E[5]*E[5];
                        }

                        integralNodes[i].level = -1;  /* remove from refined mesh */
                        newNodes += 8;
                    }
                    else{
                        maxReached = 1;
                        break;
                    }
                }
            }
        }
        if( maxReached == 0 ){
            size += newNodes;
            sublevel++;
            oldSum = sum;
        }
    }while( (newNodes > 0) && (maxReached == 0) );
    if(maxReached == 1) 
        fprintf(fp,"\n*** WARNING: MAXNODES VALUE EXCEEDED ! ***\n");
    fclose(fp);
    
    /* Save refined mesh and normalize Esq */
    b0 = vMatParam[1][0]/(vMatParam[1][1]*pi4);
    fp = fopen("refined-domain-mesh.dat","w");
    if(maxReached == 1) 
        fprintf(fp,"\n# *** WARNING: MAXNODES VALUE EXCEEDED ! ***\n");
    for(i = 0; i < size; i++)
        integralNodes[i].Esq = b0*integralNodes[i].Esq*integralNodes[i].vol;
    for(i = 0; i < size; i++){
        if(integralNodes[i].level > -1){
            fprintf(fp,"%le %le %le %le\n",integralNodes[i].x,
            integralNodes[i].y,integralNodes[i].z,integralNodes[i].Esq);
        }
    }
    fclose(fp);

    /* Calculate domain integral contribution using refined mesh */
    for(i = 0 ; i < nNodes; i++){
        Xin[0] = mNodes[i][0];
        Xin[1] = mNodes[i][1];
        Xin[2] = mNodes[i][2];
        bc = vBCType[i];
        if(bc == 0 || bc == 6) getNormal_tria6(i+1,mNodes,mElems,normal);   
        
        for(j = 0; j < size; j++){
            if(integralNodes[j].level > -1){
                dx = Xin[0] - integralNodes[j].x;
                dy = Xin[1] - integralNodes[j].y;
                dz = Xin[2] - integralNodes[j].z;
                dr = sqrt( dx*dx + dy*dy + dz*dz );
                
                if(bc == 1) vDomain[i] -= integralNodes[j].Esq/dr;
                else{
                    nr = dx*normal[0] + dy*normal[1] + dz*normal[2];
                    vDomain[i] += integralNodes[j].Esq*nr/(dr*dr*dr);
                }
            }
        }
    }

    domainSize = size;
    
    return 0;
}
