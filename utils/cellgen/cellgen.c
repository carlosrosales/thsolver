/******************************************************************************
* File		: cellgen.c                                                       *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Generates a regular 3D mesh with limits [(xmin,xmax),(ymin,ymax),(zmin,     *
* zmax)]. The output files are called "cells.dat" and "cellnodes.dat".        *
*                                                                             *
* Example:    ./cellgen 100 20 10 -50 50 -10 10 -25 25                        *
* Produces a regular mesh in the range x = (-50,50),y = (-10,10),z = (-25,25) *
* with 100 points in the x direction, 20 in the y direction and 10 in z.      *
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* cellgen is free software; you can redistribute it and/or modify             *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* cellgen is distributed in the hope that it will be useful,                  *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with cellgen; if not, write to the Free Software                      *
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "errorHandler.c"

int main(int argc, char *argv[])
{
	FILE *nodefile, *cellfile;
	unsigned int count, i, j, k, nx, ny, nz, node1, node2, node3, node4, node5, node6, node7, node8;
	double dx, dy, dz, x, xmin, xmax, y, ymin, ymax, z, zmin, zmax;

	/* EXIT IN INPUT/OUTPUT ERROR */
	if(argc == 2 && strcmp(argv[1],"-h") == 0){
		printf("\n\nCall as:\n\n\t./cellgen nx ny nz xmin xmax ymin ymax ");
		printf("zmin zmax\n\nGenerates a regular 3D mesh with");
		printf(" limits [(xmin,xmax),(ymin,ymax),(zmin,\nzmax)]. The output ");
		printf("files are called 'cells.dat' and 'cellnodes.dat'.\n\n");
		printf("Example of use:\n\n\tmesh 100 20 10 -50 50 -10 10 -25 25\n\n");
		printf("This produces a regular mesh in the range x = (-50,50), ");
		printf("y = (-10,10),\nz = (-25,25) with 100 points in the x ");
		printf("direction, 20 in the y direction\nand 10 in z direction.\n\n");
		exit(0);
	}
	else if(argc != 10){
	   printf("\nCorrect syntaxis is:\n\n\t./cellgen nx ny nz xmin xmax ymin ");
	   printf("ymax zmin zmax\n");
	   errorHandler("Type './cellgen -h' for help.\n");
	}
	else{
    	nx = atoi(argv[1]);
    	ny = atoi(argv[2]);
    	nz = atoi(argv[3]);
    	xmin = atof(argv[4]);
    	xmax = atof(argv[5]);
    	ymin = atof(argv[6]);
    	ymax = atof(argv[7]);
    	zmin = atof(argv[8]);
    	zmax = atof(argv[9]);
	}

	/* Initialize */
	dx = (xmax - xmin)/(nx - 1);
	dy = (ymax - ymin)/(ny - 1);
    dz = (zmax - zmin)/(nz - 1);
	if((nodefile = fopen("cellnodes.dat","w")) == NULL)
	    errorHandler("Error: Can't open output file cellnodes.dat");
	if((cellfile = fopen("cells.dat","w")) == NULL)
	    errorHandler("Error: Can't open output file cells.dat");
    
	/* Save node file */
	count = 0;
	for(k = 0; k < nz; k++){
		z = zmin + k*dz;
		for(j = 0; j < ny; j++){
			y = ymin + j*dy;
			for(i = 0; i < nx; i++){
    			x = xmin + i*dx;
    			fprintf(nodefile,"%d\t%le\t%le\t%le\n",count+1,x,y,z);
    			count++;
			}
		}
	}
	fclose(nodefile);
	
	/* Save connectivity file */
	count = 1;
	for(k = 0; k < (nz - 1); k++){
		for(j = 0; j < (ny - 1); j++){
			for(i = 0; i < (nx - 1); i++){
    			node1 = i + j*nx + k*nx*ny;
    			node2 = node1 + 1;
    			node3 = node2 + nx;
    			node4 = node1 + nx;
    			node5 = node1 + nx*ny;
    			node6 = node2 + nx*ny;
    			node7 = node3 + nx*ny;
    			node8 = node4 + nx*ny;
    			fprintf(cellfile,"%d\t%d\t%d\t",count,node1+1,node2+1);
    			fprintf(cellfile,"%d\t%d\t%d\t",node3+1,node4+1,node5+1);
    			fprintf(cellfile,"%d\t%d\t%d\n",node6+1,node7+1,node8+1);
    			count++;
			}
		}
	}
    fclose(cellfile);
    				
    printf("\nNumber of nodes:\t%d",nx*ny*nz);
    printf("\nNumber of cells:\t%d\n\n",(nx-1)*(ny-1)*(nz-1));
    
	return 0;
}

 
