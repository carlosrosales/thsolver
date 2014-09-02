/******************************************************************************
* File      : thPost.c                                                        *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Takes the file "temperature.dat" and the number of points where the         *
* temperature is calculated and produces the outputs "gnu-temperature.dat" to *
* be used with gnuplot and "gnu-temperature-big.dat" with an interpolated set *
* of values for the temperature for prettier plotting. Uses bilinear          *
* interpolation to obtain 2n-1 nodes instead of n.                            *
* n   : number of points in each direction.                                   *
* col : constant column index.                                                *
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* This file is part of thPost.                                                *
*                                                                             *
* thPost is free software; you can redistribute it and/or modify              *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* thPost is distributed in the hope that it will be useful,                   *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with thPost; if not, write to the Free Software                       *
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errorHandler.c"
#include "doubleMatrix.c"
#include "doubleVector.c"
#include "freeDoubleMatrix.c"

int main(int argc, char *argv[])
{
    FILE *fin, *fout;
    char cBuffer[32];
    unsigned int col, i, j, k, k1, k2, m, mx, my, n, size, sizeInterp, slow;
    double dx, dy, dz, dxi, dyi, dzi, xmax, xmin, ymax, ymin, zmax, zmin;
    double **T, **Ti, **f, **fi, **fy, *x, *y, *z, *xi, *yi, *zi;

    /* Check for correct input */
    if(argc == 2 && strcmp(argv[1],"-h") == 0){
        printf("\n\nCall as:\n\n\t./thpost n col\n");
        printf("\nTakes the file 'temperature.dat' and the number of points");
        printf(" where the\ntemperature is calculated and produces the ");
        printf("outputs 'gnu-temperature.dat' to\nbe used with gnuplot and ");
        printf("'gnu-temperature-big.dat' with an interpolated set\nof values");
        printf(" for the temperature for prettier plotting. Uses bilinear\n");
        printf("interpolation to obtain 2n-1 nodes instead of n.\n\n");
        printf("n   : number of points in each direction.\n");
        printf("col : constant column.\n\n");
        exit(0);
    }
    else if(argc != 3){
        printf("Correct syntaxis is:\n\n\t./thpost n col\n");
        errorHandler("Type './thpost -h' for help.\n");
    }
    else{
        n = atoi(argv[1]);
        col = atoi(argv[2]);
    }
    slow = col%3;
    printf("Constant Column is %d, slowest changing column is %d\n",col,slow+1);

    /* Allocate storage fot Temperature values */
    size = n*n;
    sizeInterp = 2*n - 1;
    T = doubleMatrix(size,4,0);
    f = doubleMatrix(n,n,0);
    fy = doubleMatrix(sizeInterp,sizeInterp,0);
    fi = doubleMatrix(sizeInterp,sizeInterp,0);
    
    x = doubleVector(n,0);
    y = doubleVector(n,0);
    z = doubleVector(n,0);
    xi = doubleVector(2*n-1,0);
    yi = doubleVector(2*n-1,0);
    zi = doubleVector(2*n-1,0);

    /* Read input file */
    if((fin = fopen("temperature.dat","r")) == NULL)
        errorHandler("Unable to open input file 'temperature.dat'");
    for(i =0; i < 4; i++) fscanf(fin,"%s",&cBuffer);
    for(i = 0; i < size; i++) 
        fscanf(fin,"%le %le %le %le",&T[i][0],&T[i][1],&T[i][2],&T[i][3]);

    /* Save output in gnuplot format */
    if((fout = fopen("gnu-temperature.dat","w")) == NULL) 
        errorHandler("Unable to open output file 'gnu-temperature.dat'");

    fprintf(fin,"# X\tY\tZ\tT\n");
    for(i = 0; i < size; i++){
        fprintf(fout,"%le\t%le\t%le\t%le\n",T[i][0],T[i][1],T[i][2],T[i][3]);
        if(i < (size-1)){ if(T[i][slow] != T[i+1][slow]) fprintf(fout,"\n");}
    }
    fclose(fout);
    
    /* Get maximum and minimum values of (x,y,z) */
    xmax = T[0][0];
    ymax = T[0][1];
    zmax = T[0][2];
    for(i = 1; i < size; i++){
        if(T[i][0] > xmax) xmax = T[i][0];
        if(T[i][1] > ymax) ymax = T[i][1];
        if(T[i][2] > zmax) zmax = T[i][2];
    }
    xmin = T[0][0];
    ymin = T[0][1];
    zmin = T[0][2];
    for(i = 1; i < size; i++){
        if(T[i][0] < xmin) xmin = T[i][0];
        if(T[i][1] < ymin) ymin = T[i][1];
        if(T[i][2] < zmin) zmin = T[i][2];
    }
    
    /* Create the non-constant vectors */
    dx = (xmax - xmin)/(n - 1);
    dy = (ymax - ymin)/(n - 1);
    dz = (zmax - zmin)/(n - 1);
    for(i = 0; i < n; i++){
        x[i] = xmin + dx*i;
        y[i] = ymin + dy*i;
        z[i] = zmin + dz*i;
    }
    
    /* Create interpolation vectors */
    dxi = 0.5*(xmax - xmin)/(n - 1);
    dyi = 0.5*(ymax - ymin)/(n - 1);
    dzi = 0.5*(zmax - zmin)/(n - 1);
    for(i = 0; i < sizeInterp; i++){
        xi[i] = xmin + dxi*i;
        yi[i] = ymin + dyi*i;
        zi[i] = zmin + dzi*i;
    }
    
    /* Create matrix version of T[i][3] */
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            f[i][j] = T[i*n + j][3];
        }
    }

    /* Interpolate along y */
    for(i = 0; i < n; i++){
        for(j = 0; j < (sizeInterp - 1); j++){
            if((j%2) == 0) fy[i][j] = f[i][j/2];
            else{
                k1 = (j - 1)/2;
                k2 = (j + 1)/2;
                fy[i][j] = 0.5*(f[i][k1] + f[i][k2]);
            }
        }
        fy[i][sizeInterp - 1] = f[i][n - 1];
    }
    
    /* Interpolate along x using all the interpolated points in y */
    for(j = 0; j < sizeInterp; j++){
        for(i = 0; i < (sizeInterp - 1); i++){
            if((i%2) == 0) fi[i][j] = fy[i/2][j];
            else{
                k1 = (i - 1)/2;
                k2 = (i + 1)/2;
                fi[i][j] = 0.5*(fy[k1][j] + fy[k2][j]);
            }
        }
        fi[sizeInterp - 1][j] = fy[n - 1][j];
    }   
    
    /* Save interpolatied output in gnuplot format */
    if((fout = fopen("gnu-temperature-big.dat","w")) == NULL) 
        errorHandler("Unable to open output file 'gnu-temperature-big.dat'");

    fprintf(fin,"# X\tY\tZ\tT\n");
    for(i = 0; i < sizeInterp; i++){
        for(j = 0; j < sizeInterp; j++){
            if(col == 1) fprintf(fout,"%le\t%le\t%le\t%le\n",xi[i],yi[i],zi[j],fi[i][j]);
            else if(col == 3) fprintf(fout,"%le\t%le\t%le\t%le\n",xi[i],yi[j],zi[i],fi[i][j]);
        }
        fprintf(fout,"\n");
    }
    fclose(fout);
    
    freeDoubleMatrix(T,size);
    freeDoubleMatrix(f,n);
    freeDoubleMatrix(fy,sizeInterp);
    freeDoubleMatrix(fi,sizeInterp);
    free(xi);
    free(yi);
    free(zi);
    free(x);
    free(y);
    free(z);
    
    return 1;
}
