/******************************************************************************
* File      : constants.h                                                     *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Definitions of all global constants (except gauss quadrature data) and      *
* pre-processor macros.                                                       *
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

#define MALLOC_CHECK_ 3    /* 0 -> no memory checks */
                           /* 1 -> generate error msg but do not kill program */
                           /* 2 -> do not generate error msg but kill program */
                           /* 3 -> generate error msg and kill program */

#define eps0    8.854187817E-12    /* Electric permittivity of vacuum (F/m) */
#define mu0     1.256637061E-06    /* Magnetic permittivity of vacuum (N/A^2) */
#define qe      1.602176531E-19    /* Elementary (electron) charge (C) */
#define pi      3.141592653589793E+00   /* Pi */
#define pi2     6.283185307179586E+00   /* 2*Pi */
#define pi4     1.256637061435917E+01   /* 4*Pi */
#define pio2    1.570796326794896E+00   /* Pi/2 */
#define pio4    7.853981633974483E-01   /* Pi/4 */
#define ip      3.183098861837907E-01   /* 1/Pi */
#define ip2     15915494309189534E-01   /* 1/(2Pi) */

#define T0      3.00E+02    /* Temperature at infinity (room temperature) */
#define NT      20          /* Number of points for T statistics (NTxNTxNT) */ 
#define MAXNODES    100000  /* Maximum number of nodes in the domain integral */

/* Constants for GMRES */
#define GMRES_RESTART   0   /* 0 = not activated; 1 = activated */
#if GMRES_RESTART
    #define MAX_OUTER_LOOP  10
    #define MAX_ITERS       20      
#else
    #define MAX_OUTER_LOOP  1
    #define MAX_ITERS       600     
#endif
#define TOLERANCE   1.0e-6

/* Data structure for integral evaluation points */
typedef struct EVALNODE{
    int level;
    double x, y, z, Esq, vol;
}evalNode;
