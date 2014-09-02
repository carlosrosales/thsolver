# Makefile for thSolver
#
# COPYRIGHT & LICENSE INFORMATION
# 
# Copyright 2006 Carlos Rosales Fernandez and The Institute of High
# Performance Computing (A*STAR)
#
# This file is part of thSolver.
#
# thSolver is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# thSolver is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with thSolver; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

SHELL = /bin/sh

# Compiler
CC = gcc

# Flags
LFLAGS = -lm -msse2 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -O3
CFLAGS = -msse2 -mfpmath=sse -funroll-loops -fprefetch-loop-arrays -O3

# Linking Step
thSolver: bodyForceAdapt_tria3.o bodyForceAdapt_tria6.o bodyForceCubat_tria3.o \
     bodyForceCubat_tria6.o bodyForce_tria3.o bodyForce_tria6.o comFilter.o    \
     thSolver.o dotProd.o doubleMatrix.o doublePointer.o doubleTensor.o        \
     doubleVector.o elemType.o errorHandler.o field_tria3.o field_tria6.o      \
     freeDoubleMatrix.o freeDoublePointer.o freeDoubleTensor.o freeUintMatrix.o\
     gaussBksb.o getNormal_tria3.o getNormal_tria6.o initRes.o                 \
     interp3d.o intG_tria3.o intG_tria6.o intH_tria3.o intH_tria6.o            \
     intSingularG_tria3.o intSingularG_tria6.o intSingularH_tria3.o            \
     intSingularH_tria6.o iterGMRES.o L2Norm.o matVectProd.o                   \
     shape_tria3.o shape_tria6.o solverGMRES.o temperatureAdapt_tria3.o        \
     temperatureAdapt_tria6.o temperatureCubat_tria3.o temperatureCubat_tria6.o\
     temperature_tria3.o temperature_tria6.o thermalFormA_tria3.o              \
     thermalFormA_tria6.o thermalPostProcessAdapt_tria3.o                      \
     thermalPostProcessAdapt_tria6.o thermalPostProcessColAdapt_tria3.o        \
     thermalPostProcessColAdapt_tria6.o thermalPostProcessColCubat_tria3.o     \
     thermalPostProcessColCubat_tria6.o thermalPostProcessCubat_tria3.o        \
     thermalPostProcessCubat_tria6.o thermalPostProcessCol_tria3.o             \
     thermalPostProcessCol_tria6.o thermalPostProcess_tria3.o                  \
     thermalPostProcess_tria6.o uintMatrix.o uintVector.o X2L_tria3.o          \
     X2L_tria6.o
	$(CC) $(LFLAGS) bodyForceAdapt_tria3.o bodyForceAdapt_tria6.o              \
     bodyForceCubat_tria3.o bodyForceCubat_tria6.o bodyForce_tria3.o           \
     bodyForce_tria6.o comFilter.o thSolver.o dotProd.o doubleMatrix.o         \
     doublePointer.o doubleTensor.o doubleVector.o elemType.o errorHandler.o   \
     field_tria3.o field_tria6.o freeDoubleMatrix.o freeDoublePointer.o        \
     freeDoubleTensor.o freeUintMatrix.o gaussBksb.o                           \
     getNormal_tria3.o getNormal_tria6.o initRes.o interp3d.o intG_tria3.o     \
     intG_tria6.o intH_tria3.o intH_tria6.o intSingularG_tria3.o               \
     intSingularG_tria6.o intSingularH_tria3.o intSingularH_tria6.o            \
     iterGMRES.o L2Norm.o matVectProd.o shape_tria3.o                          \
     shape_tria6.o solverGMRES.o temperatureAdapt_tria3.o                      \
     temperatureAdapt_tria6.o temperatureCubat_tria3.o temperatureCubat_tria6.o\
     temperature_tria3.o temperature_tria6.o thermalFormA_tria3.o              \
     thermalFormA_tria6.o thermalPostProcessAdapt_tria3.o                      \
     thermalPostProcessAdapt_tria6.o thermalPostProcessColAdapt_tria3.o        \
     thermalPostProcessColAdapt_tria6.o thermalPostProcessColCubat_tria3.o     \
     thermalPostProcessColCubat_tria6.o thermalPostProcessCubat_tria3.o        \
     thermalPostProcessCubat_tria6.o thermalPostProcessCol_tria3.o             \
     thermalPostProcessCol_tria6.o thermalPostProcess_tria3.o                  \
     thermalPostProcess_tria6.o uintMatrix.o uintVector.o X2L_tria3.o          \
     X2L_tria6.o -o thSolver
	mv thSolver ../bin/thSolver


# Compiling Step
bodyForceAdapt_tria3.o: bodyForceAdapt_tria3.c bodyForceAdapt_tria3.h constants.h
	$(CC) -c $(CFLAGS) bodyForceAdapt_tria3.c -o bodyForceAdapt_tria3.o

bodyForceAdapt_tria6.o: bodyForceAdapt_tria6.c bodyForceAdapt_tria6.h constants.h
	$(CC) -c $(CFLAGS) bodyForceAdapt_tria6.c -o bodyForceAdapt_tria6.o

bodyForceCubat_tria3.o: bodyForceCubat_tria3.c bodyForceCubat_tria3.h constants.h
	$(CC) -c $(CFLAGS) bodyForceCubat_tria3.c -o bodyForceCubat_tria3.o

bodyForceCubat_tria6.o: bodyForceCubat_tria6.c bodyForceCubat_tria6.h constants.h
	$(CC) -c $(CFLAGS) bodyForceCubat_tria6.c -o bodyForceCubat_tria6.o

bodyForce_tria3.o: bodyForce_tria3.c bodyForce_tria3.h constants.h
	$(CC) -c $(CFLAGS) bodyForce_tria3.c -o bodyForce_tria3.o

bodyForce_tria6.o: bodyForce_tria6.c bodyForce_tria6.h constants.h
	$(CC) -c $(CFLAGS) bodyForce_tria6.c -o bodyForce_tria6.o

comFilter.o: comFilter.c
	$(CC) -c $(CFLAGS) comFilter.c -o comFilter.o

thSolver.o: thSolver.c thSolver.h gaussData.h constants.h
	$(CC) -c $(CFLAGS) thSolver.c -o thSolver.o

dotProd.o: dotProd.c
	$(CC) -c $(CFLAGS) dotProd.c -o dotProd.o

doubleMatrix.o: doubleMatrix.c 
	$(CC) -c $(CFLAGS) doubleMatrix.c -o doubleMatrix.o

doublePointer.o: doublePointer.c
	$(CC) -c $(CFLAGS) doublePointer.c -o doublePointer.o

doubleTensor.o: doubleTensor.c
	$(CC) -c $(CFLAGS) doubleTensor.c -o doubleTensor.o

doubleVector.o: doubleVector.c
	$(CC) -c $(CFLAGS) doubleVector.c -o doubleVector.o

elemType.o: elemType.c
	$(CC) -c $(CFLAGS) elemType.c -o elemType.o

errorHandler.o: errorHandler.c
	$(CC) -c $(CFLAGS) errorHandler.c -o errorHandler.o

field_tria3.o: field_tria3.c field_tria3.h constants.h
	$(CC) -c $(CFLAGS) field_tria3.c -o field_tria3.o

field_tria6.o: field_tria6.c field_tria6.h constants.h
	$(CC) -c $(CFLAGS) field_tria6.c -o field_tria6.o

freeDoubleMatrix.o: freeDoubleMatrix.c
	$(CC) -c $(CFLAGS) freeDoubleMatrix.c -o freeDoubleMatrix.o

freeDoublePointer.o: freeDoublePointer.c
	$(CC) -c $(CFLAGS) freeDoublePointer.c -o freeDoublePointer.o

freeDoubleTensor.o: freeDoubleTensor.c
	$(CC) -c $(CFLAGS) freeDoubleTensor.c -o freeDoubleTensor.o

freeUintMatrix.o: freeUintMatrix.c
	$(CC) -c $(CFLAGS) freeUintMatrix.c -o freeUintMatrix.o

gaussBksb.o: gaussBksb.c
	$(CC) -c $(CFLAGS) gaussBksb.c -o gaussBksb.o

getNormal_tria3.o: getNormal_tria3.c
	$(CC) -c $(CFLAGS) getNormal_tria3.c -o getNormal_tria3.o

getNormal_tria6.o: getNormal_tria6.c
	$(CC) -c $(CFLAGS) getNormal_tria6.c -o getNormal_tria6.o

initRes.o: initRes.c initRes.h
	$(CC) -c $(CFLAGS) initRes.c -o initRes.o

interp3d.o: interp3d.c
	$(CC) -c $(CFLAGS) interp3d.c -o interp3d.o

intG_tria3.o: intG_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intG_tria3.c -o intG_tria3.o

intG_tria6.o: intG_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intG_tria6.c -o intG_tria6.o

intH_tria3.o: intH_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intH_tria3.c -o intH_tria3.o

intH_tria6.o: intH_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intH_tria6.c -o intH_tria6.o

intSingularG_tria3.o: intSingularG_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intSingularG_tria3.c -o intSingularG_tria3.o

intSingularG_tria6.o: intSingularG_tria6.c integral_tria6.h gaussData.h
	$(CC) -c $(CFLAGS) intSingularG_tria6.c -o intSingularG_tria6.o

intSingularH_tria3.o: intSingularH_tria3.c integral_tria3.h gaussData.h
	$(CC) -c $(CFLAGS) intSingularH_tria3.c -o intSingularH_tria3.o

intSingularH_tria6.o: intSingularH_tria6.c integral_tria6.h
	$(CC) -c $(CFLAGS) intSingularH_tria6.c -o intSingularH_tria6.o

iterGMRES.o: iterGMRES.c iterGMRES.h
	$(CC) -c $(CFLAGS) iterGMRES.c -o iterGMRES.o

L2Norm.o: L2Norm.c
	$(CC) -c $(CFLAGS) L2Norm.c -o L2Norm.o

matVectProd.o: matVectProd.c
	$(CC) -c $(CFLAGS) matVectProd.c -o matVectProd.o

shape_tria3.o: shape_tria3.c
	$(CC) -c $(CFLAGS) shape_tria3.c -o shape_tria3.o

shape_tria6.o: shape_tria6.c
	$(CC) -c $(CFLAGS) shape_tria6.c -o shape_tria6.o

solverGMRES.o: solverGMRES.c solverGMRES.h
	$(CC) -c $(CFLAGS) solverGMRES.c -o solverGMRES.o

temperatureAdapt_tria3.o: temperatureAdapt_tria3.c temperatureAdapt_tria3.h constants.h
	$(CC) -c $(CFLAGS) temperatureAdapt_tria3.c -o temperatureAdapt_tria3.o

temperatureAdapt_tria6.o: temperatureAdapt_tria6.c temperatureAdapt_tria6.h constants.h
	$(CC) -c $(CFLAGS) temperatureAdapt_tria6.c -o temperatureAdapt_tria6.o

temperatureCubat_tria3.o: temperatureCubat_tria3.c temperatureCubat_tria3.h constants.h
	$(CC) -c $(CFLAGS) temperatureCubat_tria3.c -o temperatureCubat_tria3.o

temperatureCubat_tria6.o: temperatureCubat_tria6.c temperatureCubat_tria6.h constants.h
	$(CC) -c $(CFLAGS) temperatureCubat_tria6.c -o temperatureCubat_tria6.o

temperature_tria3.o: temperature_tria3.c temperature_tria3.h constants.h
	$(CC) -c $(CFLAGS) temperature_tria3.c -o temperature_tria3.o

temperature_tria6.o: temperature_tria6.c temperature_tria6.h constants.h
	$(CC) -c $(CFLAGS) temperature_tria6.c -o temperature_tria6.o

thermalFormA_tria3.o: thermalFormA_tria3.c thermalFormA_tria3.h constants.h
	$(CC) -c $(CFLAGS) thermalFormA_tria3.c -o thermalFormA_tria3.o

thermalFormA_tria6.o: thermalFormA_tria6.c thermalFormA_tria6.h constants.h
	$(CC) -c $(CFLAGS) thermalFormA_tria6.c -o thermalFormA_tria6.o

thermalPostProcessAdapt_tria3.o: thermalPostProcessAdapt_tria3.c thermalPostProcessAdapt_tria3.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcessAdapt_tria3.c -o thermalPostProcessAdapt_tria3.o

thermalPostProcessAdapt_tria6.o: thermalPostProcessAdapt_tria6.c thermalPostProcessAdapt_tria6.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcessAdapt_tria6.c -o thermalPostProcessAdapt_tria6.o

thermalPostProcessColAdapt_tria3.o: thermalPostProcessColAdapt_tria3.c thermalPostProcessColAdapt_tria3.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcessColAdapt_tria3.c -o thermalPostProcessColAdapt_tria3.o

thermalPostProcessColAdapt_tria6.o: thermalPostProcessColAdapt_tria6.c thermalPostProcessColAdapt_tria6.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcessColAdapt_tria6.c -o thermalPostProcessColAdapt_tria6.o

thermalPostProcessColCubat_tria3.o: thermalPostProcessColCubat_tria3.c thermalPostProcessColCubat_tria3.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcessColCubat_tria3.c -o thermalPostProcessColCubat_tria3.o

thermalPostProcessColCubat_tria6.o: thermalPostProcessColCubat_tria6.c thermalPostProcessColCubat_tria6.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcessColCubat_tria6.c -o thermalPostProcessColCubat_tria6.o

thermalPostProcessCubat_tria3.o: thermalPostProcessCubat_tria3.c thermalPostProcessCubat_tria3.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcessCubat_tria3.c -o thermalPostProcessCubat_tria3.o

thermalPostProcessCubat_tria6.o: thermalPostProcessCubat_tria6.c thermalPostProcessCubat_tria6.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcessCubat_tria6.c -o thermalPostProcessCubat_tria6.o

thermalPostProcessCol_tria3.o: thermalPostProcessCol_tria3.c thermalPostProcessCol_tria3.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcessCol_tria3.c -o thermalPostProcessCol_tria3.o

thermalPostProcessCol_tria6.o: thermalPostProcessCol_tria6.c thermalPostProcessCol_tria6.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcessCol_tria6.c -o thermalPostProcessCol_tria6.o

thermalPostProcess_tria3.o: thermalPostProcess_tria3.c thermalPostProcess_tria3.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcess_tria3.c -o thermalPostProcess_tria3.o

thermalPostProcess_tria6.o: thermalPostProcess_tria6.c thermalPostProcess_tria6.h constants.h
	$(CC) -c $(CFLAGS) thermalPostProcess_tria6.c -o thermalPostProcess_tria6.o

uintMatrix.o: uintMatrix.c
	$(CC) -c $(CFLAGS) uintMatrix.c -o uintMatrix.o

uintVector.o: uintVector.c
	$(CC) -c $(CFLAGS) uintVector.c -o uintVector.o

X2L_line2.o: X2L_line2.c
	$(CC) -c $(CFLAGS) X2L_line2.c -o X2L_line2.o

X2L_tria3.o: X2L_tria3.c
	$(CC) -c $(CFLAGS) X2L_tria3.c -o X2L_tria3.o

X2L_tria6.o: X2L_tria6.c
	$(CC) -c $(CFLAGS) X2L_tria6.c -o X2L_tria6.o

# Eliminate all object files
clean:
	rm -f *.o

