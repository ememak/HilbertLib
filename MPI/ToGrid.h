#ifndef TOGRID_H_INCLUDED
#define TOGRID_H_INCLUDED
#include<stdio.h>
#include<mpi.h>
#include"AxesTranspose.h"

void toGrid(double* in, coord_t* out, double* MIN, double* MAX, int Dimensions, int BitsPrecision);
void toPos(coord_t* in, double* out, double* MIN, double* MAX, int Dimensions, int BitsPrecision);
void getMinMax(double* data, int DataCount, double** mind, double** maxd, int Dimensions);
void allToGrid(double* in, double* mind, double* maxd, coord_t** out, int PointsCount, int Dimensions, int BitsPrecision);
#endif
