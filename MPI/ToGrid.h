#ifndef TOGRID_H_INCLUDED
#define TOGRID_H_INCLUDED
#include<stdio.h>
#include<mpi.h>
#include"AxesTranspose.h"

// function takes double coordinates and outputs position on grid (centers) (2^BitsPrecision)^Dimensions
void toGrid(double* in, // one point double coordinates (input)
	    coord_t* out, // one point coordinates on grid (output)
	    double* mind, // global min on each dimension (input)
	    double* maxd, // global max on each dimension (input)
	    int Dimensions, // |dimensions| (input)
	    int BitsPrecision // curve precision (input)
	    );
// function takes grid coordinates on grid and outputs its double equivalent
void toPos  (coord_t* in, // one point coordinates on grid (input)
	    double* out, //  one point double coordinates (output)
	    double* mind, // global min on each dimension (input)
	    double* maxd, // global max on each dimension (input)
	    int Dimensions, // |dimensions| (input)
	    int BitsPrecision // curve precision (input)
	    );
// function takes double coordinates and outputs global min and max, used in other grid functions
void getMinMax	(double* in, // local points represented as Dimensions doubles (input)
		int PointsCount, //number of my points to get min and max on coordinates (input)
		double** mind, // global min on each dimension (output)
		double** maxd, // global max on each dimension (output)
		int Dimensions // |dimensions| (input)
		);
// function takes many points and outputs grid positions
void allToGrid	(double* in, // points represented as Dimensions doubles (input)
		double* mind, // global min on each dimension, could be nullptr or output of getMinMax (input)
		double* maxd, // global max on each dimension, could be nullptr or output of getMinMax (input)
		coord_t** out, // grid coordinates represented as Dimensions+1 ints, first is local idx, next are coordinates (output)
		int PointsCount, // number of points to convert to grid (input)
		int Dimensions, // |dimensions| (input)
		int BitsPrecision // curve precision (input)
		);
#endif
