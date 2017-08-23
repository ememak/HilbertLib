#ifndef MYTREEINCLUDED
#define MYTREEINCLUDED

#include "AxesTranspose.h"
#include "PtrVector.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MYTREEMINMAX

struct MTNode
{
	int dim;
	coord_t val;		// or size
	void *left, *right;	// if right is NULL : left = pointer to the same(coordinates) object, val is their size
#ifdef MYTREEMINMAX
	coord_t min, max;
#endif
};

typedef struct MTNode MTNode;

void makeMTNode (MTNode * foo, int dimdiv,
		 coord_t val);

void MTmake (MTNode * Node,
	     coord_t * *Data,
	     int DataSize,
	     int Dimensions, int ActDim);
void MTQuery (MTNode * Node,
	      coord_t * LD,
	      coord_t * RD,
	      coord_t ** *Res,
	      int *ResSize, int Dimensions);
void MTQuery2 (MTNode * Node,
	      coord_t * LD,
	      coord_t * RD,
	      int *Res,
	      int Dimensions);
void MTDelete (MTNode * Node);
#endif
