#ifndef HILBERTLIBDEFINED
#define HILBERTLIBDEFINED
#include "MyTree.h"

void HilbertLibPartition (coord_t * MyPoints,
			  int MyPointsCount,
			  int RootRank,
			  int Dimensions,
			  int BitsPrecision,
			  int rank,
			  int size,
			  double * oldPoints,
			  int * oldIdx,
			  double * *NewDataPtr,
			  int * *NewDataIdx,
			  int *NewDataSize);

MTNode *HilbertLibPrepareNodeForQueries (
					coord_t * Data,
					int DataSize,
					int Dimensions);

void answerPointQuery (int QuerySender,
		    int Dimensions,
		    coord_t * Data,
		    int DataSize,
		    MTNode * Root,
		    int MyRank,
		    double* data,
		    int * dataIdx,
		    coord_t ** *SelfQueryResult,
		    int *SelfQueryResultCount,
		    coord_t* LD,
		    coord_t* RD);

void answerProcessQuery (int QuerySender,
		    int Dimensions,
		    MTNode * Root,
		    int MyRank,
		    int ** Processes,
		    int * ProcessesCount,
		    coord_t* LD,
		    coord_t* RD);

void recvPointQuery(double * *NewNeighbours,
		int * *NewNeighboursIdx,
		int *NewNeighboursSize,
		double * *Results,
		int * *ResultsIdx,
		int *ResultsSize,
		int Dimensions,
		int ProcessCount,
		double * data,
		int * dataIdx,
		coord_t ** SelfQueryResult,
		int SelfQueryResultCount,
		int MyRank);
#endif
