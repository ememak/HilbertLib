#ifndef HILBERTLIBDEFINED
#define HILBERTLIBDEFINED

void HilbertLibPartition (coord_t * MyPoints,
			  int MyPointsCount,
			  int RootRank,
			  int Dimensions,
			  int BitsPrecision,
			  int rank,
			  int size,
			  coord_t * *NewDataPtr,
			  coord_t * *NewDataIdx,
			  int *NewDataSize);
/*
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

void recvPointQuery (coord_t * *NewNeighbours,
		  int *NewNeighboursSize,
		  coord_t *** Results,
		  int * ResultsSize,
		  int Dimensions,
		  int ProcessCount,
		  coord_t ** SelfQueryResult,
		  int SelfQueryResultCount,
		  int MyRank);*/
#endif
