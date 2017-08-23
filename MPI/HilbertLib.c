#include <assert.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "AxesTranspose.h"
#include "HilbertLib.h"
#define calloc(a,b) (a==0 ? NULL : calloc(a,b))
#define malloc(a) (a==0 ? NULL : malloc(a))
#define QUERY_TAG 1000
#define ANSWER_TAG 10000

hilpos_t *HilbertPos;

int
HilbertLibCurveSortComparator (const void *_elem1, const void *_elem2) //comparing positions on curve, HilbertPos[x][0] is local idx of point x
{
	coord_t *__elem1 = (coord_t *) _elem1;
	hilpos_t pos = HilbertPos[(__elem1)[0]];

	coord_t *__elem2 = (coord_t *) _elem2;
	hilpos_t pos2 = HilbertPos[(__elem2)[0]];
	return (pos<pos2) ? -1 : (pos>pos2);
}

// X - data, Datasize - |data|, b - bits, n - |dimensions| [Node]
void
HilbertLibNodeCurveSort (coord_t * X,	// particles represented as MDPoints (input & output), sort in place
			 hilpos_t * *HCoordinates,	// Hilbert coordinates to return (output)
			 int Datasize,	// |data| - number of particles in this node (input)
			 int b,	// precision (input)
			 int n	// dimensions (input)
	)
{
	if (Datasize == 0)
	  {
		  (*HCoordinates) = NULL;
		  return;
	  }
	int i = 0;
	coord_t *tmp = malloc (n * sizeof (coord_t));
	coord_t *tmp2 = malloc (n * sizeof (coord_t));
	hilpos_t *first_elem = malloc (Datasize * sizeof (hilpos_t));
	for (i = 0; i < Datasize; i++)
	  {
		  first_elem[i] =
			  GetHCoordinate (&(X[(n+1) * i + 1]),
					  tmp2, tmp, b, n);
	  }
	free (tmp);
	free (tmp2);

	HilbertPos = first_elem;
	qsort (X, Datasize, sizeof (coord_t)*(n+1),
	       HilbertLibCurveSortComparator);

	(*HCoordinates) = malloc (Datasize * sizeof (hilpos_t));
	for (i = 0; i < Datasize; i++)
	    (*HCoordinates)[i] = first_elem[X[(n+1)*i]];
	free (first_elem);
}

// how many particles have hcoordinates <= Right [Node]
int
HilbertLibNodeBinSearch (hilpos_t * HCoordinates,
			 int Datasize, hilpos_t Right)
{
	// binary search left,right,middle
	int bsleft = 0, bsright = Datasize, bsmiddle;
	// properly it is binsearching the first node which have hcoordinate > Right
	while (bsleft < bsright)
	  {
		  bsmiddle = (bsleft + bsright) / 2;
		  if (HCoordinates[bsmiddle] <= Right)
		    {
			    bsleft = bsmiddle + 1;
		    }
		  else
		    {
			    bsright = bsmiddle;
		    }
	  }
	return bsleft;
}

// counting all nodes(hc) which satisfy : Left < hc <= Right [Node]
int
HilbertLibNodeHowMany (hilpos_t * HCoordinates,
		       int Datasize, hilpos_t Left, hilpos_t Right)
{
	if (Datasize == 0)
		return 0;
	return HilbertLibNodeBinSearch
		(HCoordinates, Datasize,
		 Right) -
		HilbertLibNodeBinSearch
		(HCoordinates, Datasize, Left);
}

// computing min and max on local points
void
HilbertLibNodeGetMINMAX (hilpos_t * HCoordinates,
			 int Datasize, hilpos_t * MIN, hilpos_t * MAX)
{
	int i;
	(*MIN) = HILPOS_T_MAX_VALUE;
	(*MAX) = 0;
	if (Datasize == 0)
		return;
	for (i = 0; i < Datasize; i++){
	    if (HCoordinates[i] < *MIN)
		*MIN = HCoordinates[i];
	}
	for(i = 0; i < Datasize; i++){
	    if(HCoordinates[i] > *MAX)
		*MAX = HCoordinates[i];
	}
}

// divide points into bins[Root] 
int
HilbertLibGetNOfParticles (int ProcessCount, // communication size (input)
			   int PointCount, // local number of points (input)
			   int RootRank, // rank of root process (input)
			   hilpos_t * MIN, // global min on hilbert coordinate (output)
			   hilpos_t * MAX, // global max on hilbert coordinate (output)
			   hilpos_t * HCoordinates // point hilbert coordinates, output of HilbertLibCurveSort (input)
			   )
{
	int i;
	int sum = 0;
	MPI_Reduce(&PointCount, &sum, 1, MPI_INT, MPI_SUM, RootRank, MPI_COMM_WORLD); // get total number of points and return it

	HilbertLibNodeGetMINMAX (HCoordinates, PointCount, MIN, MAX);
	hilpos_t *sendBuf = calloc (2, sizeof (hilpos_t));
	sendBuf[0] = *MIN;
	sendBuf[1] = *MAX;
	hilpos_t *recvBuf = calloc (2 * ProcessCount, sizeof (hilpos_t));

	MPI_Gather (sendBuf, 2, MPI_HILPOS_T,
		    recvBuf, 2, MPI_HILPOS_T,
		    RootRank, MPI_COMM_WORLD);

	//gather Min and Max from other processes
	hilpos_t MinRes = HILPOS_T_MAX_VALUE, MaxRes = 0;
	for (i = 0; i < ProcessCount * 2; i += 2)
	  {
		  if (recvBuf[i] < MinRes)
			  MinRes = recvBuf[i];
		  if (recvBuf[i + 1] > MaxRes)
			  MaxRes = recvBuf[i + 1];
	  }
	*MIN = MinRes;
	*MAX = MaxRes;
	free (sendBuf);
	free (recvBuf);
	return sum;
}

// Calculates the next boundary [Root]
hilpos_t
HilbertLibCalculateNextBoundary (hilpos_t a, //left border of not assigned space in hilbert coordinate (input)
				 hilpos_t b, //right border of not assigned space in hilbert coordinate (input)
				 hilpos_t * MyHCoordinates, //hilbert coordinates of my points, output of HilbertLibCourveSort (input)
				 int MyPointsCount, // local points count (input)
				 int particlesRate, // approximation of how many points should next process have (input)
				 int RootRank, // rank of root process (input)
				 int *how_many_used // how many point are assigned to next process (output)
				 )
{
	hilpos_t bsleft = a, bsright = b, bsmiddle;
	int i = 0;
	int suma = 0;
	hilpos_t *sendBuff = calloc (2, sizeof (hilpos_t));
	while (bsright > bsleft * (1 + 2.0*HILPOS_EPS)+HILPOS_EPS)
	{
		  bsmiddle = (bsleft + bsright +
			   HILPOS_BS_1) / HILPOS_BS_2;

		  sendBuff[0] = a;
		  sendBuff[1] = bsmiddle;

		  MPI_Bcast (sendBuff, 2,
			     MPI_HILPOS_T, RootRank, MPI_COMM_WORLD);
		  int singlebuff =
			  HilbertLibNodeHowMany (MyHCoordinates,
						 MyPointsCount,
						 sendBuff[0],
						 sendBuff[1]);

		  MPI_Reduce (&singlebuff, &suma,
			      1, MPI_INT, MPI_SUM,
			      RootRank, MPI_COMM_WORLD);
		  if (suma > particlesRate)
		    {
			    bsright = bsmiddle*(1.0-2.0*HILPOS_EPS);
		    }
		  else
		    {
			    bsleft = bsmiddle;
		    }
		  suma=0;
	  }
	// Getting the information about how many points are in (a,bsleft>
	sendBuff[0] = a;
	sendBuff[1] = bsleft;
	MPI_Bcast (sendBuff, 2, MPI_HILPOS_T,
		   RootRank, MPI_COMM_WORLD);

	int singlebuff = HilbertLibNodeHowMany (MyHCoordinates,
						MyPointsCount,
						sendBuff[0],
						sendBuff[1]);

	MPI_Reduce (&singlebuff, &suma, 1,
		    MPI_INT, MPI_SUM, RootRank, MPI_COMM_WORLD);
	(*how_many_used) = suma;
	free (sendBuff);
	return bsright;
}

// Make Bins [Root]
hilpos_t *
HilbertLibRootMakeBins (int RootRank,	// rank of root process (input)
			int size, // size of communication (input)
			hilpos_t * MyHCoordinates,	// Hilbert Coordinates of root points (input)
			int MyPointsCount,	// Number of points assigned to root (input)
			int dimensions	// |dimensions|
	)
{
	hilpos_t MaxCoord, MinCoord;
	int allPointsCount = HilbertLibGetNOfParticles (size,
							MyPointsCount,
							RootRank,
							&MinCoord,
							&MaxCoord,
							MyHCoordinates);

	int particlesRate = allPointsCount / size;
	int i = 0;
	// in binsBoundaries[i] there will be last Hilbert Coordinate for process i
	hilpos_t *binsBoundaries = malloc (size * sizeof (hilpos_t));
	hilpos_t lastUsed = MinCoord*(1.0-2.0*HILPOS_EPS);
	for (i = 0; i < size; i++)
	  {
		  int how_many_used = 0;

		  binsBoundaries[i] =
			  HilbertLibCalculateNextBoundary
			  (lastUsed, MaxCoord,
			   MyHCoordinates,
			   MyPointsCount,
			   particlesRate,
			   RootRank, &how_many_used);

		  allPointsCount -= how_many_used;
		  lastUsed = binsBoundaries[i];
		  if (i != (size - 1))
			  particlesRate =
				  allPointsCount /
				  (size - i - 1);
	  }
	// last bcast to indicate the end of queries
	hilpos_t *sendbuf = calloc (2, sizeof (hilpos_t));
	sendbuf[0] = 1;
	sendbuf[1] = 0;
	MPI_Bcast (sendbuf, 2, MPI_HILPOS_T,
		   RootRank, MPI_COMM_WORLD);
	free (sendbuf);
	return binsBoundaries;
}

//Sends # of particles to Root [Node]
void
HilbertLibSendNOfParticles (int DataSize,
			    int RootRank, hilpos_t * HCoordinates)
{
	int suma;
	MPI_Reduce(&DataSize, &suma, 1, MPI_INT, MPI_SUM, RootRank, MPI_COMM_WORLD);
	hilpos_t *sendBuff = calloc (2, sizeof (hilpos_t));
	hilpos_t MIN, MAX;
	HilbertLibNodeGetMINMAX (HCoordinates, DataSize, &MIN, &MAX);
	sendBuff[0] = MIN;
	sendBuff[1] = MAX;
	MPI_Gather (sendBuff, 2, MPI_HILPOS_T,
		    NULL, 0, MPI_INT, RootRank, MPI_COMM_WORLD);
	free (sendBuff);
}

//divide points into bins[Node]
void
HilbertLibNodeMakeBins (hilpos_t * MyHCoordinates,
			size_t MyParticlesCount, int RootRank)
{
	//Send number of particles i have
	int suma;
	HilbertLibSendNOfParticles (MyParticlesCount, RootRank, MyHCoordinates);
	hilpos_t *recvBuff = calloc (2, sizeof (hilpos_t));

	while (true)
	  {
		  suma=0;
		  MPI_Bcast (recvBuff, 2,
			     MPI_HILPOS_T, RootRank, MPI_COMM_WORLD);
		  if (recvBuff[1] < recvBuff[0])
		    {
			    break;
		    }

		  int sendbuf = HilbertLibNodeHowMany (MyHCoordinates,
						       MyParticlesCount,
						       recvBuff[0],
						       recvBuff[1]);

		  MPI_Reduce (&sendbuf,
			      &suma,
			      1, MPI_INT, MPI_SUM,
			      RootRank, MPI_COMM_WORLD);
	  }
	free (recvBuff);
}

// Relocate the points, according to Hilbert Curve
#define SENDING_TAG1  111
#define SENDING_TAG2  222
void
HilbertLibRelocate (coord_t * Data, // old points on grid as Dimensions+1 ints, first is local idx, next are coordinates on each dimension (input)
		    hilpos_t * HCoordinates, // hilbert coordinates of each point (input)
		    hilpos_t * Boundaries, // hilbert coordinate of last point of each process (input)
		    int MyPointsCount, // my points count (input)
		    int ProcessCount, // size of communication (input)
		    int Dimensions, // |dimensions| (input)
		    double * oldPoints, // old points precise coordinates (input)
		    int * oldIdx, // old points global idx (input)
		    double * *NewData, // new points precise coordinates (output)
		    int * *NewIdx, // new points global idx (output)
		    int *NewDataCount // number of new points (output)
		    )
{
	//first send the amounts of particles to send later
	int *sendAmounts = calloc (ProcessCount, sizeof (int));
	int i, j;
	int wsk = 0, wskBuf = 0;
	double *sendBuf = malloc (MyPointsCount * Dimensions * sizeof(double));
	int *idxSendBuf = malloc (MyPointsCount * sizeof (coord_t));
	for (i = 0; i < MyPointsCount; i++)
	  {
		  memcpy(sendBuf+i*Dimensions,
			&oldPoints[Dimensions * Data[i*(Dimensions+1)]],
			sizeof(double)*Dimensions);
		  idxSendBuf[i] = oldIdx[Data[i * (Dimensions+1)]];
		  while (Boundaries[wsk] < HCoordinates[i]*(1.0-HILPOS_EPS))
		    {
			    wsk++;
			    assert (wsk < ProcessCount);
		    }
		  sendAmounts[wsk]++;
	  }
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int *recvAmounts = malloc (ProcessCount * sizeof (int));
	MPI_Alltoall (sendAmounts, 1, MPI_INT,
		      recvAmounts, 1, MPI_INT, MPI_COMM_WORLD);
	// All processes now know who will give them what amount of particles
	int pref = 0;
	// iSending
	MPI_Request *requestList, requestNull;

	for (i = 0; i < ProcessCount; i++)
	  {
		  if (sendAmounts[i] != 0){
		  assert (pref < MyPointsCount * (Dimensions+1));
		  assert (pref +
			  sendAmounts[i] *
			  (Dimensions+1) <= MyPointsCount * (Dimensions+1));
		  MPI_Isend (sendBuf + pref,
			     sendAmounts[i] * Dimensions,
			     MPI_DOUBLE, i,
			     SENDING_TAG1,
			     MPI_COMM_WORLD, &requestNull);
		  MPI_Request_free (&requestNull);
		  MPI_Isend (idxSendBuf + (pref/Dimensions),
			    sendAmounts[i],
			    MPI_INT, i,
			    SENDING_TAG2,
			    MPI_COMM_WORLD, &requestNull);
		  MPI_Request_free (&requestNull);

		  pref += sendAmounts[i] * Dimensions;
		  }
	  }
	free (sendAmounts);
	// iRecving
	int myNewPointsSize = 0;
	int zeros = 0;
	for (i = 0; i < ProcessCount; i++)
	  {
		  myNewPointsSize += recvAmounts[i];
		  if (recvAmounts[i] == 0)
			  zeros += 1;
	  }
	int actual_size = 2*(ProcessCount - zeros);
	if (actual_size > 0)
		requestList =
			calloc (actual_size, sizeof (MPI_Request));
	else
		requestList = NULL;

	*NewDataCount = myNewPointsSize;
	(*NewData) = malloc (myNewPointsSize * sizeof (double)*Dimensions);
	(*NewIdx) = malloc (myNewPointsSize * sizeof(int));
	pref = 0;
	unsigned int reqptr = 0;
	for (i = 0; i < ProcessCount; i++)
	  {
		  if (recvAmounts[i] == 0)
			  continue;
		  MPI_Irecv((*NewData) + pref,
			    recvAmounts[i] *
			    Dimensions,
			    MPI_DOUBLE, i,
			    SENDING_TAG1,
			    MPI_COMM_WORLD,
			    &requestList[reqptr*2]);
		  MPI_Irecv ((*NewIdx) + (pref/Dimensions),
			    recvAmounts[i],
			    MPI_INT, i,
			    SENDING_TAG2,
			    MPI_COMM_WORLD,
			    &requestList[reqptr*2+1]);
		  reqptr++;

		  pref += recvAmounts[i] * Dimensions;
	  }
	// Waiting
	if (ProcessCount - zeros != 0)
		MPI_Waitall (actual_size,
			     requestList, MPI_STATUSES_IGNORE);
	MPI_Barrier (MPI_COMM_WORLD);
	free (recvAmounts);
	free (requestList);
	free (sendBuf);
	free (idxSendBuf);
}

void
HilbertLibPartition (coord_t * MyPoints,
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
		     int *NewDataSize)
{
	int i, j;
	hilpos_t *HCoordinates;
	HilbertLibNodeCurveSort (MyPoints,
				 &HCoordinates,
				 MyPointsCount,
				 BitsPrecision, Dimensions);

	hilpos_t *boundaries = NULL;
	if (rank == RootRank)
	  {
		  boundaries = HilbertLibRootMakeBins (RootRank, size, HCoordinates, MyPointsCount, Dimensions);
		  MPI_Bcast(boundaries, size, MPI_HILPOS_T, RootRank, MPI_COMM_WORLD);
	  }
	else
	  {
		  HilbertLibNodeMakeBins (HCoordinates, MyPointsCount, RootRank);
		  boundaries = malloc (size * sizeof(hilpos_t));
		  MPI_Bcast(boundaries, size, MPI_HILPOS_T, RootRank, MPI_COMM_WORLD);
	  }

//for tests with no migration, this HilbertLibRelocate call was commented
	double *NewData;
	int *NewIdx;
	int NewDataCount;
	HilbertLibRelocate (MyPoints,
			    HCoordinates,
			    boundaries,
			    MyPointsCount,
			    size, Dimensions,
			    oldPoints, oldIdx,
			    &NewData, &NewIdx,
			    &NewDataCount);

	(*NewDataPtr)=NewData;
	(*NewDataIdx)=NewIdx;
	(*NewDataSize) = NewDataCount;
	free (MyPoints);
	free (HCoordinates);
	free (boundaries);
	free (oldPoints);
	free (oldIdx);
}

MTNode *
HilbertLibPrepareNodeForQueries (coord_t * Data,
				 int DataSize, int Dimensions)
{
	MTNode *root = calloc (1, sizeof (MTNode));
	makeMTNode (root, 0, 0);
	coord_t **temp = calloc (DataSize,
				 sizeof (coord_t *));
	int i;
	for (i = 0; i < DataSize; i++)
	  {
		  temp[i] = &(Data[(Dimensions+1)*i]);
	  }
	MTmake (root, temp, DataSize, Dimensions, 0);
	//free (temp);
	return root;
}

void
answerPointQuery    (int QuerySender,
	       int Dimensions,
	       coord_t * Data,
	       int DataSize,
	       MTNode * Root,
	       int MyRank,
	       double* data,
	       int * dataIdx,
	       coord_t ** *SelfQueryResult,
	       int *SelfQueryResultCount,
	       coord_t* LD, coord_t* RD)
{
	int i;
	int number_of_bytes = sizeof (coord_t) * 2 * Dimensions;
	coord_t **Res;
	int resSize=0;
	unsigned char* buff = malloc(number_of_bytes);
	if( QuerySender == MyRank ) {
	    memcpy(buff, LD, sizeof(coord_t)*Dimensions);
	    memcpy(buff + sizeof(coord_t)*Dimensions, RD, sizeof(coord_t)*Dimensions);
	    MPI_Bcast(buff, number_of_bytes, MPI_BYTE, QuerySender, MPI_COMM_WORLD);
	} 
	else {
	    MPI_Bcast(buff, number_of_bytes, MPI_BYTE, QuerySender, MPI_COMM_WORLD);
	    LD = (coord_t*) buff;
	    RD = (coord_t*) (buff + Dimensions* sizeof (coord_t));
	}
	MTQuery (Root, LD, RD,
	        &Res,
	        &resSize, Dimensions);
	MPI_Barrier(MPI_COMM_WORLD);
	free(buff);
	if (MyRank == QuerySender)
	{
	    (*SelfQueryResult) = Res;
	    (*SelfQueryResultCount) = resSize;
	    return;
	}
	MPI_Request req_null;
	MPI_Isend (&resSize, 1,
		MPI_INT, QuerySender, ANSWER_TAG-1, MPI_COMM_WORLD, &req_null);
	MPI_Request_free(&req_null);

	int real_buff_size = resSize * Dimensions;
	double *real_buff = malloc (real_buff_size * sizeof(double));
	double *real_buff_ptr = real_buff;
	int * idx_buff = malloc(resSize * sizeof(int));
	for (i = 0; i < resSize; i++)
	{
	    memcpy(real_buff_ptr, &data[Dimensions * Res[i][0]], sizeof(coord_t));
	    real_buff_ptr += Dimensions;
	    idx_buff[i] = dataIdx[Res[i][0]];
	}
	MPI_Ssend (real_buff,
		real_buff_size,
		MPI_DOUBLE, QuerySender, ANSWER_TAG+2*MyRank-1, MPI_COMM_WORLD);
	MPI_Ssend (idx_buff,
		resSize,
		MPI_INT, QuerySender, ANSWER_TAG+2*MyRank, MPI_COMM_WORLD);
	free(Res);
	free(real_buff);
	free(idx_buff);
}

void
answerProcessQuery(int QuerySender,
	       int Dimensions,
	       MTNode * Root,
	       int MyRank,
	       int ** Processes,
	       int * ProcessesCount,
	       coord_t* LD, coord_t* RD)
{
	int number_of_bytes = sizeof (coord_t) * 2 * Dimensions;
	int Res = 0, i;
	unsigned char* buff = malloc(number_of_bytes);

	if( QuerySender == MyRank ) {
	    memcpy(buff, LD, sizeof(coord_t)*Dimensions);
	    memcpy(buff + sizeof(coord_t)*Dimensions, RD, sizeof(coord_t)*Dimensions);
	    MPI_Bcast(buff, number_of_bytes, MPI_BYTE, QuerySender, MPI_COMM_WORLD);
	} 
	else {
	    MPI_Bcast(buff, number_of_bytes, MPI_BYTE, QuerySender, MPI_COMM_WORLD);
	    LD = (coord_t*) buff;
	    RD = (coord_t*) (buff + Dimensions* sizeof (coord_t));
	}
	free (buff);

	MTQuery2 (Root, LD, RD,
	        &Res, Dimensions);

	(*ProcessesCount)=0;
	MPI_Reduce(&Res, ProcessesCount, 1, MPI_INT, MPI_SUM, QuerySender, MPI_COMM_WORLD);
	if(Res==0 && MyRank!=QuerySender)
	return;
	(*Processes)=calloc(*ProcessesCount, sizeof(int));
	if(MyRank==QuerySender){
		for(i=0; i<((*ProcessesCount)-Res); i++){
			MPI_Recv(&((*Processes)[i]), 1, MPI_INT, MPI_ANY_SOURCE, ANSWER_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	}
	else{
		MPI_Send(&MyRank, 1, MPI_INT, QuerySender, ANSWER_TAG, MPI_COMM_WORLD);
	}
}

void
recvPointQuery	(double * *NewNeighbours,
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
		int MyRank)
{
	int i, j;
	int *cntBuffers = calloc (ProcessCount, sizeof (int));
	MPI_Request* req=calloc(ProcessCount-1, sizeof(MPI_Request));
	for (i = 0; i < ProcessCount; i++)
	  {
		if(i==MyRank)
		continue;
		MPI_Irecv(&cntBuffers[i], 1, MPI_INT, i, ANSWER_TAG-1, MPI_COMM_WORLD, &req[i-(i>MyRank)]);
	  }
	for(i=0; i<ProcessCount-1; i++){
		MPI_Waitany(ProcessCount-1, req, &j, MPI_STATUS_IGNORE);
		(*NewNeighboursSize) += cntBuffers[j];
	}
	(*NewNeighbours) = malloc((*NewNeighboursSize) * Dimensions * sizeof(double));
	(*NewNeighboursIdx) = malloc((*NewNeighboursSize) * sizeof(int));
	int newNeighboursPtr = 0;
	for (i = 0; i < ProcessCount; i++)
	  {
		  if (i == MyRank)
			  continue;
		  int buff_size = cntBuffers[i] * Dimensions;

		  MPI_Recv(&((*NewNeighbours)[newNeighboursPtr]), buff_size, MPI_DOUBLE, i, ANSWER_TAG+2*i-1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  newNeighboursPtr += buff_size;
		  MPI_Recv(&((*NewNeighboursIdx)[newNeighboursPtr/Dimensions]), cntBuffers[i], MPI_INT, i, ANSWER_TAG+2*i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		  // next points from input 
	  }
	free (cntBuffers);
	free (req);
	(*ResultsSize) = (SelfQueryResultCount + (*NewNeighboursSize));
	(*Results) = malloc ((*ResultsSize) * Dimensions * sizeof(double));
	(*ResultsIdx) = malloc ((*ResultsSize) * sizeof(int));

	memcpy((*Results), (*NewNeighbours), (*NewNeighboursSize) * Dimensions * sizeof(double));
	memcpy((*ResultsIdx), (*NewNeighboursIdx), (*NewNeighboursSize) * sizeof(int));

	for(i=0; i<SelfQueryResultCount; i++){
	    memcpy ((*Results)+(*NewNeighboursSize)+ i * Dimensions,
		    &data[SelfQueryResult[i][0]],
		    Dimensions * sizeof (double));
	    (*ResultsIdx)[i + (*NewNeighboursSize)] = dataIdx[SelfQueryResult[i][0]];
	}
}
