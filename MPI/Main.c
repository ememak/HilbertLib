#include <assert.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "AxesTranspose.h"
#include "HilbertLib.h"
#include "MyTree.h"
#include "ToGrid.h"
#define calloc(a,b) (a==0 ? NULL : calloc(a,b))

#define ROOT 0
#define DIMENSIONS 3
#define BITS_PRECISION 25
#define NR_OF_QUERIES 1


coord_t
rand_coord_t ()
{
	return rand () % (1ll << BITS_PRECISION - 1);
}

// Code used for tests, for testing query loop that is calculated tests times should be deleted (and freeing variables should be done later)
int
main (int argc, char *argv[])
{
	int xxx, suf, seed, tests;
	//int xxx=10000;
	// Initialization
	MPI_Init (&argc, &argv);

	int rank, size;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	
	if(argc<5){
	    if(rank==0)
		fprintf(stderr, "Parameters are missing.\n");
	    return 1;
	}
	else{
	    xxx = atoi(argv[1]);
	    seed = atoi(argv[2]) * (rank+1);
	    tests = atoi(argv[3]);
	    suf = atoi(argv[4]);
	}
	
	FILE *File;
	char *arr = calloc (100, sizeof (char));;
	sprintf (arr, "./out/MainOutput%d", suf);
	File = fopen (arr, "w");
	free (arr);
	
	coord_t *MyPoints=NULL;
	int MyPointsCount=0;
	int i, j, k;
	double * inp=NULL;
	double *mind=NULL, *maxd=NULL;
	int* idx=NULL;
	double time = 0.0;
	//int seed = time(0)+rank+size;
	// Random Input Generation
	for(k=0; k<tests; k++){
	srand(seed);
	MPI_Barrier(MPI_COMM_WORLD);
	double begin = MPI_Wtime();
	MyPointsCount = xxx;
	inp = malloc(DIMENSIONS * MyPointsCount * sizeof(double));
	idx = malloc(MyPointsCount * sizeof(coord_t));
	for(i=0; i<MyPointsCount; i++){
	    inp[DIMENSIONS * i] = (double)(rand())/(double)(RAND_MAX);
	    inp[DIMENSIONS * i + 1] = (double)(rand())/(double)(RAND_MAX);
	    inp[DIMENSIONS * i + 2] = (double)(rand())/(double)(RAND_MAX);
	    idx[i] = (rank*MyPointsCount+i);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
	    fprintf(File, "End of point generation, time elapsed: %lf\n", MPI_Wtime()-begin);

	getMinMax(inp, MyPointsCount, &mind, &maxd, DIMENSIONS);
	allToGrid(inp, mind, maxd, &MyPoints, MyPointsCount, DIMENSIONS, BITS_PRECISION);

	MPI_Barrier(MPI_COMM_WORLD);
	double postbegin=MPI_Wtime();
	if(rank==0)
	    fprintf(File, "End of point conversion, time elapsed: %lf\n", postbegin-begin);

	//Printing genereated points
	/*for (i = 0; i < MyPointsCount; i++)
	  {
		  fprintf (File, "Punkt #%d : ", i);
		  fprintf (File, "%d %d %d ", MyPoints[(DIMENSIONS+1)*i+1], MyPoints[(DIMENSIONS+1)*i+2], MyPoints[(DIMENSIONS+1)*i+3]);
		  fprintf (File, "\n");
	  }
	fflush(File);*/

	double premid, postmid;
	double *NewData = NULL;
	int *NewIdx = NULL;
	int NewDataCount = 0;
	premid = MPI_Wtime();
	HilbertLibPartition (// MyPoints is freed
			    MyPoints,
			    MyPointsCount,
			    ROOT,
			    DIMENSIONS,
			    BITS_PRECISION,
			    rank, size,
			    inp, idx,
			    &NewData, &NewIdx,
			    &NewDataCount);

	allToGrid(NewData, mind, maxd, &MyPoints, NewDataCount, DIMENSIONS, BITS_PRECISION); // conversion of new points to grid

	MPI_Barrier(MPI_COMM_WORLD);
	postmid=MPI_Wtime();
	if(rank==0)
	    fprintf(File, "Partition end, time elapsed: %lf\n", postmid-premid);
	time+=postmid-premid;
	free(MyPoints);
	free(NewData);
	free(mind);
	free(maxd);
	}

	/*fprintf (File, "NewDataCount[%d] = %d\n", rank, NewDataCount);
	for (i = 0; i < NewDataCount; i++)
	  {
		  fprintf (File, "rank:#%d point:#%d @@@   ", rank, i);
		  fprintf (File, "%lf %lf %lf %lf\n", NewData[DIMENSIONS*i], NewData[DIMENSIONS*i+1], NewData[DIMENSIONS*i+2], GetHCoordinate(&MyPoints[(DIMENSIONS+1)*i+1], X, Y, BITS_PRECISION, DIMENSIONS));
	  }*/

	MPI_Barrier(MPI_COMM_WORLD);
	double end = MPI_Wtime();
	if(rank == 0)
		fprintf(File, "AVG TIME: %lf\n", time/tests);

	/*MTNode *MyTreeRoot = HilbertLibPrepareNodeForQueries (MyPoints,
							      NewDataCount,
							      DIMENSIONS);

	coord_t **SelfQueryResult=NULL;
	int SelfQueryResultCount=0;
	coord_t *LD = calloc (DIMENSIONS, sizeof (coord_t));
	coord_t *RD = calloc (DIMENSIONS, sizeof (coord_t));

	if (rank == 0)
	  {
		double *LDD = calloc(DIMENSIONS, sizeof(double));
		double *RDD = calloc(DIMENSIONS, sizeof(double));
		    for (i = 0; i < DIMENSIONS; i++)
		    {
			LD[i] = rand_coord_t ();
			RD[i] = rand_coord_t ();
			if (LD[i] > RD[i])
			    {
				LD[i] ^= RD[i];
				RD[i] ^= LD[i];
				LD[i] ^= RD[i];
			    }
			assert(LD[i]<=RD[i]);
			}
			toPos(LD, LDD, mind, maxd, DIMENSIONS, BITS_PRECISION);
			toPos(RD, RDD, mind, maxd, DIMENSIONS, BITS_PRECISION);
			printf("Query from process: %d\n", rank);
			for(i=0;i<DIMENSIONS;i++) {
			    printf("%d LDD = %lf, RDD = %lf\n",i,LDD[i],RDD[i]);
			}
			free(LDD);
			free(RDD);
		fflush(stdout);
	}
	double befpointquery=MPI_Wtime();
	int * procsRes=NULL;
	int procsResCount=0;
	/*answerPointQuery(0,
			DIMENSIONS,
			MyPoints,
			NewDataCount,
			MyTreeRoot,
			rank,
			NewData,
			NewIdx,
			&SelfQueryResult,
			&SelfQueryResultCount,
			LD, RD);
	double *NewNeighbours=NULL;
	int * NewNeighboursIdx=NULL;
	int NewNeighboursSize = 0;
	double *Results=NULL;
	int * ResultsIdx=NULL;
	int ResultsSize=0;

	if (rank == 0)
	  {
		  recvPointQuery(&NewNeighbours,
				&NewNeighboursIdx,
				&NewNeighboursSize,
				&Results, &ResultsIdx, &ResultsSize,
				DIMENSIONS, size,
				NewData, NewIdx,
				SelfQueryResult,
				SelfQueryResultCount,
				rank);
		  printf("PointQuery end, time elapsed: %lf\n", MPI_Wtime()-befpointquery);
		  printf("ResultsSize: %d\n", ResultsSize);
		  /*for(i=0; i<ResultsSize; i++){
		    printf("Point %d: ", ResultsIdx[i]);
		    for(j=0; j<DIMENSIONS; j++){
			printf("%lf ", Results[i*DIMENSIONS+j]);
		    }
		    printf("\n");
		  }*/
	  //}
	/*double befquery = MPI_Wtime();
	answerProcessQuery(0,
			DIMENSIONS,
			MyTreeRoot,
			rank,
			&procsRes,
			&procsResCount,
			LD, RD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
	printf("Query Time: %lf\n", MPI_Wtime()-befquery);
	if(rank==0){
	    printf("ProcsRes:\n");
	    for(i=0; i<procsResCount; i++){
		printf("%d ", procsRes[i]);
	    }
	    printf("\n");
	}*/
	//free(MyPoints);
	//free(procsRes);
	//free(mind);
	//free(maxd);
	//free(LD);
	//free(RD);
	//free (NewData);
	//free (NewNeighbours);
	//MTDelete (MyTreeRoot);
	//free (SelfQueryResult);
	//free(Results);
	//free(ResultsIdx);
	//fclose (File);
	MPI_Finalize ();

	return 0;
}

