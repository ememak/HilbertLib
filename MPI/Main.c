#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "AxesTranspose.h"
#include "HilbertLib.h"
#define calloc(a,b) (a==0 ? NULL : calloc(a,b))

#define ROOT 0
#define DIMENSIONS 3
#define BITS_PRECISION 25
#define NR_OF_QUERIES 1
const int SWITCH=1;


coord_t
rand_coord_t ()
{
	return (RAND_MAX*rand ()+rand()) % (1ll << BITS_PRECISION - 1);
}

int
main (int argc, char *argv[])
{
	FILE *File = fopen ("./out/MainOutput0", "w");
	
	coord_t *MyPoints=NULL;
	int MyPointsCount=0;
	int i, j, k;
	double * inp=NULL;
	double *mind=NULL, *maxd=NULL;
	coord_t* idx=NULL;
	//int seed = time(0)+rank+size;
	// Random Input Generation
	srand(time(0)*516056890);
	MyPointsCount = 1000000;
	for(k=0; k<100; k++){
	MyPoints = malloc((DIMENSIONS+1) * MyPointsCount * sizeof(coord_t));
	for(i=0; i<MyPointsCount; i++){
	    MyPoints[i*(DIMENSIONS+1)] = i;
	    MyPoints[i*(DIMENSIONS+1)+1] = rand_coord_t();
	    MyPoints[i*(DIMENSIONS+1)+2] = rand_coord_t();
	    MyPoints[i*(DIMENSIONS+1)+3] = rand_coord_t();
	    //MyPoints[i*(DIMENSIONS+1)+4] = rand_coord_t();
	    //MyPoints[i*(DIMENSIONS+1)+5] = rand_coord_t();
	}
	//Printing genereated points
	/*for (i = 0; i < MyPointsCount; i++)
	  {
		  fprintf (File, "Punkt #%d : ", i);
		  fprintf (File, "%d %d %d ", MyPoints[(DIMENSIONS+1)*i+1], MyPoints[(DIMENSIONS+1)*i+2], MyPoints[(DIMENSIONS+1)*i+3]);
		  fprintf (File, "\n");
	  }
	fflush(File);*/
	coord_t *NewData = NULL;
	coord_t *NewIdx = NULL;
	int NewDataCount = 0;
		HilbertLibPartition (// MyPoints is freed
				    MyPoints,
				    MyPointsCount,
				    ROOT,
				    DIMENSIONS,
				    BITS_PRECISION,
				    0, 1,
				    &NewData, &NewIdx,
				    &NewDataCount);
	/*coord_t *X=calloc(DIMENSIONS, sizeof(coord_t)), *Y=calloc(DIMENSIONS, sizeof(coord_t));
	for (i = 0; i < MyPointsCount; i++)
	  {
		  fprintf (File, "point:#%d @@@   ", i);
		  fprintf (File, "%d %d %d %lf\n", NewData[(DIMENSIONS+1)*i+1], NewData[(DIMENSIONS+1)*i+2], NewData[(DIMENSIONS+1)*i+3], GetHCoordinate(&NewData[(DIMENSIONS+1)*i+1], X, Y, BITS_PRECISION, DIMENSIONS));
	  }
	free(X);
	free(Y);*/
	free(NewData);
	}
	fclose(File);
	return 0;
	/*MTNode *MyTreeRoot = HilbertLibPrepareNodeForQueries (NewData,
							      NewDataCount,
							      DIMENSIONS);
	MDPoint **SelfQueryResult=NULL;
	int SelfQueryResultCount=0;
	coord_t *LD = calloc (DIMENSIONS, sizeof (coord_t));
	coord_t *RD = calloc (DIMENSIONS, sizeof (coord_t));
		
	if (rank == 0)
	  {
		double *LDD=NULL;
		double *RDD=NULL;
		if(SWITCH==2 || SWITCH==3){
		    LDD=calloc(DIMENSIONS, sizeof(double));
		    RDD=calloc(DIMENSIONS, sizeof(double));
		}
		for(j=0; j<NR_OF_QUERIES; j++){
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
		    if(SWITCH==2 || SWITCH==3){
			toPos(LD, LDD, mind, maxd, DIMENSIONS, BITS_PRECISION);
			toPos(RD, RDD, mind, maxd, DIMENSIONS, BITS_PRECISION);
			printf("Query: %d from process: %d\n", j, rank);
			for(i=0;i<DIMENSIONS;i++) {
			    printf("%d LDD = %lf, RDD = %lf\n",i,LDD[i],RDD[i]);
			}
			free(LDD);
			free(RDD);
		    }
		    else{
			printf("Query: %d from process: %d\n", j, rank);
			for(i=0;i<DIMENSIONS;i++) {
			    printf("%d LD = %d, RD = %d\n",i,LD[i],RD[i]);
			}
		    }
		}
		fflush(stdout);
	}
	int * procsRes=NULL;
	int procsResCount=0;
	answerPointQuery(0,
			DIMENSIONS,
			NewData,
			NewDataCount,
			MyTreeRoot,
			rank,
			&SelfQueryResult,
			&SelfQueryResultCount,
			LD, RD);
	MDPoint *NewNeighbours=NULL;
	int NewNeighboursSize = 0;
	MDPoint **Results=NULL;
	int ResultsSize=0;

	if (rank == 0)
	  {
		  recvPointQuery(&NewNeighbours,
				&NewNeighboursSize,
				&Results, &ResultsSize, 
				DIMENSIONS, size, 
				SelfQueryResult,
				SelfQueryResultCount,
				rank);
		  if(SWITCH==1){
			printf("ResultsSize: %d\nResults:\n", ResultsSize);
			for(i=0; i<ResultsSize; i++){
				for(j=0; j<DIMENSIONS; j++) {
					printf("%d ", Results[i]->coordinates[j]);
				}
				printf("\n");
			}
		  }
		  if(SWITCH==2 || SWITCH ==3){
			printf("ResultsSize: %d\nResults:\n", ResultsSize);
			for(i=0; i<ResultsSize; i++) {
			    for(j=0; j<DIMENSIONS; j++) {
				printf("%lf ", inp[DIMENSIONS * Results[i]->own_data_id + j]);
			    }
			printf("\n");
			}
		    free(inp);
		    free(idx);
		  }
	  }
	answerProcessQuery(0,
			DIMENSIONS,
			MyTreeRoot,
			rank,
			&procsRes,
			&procsResCount,
			LD, RD);
	if(rank==0){
	    printf("ProcsRes:\n");
	    for(i=0; i<procsResCount; i++){
		printf("%d ", procsRes[i]);
	    }
	    printf("\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(SWITCH ==2 || SWITCH==3){
	    free(mind);
	    free(maxd);
	}
	free(LD);
	free(RD);
	for (i = 0; i < NewDataCount; i++)
	  {
		MDPointRemove (&NewData[i]);
	  }
	free (NewData);
	for (i=0;i<NewNeighboursSize;i++) {
		MDPointRemove (&NewNeighbours[i]);
	}
	free (NewNeighbours);
	MTDelete (MyTreeRoot);
	for(i=0; i<SelfQueryResultCount; i++){
		MDPointRemove(SelfQueryResult[i]);
	}
	free (SelfQueryResult);
	if(rank==0){
		for(i=0; i<ResultsSize; i++){
		    MDPointRemove(Results[i]);
		}
	}
	free(Results);
	//fclose (File);
	MPI_Finalize ();

	return 0;
	/*int *newPointsCount = calloc(size,sizeof(int));
	   MPI_AlltoAll(&NewDataCount, 1, MPI_INT, newPointsCount, 1, MPI_INT, MPI_COMM_WORLD);

	   //fprintf("%d\n",NewDataCount);
	   MPI_File *fh;
	   MPI_File_open(MPI_COMM_WORLD, "mainResult.dat", MPI_MODE_RDWR, MPI_INFO_NULL, fh);
	   MPI_File_set_size(fh,0);

	   int myoffset = 0;
	   for(i=0;i<rank;i++)
	   myoffset += newPointsCount[i];
	   MPI_File_seek(offset,1,MPI_SEEK_SET);
	   char *line = calloc(128,sizeof(char));

	   for(i=0;i<NewDataCount;i++) {

	   for(j=0;j<DIMENSIONS;j++) {
	   sfprintf(line, "%d ", NewData[i].coordinates[j]);
	   MPI_File_write(
	   fh,


	   }
	   //fprintf("data_id = %d",NewData[i].own_data_id);
	   sfprintf(line,"%d", rank);
	   sfprintf(line,"\n");
	   fprintf("%s",line);
	   } */

	//free(line);*/
	//MPI_File_close(fh);
}

