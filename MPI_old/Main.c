#include <assert.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "AxesTranspose.h"
#include "HilbertLib.h"
#include "MDPoint.h"
#include "MyTree.h"
#include "Pair.h"
#include "ToGrid.h"
#define calloc(a,b) (a==0 ? NULL : calloc(a,b))

#define ROOT 0
#define DIMENSIONS 3
#define BITS_PRECISION 15
#define NR_OF_QUERIES 1
const int SWITCH=1;


coord_t
rand_coord_t ()
{
	return rand () % (1ll << BITS_PRECISION - 1);
}

int
main (int argc, char *argv[])
{
//	int xxx = atoi(argv[1]);
	int xxx=1000;
	// Initialization
	MPI_Init (&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);
	double begin = MPI_Wtime();

	int rank, size;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	FILE *File, *in;
	char *arr = calloc (100, sizeof (char));;
	sprintf (arr, "./out/MainOutput%d", rank);
	File = fopen (arr, "w");
	free (arr);
	MPI_Comm_size (MPI_COMM_WORLD, &size);
	MDPoint *MyPoints=NULL;
	int MyPointsCount=0;
	int i, j, q;
	double ** inp=NULL;
	double *mind=NULL, *maxd=NULL;
	tag_t* idx=NULL;
	// Random Input Generation
	if (SWITCH==1)
	  {
		  MyPointsCount = xxx;
		  //srand (time (0) + rank + size);
		  srand(50000 + rank + size);
		  MyPoints = calloc (MyPointsCount, sizeof (MDPoint));
		  for (i = 0; i < MyPointsCount; i++)
		    {
			    make_MDPoint (&MyPoints[i], DIMENSIONS);
			    for (j = 0; j < DIMENSIONS; j++)
			      {
				      MyPoints[i].coordinates
					      [j] = rand_coord_t ();
			      }
			    MyPoints[i].own_data_id = MyPoints[i].coordinates[0] ^ 17;
		    }
	  }
	if(SWITCH==2)
	  {
	  if(rank==0){
		  MyPointsCount = 11659;
		  in=fopen("./tests/vessel.txt", "r");
		  inp=calloc(MyPointsCount, sizeof(double*));
		  idx=calloc(MyPointsCount, sizeof(int));
		  MyPoints = calloc (MyPointsCount, sizeof (MDPoint));
		  for (i = 0; i < MyPointsCount; i++)
		    {
			inp[i]=calloc(DIMENSIONS, sizeof(double));
			for (j = 0; j < DIMENSIONS; j++)
			  {
			    fscanf (in, "%lf ", &inp[i][j]);
			  }
			idx[i]=(tag_t)i;
		    }
		  fclose(in);
		  }
		getMinMax(inp, MyPointsCount, &mind, &maxd, DIMENSIONS);
		allToGrid(inp, idx, mind, maxd, MyPointsCount, &MyPoints, DIMENSIONS, BITS_PRECISION);
	  }
	if(SWITCH==3){
	    in=fopen("./tests/dane.bin", "rb");
	    fseek(in, 0, SEEK_END);
	    int filesize=ftell(in);
	    
	    if(sizeof(double) * DIMENSIONS * rank * 300 < filesize) //300 - how many points one process will have at the begining
	    {
		if(filesize < sizeof(double) * rank * 301 * DIMENSIONS){
		    MyPointsCount = filesize - sizeof(double) * rank * 300 * DIMENSIONS;
		}
		else{
		    MyPointsCount = 300;
		}
		fseek(in, sizeof(double)*rank*DIMENSIONS*300, SEEK_SET);
		inp=calloc(MyPointsCount, sizeof(double *));
		idx=calloc(MyPointsCount, sizeof(int));

		for(i=0; i<MyPointsCount; i++){
		    inp[i] = calloc(DIMENSIONS, sizeof(double));
		    idx[i] = (tag_t)(i+rank*300);
		    fread(inp[i], sizeof(double), 3, in);
		}
	    }
	    else{
		MyPointsCount = 0;
	    }
	    MyPoints=calloc(MyPointsCount, sizeof(MDPoint));
	    fclose(in);
	    getMinMax(inp, MyPointsCount, &mind, &maxd, DIMENSIONS);
	    allToGrid(inp, idx, mind, maxd, MyPointsCount, &MyPoints, DIMENSIONS, BITS_PRECISION);
	}

	//Printing genereated points
	for (i = 0; i < MyPointsCount; i++)
	  {
		  fprintf (File, "Punkt #%d : ", i);
		  for (j = 0; j < DIMENSIONS; j++)
		    {
			    fprintf (File, "%u ",
				     MyPoints[i].coordinates[j]);
		    }
		  fprintf (File, "\n");
	  }
	fflush(File);
	MDPoint *NewData = NULL;
	int NewDataCount = 0;
	HilbertLibPartition (	// MyPoints is freed
				    MyPoints,
				    MyPointsCount,
				    ROOT,
				    DIMENSIONS,
				    BITS_PRECISION,
				    rank,
				    size, &NewData, &NewDataCount);
	fprintf (File, "NewDataCount[%d] = %d\n", rank, NewDataCount);
	for (i = 0; i < NewDataCount; i++)
	  {

		  fprintf (File, "rank:#%d point:#%d @@@   ", rank, i);
		  for (j = 0; j < DIMENSIONS; j++)
		    {
			    fprintf (File, "%d ", NewData[i].coordinates[j]);
		    }
		  fprintf (File, "\n");
	  }
	MPI_Barrier(MPI_COMM_WORLD);
	double end = MPI_Wtime();
	if(rank == 0)
		printf("OVERALL TIME : %f\n",end-begin);
	//MPI_Barrier(MPI_COMM_WORLD);// do tego momentu kod jest chyba sprawdzony, dalej czesc wyglada dziwnie
	MTNode *MyTreeRoot = HilbertLibPrepareNodeForQueries (NewData,
							      NewDataCount,
							      DIMENSIONS);
	int *number_of_queries;
	MPI_Request *req1=calloc(size, sizeof(MPI_Request));
	unsigned char *BigBuff = NULL;
	MDPoint ***SelfQueriesResult=NULL;
	int *SelfQueriesResultCount=NULL;
	int *SelfQueriesRank=NULL;
	int SelfQueriesCount=0;
	int QueryCounter=0; // +1 in every sendQuery call

	//if (rank == 0 || rank==1)
	  //{
		exchangeNumberOfQueries (&number_of_queries, size, NR_OF_QUERIES);
		coord_t *LD = calloc (DIMENSIONS, sizeof (coord_t));
		coord_t *RD = calloc (DIMENSIONS, sizeof (coord_t));
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
		     sendQuery 	(LD,
				RD,
			        size,
			        DIMENSIONS,
			        rank, &QueryCounter, /*&req1, */&BigBuff, NR_OF_QUERIES);
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
			/*printf("Query: %d from process: %d\n", j, rank);
			for(i=0;i<DIMENSIONS;i++) {
			    printf("%d LD = %d, RD = %d\n",i,LD[i],RD[i]);
			}*/
		    }
		}
		free(LD);
		free(RD);
		printf("free %p\n",&req1);
		fflush(stdout);
	  /*}
	else
	  {
		  exchangeNumberOfQueries (&number_of_queries, size, 0);
	  }*/
	answerQueries	(size,
			DIMENSIONS,
			NewData,
			NewDataCount,
			MyTreeRoot,
			number_of_queries,
			rank,
			&SelfQueriesResult,
			&SelfQueriesResultCount,
			&SelfQueriesRank, &SelfQueriesCount,
			BigBuff, req1);
	//MPI_Barrier(MPI_COMM_WORLD);
	MDPoint *NewNeighbours=NULL;
	int NewNeighboursSize = 0;
	MDPoint ***Results=NULL;
	int *ResultsSize=NULL;

	/*if (rank == 0)
	  {*/
		  recvQueries (&NewNeighbours,
			       &NewNeighboursSize,
			       &Results, &ResultsSize, 
				DIMENSIONS, size, NR_OF_QUERIES,
				number_of_queries[rank],
				SelfQueriesResult,
				SelfQueriesResultCount,
				SelfQueriesRank,
				rank);
		  if(SWITCH==1){
		    /*for(q=0; q<NR_OF_QUERIES; q++) {
			  printf("Wyniki %d, jest ich %d\n", q, ResultsSize[q]);
			  for(i=0; i<ResultsSize[q]; i++) {
			    if(i%10000==0){
				for(j=0; j<DIMENSIONS; j++) {
					printf("%d ", Results[q][i]->coordinates[j]);
				}
				printf("\n");
				}
			    }
		    }*/
		  }
		  if(SWITCH==2 || SWITCH ==3){
		    for(q=0; q<NR_OF_QUERIES; q++) {
			printf("Wyniki %d, jest ich %d\n", q, ResultsSize[q]);
			for(i=0; i<ResultsSize[q]; i++) {
			    for(j=0; j<DIMENSIONS; j++) {
				printf("%lf ", inp[Results[q][i]->own_data_id][j]);
			    }
			printf("\n");
			}
		    }
		    for(i=0; i<MyPointsCount; i++){
			free(inp[i]);
		    }
		    free(inp);
		    free(idx);
		  }
	  //}
	for(i=0, j=0; i<size; i++){
	    if(number_of_queries[i]!=0 && i!=rank)
	    j++;
	}
	MPI_Request *req2=calloc(j, sizeof(MPI_Request));
	for(i=0, j=0; i<size; i++){
	    if(number_of_queries[i]!=0 && i!=rank)
	    req2[j++]=req1[i];
	}
	MPI_Waitall(j, req2, MPI_STATUSES_IGNORE);
	MPI_Barrier(MPI_COMM_WORLD);
	if(SWITCH ==2 || SWITCH==3){
	    free(mind);
	    free(maxd);
	}
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
	free (BigBuff);
	for(j=0; j<number_of_queries[rank]; j++){
	    for(int k=0; k<SelfQueriesResultCount[j]; k++){
		MDPointRemove(SelfQueriesResult[j][k]);
	    }
	    free(SelfQueriesResult[j]);
	}
	free (SelfQueriesResult);
	free (SelfQueriesResultCount);
	free (SelfQueriesRank);
	free (number_of_queries);
	//if(rank==0){	//probably bad way to clean it
	    for(q=0; q<NR_OF_QUERIES; q++){
		for(i=0; i<ResultsSize[q]; i++){
		    MDPointRemove(Results[q][i]);
		}
		free(Results[q]);
	    }
	//}
	free(Results);
	free (ResultsSize);
	fclose (File);
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

