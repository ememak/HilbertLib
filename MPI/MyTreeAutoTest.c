#include "MyTree.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#define DIMENSIONS 5
#define DATASIZE 1000000
#define TESTSCNTMAX 100

int jest[DATASIZE+5];

int main (int argc, char *argv[])
{
	srand (time (NULL));
	MTNode Root;
	makeMTNode (&Root, 0, 0);
	coord_t *Data = malloc (DATASIZE*(DIMENSIONS+1)*sizeof (coord_t));
	coord_t **DataPtr = malloc (DATASIZE*sizeof (coord_t *));
	int i,j;
	printf("Points-Start\n");
	for (i = 0; i < DATASIZE; i++)
	  {
		  //printf("Data[%d] = ",i);
		  Data[(DIMENSIONS+1)*i]=i;
		  for (j = 0; j < DIMENSIONS; j++)
		    {
			    Data[(DIMENSIONS+1)*i+j+1] = rand () % COORD_T_MAX;
		//	    printf ("%d ", Data[i].coordinates[j]);
		    }
		  DataPtr[i] = &Data[(DIMENSIONS+1)*i];
		  //printf ("\n");
	  }
	printf("Points-End\n");
	MTmake (&Root, DataPtr, DATASIZE, DIMENSIONS, 0);
	int tests = TESTSCNTMAX;
	int it;
	for(it=0;it<tests;it++) {
		for(i=0; i<DATASIZE; i++)
		jest[i]=0;
		printf("TEST %d\n",it);
		coord_t *LD = calloc (DIMENSIONS, sizeof (coord_t));
		coord_t *RD = calloc (DIMENSIONS, sizeof (coord_t));
		//printf ("OKREŚL OBSZAR POSZUKIWAN (l1,r1), (l2,r2) ..., (ln,rn)\n");
		for (i = 0; i < DIMENSIONS; i++)
		  {
			  LD[i] = rand()%COORD_T_MAX;
			  RD[i] = rand()%COORD_T_MAX;
			  if(RD[i] < LD[i]) {
				RD[i] = LD[i]^RD[i];
				LD[i] = LD[i]^RD[i];
				RD[i] = LD[i]^RD[i];
			  }
			
		  }
		coord_t **Res = NULL;
		int ResSize = 0;
		MTQuery (&Root, LD, RD, &Res, &ResSize, DIMENSIONS);
		for(i = 0; i < DATASIZE; i++) {
			int good = 1;
			for(j=0;j<DIMENSIONS;j++) {
				if((LD[j] > Data[(DIMENSIONS+1)*i+j+1]) || (RD[j] < Data[(DIMENSIONS+1)*i+j+1])) {
					good = 0;
					break;
				}
			}
			if(good == 1)
				jest[i] = 1;
			else
				jest[i] = 0;
		}
		int error = 0;
		int how_many_points = 0;
		for(i=0;i<ResSize;i++) {
			jest[Res[i][0]]+=2;
			how_many_points += 1;
		}
		int bad_point;
		for(i=0;i<ResSize;i++) {
			if(jest[i] == 1) {
				error = 1;
				printf("Brute-force solution found some other point(%d) lying in queried space.\n", i);
				bad_point = i;
				break;
			} else if(jest[i] == 2) {
				error = 1;
				printf("Tree solution found some other point(%d) lying in queried space.\n",i);
				bad_point = i;
				break;
			}
		}
		if(error == 1) {
			printf("WRONG ANSWER !!!\n");
			for(i=0;i<DIMENSIONS;i++) {
			  printf ("Dimension %d : ", i);
			  printf ("%d %d %d\n", LD[i], RD[i], Data[bad_point*(DIMENSIONS+1)+i+1]);
			  if(!(LD[i] <= Data[bad_point*(DIMENSIONS+1)+i+1] && RD[i] >= Data[bad_point*(DIMENSIONS+1)+i+1]))
				  printf("%dth dimension requirement is not satisfied\n\n\n",i);
			}
			printf("\n");
		} else {
			printf("OK %d points found\n",how_many_points);
		}
		free (Res);
		free (LD);
		free (RD);
	}
	free (Data);
	MTDelete (&Root);
	return 0;
}
