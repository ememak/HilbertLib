#include "MyTree.h"
#include <assert.h>
#include <stdio.h>

void
makeMTNode (MTNode * foo, int dimdiv, coord_t val)
{
	foo->dim = dimdiv;
	foo->val = val;
	foo->left = NULL, foo->right = NULL;
#ifdef MYTREEMINMAX
	foo->min = 0;
	foo->max = COORD_T_MAX;
#endif
}

void
coordinatesMINMAX (coord_t * *Data, int DataSize,
		   int Dim, coord_t * MIN, coord_t * MAX)
{
	int i;
	*MIN = Data[0][Dim+1];
	*MAX = Data[0][Dim+1];
	for (i = 1; i < DataSize; i++)
	  {
		  if (Data[i][Dim+1] > *MAX)
		    {
			    *MAX = Data[i][Dim+1];
		    }
		  else if (Data[i][Dim+1] < *MIN)
		    {
			    *MIN = Data[i][Dim+1];
		    }
	  }
}

void
MTmake (MTNode * Node,
	coord_t * *Data, int DataSize, int Dimensions, int ActDim)
{
	//printf("Wchodze do %p\n",Node);
	int i;
	/*printf("Data\n");
	   for(i=0;i<DataSize;i++) {
	   printf("%d:  ",i);
	   int j;
	   for(j=0;j<Dimensions;j++) {
	   printf("%d ", Data[i]->coordinates[j]);
	   }
	   printf("\n");
	   } */
	if (DataSize == 1)
	  {
		  Node->right = NULL;
		  Node->left = calloc (DataSize, sizeof (coord_t *));
		  int i;
		  for (i = 0; i < DataSize; i++)
			  ((coord_t **) (Node->left))[i] = Data[i];
		  Node->val = DataSize;

	  }
	if (DataSize <= 1)
	  {
		  free (Data);
		  return;
	  }


	coord_t **leftData = NULL;
	coord_t **rightData = NULL;
	coord_t pivot;
	int countSmaller;
	coord_t MIN, MAX;
	int started = 0;
	while (1)
	  {
		  if (started == Dimensions)
		    {		// Everything is the same
			    Node->right = NULL;
			    Node->left = calloc(DataSize, sizeof (coord_t *));
			    int i;
			    for (i = 0; i < DataSize; i++)
				    ((coord_t **)(Node->left))[i] = Data[i];
			    Node->val = DataSize;
			    free (Data);
			    return;
		    }
		  pivot = Data[rand()%DataSize][ActDim+1];
		  countSmaller = 0;
		  for (i = 0; i < DataSize; i++)
		    {
			    if (Data[i][ActDim+1] <= pivot)
			      {
				      countSmaller++;
			      }
		    }
		  if (countSmaller == 0
		      || DataSize - countSmaller == 0)
		    {
			    coordinatesMINMAX
				    (Data, DataSize, ActDim, &MIN, &MAX);
			    if (MIN != MAX)
			      {
				      pivot = (MIN + MAX) / 2;
				      break;
			      }
			    started++;
			    ActDim += 1;
			    ActDim %= Dimensions;
		    }
		  else
		    {
			    break;
		    }
	  }

	//printf("pivot : %d Node->dim = %d\n",pivot,ActDim);
	Node->dim = ActDim;
	Node->val = pivot;


	countSmaller = 0;
#ifdef MYTREEMINMAX
	Node->min = Data[0][ActDim+1];
	Node->max = Data[0][ActDim+1];
#endif
	for (i = 0; i < DataSize; i++)
	  {
		  if (Data[i][ActDim+1] <= pivot)
		    {
			    countSmaller++;
		    }
#ifdef MYTREEMINMAX
		  if (Data[i][ActDim+1] < Node->min)
		    {
			    Node->min = Data[i][ActDim+1];
		    }
		  if (Data[i][ActDim+1] > Node->max)
		    {
			    Node->max = Data[i][ActDim+1];
		    }
#endif

	  }
	assert (countSmaller != 0);
	assert ((DataSize - countSmaller) != 0);
	leftData = calloc (countSmaller, sizeof (coord_t *));
	rightData =
		calloc (DataSize - countSmaller, sizeof (coord_t *));
	int leftDataPtr = 0, rightDataPtr = 0;
	for (i = 0; i < DataSize; i++)
	  {
		  if (Data[i][ActDim+1] <= pivot)
			{
				leftData[leftDataPtr] = Data[i];
				leftDataPtr += 1;
			}
		  else
		    {
			    rightData[rightDataPtr] = Data[i];
			    rightDataPtr += 1;
		    }
	  }
	Node->left = calloc (1, sizeof (MTNode));
	Node->right = calloc (1, sizeof (MTNode));
	makeMTNode (Node->left, 0, 0);
	makeMTNode (Node->right, 0, 0);
	free (Data);
	MTmake (Node->left, leftData,
		countSmaller, Dimensions, (ActDim + 1) % Dimensions);
	MTmake (Node->right, rightData,
		DataSize - countSmaller,
		Dimensions, (ActDim + 1) % Dimensions);
}

void
MTDelete (MTNode * Node)
{
	if (Node->left != NULL && Node->right != NULL)
	  {
		  MTDelete (Node->left);
		  MTDelete (Node->right);
		  free (Node->left);
		  free (Node->right);
	  }
	else
	  {
		  free ((coord_t **) Node->left);
	  }
}

void
MTQueryLocal (MTNode * Node,
	      coord_t * LD,
	      coord_t * RD, PtrVector * vec, int Dimensions)
{
	if (Node->right == NULL)
	  {
		  if (Node->left == NULL)
			  return;
		  int i;
		  for (i = 0; i < Dimensions; i++)
		    {
			    coord_t x =(((coord_t **)
					 (Node->left))[0])[i+1];
			    if ((x < LD[i]) || (x > RD[i]))
				    return;
		    }
		  for (i = 0; i < (Node->val); i++)
		    {
			    PtrVectorPB (vec,
					 ((coord_t **)
					  Node->left)[i]);
		    }
		  return;
	  }
#ifdef MYTREEMINMAX
	if (Node->max < LD[Node->dim] || Node->min > RD[Node->dim])
		return;
#endif
	if (LD[Node->dim] <= Node->val)
	  {
		  MTQueryLocal (Node->left, LD, RD, vec, Dimensions);
	  }
	if (RD[Node->dim] > Node->val)
	  {
		  MTQueryLocal (Node->right, LD, RD, vec, Dimensions);
	  }
}

void
MTQueryLocal2 (MTNode * Node,
	      coord_t * LD,
	      coord_t * RD, int * r, int Dimensions)
{
	if (Node->right == NULL)
	  {
		  if (Node->left == NULL)
		      return;
		  int i;
		  for (i = 0; i < Dimensions; i++)
		    {
		        coord_t x =(((coord_t **)
				 (Node->left))[0])[i+1];
			    if ((x < LD[i]) || (x > RD[i]))
			    return;
		    }
		    (*r)=1;
		  return;
	  }
#ifdef MYTREEMINMAX
	if (Node->max < LD[Node->dim] || Node->min > RD[Node->dim])
		return;
#endif
	if (LD[Node->dim] <= Node->val)
	  {
		  MTQueryLocal2 (Node->left, LD, RD, r, Dimensions);
	  }
	if (RD[Node->dim] > Node->val)
	  {
		  MTQueryLocal2 (Node->right, LD, RD, r, Dimensions);
	  }
}

void
MTQuery (MTNode * Node,
	 coord_t * LD,
	 coord_t * RD, coord_t ** *Res, int *ResSize, int Dimensions)
{
	PtrVector vec;
	makePtrVector (&vec);
	MTQueryLocal (Node, LD, RD, &vec, Dimensions);
	*Res = calloc (vec.size, sizeof (coord_t *));
	*ResSize = vec.size;
	memcpy (*Res, vec.arr, (*ResSize) * sizeof (coord_t *));
	PtrVectorDeallocate (&vec);
}

void
MTQuery2 (MTNode * Node,
	 coord_t * LD,  coord_t * RD, int * Res, int Dimensions)
{
	int w = 0;
	MTQueryLocal2 (Node, LD, RD, &w, Dimensions);
	(*Res)=w;
}
