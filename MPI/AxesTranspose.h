#ifndef HLAXESTRANSPOSE_H_INCLUDED
#define HLAXESTRANSPOSE_H_INCLUDED
#include<limits.h>
#include<stdlib.h>
#include<string.h>

typedef unsigned int coord_t;
#define COORD_T_MAX UINT_MAX
#define MPI_COORD_T MPI_UNSIGNED

//unsigned int hilpos_t
/*
typedef unsigned int hilpos_t;
#define MPI_HILPOS_T MPI_UNSIGNED
#define HILPOS_T_MAX_VALUE ((hilpos_t)1<<31)
#define HILPOS_EPS 1
//change to 2.0
#define HILPOS_BS_2 2
//change to 1.0
#define HILPOS_BS_1 1 
//change to HILPOS_EPS
#define HILPOS_MARGIN 0*/

//double hilpos_t
typedef double hilpos_t;
#define MPI_HILPOS_T MPI_DOUBLE
#define HILPOS_T_MAX_VALUE ((hilpos_t)1e200)
#define HILPOS_EPS ((hilpos_t)(1.0e-10))
#define HILPOS_BS_1 ((hilpos_t)0)
#define HILPOS_BS_2 ((hilpos_t)2)


static inline hilpos_t
GetHCoordinate (coord_t * Z, coord_t * X, coord_t * Y, int b, int n)
{				// position,bits,dimensions, Z,Y is used
	memcpy (X, Z, sizeof (coord_t) * n);
	coord_t M = (1 << (b - 1)), P, Q, t;
	int i, j, tmp;
	for (Q = M, j = b - 1; Q > 1; Q >>= 1, j--)
	  {
		  P = Q - 1;
		  for (i = 0; i < n; i++)
			  /*if ((X[i] & Q) > 0)
				  X[0] ^= P;
			  else
			    {
				    t = ((X[0] ^ X[i]) & P);
				    X[0] ^= t;
				    X[i] ^= t;
			    }*/ //old, more readable (but slower) version of code below
			  {
				    tmp = (X[i] & Q) >> j;
				    t = ((1 ^ tmp) * ((X[0] ^ X[i]) & P));
				    X[0] ^= ((tmp * P) ^ t);
				    X[i] ^= t;
			  }
	  }
	for (i = 1; i < n; i++)
		X[i] ^= X[i - 1];
	t = 0;
	for (j = b-1; j > 0;j--)
		{
		    t ^= ((X[n-1] & (1 << j)) >> j)*((1 << j) - 1);
		}
	/*for(Q=M; Q>1; Q>>=1){ // old, more readable version of loob above
	    if ((X[n - 1] & Q) > 0)
		t ^= (Q - 1);
	}*/
	for (i = 0; i < n; i++)
		X[i] ^= t;
	int wsk = 0;
	int wskbits = 0;
	memset(Y, 0, n * sizeof(coord_t));
	for (i = b - 1; i >= 0; i--)
	  {
		  for (j = 0; j < n; j++)
		    {
			    //Y[wsk] <<= 1; //use this without this last bitshift from line below
			    Y[wsk] += ((X[j] >> i) & 1) << (b - 1 - wskbits);
			    wskbits++;
			    if (wskbits == b)
			      {
				      wskbits = 0;
				      wsk++;
			      }
		    }
	  }

	hilpos_t res = 0;
	hilpos_t akt_val = 1;
	for (i = n - 1; i >= 0; i--)
	  {
		  res += akt_val * Y[i];
		  akt_val *= (hilpos_t) (1 << b);
	  }

	return res;
}

#endif
