#include <stdio.h>
#include <stdlib.h>
#include <particles.h>
#include <mpi.h>

/* This function generates particle data */
int generateParticles(int np, int seed) {

  int p;

  srand(seed);

  lnp=np;

  particles = (struct particle_data*) malloc(lnp*2*sizeof(struct particle_data));

  for(p=0;p<lnp;p++) {

    particles[p].gid=(ZOLTAN_ID_TYPE)(rank*lnp+p);

    particles[p].position.x=(double)(rand())/(double)(RAND_MAX);
    particles[p].position.y=(double)(rand())/(double)(RAND_MAX);
    particles[p].position.z=(double)(rand())/(double)(RAND_MAX);

  }

  return 0;

}


