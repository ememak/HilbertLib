#include <zoltan.h>

struct vector3d{
  float x;
  float y;
  float z;
};

struct particle_data{
  ZOLTAN_ID_TYPE gid;
  struct vector3d position;
};

struct particle_data *particles;

struct export_list_data{
  int particle;
  int proc;
};

int rank;
int size;

int niter;
int np;	// global number of particles
int lnp; // local (per process) number of particles

float boxsize;
float epsilon;

#define G 6.674
