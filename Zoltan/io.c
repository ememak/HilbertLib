#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <io.h>
#include <particles.h>

int write_particles(int state) {

  int i,j,k;
  struct vector3d vector_field[lnp];
  float scalar_float_field[lnp];
  int scalar_int_field[lnp];
  int tlnp[size];
  char fname[64];
  struct ioheader head;
  int dummy;
  MPI_File fh;
  MPI_Info info;
  MPI_Offset loffset,goffset;

  // File name
  sprintf(fname,"results/step%04d",state);

  dummy=sizeof(struct ioheader);

  // Open file 
  MPI_File_open(MPI_COMM_WORLD,fname,MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&fh);

  // Gather information on number of particles
  MPI_Allgather(&lnp,1,MPI_INT,tlnp,1,MPI_INT,MPI_COMM_WORLD);

  // Prepare file header
  head.npart[0]=np; 
  head.nall[0]=np;
  for(i=1;i<6;i++) { head.npart[i]=0; head.nall[i]=0;}
  for(i=0;i<6;i++) { head.massarr[i]=0.0; }

  // Write file header (main process)
  if(rank==0) {
    MPI_File_write(fh,&dummy,1,MPI_INT,MPI_STATUS_IGNORE);
    MPI_File_write(fh,&head,sizeof(struct ioheader),MPI_BYTE,MPI_STATUS_IGNORE);
    MPI_File_write(fh,&dummy,1,MPI_INT,MPI_STATUS_IGNORE);
  }

  goffset+=2*sizeof(int)+sizeof(struct ioheader);

  // Write position
  goffset+=sizeof(int); 
  loffset=goffset; 
  for(i=0;i<rank;i++) loffset=loffset+tlnp[i]*sizeof(struct vector3d);
  MPI_File_seek(fh,loffset,MPI_SEEK_SET);
  for(i=0;i<lnp;i++) vector_field[i]=particles[i].position;
  MPI_File_write(fh,vector_field,3*lnp,MPI_FLOAT,MPI_STATUS_IGNORE);
  goffset+=np*sizeof(struct vector3d);
  goffset+=sizeof(int);

  // Write IDs
  goffset+=sizeof(int);
  loffset=goffset;
  for(i=0;i<rank;i++) loffset=loffset+tlnp[i]*sizeof(int);
  MPI_File_seek(fh,loffset,MPI_SEEK_SET);
  for(i=0;i<lnp;i++) scalar_int_field[i]=(int)particles[i].gid;
  MPI_File_write(fh,scalar_int_field,lnp,MPI_INT,MPI_STATUS_IGNORE);
  goffset+=np*sizeof(int);
  goffset+=sizeof(int);

  // Write ranks
  goffset+=sizeof(int);
  loffset=goffset;
  for(i=0;i<rank;i++) loffset=loffset+tlnp[i]*sizeof(float);
  MPI_File_seek(fh,loffset,MPI_SEEK_SET);
  for(i=0;i<lnp;i++) scalar_float_field[i]=(float)rank; 
  MPI_File_write(fh,scalar_float_field,lnp,MPI_FLOAT,MPI_STATUS_IGNORE);
  goffset+=np*sizeof(float);
  goffset+=sizeof(int);

  // Skip
  goffset+=np*sizeof(float)+2*sizeof(int);

  // Write density
  goffset+=sizeof(int);
  loffset=goffset;
  for(i=0;i<rank;i++) loffset=loffset+tlnp[i]*sizeof(float);
  MPI_File_seek(fh,loffset,MPI_SEEK_SET);
  for(i=0;i<lnp;i++) scalar_float_field[i]=(float)rank;
  MPI_File_write(fh,scalar_float_field,lnp,MPI_FLOAT,MPI_STATUS_IGNORE);
  goffset+=np*sizeof(float);
  goffset+=sizeof(int);


  // Write file footer
  if(rank==0) {
    MPI_File_seek(fh,goffset,MPI_SEEK_SET);
    MPI_File_write(fh,&dummy,1,MPI_INT,MPI_STATUS_IGNORE);
  }

  // Close file
  MPI_File_close(&fh);

}