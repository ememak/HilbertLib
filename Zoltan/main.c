#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<mpi.h>
#include<particles.h>

int main(int argc,char **argv) {

  int i, k;
  int iter, suf;
  int seed, tests;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if(argc<5) {
    if(rank==0) {
      fprintf(stderr,"Parameters are missing.\n");
    }
    exit(1); 
  } else {
    np = atoi(argv[1]);
    seed = atoi(argv[2])+rank;
    niter = atoi(argv[3]);
    suf = atoi(argv[4]);
  }
    FILE *File;
    char *arr = calloc (100, sizeof (char));;
    sprintf (arr, "./results/MainOutput%d", suf);
    File = fopen (arr, "w");
    free (arr);

    double time=0.0;
    for(iter = 0; iter<niter; iter++){
    double begin=MPI_Wtime();
    generateParticles(np, seed);

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0) fprintf(File, "End of points generation, elapsed time: %lf\n", MPI_Wtime()-begin);
    double premid, postmid;
    /*for (i = 0; i < np/size; i++)
      {
	    fprintf (File, "Punkt #%d : ", i);
	    fprintf (File, "%lf %lf %lf\n",
	         particles[i].position.x, particles[i].position.y, particles[i].position.z);
      }
    fflush(File);*/
    decompositionInit(argc,argv);

	MPI_Barrier(MPI_COMM_WORLD);
	premid=MPI_Wtime();

	decompose();

	MPI_Barrier(MPI_COMM_WORLD);
	postmid=MPI_Wtime();
	if(rank==0)
	fprintf(File, "Iteration %d end, elapsed time: %lf\n", iter, postmid-premid);
	time+=postmid-premid;
	free(particles);
    }
  /*fprintf (File, "NewData:\n");
    for (i = 0; i < np/size; i++)
      {
	    fprintf (File, "rank:#%d point:#%d @@@   ", rank, i);
	    fprintf (File, "%lf %lf %lf\n", particles[i].position.x, particles[i].position.y, particles[i].position.z);
      }*/

    /*double befquery=MPI_Wtime();
    int procs[size];
    int numprocs;
    query(0.0, 0.0, 0.0, 0.5, 0.5, 0.5, procs, &numprocs);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
	printf("Query Time: %lf\n", MPI_Wtime()-befquery);
	printf("procsRes:\n");
	for(i=0; i<numprocs; i++)
	printf("%d ", procs[i]);
	printf("\n");
    }*/

    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
    fprintf(File, "AVG Time: %lf\n", time/niter);

    fclose(File);
    MPI_Finalize();

    return 0;
}
