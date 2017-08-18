#include"ToGrid.h"

void toGrid(double* in, coord_t* out, double* MIN, double* MAX, int Dimensions, int BitsPrecision){//Points are in space <min[i], max[i])
    int i;
    double len=(double)(1<<BitsPrecision);
    for(i=0; i<Dimensions; i++){
	out[i]=(int)(((in[i]-MIN[i])/MAX[i])*len);
    }
    return;
}

void toPos(coord_t* in, double* out, double* MIN, double* MAX, int Dimensions, int BitsPrecision){//C.A. reverse of toGrid
    int i;
    double len=1.0/(double)(1<<BitsPrecision);
    for(i=0; i<Dimensions; i++){
	out[i]=((double)in[i])*len*MAX[i]+MIN[i];
    }
    return;
}

void getMinMax(double* data, int DataCount, double** mind, double** maxd, int Dimensions){
    int i, j;
    (*mind) = malloc(Dimensions * sizeof(double));
    (*maxd) = malloc(Dimensions * sizeof(double));
    double *locmind=malloc(Dimensions * sizeof(double));
    double *locmaxd=malloc(Dimensions * sizeof(double));
    if(DataCount>0){
	for(i=0; i<Dimensions; i++){
	    locmind[i]=data[i];
	    locmaxd[i]=data[i];
	}
    }
    else{
	for(i=0; i<Dimensions; i++){
	    locmind[i]=1.0e300;
	    locmaxd[i]=-1.0e300;
	}
    }
    for(i=1; i<DataCount; i++){
        for(j=0; j<Dimensions; j++){
	    if(data[Dimensions * i + j] < locmind[j])
		locmind[j]=data[Dimensions * i + j];
	    if(data[Dimensions * i + j]>locmaxd[j])
		locmaxd[j]=data[Dimensions * i + j];
        }
    }
    MPI_Allreduce(locmind, (*mind), Dimensions, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(locmaxd, (*maxd), Dimensions, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    free(locmind);
    free(locmaxd);
}

//use to change double coords to ints, mind and maxd should be results of getMinMax or NULL
void allToGrid(double* in, double* mind, double* maxd, coord_t** out, int PointsCount, int Dimensions, int BitsPrecision){
    int i;
    (*out) = malloc(PointsCount * sizeof(coord_t) * (Dimensions+1));
    if(mind==NULL || maxd==NULL)
    getMinMax(in, PointsCount, &mind, &maxd, Dimensions);
    for(i=0; i<PointsCount; i++){
		(*out)[(Dimensions+1) * i]=i;
		toGrid(&(in[i * Dimensions]), &((*out)[(Dimensions+1)*i + 1]), mind, maxd, Dimensions, BitsPrecision);
    }
    return;
}
