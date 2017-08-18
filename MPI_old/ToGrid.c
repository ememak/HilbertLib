#include"ToGrid.h"

void toGrid(double* in, coord_t* out, double* MIN, double* MAX, int Dimensions, int BitsPrecision){//Points are in space <min[i], max[i])
    double len=(double)(1<<BitsPrecision);
    for(int i=0; i<Dimensions; i++){
	out[i]=(int)(((in[i]-MIN[i])/MAX[i])*len);
    }
    return;
}

void toPos(coord_t* in, double* out, double* MIN, double* MAX, int Dimensions, int BitsPrecision){//C.A. reverse of toGrid
    double len=1.0/(double)(1<<BitsPrecision);
    for(int i=0; i<Dimensions; i++){
	out[i]=((double)in[i])*len*MAX[i]+MIN[i];
    }
    return;
}

void getMinMax(double** data, int DataCount, double** mind, double** maxd, int Dimensions){
    (*mind) = calloc(Dimensions, sizeof(double));
    (*maxd) = calloc(Dimensions, sizeof(double));
    double *locmind=calloc(Dimensions, sizeof(double));
    double *locmaxd=calloc(Dimensions, sizeof(double));
    if(DataCount>0){
	for(int i=0; i<Dimensions; i++){
	    locmind[i]=data[0][i];
	    locmaxd[i]=data[0][i];
	}
    }
    else{
	for(int i=0; i<Dimensions; i++){
	    locmind[i]=1.0e300;
	    locmaxd[i]=-1.0e300;
	}
    }
    for(int i=1; i<DataCount; i++){
        for(int j=0; j<Dimensions; j++){
	    if(data[i][j]<locmind[j])
		locmind[j]=data[i][j];
	    if(data[i][j]>locmaxd[j])
		locmaxd[j]=data[i][j];
        }
    }
    MPI_Allreduce(locmind, (*mind), Dimensions, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(locmaxd, (*maxd), Dimensions, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
}

//use to change double coords to ints, mind and maxd should be results of getMinMax or NULL
void allToGrid(double** in, tag_t* idx, double* mind, double* maxd, int PointsCount, MDPoint** out, int Dimensions, int BitsPrecision){
    (*out) = calloc(PointsCount, sizeof(MDPoint));
    if(mind==NULL || maxd==NULL)
    getMinMax(in, PointsCount, &mind, &maxd, Dimensions);
    for(int i=0; i<PointsCount; i++){
	make_MDPoint(&(*out)[i], Dimensions);
	(*out)[i].own_data_id=idx[i];
	toGrid(in[i], (*out)[i].coordinates, mind, maxd, Dimensions, BitsPrecision);
    }
    return;
}
