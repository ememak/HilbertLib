#include <stdio.h>
#include <stdlib.h>
#include <particles.h>
#include <mpi.h>
#include <zoltan.h>
#include <decompose.h>


/* This function performs the initialization of Zoltan decomposition */
int decompositionInit(int argc,char **argv) {
 
  float version;

  Zoltan_Initialize(argc,argv,&version);
  if(rank==0) printf("Zoltan Version %.3f. Initialized.\n",version);

  ztn=Zoltan_Create(MPI_COMM_WORLD);

  /* Hilbert Space-Filling Curve Partitioning */
  Zoltan_Set_Param(ztn,"LB_METHOD","HSFC"); 
  /* Global ID is 1 integer */
  Zoltan_Set_Param(ztn, "NUM_GID_ENTRIES", "1"); 
  /* Local ID is 1 integer */
  Zoltan_Set_Param(ztn, "NUM_LID_ENTRIES", "1");
  /* Enable object weights */
  Zoltan_Set_Param(ztn, "OBJ_WEIGHT_DIM", "1");
  /* Quiet mode; no output unless an error or warning is produced */
  Zoltan_Set_Param(ztn, "DEBUG_LEVEL","0");
  /* Information about cuts and bounding box */
  Zoltan_Set_Param(ztn, "KEEP_CUTS","1");
  /* Automatic migration turned on */
  Zoltan_Set_Param(ztn, "AUTO_MIGRATE", "0");

  /* Returns the number of values needed to express the geometry of an object */
  Zoltan_Set_Fn(ztn, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) ztn_return_dimension, particles); 	
  /* Returns a vector of geometry values for a given object */
  Zoltan_Set_Fn(ztn, ZOLTAN_GEOM_FN_TYPE, (void (*)()) ztn_return_coords, particles); 
  /* Returns the number of objects that are currently assigned to the processor */  
  Zoltan_Set_Fn(ztn, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) ztn_return_num_node, particles);
  /* Fills arrays with information about the objects currently assigned to the processor */
  Zoltan_Set_Fn(ztn, ZOLTAN_OBJ_LIST_FN_TYPE, (void (*)()) ztn_return_owned_nodes, particles);
  /* Returns the size (in bytes) of the data buffer that is needed to pack all of a single object's data */
  Zoltan_Set_Fn(ztn, ZOLTAN_OBJ_SIZE_FN_TYPE,(void (*)()) ztn_return_particle_data_size, particles);

  /* Allows the application to specify how to copy all needed data for a given object into a communication buffer */
  Zoltan_Set_Fn(ztn, ZOLTAN_PACK_OBJ_FN_TYPE,(void (*)()) ztn_pack, particles);  			
  /* Allows the application to specify how to copy all needed data for a given object from a communication buffer */
  Zoltan_Set_Fn(ztn, ZOLTAN_UNPACK_OBJ_FN_TYPE,(void (*)()) ztn_unpack, particles);  
  /* Performs any pre-processing desired by the application */
  Zoltan_Set_Fn(ztn, ZOLTAN_PRE_MIGRATE_PP_FN_TYPE,(void (*)()) ztn_pre, particles);
  /* Performs any processing desired by the application between the packing and unpacking of objects being migrated */
  Zoltan_Set_Fn(ztn, ZOLTAN_MID_MIGRATE_PP_FN_TYPE,(void (*)()) ztn_mid, particles);
  /* Performs any post-processing desired by the application */
  Zoltan_Set_Fn(ztn, ZOLTAN_POST_MIGRATE_PP_FN_TYPE,(void (*)()) ztn_post, particles);

  return 0;

}

/* This function performes the decomposition */
int decompose() {

  int rc;
  int i;

  rc = Zoltan_LB_Partition(ztn, /* input (all remaining fields are output) */
        &changes,        	/* 1 if partitioning was changed, 0 otherwise */
        &numGidEntries,  	/* Number of integers used for a global ID */
        &numLidEntries,  	/* Number of integers used for a local ID */
        &numImport,      	/* Number of objects to be sent to me */
        &importGlobalGids,  	/* Global IDs of objects to be sent to me */
        &importLocalGids,   	/* Local IDs of objects to be sent to me */
        &importProcs,    	/* Process rank for source of each incoming object */
        &importToPart,   	/* New partition for each incoming object */
        &numExport,      	/* Number of objects I must send to other processes*/
        &exportGlobalGids,  	/* Global IDs of the objects I must send */
        &exportLocalGids,   	/* Local IDs of the objects I must send */
        &exportProcs,    	/* Process to which I send each of the objects */
        &exportToPart);  	/* Partition to which each object will belong */

  if (rc != ZOLTAN_OK)
    printf("Error in Zoltan library\n");

  // Free the arrays allocated by Zoltan_LB_Partiotion
  Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids,&importProcs, &importToPart);
  Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids,&exportProcs, &exportToPart);

  return 0;

}

int query(double xmin, double ymin, double zmin, double xmax, double ymax, double zmax, int * procs, int * numprocs){
    Zoltan_LB_Box_Assign(ztn, xmin, ymin, zmin, xmax, ymax, zmax, procs, numprocs);
}

int ztn_return_dimension(void *data, int *ierr) {
  return 3;
}

void ztn_return_coords(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr) {
  struct particle_data *p = (struct particle_data*) data;
  geom_vec[0]=p[(int)(*local_id)].position.x;
  geom_vec[1]=p[(int)(*local_id)].position.y;
  geom_vec[2]=p[(int)(*local_id)].position.z;
}

int ztn_return_num_node(void *data, int *ierr) {
  return lnp;
}

void ztn_return_owned_nodes(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr) {
   int i;
   struct particle_data *p = (struct particle_data*) data;
   for(i=0;i<lnp;i++) {
     global_ids[i*num_gid_entries]=p[i].gid;
     local_ids[i*num_lid_entries]=i;
     obj_wgts[i]=1.0;
   }
}

int ztn_return_particle_data_size(void *data,int num_gid_entries,int num_lid_entries,ZOLTAN_ID_PTR global_id,ZOLTAN_ID_PTR local_id,int *ierr) {
   return sizeof(struct particle_data);
}

void ztn_pack(void *data,int num_gid_entries,int num_lid_entries,ZOLTAN_ID_PTR global_id,ZOLTAN_ID_PTR local_id,int dest,int size,char *buf,int *ierr) {
   struct particle_data *p = (struct particle_data*) data;
   memcpy(buf,&(p[(int)(*local_id)]),sizeof(struct particle_data));
   p[(int)(*local_id)].gid=-1; // Mark local particle as exported
}

void ztn_pre(void *data, int num_gid_entries, int num_lid_entries, int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr) {
  // Any pre communication operations should go here
  // Example: print decomposition statistics
  //printf("Process: %d Exports: %d Imports: %d\n",rank,num_export,num_import);
}

void ztn_mid(void *data, int num_gid_entries, int num_lid_entries, int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr) {
  int pos,i;
  struct particle_data *p = (struct particle_data*) data;
  pos=0;
  for(i=0;i<lnp;i++) {
    if(i!=pos && p[i].gid!=-1) { p[pos]=p[i]; p[pos]=p[i]; }
    if(p[i].gid!=-1) pos++;
  }
  lnp=lnp-num_export;
}

void ztn_post(void *data, int num_gid_entries, int num_lid_entries, int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr) {
  // Any post communication operations should go here
}

void ztn_unpack(void *data,int num_gid_entries,ZOLTAN_ID_PTR global_id,int size,char *buf,int *ierr) {
   struct particle_data *p = (struct particle_data*) data;
   memcpy(&p[lnp],buf,sizeof(struct particle_data));
   lnp++;
}

