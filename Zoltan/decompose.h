#include <zoltan.h>

/* Zoltan arrays and variables */
struct Zoltan_Struct *ztn;
int changes; 				// 1 if partitioning was changed, 0 otherwise
int numGidEntries; 			// Number of integers used for a global ID 
int numLidEntries; 			// Number of integers used for a local ID 
int numImport; 				// Number of objects to be sent to me 
ZOLTAN_ID_PTR importGlobalGids; 	// Global IDs of objects to be sent to me 
ZOLTAN_ID_PTR importLocalGids; 		// Local IDs of objects to be sent to me 
int *importProcs; 			// Process rank for source of each incoming object 
int *importToPart; 			// New partition for each incoming object 
int numExport; 				// Number of objects I must send to other processes
ZOLTAN_ID_PTR exportGlobalGids; 	// Global IDs of the objects I must send 
ZOLTAN_ID_PTR exportLocalGids; 		// Local IDs of the objects I must send 
int *exportProcs; 			// Process to which I send each of the objects 
int *exportToPart; 			// Partition to which each object will belong 

/* Zoltan query functions */
int ztn_return_dimension(void *data, int *ierr);
void ztn_return_coords(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr);
int ztn_return_num_node(void *data, int *ierr);
void ztn_return_owned_nodes(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr);
int ztn_return_particle_data_size(void *data,int num_gid_entries,int num_lid_entries,ZOLTAN_ID_PTR global_id,ZOLTAN_ID_PTR local_id,int *ierr);
void ztn_pack(void *data,int num_gid_entries,int num_lid_entries,ZOLTAN_ID_PTR global_id,ZOLTAN_ID_PTR local_id,int dest,int size,char *buf,int *ierr);
void ztn_unpack(void *data,int num_gid_entries,ZOLTAN_ID_PTR global_id,int size,char *buf,int *ierr);
void ztn_pre(void *data, int num_gid_entries, int num_lid_entries, int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr);
void ztn_mid(void *data, int num_gid_entries, int num_lid_entries, int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr);
void ztn_post(void *data, int num_gid_entries, int num_lid_entries, int num_import, ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids, int *import_procs, int *import_to_part, int num_export, ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids, int *export_procs, int *export_to_part, int *ierr);
