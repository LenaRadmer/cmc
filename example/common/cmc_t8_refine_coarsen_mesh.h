#ifndef CMC_T8_REFINE_COARSEN_MESH_H /*wICHTIG!*/
#define CMC_T8_REFINE_COARSEN_MESH_H

#include "t8_cmesh.h"                           /*cmesh definition and basic interface. */
#include <t8_forest/t8_forest_general.h>        /* forest definition and general interface. */

t8_cmesh_t
t8_build_cube_coarse_mesh (sc_MPI_Comm comm);

t8_forest_t
t8_build_uniform_forest (sc_MPI_Comm comm, t8_cmesh_t cmesh, int level);

void t8_write_forest_vtk (t8_forest_t forest, const char *prefix);

void t8_destroy_forest (t8_forest_t forest);

void t8_print_forest_information (t8_forest_t forest);

t8_forest_t
t8_adapt_forest (t8_forest_t forest);

t8_forest_t
t8_decompressed_forest (t8_forest_t forest, const int level);


#endif              /* !CMC_T8_REFINE_COARSEN_MESH_H*/