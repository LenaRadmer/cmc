#include <../example/common/cmc_t8_refine_coarsen_mesh.h> 
#include <iostream>
#include <t8.h>                                 /*General t8code header. */
#include <t8_cmesh_vtk_writer.h>                /*cmesh-writer interface. */
#include "t8_cmesh/t8_cmesh_examples.h"         /*collection of exemplary cmeshes. */
#include <t8_forest/t8_forest_io.h>             /* forest io interface. */
#include <t8_schemes/t8_default/t8_default_cxx.hxx> /* default refinement scheme. */
#include <t8_forest/t8_forest_geometrical.h>    /*geometrical information of the forest*/
#include <t8_vec.h>                             /*Basic operations on 3D vectors*/
#include <t8_element_c_interface.h>
#include <t8_element.h>

struct my_adapt_data {
    double midpoint[3];
    double coarsen_outside_radius;
};

struct my_refinement_data {
    int initial_uniform_level;
};

#if 0
// 2) baue cmesh (siehe t8_cmesh.h und cmesh/t8_cmesh_examples.h)
/*static t8_cmesh_t
t8_build_coarse_cmesh (sc_MPI_Comm comm) {
    t8_cmesh_t cmesh;*/
    /*build a unit square with 1 quadrialateral tree. Second argument is the MPI communicator to use for this cmesh. */
 /*   cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);

    return cmesh;
}*/
/*Write vtk file of the cmesh. */
/*static void t8_write_cmesh_vtk (t8_cmesh_t cmesh, const char *prefix) {
    t8_cmesh_vtk_write_file (cmesh, prefix, 1.0);
}
*/
/*Destroy a cmesh (free all allocated memory)*/
/*static void t8_destroyCmesh (t8_cmesh_t cmesh) {
    t8_cmesh_destroy (&cmesh);
}*/
#endif

// 3) baue einen forest basierend auf dem cmesh (siehe forest/t8_forest.h)
// 4) verfeinere den Forest uniform auf ein bestimmtes initiales "Refinement-Level"
t8_cmesh_t
t8_build_cube_coarse_mesh (sc_MPI_Comm comm) {
    t8_cmesh_t cmesh;

    /* Build a coarse mesh. */
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
    t8_global_productionf (" [step0] Constructed coarse mesh. \n");

    return cmesh;
}

void t8_write_forest_vtk (t8_forest_t forest, const char *prefix) {
    t8_forest_write_vtk (forest, prefix);
}

void t8_destroy_forest (t8_forest_t forest) {
    t8_forest_unref (&forest);
}

t8_forest_t 
t8_build_uniform_forest (sc_MPI_Comm comm, t8_cmesh_t cmesh, int level) {
    t8_forest_t forest;
    t8_scheme_cxx_t *scheme;

    /* Create the refinement scheme. */
    scheme = t8_scheme_new_default_cxx ();
    /* Create the uniform forest. */
    forest = t8_forest_new_uniform (cmesh, scheme, level, 0, comm);

    return forest;
}

// 5) vergröbere den Forest mehrmals/iterativ adaptiv (mit einem beliebigen Adapt-Kriterium, evtl. siehe Tutorial 3)

static int t8_adapt_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id, 
                        t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[]) {
     /*Hold element for adapt criterion. */
    double centroid[3];    
    /*User pointer to forest*/                 
    const struct my_adapt_data *adapt_data = (const struct my_adapt_data *) t8_forest_get_user_data (forest);    
    double dist;                    /*Hold distance from centroid and element´s midpoint. */

    /* Want to make sure that we actually did set a user pointerto forest and didn´t get a NULL pointer. */
    T8_ASSERT (adapt_data != NULL);
    
    /* Compute the element's centroid coordinates*/
    t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

    /* Compute the distance from centroid and element´s midpoint*/
    dist = t8_vec_dist (centroid, adapt_data->midpoint);
    
    /* Coarsening criterion*/
    if (is_family && dist > adapt_data->coarsen_outside_radius) {
        /* Coarsen this family. First we need to check if we have a family that we can coarsen*/
        return -1;
    }
    /* Do not refine or coarsen this element*/
    return 0;
}
/* Adapt a forest according to our t8_adapt_callback function. We will create a new forest and return it*/
t8_forest_t
t8_adapt_forest (t8_forest_t forest) {
    t8_forest_t forest_adapt;
    /*Variable for holding number of global elements of the previous forest for our while criterion*/
    t8_gloidx_t num_elements_former_forest;
    num_elements_former_forest = 0;

    struct my_adapt_data adapt_data;
    adapt_data.midpoint[0] = 0.3;
    adapt_data.midpoint[1] = 0.3;
    adapt_data.midpoint[2] = 0.0;
    adapt_data.coarsen_outside_radius = 0.3;

    /* Apply the coarsening as often as possible*/
    while (num_elements_former_forest != t8_forest_get_global_num_elements (forest)) {      /*Eventuell Level !=0 hinzufügen!*/
        /*Update num_elements_former_forest*/
        num_elements_former_forest = t8_forest_get_global_num_elements (forest);

        forest_adapt = t8_forest_new_adapt (forest, t8_adapt_callback, 0, 0, &adapt_data);
        /*Update the forest*/
        forest = forest_adapt;
    }


    return forest;
}

/*Print the local and global number of elements of a forest*/
void t8_print_forest_information (t8_forest_t forest) {
    t8_locidx_t local_num_elements;
    t8_gloidx_t global_num_elements;

    /*Check that forest is a committed, that is valid and usable, forest*/
    T8_ASSERT (t8_forest_is_committed (forest));

    /* Get the local number of elements. */
    local_num_elements = t8_forest_get_local_num_elements (forest);
    /* Get the global number of elements. */
    global_num_elements = t8_forest_get_global_num_elements (forest);

    t8_global_productionf (" Local number of elements:\t\t%i\n", local_num_elements);
    t8_global_productionf (" Global number of elements:\t%li\n", global_num_elements);

}
// 6) verfeinere den Forest (iterativ) auf das initiale level

static int t8_decompression_callback(t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id, 
                                t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[]) {
    int elem_level;
    /*Compute element id level*/
    elem_level = t8_element_level(ts, elements[0]);
    /*User pointer to forest*/
    const struct my_refinement_data *refinement_data = (const struct my_refinement_data *) t8_forest_get_user_data (forest);
    /*Refinement criterion*/
    if (elem_level < refinement_data->initial_uniform_level) {
        /*Refine element if the level is lower than the initial uniform level*/
        return 1;
    }
    return 0;       /*Level is the same as initial uniform level*/

}
/*Decompression of adapted forest back to uniform forest*/
t8_forest_t
t8_decompressed_forest (t8_forest_t forest, const int level) {
    /*New decrompressed forest*/
    t8_forest_t forest_decompressed;
    /*Variable for holding number of global elements of the previous forest for our while criterion*/
    t8_gloidx_t num_elements_former_forest;
    num_elements_former_forest = 0;

    struct my_refinement_data refinement_data;
    refinement_data.initial_uniform_level = level;

    /*Apply the refinement as often as possible*/
    while (num_elements_former_forest != t8_forest_get_global_num_elements (forest)) {      /*Eventuell Level !=0 hinzufügen!*/
        /*Update num_elements_former_forest*/
        num_elements_former_forest = t8_forest_get_global_num_elements (forest);

        forest_decompressed = t8_forest_new_adapt (forest, t8_decompression_callback, 0, 0, &refinement_data);
        /*Update the forest*/
        forest = forest_decompressed;
    }
    return forest;
}