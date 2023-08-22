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

int
main(int argc, char **argv)
{
    
    int mpiret;
    sc_MPI_Comm comm;
    t8_cmesh_t cmesh;
    t8_forest_t forest;
    /*Prefix for our output files.*/
    const char *prefix_uniform = "t8_my_first_uniform_forest";
    const char *prefix_adapt = "t8_my_first_adapted_forest";
    const char *prefix_decompressed = "t8_my_first_decompressed_forest";
    const int level = 5;                            /*Besser 6 nicht ueberschreiten, sonst wird das Ganze ganz schoen riesig*/
    //t8_locidx_t local_num_trees;                  /*Frage: Warum war mein Hypercube nicht 3D? Antwort: Liegt an  t8_eclass Element. Neue Frage: Welche Klasse muss ich denn dann nehmen? Antwort: T8_CLASS_HEX!*/
    //t8_gloidx_t global_num_trees;
    //t8_locidx_t local_num_elements;
    //t8_gloidx_t global_num_elements;

    //std::cout << "Hello!" << std::endl;

    // 1) initialize t8code

    /*Order: MPI, libsc, t8code*/

    /*Initialize MPI, before we can initialize sc and t8code. */
    mpiret = sc_MPI_Init (&argc, &argv);
    /* Error check the MPI return value. */
    SC_CHECK_MPI (mpiret);

    /* Initialize the sc library, before we can initialize t8code too. */
    sc_init(sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
    /*Initialize t8code with log level SC_LP_PRODUCTION. */
    t8_init (SC_LP_PRODUCTION);

    /*Print message*/
    t8_global_productionf (" \n");
    t8_global_productionf (" Hello, this is my first t8code example :) \n");
    t8_global_productionf (" In this example i build my first uniform forest, coarsen it and output it to vtu files. \n");
    t8_global_productionf (" \n");

    /* We use MPI_COMM_WORLD as a communicator. */
    comm = sc_MPI_COMM_WORLD;

    /* Create cmesh from first step. */
    cmesh = t8_build_cube_coarse_mesh (comm);

    /* Build the uniform forest, it is automatically partitioned among the processes and print informations. */
    forest = t8_build_uniform_forest (comm, cmesh, level);
    t8_global_productionf (" Created uniform forest. \n");
    t8_global_productionf (" Refinement level: \t%i\n", level);
    t8_print_forest_information (forest);    
    
    /* Write forest to vtu files. */
    t8_write_forest_vtk (forest, prefix_uniform);
    t8_global_productionf (" Wrote forest to vtu files: %s*\n", prefix_uniform);


    /*Adapt the forest. Reuse the forest variable, but notice, that the adapted forest is a new forest*/
    forest = t8_adapt_forest (forest);


    /*Print information of our new forest*/
    t8_global_productionf (" Adapted forest. \n");
    t8_print_forest_information (forest);

    /* Write forest to vtu files. */
    t8_write_forest_vtk (forest, prefix_adapt);
    t8_global_productionf (" Wrote forest to vtu files: %s*\n", prefix_adapt);

    /*Decompress the forest. Reuse the forest variable, but notice, that the decompressed forest is a new forest*/
    forest = t8_decompressed_forest (forest, level);

    /*Print information of our new forest*/
    t8_global_productionf (" Decompressed forest. \n");
    t8_print_forest_information (forest);

    /* Write forest to vtu files. */
    t8_write_forest_vtk (forest, prefix_decompressed);
    t8_global_productionf (" Wrote forest to vtu files: %s*\n", prefix_decompressed);
    

    /* Destroy the forest. */
    t8_destroy_forest (forest);
    t8_global_productionf (" Destroyed forest. \n");
    // 2) baue cmesh (siehe t8_cmesh.h und cmesh/t8_cmesh_examples.h)
    
    /*Build the actual coarse mesh. */
    //cmesh = t8_build_coarse_cmesh (sc_MPI_COMM_WORLD);
    /* Compute local and global number of trees. */
    //local_num_trees = t8_cmesh_get_num_local_trees (cmesh);
    //global_num_trees = t8_cmesh_get_num_trees (cmesh);
    //t8_global_productionf (" Created coarse mesh, \n");
    //t8_global_productionf (" Local number of trees:\t%i\n", 
    //                        local_num_trees);
    //t8_global_productionf (" Global number of trees:\t%li\n",
    //                        global_num_trees);
    //t8_write_cmesh_vtk (cmesh, prefix);
    //t8_global_productionf (" Wrote coarse mesh to vtu files: %s*\n",
    //                        prefix);
    //t8_destroyCmesh (cmesh);
    //t8_global_productionf (" Destroyed coarse mesh. \n");

    // 3) baue einen forest basierend auf dem cmesh (siehe forest/t8_forest.h)

    // 4) verfeinere den Forest uniform auf ein bestimmtes initiales "Refinement-Level"

    // 5) vergröbere den Forest mehrmals/iterativ adaptiv (mit einem beliebigen Adapt-Kriterium, evtl. siehe Tutorial 3)

    // 6) verfeinere den Forest (iterativ) auf das initiale level

    // 7) finalize t8code

    sc_finalize ();
    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);

    return 0;
}