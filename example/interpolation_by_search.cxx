#include "cmc.h"
#include "utilities/cmc_util.h"
#include "utilities/cmc_log_functions.h"
#include "t8code/cmc_t8code.h"
#include "t8code/cmc_t8code_data.h"
#include "t8_vec.h" 
#include "t8_geometry/t8_geometry_base.h"
#include <chrono>
#include <numeric>

/************************************************************************/
/**** General functions (without featuring coarsening or refinement) ****/
/**
 * @brief Write a forest out to .vtu files
 * 
 * @param forest The forest to write out
 * @param file_prefix The name the .vtu files will have
 */
void
cmc_t8_forest_write_vtk(t8_forest_t forest, const char* file_prefix)
{
    /* Write out the forest */
    int vtk_err = t8_forest_vtk_write_file (forest, file_prefix, 1, 1, 1, 1, 0, 0, NULL);
    /* Check if an error occured during the write operation */
    if (vtk_err == 0) { cmc_err_msg("An error occrued during the creation of the t8code-forest vtk-file.");}
}

/**
 * @brief A function overload for writing out .vtu files of the forest. This function writes the forest as well as a single
 * data variable which is defined on the forest.
 * 
 * @param forest The forest to write out
 * @param file_prefix The name of the .vtu files
 * @param data_variable The data variable which will be written to the .vtu files as well (#data_points == #elements_in_the_forest)
 */
void
cmc_t8_forest_write_vtk(t8_forest_t forest, const char* file_prefix, std::vector<double>& data_variable)
{
    /* Check if the vector contains exactly the same amount of data points as elements are contained in the forest */
    if(data_variable.size() != static_cast<size_t>(t8_forest_get_local_num_elements(forest))){ cmc_warn_msg("There are not exactly as many values in the vector as there are element in the forest (this may be caused by not using the push_back() function of the vector)");}
    
    /* For each user defined data field we need one t8_vtk_data_field_t variable */
    t8_vtk_data_field_t vtk_data;

    /* Set the type of this variable. Since we have one value per element, we pick T8_VTK_SCALAR */
    vtk_data.type = T8_VTK_SCALAR;

    /* The name of the field as should be written to the file. */
    strcpy (vtk_data.description, "Example_Data_Var");
    vtk_data.data = data_variable.data();

    /* Write out the forest as well as the data variable which is defined on the forest */
    t8_forest_write_vtk_ext (forest, file_prefix, 0, 0, 1, 1, 0, 0, 0, 1, &vtk_data);
}

/**
 * @brief Create some example data which can be associated with the forest. Each element will be associated with exactly 
 * one data point. The data points are ordered just the same as the forest elements, e.g. the element with index=12 (in the forest) will
 * associated with the data point contained in the 12th position of the element_data vector. Therefore, the element_data
 * vector will contain exactly "t8_forest_get_global_num_elements(forest)" data points.
 * 
 * @param forest The forest with whom the data array will be associated
 * @param element_data A vector which will be filled by this function contaning the associated data
 */
void
cmc_create_example_data(t8_forest_t forest, std::vector<double>& element_data)
{
    /* Allocate memory for the vector holding the element data */
    element_data.reserve(t8_forest_get_local_num_elements(forest));

    /* A vector for the midpoint coordinates of the elements */
    double elem_coords[3] = {0.0, 0.0, 0.0};

    /* Get the dimension of the forest */
    const int dimension_of_forest = t8_geom_get_dimension(t8_cmesh_get_tree_geometry(t8_forest_get_cmesh(forest), 0));

    /* Declare a variable which stores the pointer to the current element */
    t8_element_t* element;

    /* Decalare a scalar indicating the frequency of of the cosine functions */
    const double freq = 4.0;

    /* Iterate over all elements of the forest */
    for (int element_id = 0; element_id < t8_forest_get_local_num_elements(forest); ++element_id)
    {
        /* Retrieve the actual element corresponding to the element id */
        element = t8_forest_get_element_in_tree(forest, 0, element_id);

        /* Compute the element's centroid coordinates. */
        t8_forest_element_centroid (forest, 0, element, elem_coords);

        /* Based on the dimension of the forest calculate an example data point for this variable */
        if (dimension_of_forest == 2)
        {
            /** A two dimensional forest **/
            /* Save the computed data point for this element */
            element_data.push_back(std::cos(freq * elem_coords[0]) * std::cos(freq * elem_coords[1]));
        } else
        {
            /** A three dimensional forest **/
            /* Save the computed data point for this element */
            element_data.push_back(std::cos(freq * elem_coords[0]) * std::cos(freq * elem_coords[1]) * std::cos(freq * elem_coords[2]));
        }
    }

}

/**
 * @brief The function checks whether both vectors (decompressed by different methods) contain the same (decompressed) data
 * 
 * @param element_data_decompr1 Element_data vector from the first decompresison approach
 * @param element_data_decompr2 Element_data vector from the second decompression approach
 */
void
cmc_check_decompression_results_for_equality(std::vector<double>& element_data_decompr1, std::vector<double>& element_data_decompr2)
{
    /* Check if the vectors have the same length */
    (element_data_decompr1.size() != element_data_decompr2.size() ? cmc_warn_msg("The vectors do not have the same size") : cmc_debug_msg("The vectors which will be compared have the same size"));

    /* The maximum length of the vectors */
    const size_t vec_length = std::max(element_data_decompr1.size(), element_data_decompr2.size());

    /* Try to compare all entries within the vectors */
    try
    {
        /* Get the iterator to the first entry of the first element vector */
        auto iter1 = element_data_decompr1.begin();
        /* Get the iterator to the first entry of the second element vector */
        auto iter2 = element_data_decompr2.begin();

        /* Iterate over all data points */
        for (size_t iteration_count = 0; iteration_count < vec_length; ++iteration_count)
        {
            if (std::abs(*(iter1 + iteration_count) - *(iter2 + iteration_count)) > 8.0 * std::numeric_limits<double>::epsilon())
            {
                cmc_warn_msg("The data within the vectors is not equal, they start differing at index: ", iteration_count);
                break;
            }
            
        }
    }
    catch (const std::exception& e)
    {
        /* Print out the thrown exception */
        std::cout << e.what() << std::endl;
    }
}

/**
 * @brief This struct resembles the user's data which holds the variable data which will interpolated between (succeding adaptive) forests
 */
struct cmc_variable_interpolation_data
{
    std::vector<double> data;
    std::vector<double> data_new;
};

/* A typedef declaring a pointer to the struct */
typedef struct cmc_variable_interpolation_data* cmc_variable_interpolation_data_t;

/************************************************************************/


/************************************************************************/
/* The funcitons below are featuring the coarsening/compression of the forest and the data */
/**
 * @brief A 't8code_replace_callback_function' which performs the interpolation of data between different forests (only differing in at maximum one refinement level per element)
 * This function performs a standard mean calucation during the coarsening process.
 */
void
cmc_t8_interpolation_function_for_coarsening(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                             t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                             int num_incoming, t8_locidx_t first_incoming)
{
    /* Retrieve the user data which contains the variable's data vectors */
    cmc_variable_interpolation_data_t interpolation_data = static_cast<cmc_variable_interpolation_data_t>(t8_forest_get_user_data(forest_new));

    if (refine == -1)
    {
        /** If a family has been coarsened **/
        /* Calculate the standard mean of the outgoing elements' value */
        /* This is done by accumulating the values of all outgoing elements and afterwards a division by the number of outgoing elements */
        const double standard_mean = std::accumulate(interpolation_data->data.begin() + first_outgoing, interpolation_data->data.begin() + first_outgoing + num_outgoing, 0.0) / num_outgoing;
        /* Since the iteration of the elements is linear (and we are focusing on a serial execution environment), we can just push_back the values to the vector holding the interpolated data (instead of assigning the value at an specific index) */
        interpolation_data->data_new.push_back(standard_mean); //This works only in the current case (serial execution, linearily ordered, only one tree, etc.), a more general assignemnt would be: interpolation_data->data_new[first_incoming] = standard_mean; or even better use: interpolation_data->data_new.insert(interpolation_data->data_new.begin() + first_incoming, standard_mean);
    } else if (refine == 0)
    {
        /** The element stays the same **/
        /* When the element stays the same, we can just copy the value over */
        /* Since the iteration of the elements is linear (and we are focusing on a serial execution environment), we can just push_back the values to the vector holding the interpolated data (instead of assigning the value at an specific index) */
        interpolation_data->data_new.push_back(interpolation_data->data[first_outgoing]); //This works only in the current case (serial execution, linearily ordered, only one tree, etc.), a more general assignemnt would be: interpolation_data->data_new[first_incoming] = interpolation_data->data[first_outgoing]; or even better use: interpolation_data->data_new.insert(interpolation_data->data_new.begin() + first_incoming, interpolation_data->data[first_outgoing]);
    
    } else
    {
        /** The element has been refined **/
        cmc_err_msg("A refinement during the coarsening/compression is illegal.");
    }
}

/**
 * @brief This is just an exemplary adapt function which coarsens some elements. This function
 * is similar to the the adapt function from the step3 of the tutorial section
 */
int
cmc_t8_arbitrary_coarsen_callback (t8_forest_t forest,
                                   t8_forest_t forest_from,
                                   t8_locidx_t which_tree,
                                   t8_locidx_t lelement_id,
                                   t8_eclass_scheme_c *ts,
                                   const int is_family,
                                   const int num_elements, t8_element_t *elements[])
{
    /* Declare an double array (capable oif holding an element's midpoint) */
    double centroid[3] = {0.0, 0.0, 0.0};
    
    /* Compute the element's centroid coordinates. */
    t8_forest_element_centroid (forest_from, which_tree, elements[0], centroid);

    /* We define a refernece point from which the distance to the element's midpoint is measured */
    const double reference_point[] = {0.6, 0.6, 0.0};

    /* Compute the distance from the element's midpoint to our refernece point. */
    const double dist = t8_vec_dist (centroid, reference_point);
    
    /* If a family is passed to the callback and the family's first element's midpoint is in a certain area, we will coarsen this family */
    if (is_family && (dist <= 0.2 || (dist >= 0.4 && dist <= 0.65)))
    {
        /* Coarsen this family */
        return -1;
    }

    /* Otherwise, the element stays the same */
    return 0;
}

/**
 * @brief This function applies an arbitrary adaptation criterion which will coarsen the uniform forest.
 * Applying this coarsening step resembles the lossy compression step.
 * 
 * @param forest A forest which was uniformly refined
 * @return t8_forest_t The function returns the coarsened forest
 */
t8_forest_t
cmc_coarsen_the_initial_forest_arbitrarily(t8_forest_t forest, std::vector<double>& variable_data)
{
    /** Within this function, we coarsen the forest iteratively several times **/
    /* Declare a new forest variable */
    t8_forest_t forest_adapt;

    /* A variable holding the name of the file in which the forest will be written */
    std::string file_prefix = std::string("cmc_t8_forest_coarsening_step");

    /* Create a new struct for the interpolation data (holding the variable data) */
    cmc_variable_interpolation_data interpolation_data;

    /* Save the variable data within the interpolation data struct (by swapping the empty vector in the struct with the variable data which was passed to this function )*/
    std::swap(interpolation_data.data, variable_data); //The data of the variable is from now on only accessible via "interpolation_data.data" and no longer via variable_data
    
    /* Perform four coarsening steps */
    for (int coarsening_step = 1; coarsening_step <= 4; ++coarsening_step)
    {
        cmc_debug_msg("Coarsening step: ", coarsening_step);

        /* Initialize the new forest */
        t8_forest_init(&forest_adapt);

        /* We need to reference the 'old' forest, otherwise it will be deleted by the "set_adapt_call" on commit. This has to be done because, we need both - the old and the new forest - in order to perform the interpolation */
        t8_forest_ref(forest);

        /* Adapt the forest accordingly to the callback function */
        t8_forest_set_adapt(forest_adapt, forest, cmc_t8_arbitrary_coarsen_callback, 0);

        /* Commit the adapted forest (perform the adaption step) */
        t8_forest_commit(forest_adapt);

        /* Set the interpolation data accordingly */
        t8_forest_set_user_data(forest_adapt, static_cast<void*>(&interpolation_data));
        
        /* Allocate memory for the interpolated data (which will be transferred from 'forest' to 'forest_adapt') equal to the amount of elements of 'forest_adapt' */
        interpolation_data.data_new.reserve(t8_forest_get_local_num_elements(forest_adapt));

        /* Interpolate the element data onto the new coarsened forest */
        t8_forest_iterate_replace(forest_adapt, forest, cmc_t8_interpolation_function_for_coarsening);
    
        /* Swap the 'old' variable data which was defined on 'forest' with the 'new' variable data which was interpolated onto 'forest_adapt' */
        std::swap(interpolation_data.data, interpolation_data.data_new);

        /* Delete the 'old' variable data (which was defined on 'forest') */
        interpolation_data.data_new.clear();

        /* After the interpolation of the variable data, we do not need the 'old' forest anymore, therefore, it can be dereferenced/deallocated */
        t8_forest_unref(&forest);

        /* Save the coarsened forest (i.e. 'swap' the 'old' forest with the 'new' forest) */
        forest = forest_adapt;

        /* Write out the forest to .vtu files and give them a reasonable name */
        cmc_t8_forest_write_vtk(forest, (file_prefix + std::to_string(coarsening_step)).c_str(), interpolation_data.data);
    }

    /* We swap the vectors holding the variable' data back, such that the coarsened/compressed variable data now resides again in the vector 'variable_data' */
    /* Such that the compressed data can be accessed by the calling function */
    std::swap(interpolation_data.data, variable_data);

    /* Return the coarsened forest */
    return forest;
}
/************************************************************************/


/************************************************************************/
/* The functions below are featuring the refinement/decompression of the forest and the data */
/** But the decompression is applied iteratively (only one level finer per adaptation step) **/
/**
 * @brief This struct resembles the user's adapt data which is needed for the refinement back onto the initial level
 */
struct cmc_refine_adapt_data
{
    int initial_refinement_level{0};
    int coarsest_element_level_present_in_the_forest{0};
};

/* A typedef declaring a pointer to the struct */
typedef struct cmc_refine_adapt_data* cmc_refine_adapt_data_t;

/**
 * @brief A 't8code_replace_callback_function' which performs the interpolation of data between different forests (only differing in at maximum one refinement level per element)
 * This function performs a simple copying mechanism during the refinement process (i.e. copy the value from the coarse element to all refined child elements).
 */
void
cmc_t8_interpolation_function_for_refinement(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                             t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                             int num_incoming, t8_locidx_t first_incoming)
{
    /* Retrieve the user data which contains the variable's data vectors */
    cmc_variable_interpolation_data_t interpolation_data = static_cast<cmc_variable_interpolation_data_t>(t8_forest_get_user_data(forest_new));

    if (refine == 1)
    {
        /** If the element is refined **/
        /* Copy the element's value from the coarse element to all finer elements which replace the coarser element */
        for (int child_id = 0; child_id < num_incoming; ++child_id)
        {
            /* Since the iteration of the elements is linear (and we are focusing on a serial execution environment), we can just push_back the values to the vector holding the interpolated data (instead of assigning the value at an specific index) */
            interpolation_data->data_new.push_back(interpolation_data->data[first_outgoing]); //This works only in the current case, a more general assignemnt would be: interpolation_data->data_new[first_incoming + child_id] = interpolation_data->data[first_outgoing]; Actually, even better would be to replace the whole for-loop by: interpolation_data->data_new.insert(interpolation_data->data_new.begin() + first_incoming, num_incoming, interpolation_data->data[first_outgoing]);
        }
    } else if (refine == 0)
    {
        /** The element stays the same **/
        /* When the element stays the same, we can just copy the value over */
        /* Since the iteration of the elements is linear (and we are focusing on a serial execution environment), we can just push_back the values to the vector holding the interpolated data (instead of assigning the value at an specific index) */
        interpolation_data->data_new.push_back(interpolation_data->data[first_outgoing]); //This works only in the current case (serial execution, linearily ordered, only one tree, etc.), a more general assignemnt would be: interpolation_data->data_new[first_incoming] = interpolation_data->data[first_outgoing]; or even better, use: interpolation_data->data_new.insert(interpolation_data->data_new.begin() + first_incoming, interpolation_data->data[first_outgoing]);
    } else
    {
        /** The element has been refined **/
        cmc_err_msg("A coarsening during the refinement/decompression is illegal.");
    }
}

/**
 * @brief This function refines an element if it's level is below (i.e. coarser) than the initial refinement level (saved within the adapt data)
 */
int
cmc_t8_adapt_callback_refine_to_initial_lvl (t8_forest_t forest,
                                             t8_forest_t forest_from,
                                             int which_tree,
                                             int lelement_id,
                                             t8_eclass_scheme_c * ts,
                                             const int is_family,
                                             const int num_elements,
                                             t8_element_t * elements[]) 
{
    /* Get the adapt data from the forest */
    cmc_refine_adapt_data_t adapt_data = static_cast<cmc_refine_adapt_data_t>(t8_forest_get_user_data(forest));
    
    /* Check if the element is already on the initial refinement level (or if it is still coarser) */
    if(t8_element_level(ts, elements[0]) < adapt_data->initial_refinement_level)
    {
        /* Eventually, update the coarsest element level which is currently present within the mesh */
        if (t8_element_level(ts, elements[0]) < adapt_data->coarsest_element_level_present_in_the_forest)
        {
            /* Save the coarsest level; this indicates how many iterations we still need to perform in order to reach the initial_refienment_level for all elements */
            adapt_data->coarsest_element_level_present_in_the_forest = t8_element_level(ts, elements[0]) + 1; // We add one because after the callback, the element is refined and therefore already one level finer than right now
        }

        /* If the element's level is still coarser than the initial refinement level, we can refine the element */
        return 1;
    } else
    {
        /* If we have already reached the initial refinement level, the element stays the same */
        return 0;
    }
}

/**
 * @brief This function refines the coarse forest iteratively (back) on to the supplied initial refinement level.
 * Applying this refinement resembles the decompression step
 * 
 * @param forest The coarse forest which will be refined
 * @param initial_refinement_level The initial refinement level, onto which all elements will be refined
 * @return t8_forest_t The function returns the refined (uniform) forest
 */
t8_forest_t
cmc_refine_the_forest_iteratively_to_the_initial_level(t8_forest_t forest, const int initial_refinement_level, std::vector<double>& variable_data)
{
    /** Within this functionn we refine the forest back onto the initial refinement level **/
    /* Declare a new forest variable */
    t8_forest_t forest_adapt;

    /* A variable holding the name of the file in which the forest will be written */
    std::string file_prefix = std::string("cmc_t8_forest_refinement_step");

    /* Define a struct for the adapt data which the forest adapt callback needs in order to "get to know" the initial refinement level */
    cmc_refine_adapt_data adapt_data;

    /* Save the intial refinement level within the adapt data */
    adapt_data.initial_refinement_level = initial_refinement_level;
    adapt_data.coarsest_element_level_present_in_the_forest = 0;

    /* Create a new struct for the interpolation data (holding the variable data) */
    cmc_variable_interpolation_data interpolation_data;

    /* Save the variable data within the interpolation data struct (by swapping the empty vector in the struct with the variable data which was passed to this function )*/
    std::swap(interpolation_data.data, variable_data); //The data of the variable is from now on only accessible via "interpolation_data.data" and no longer via variable_data
    
    /* A counter for the refinement steps */
    int refinement_step = 1;

    /* We refine until we have reached the intial level for each element (this is the case
     * when the elements of the former and the current forest do not change,
     * but at latest when we have perfomed an amount of refinement steps equal to the intial refinement level) */
    while (adapt_data.coarsest_element_level_present_in_the_forest < initial_refinement_level && refinement_step <= initial_refinement_level)
    {
        cmc_debug_msg("Refinement step: ", refinement_step);
        
        /* Reset the coarsest present element level within the forest */
        adapt_data.coarsest_element_level_present_in_the_forest = initial_refinement_level;

        /* Initialize the new forest */
        t8_forest_init(&forest_adapt);

        /* Set the user data (needed for the adaption step) */
        t8_forest_set_user_data(forest_adapt, static_cast<void*>(&adapt_data));

        /* We need to reference the 'old' forest, otherwise it will be deleted by the "set_adapt_call" on commit. This has to be done because, we need both - the old and the new forest - in order to perform the interpolation */
        t8_forest_ref(forest);

        /* Adapt the forest accordingly to the callback function */
        t8_forest_set_adapt(forest_adapt, forest, cmc_t8_adapt_callback_refine_to_initial_lvl, 0);

        /* Commit the adapted forest (perform the adaption step) */
        t8_forest_commit(forest_adapt);

        /* Set the interpolation data accordingly */
        t8_forest_set_user_data(forest_adapt, static_cast<void*>(&interpolation_data));
        
        /* Allocate memory for the interpolated data (which will be transferred from 'forest' to 'forest_adapt') equal to the amount of elements of 'forest_adapt' */
        interpolation_data.data_new.reserve(t8_forest_get_local_num_elements(forest_adapt));

        /* Interpolate the element data onto the new refined forest */
        t8_forest_iterate_replace(forest_adapt, forest, cmc_t8_interpolation_function_for_refinement);
    
        /* Swap the 'old' variable data which was defined on 'forest' with the 'new' variable data which was interpolated onto 'forest_adapt' */
        std::swap(interpolation_data.data, interpolation_data.data_new);

        /* Delete the 'old' variable data (which was defined on 'forest') */
        interpolation_data.data_new.clear();

        /* After the interpolation of the variable data, we do not need the 'old' forest anymore, therefore, it can be dereferenced/deallocated */
        t8_forest_unref(&forest);

        /* Save the refined forest */
        forest = forest_adapt;

        /* We only want to write out the forest if something has changed in comparison to the former iteration */
        cmc_t8_forest_write_vtk(forest, (file_prefix + std::to_string(refinement_step)).c_str(), interpolation_data.data);

        /* Update iteration count variable */
        ++refinement_step;
    }

    /* We swap the vectors holding the variable' data back, such that the refined/decompressed variable data now resides again in the vector 'variable_data' */
    /* Such that the compressed data can be accessed by the calling function */
    std::swap(interpolation_data.data, variable_data);

    /* Return the forest which was refined back onto the initial refinement level */
    return forest;
}
/************************************************************************/


struct my_element {
    std::array<double, 3> midpoint{0.0, 0.0, 0.0}; 
    //int    is_inside_partition; /*Frage: Brauche ich das wirklich? Vermutlich schon, aber ...*/
    //double value; /*Contains the value that we want to copy in our uniform forest. */
    t8_locidx_t corresponding_coarse_element_id; 
};

static int t8_search_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf, 
                                t8_element_array_t *leaf_elements, t8_locidx_t tree_leaf_index, void *query, size_t query_index)
{
    T8_ASSERT (query == NULL);         /*Frage: Nicht !=NULL?*/     
    
    return 1;
}

static int t8_search_query_callback (t8_forest_t forest, t8_locidx_t ltreeid, const t8_element_t *element, const int is_leaf,
                                    t8_element_array_t *leaf_elements, t8_locidx_t tree_leaf_index, void *query, size_t query_index)
{
    int particle_is_inside_element;

    T8_ASSERT (query != NULL);

    /* Cast the query pointer to an element pointer. */
    my_element *particle = (my_element *) query;

    /*Numerical tolerance for the is inside element check. */
    const double tolerance = 1e-8;

    //T8_ASSERT (query != NULL);

    /*Beachte src/t8_forest/t8_forest_general.h Zeile 764. */
    /*Test whether this particle is inside this element. */
    particle_is_inside_element = t8_forest_element_point_inside (forest, 0, element, particle->midpoint.data(), tolerance);

    if (particle_is_inside_element) {
        if (is_leaf) {
            /* The particle is inside the element and is a leaf element. 
            * We mark the particle for being inside the partition. */
            /* In order to find the index of the element inside the array, we compute the 
            * index of the element within the tree. */
            /*Beachte src/t8_forest/t8_forest_general.h Zeile 662. */

            /*Nochmal schauen, ob das wirklich die richtige Funktion ist. */
            /*Die ID würde ich jetzt gerne in  t8_locidx_t corresponding_coarse_element_id abspeichern, aber 
            * bin mir noch nicht wirklich sicher, wie das klappen soll. Muss ich mir genauer anschauen.; */
            particle->corresponding_coarse_element_id = tree_leaf_index;
            
        }
        /* The particle is inside the element. The query should remain active.
        * If this element is not a leaf the search will continue with its children. */
        return 1;
    }
    /*The particle is not inside the element. Deactivate the query. 
    * If no active queries are left, the search will stop for this element and its children. */
    return 0;

}
/************************************************************************/
/* The function below applies the (complete) decompression in just a single adaptation step */
/********* This is done with the help of the search-functions implemented in t8code *********/
/**
 * @brief This function refines the coarse forest (i.e. decompresses it) with the help of the search algorithm back onto the
 * supplied initial refinement level.
 * The search algorithm is used to perform the adaptation/refinement from the coarse resolution onto the initial resolution (given by the
 * initial refinement level) in just one step, unlike the other approach which will perform the decompression iteratively (getting only one level
 * finer per iteration until the initial refinement level is finally reached).
 *
 * @param coarse_forest The coarse forest which will be refined
 * @param initial_refinement_level The initial refinement level onto which the elements will be refined
 * @return t8_forest_t This functions returns the refined/decompressed forest
 */
t8_forest_t
cmc_refine_the_forest_by_search_algorithm_to_the_initial_level(t8_forest_t coarse_forest, const int initial_refinement_level, std::vector<double>& variable_data)
{

    // Here continues your part
    // TODO: implement
    
    // Some ideas:
    //     1. Create a uniform forest on the initial_refinement_level
    //     2. Compute for each element of the uniform forest the midpoint (e.g. see function "t8_forest_element_centroid(...)" (within src/t8_forest/t8_forest_geometrical.h))
    //     3. Search the midpoints of the (uniform forest's) elements in the coarse_forest -> find the (leaf) element in the coarse_forest containing the (uniform forest's) elements' midpoints
    //         3b. Take a look at the function "t8_forest_element_point_inside(...)" (within src/t8_forest/t8_forest_general.h) (this checks whether a point is inside an element or not, i.e. check if the midpoints of the uniform forest's elements lie within the corse forest's elements)
    //     4. Copy the element data from the coarse_forest's element to the uniform forest's element
    //
    // Note: Since the level of the elements within forests differ by more than one level generally, we cannot use "t8_forest_iterate_replace(...)"
    //       We need to copy the data manually to new vector (but this should not be too much overhead and should not be too complicated).


    /* Create a cmesh (based only on one quadrilateral or one hexahedron) */    
    t8_cmesh_t cmesh;

    const int dimension_of_forest = t8_geom_get_dimension(t8_cmesh_get_tree_geometry(t8_forest_get_cmesh(coarse_forest), 0));

    if (dimension_of_forest == 2)
    {   
        cmesh = t8_cmesh_new_hypercube(T8_ECLASS_QUAD, MPI_COMM_WORLD, 0, 0, 1);
    
    } else 
    {
        cmesh = t8_cmesh_new_hypercube(T8_ECLASS_HEX, MPI_COMM_WORLD, 0, 0, 1);
    }

    /* Declare a new forest which is uniform and on the initial refinement level. */
    t8_forest_t forest_initial_level = t8_forest_new_uniform(cmesh, t8_scheme_new_default_cxx(), initial_refinement_level, 0, MPI_COMM_WORLD);

    sc_array_t          *list_of_elements;

    /* Create an array for num_particles many particles. */
    list_of_elements = sc_array_new_count (sizeof (my_element), t8_forest_get_global_num_elements(forest_initial_level));

    t8_element_t* element;

    for (t8_locidx_t element_id = 0; element_id < t8_forest_get_global_num_elements(forest_initial_level); element_id++) 
    {
        /* Get this particle's pointer. */
        my_element* my_elem = (my_element *) sc_array_index_int (list_of_elements, element_id);

        element = t8_forest_get_element_in_tree (forest_initial_level, 0, element_id);

        /* Compute the coordinates of the centroid of the current element of forest_initial_level. */
        t8_forest_element_centroid (forest_initial_level, 0, element, my_elem->midpoint.data());
    }

    /* Perform the search of the forest. The second argument is the search callback function,
    * then the query callback function and the last argument is the array of queries. */
    t8_forest_search (coarse_forest, t8_search_callback,
                    t8_search_query_callback, list_of_elements);

                    /*Frage: Wie werden nun die benötigten Variablen an die Funktionen t8_search_query_callback und t8_search_callback 
                    * übergeben? */
    /*Frage (mehr an mich): haben wir jetzt eine Liste, in der steht, wo sich welche Daten befinden? Sind diese sortiert? Damit wir einfach 
    * wieder durch iterieren können?*/

    /*Mappen von corresponding_coarse_element_id auf uniform forest und daraus vtu datei machen mit daten drauf: 
    *Double Array erstellen, in größe von uniformen forest und da corresponding_coarse_element_id als double speichern (casten).
    *forest mit double array rausschreiben.
    *Beispiel: zeile 235 - 251. 
    *Dann schauen wir, obs passt.--- Das sieht so aus, als ob es passt. */    
    #if 0
    /*Vector to contain the corresponding_coarse_element_id. */
    std::vector<double> corresponding_coarse_element_id_to_vtu;
    corresponding_coarse_element_id_to_vtu.reserve(t8_forest_get_global_num_elements(forest_initial_level));
    for (t8_locidx_t i = 0; i < t8_forest_get_global_num_elements(forest_initial_level); ++i) {
        
        double casted_element_id = static_cast<double>(((my_element *) t8_sc_array_index_locidx (list_of_elements, i))->corresponding_coarse_element_id);
        corresponding_coarse_element_id_to_vtu.push_back(casted_element_id);
    }

    /*for - schleife zum kopieren. neuen double (data_new) vector erzeugen, in der länge von initial forest, 
    * und dann daten von variable_data kopieren.
    * daten aus neuem vector in variable_data kopieren.*/
    
    /*Vector to contain the data for the uniform forest. */
    /**/
    std::vector<double> data_new;
    data_new.reserve(t8_forest_get_global_num_elements(forest_initial_level));
    for (t8_locidx_t i = 0; i < t8_forest_get_global_num_elements(forest_initial_level); ++i) {
        /*Get the element, which contains the element of the uniform forest, from the coarse forest 
        * and copy the data according to that. */
        double element_to_search_for = corresponding_coarse_element_id_to_vtu[i];
        data_new.push_back(variable_data[element_to_search_for]);
    }
    #endif
    /*Alternative way without corresponding_coarse_element_id_to_vtu*/
    std::vector<double> data_new;
    data_new.reserve(t8_forest_get_global_num_elements(forest_initial_level));
    for (t8_locidx_t i = 0; i < t8_forest_get_global_num_elements(forest_initial_level); ++i) {
        /*Get the element, which contains the element of the uniform forest, from the coarse forest 
        * and copy the data according to that. */
        double element_to_search_for = static_cast<double>(((my_element *) t8_sc_array_index_locidx (list_of_elements, i))->corresponding_coarse_element_id);
        data_new.push_back(variable_data[element_to_search_for]);

    }

    t8_vtk_data_field_t enhanced_decompression_data;

    enhanced_decompression_data.data = (double *) data_new.data();
    strcpy (enhanced_decompression_data.description, "Example data new");
    enhanced_decompression_data.type = T8_VTK_SCALAR;

    t8_forest_write_vtk_ext (forest_initial_level, "Enhanced_decompressed_forest", 1, 1, 1, 1, 0, 0, 0, 1, &enhanced_decompression_data);

    /* We swap the vectors holding the variable' data back, such that the refined/decompressed variable data now resides again in the vector 'variable_data' */
    /* Such that the decompressed/interpolated data can be accessed by the calling function */
    std::swap(data_new, variable_data);
    sc_array_destroy(list_of_elements);
    t8_forest_unref(&coarse_forest);
    return forest_initial_level;

}
/************************************************************************/


/**
 * @brief The main function of this example program. It creates a uniform forest.
 * Afterwards the forest will be arbitrarily coarsened and finally refined back to the initial refinement level 
 */
int
main(int argc, char* argv[])
{
    /* Initialize cmc (the SC-library, t8code and eventually MPI will be initialized by this function call) */
    cmc_initialize();

    /* Create a cmesh (based only on one quadrilateral or one hexahedron) */
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube(T8_ECLASS_QUAD, MPI_COMM_WORLD, 0, 0, 1);

    cmc_debug_msg("A new cmesh based on a quadrilateral has been created.");

    /* Define an arbitrary initial refinement level */
    const int initial_refinement_level = 7;

    /* Create a new forest from the cmesh and uniformly refine the forest onto level 7 */
    t8_forest_t forest = t8_forest_new_uniform(cmesh, t8_scheme_new_default_cxx(), initial_refinement_level, 0, MPI_COMM_WORLD);

    cmc_debug_msg("A uniform forest of refinement level ", initial_refinement_level, " has been created.");

    /* We want to associate some data with the initial forest */
    std::vector<double> element_data; /* Declare a vector which will conatain the element data */
    cmc_create_example_data(forest, element_data); /* Fill the vector with some example data */

    /* We write out the intial forest to .vtu files */
    cmc_t8_forest_write_vtk(forest, "cmc_t8_initial_forest", element_data);

    /* Coarsen the forest */
    forest = cmc_coarsen_the_initial_forest_arbitrarily(forest, element_data); //The element data vector will contain the variable's data which was interpolated onto the coarsened forest

    cmc_debug_msg("The initial uniform forest has been coarsened.");

    /****************************************************************/
    /* We reference the coarse forest once more, which means keep the forest after the first decompression/refinement 
     * in order to refine/decompress the coarse forest again by the new approach (refinement by search)
     */
    t8_forest_ref(forest);

    /* We copy the variable's data in order to refine/decompress the coarse forest with the coarse variable's data again by the new approach (refinement by search) */
    std::vector<double> element_data_copy = element_data;

    /******** Here starts the interesting decompression part ********/
    
    /* Declare a new forest variable */
    t8_forest_t decompressed_forest;

    /* Get the current time point before the iterative decompression */
    auto start_time_iterative = std::chrono::high_resolution_clock::now();

    /* Refine the forest (iteratively) back onto the intial refinement level */
    decompressed_forest = cmc_refine_the_forest_iteratively_to_the_initial_level(forest, initial_refinement_level, element_data_copy);

    /* Get the time point after the iterative decompression */
    auto end_time_iterative = std::chrono::high_resolution_clock::now();

    /* Calculate the duration of the iterative decompression */
    auto duration_iterative = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_iterative - start_time_iterative);

    cmc_debug_msg("The duration of the iterative decompression took: ", duration_iterative.count(), "ms"); //Eventually change milliseconds to milliseconds or something else

    /* Deallocate the iteratively decompressed forest */
    t8_forest_unref(&decompressed_forest);

    cmc_debug_msg("The iteratively decompressed forest has been deallocated");

    // Here starts your part

    /* Get the current time point before the decompression by search */
    auto start_time_by_search = std::chrono::high_resolution_clock::now();

    /* Refine the forest (interpolation by search) back onto the initial refinement level */
    decompressed_forest = cmc_refine_the_forest_by_search_algorithm_to_the_initial_level(forest, initial_refinement_level, element_data);

    /* Get the time point after the decompression by search */
    auto end_time_by_search = std::chrono::high_resolution_clock::now();

    /* Calculate the duration of the iterative decompression by search */
    auto duration_by_search = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_by_search - start_time_by_search);

    cmc_debug_msg("The duration of the decompression by search took: ", duration_by_search.count(), "ms"); //Eventually change milliseconds to milliseconds or something else

    /* Deallocate the iteratively decompressed forest */
    t8_forest_unref(&decompressed_forest);

    cmc_debug_msg("The decompressed forest with the help of the search algorithm has been deallocated");

    /****************************************************************/

    /* Let's check whether the data has been decompressed equally. Therefore, both element_data vectors should contain the same values at the same index */
    cmc_check_decompression_results_for_equality(element_data, element_data_copy);

    /* Finalize cmc (the SC-library, t8code and eventually MPI will be finalized by this function call) */
    cmc_finalize();
}
