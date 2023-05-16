#include "cmc.h"
#include "mpi/cmc_mpi.h"
#include "utilities/cmc_log_functions.h"
#include "netcdf/cmc_netcdf.h"
#include "netcdf/cmc_netcdf.hxx"
#include "t8code/cmc_t8code.h"
#include "t8code/cmc_t8code_data.hxx"
#include "t8code/cmc_t8code_geo_mesh.h"
#include "t8code/cmc_t8code_geo_data.h"
#include "component_interfaces/cmc_t8_nc.h"
#include "lossy/cmc_amr_compressor.h"


int
main(int argc, char* argv[])
{
  /* Initialize cmc */
  cmc_initialize();

  {
    /* Open a netCDF file with geo-spatial data */
    int ncid = cmc_nc_open("../../data/ECMWF_ERA-40_subset.nc");

    /* Create a class holding the data from the netCDF-File */
    cmc_nc_data_t nc_data = cmc_nc_create(ncid);

    /* Define a hyperslab of the data */
    /* This hyperslab will be read in and used for the compression */
    const size_t start_ptr[3] = {0,0,0};  //Example netCDF File
    const size_t count_ptr[3] = {1,73,144}; //Example netCDF File

    /* Inquire given data/variables which are defined on a geo-spatial domain (latitude, longitude, height, (time)) */
    cmc_nc_inquire_vars(nc_data, start_ptr, count_ptr, "p2t");

    /* Close the netCDF file */
    cmc_nc_close(ncid);

    /* Define data classes holding the variable data and forests as well as additional information during the compression process */
    cmc_amr_data_t amr_data;

    /* Create/Allocate a new 'AMR_DATA' class for the compression of netCDF inquired data */
    amr_data = cmc_create_amr_compression_data(nc_data, MPI_COMM_WORLD);

    /* Set a compression criterium - e.g. error threshold woth a predefined tolerance */
    cmc_amr_pre_setup_set_compression_criterium_error_threshold(amr_data, 0.02);

    /* Setup the compression for a given 'compression mode' */
    cmc_amr_setup_compression(amr_data, CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE_2D);

    /* Execute the adaptation/compression */
    cmc_amr_compress(amr_data);

    /* Write a vtk file of the compressed data */
    cmc_amr_write_vtk_file(amr_data, "compressed_data");

    /* Decompress the data */
    cmc_amr_decompress(amr_data);

  	/* Write a vtk file of the decompressed data */
    cmc_amr_write_vtk_file(amr_data, "decompressed_data");

    /* Deallocate the netCDF data */
    cmc_nc_destroy(nc_data);

    /* Deallocate the Lossy AMR Compression data */
    cmc_amr_destroy(amr_data);
  }

  /* Finalize cmc */
  cmc_finalize();

  return 0;
}
