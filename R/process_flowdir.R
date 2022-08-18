
<<<<<<< HEAD
#' Title
#'
#' @param dem
#' @param threshold
#' @param output_filename
#' @param return_products
#' @param temp_dir
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
#'
=======
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d

process_flowdir<-function(
    dem,
    threshold,
<<<<<<< HEAD
    output_filename,
=======
    save_dir,
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
    return_products=F,
    temp_dir=NULL,
    verbose=F
) {

  require(sf)
  require(terra)
  require(whitebox)
  require(tidyverse)

  if (!is.integer(threshold)) stop("'threshold' must be an integer value")
<<<<<<< HEAD

=======
  if (!dir.exists(save_dir)) dir.create(save_dir)
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
<<<<<<< HEAD
  temp_dir<-normalizePath(temp_dir)
  output_filename<-normalizePath(output_filename,mustWork =F)
  if (!grepl("\\.zip$",output_filename)) stop("output_filename must be a character ending in '.zip'")
=======
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d

  wbt_options(exe_path=wbt_exe_path(),
              verbose=verbose,
              wd=temp_dir)

<<<<<<< HEAD
  terra::terraOptions(verbose = verbose,
                      tempdir = temp_dir
  )

=======
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
  dem<-hydroweight::process_input(dem,input_name="dem")
  if (!inherits(dem,"SpatRaster")) stop("dem must be a class 'SpatRaster'")
  target_crs<-crs(dem)

<<<<<<< HEAD
  writeRaster(dem,file.path(temp_dir,"dem_final.tif"),overwrite=T,gdal="COMPRESS=NONE")
=======
  writeRaster(dem,file.path(temp_dir,"dem_final.tif"),overwrite=T)
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d


  # Process flow dirrection -------------------------------------------------
  ## Generate d8 pointer
  if (verbose) print("Generating d8 pointer")
  wbt_d8_pointer(
    dem="dem_final.tif",
    output="dem_d8.tif"
  )

  ## Generate d8 flow accumulation in units of cells
  if (verbose) print("Generating d8 flow accumulation")
  wbt_d8_flow_accumulation(
    input = "dem_d8.tif",
    output = "dem_accum_d8.tif",
    out_type = "cells",
    pntr=T
  )
  wbt_d8_flow_accumulation(
    input = "dem_d8.tif",
<<<<<<< HEAD
    output = "dem_accum_d8_sca.tif",
=======
    output = "dem_accum_d8.tif",
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
    out_type = "sca",
    pntr=T
  )


  ## Generate streams with a stream initiation threshold
  if (verbose) print("Extracting Streams")
  wbt_extract_streams(
    flow_accum = "dem_accum_d8.tif",
    output = "dem_streams_d8.tif",
    threshold = threshold
  )

<<<<<<< HEAD
  # Prepare Output ----------------------------------------------------------

=======
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
  dist_list_out<-list(
    "dem_final.tif",
    "dem_d8.tif",
    "dem_accum_d8.tif",
<<<<<<< HEAD
    "dem_accum_d8_sca.tif",
    "dem_streams_d8.tif"
=======
    "dem_streams_d8.tif",
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
  )

  dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

<<<<<<< HEAD
  out_file<-output_filename
  out_dir<-gsub(basename(out_file),"",out_file)
  if (!dir.exists(out_dir)) dir.create(out_dir)

  if (file.exists(out_file)) {
    out_file<-file.path(out_dir,gsub("\\.zip","",basename(out_file)),"_",paste0(basename(tempfile()), ".zip"))
=======
  out_file<-file.path(save_dir,paste0("Processed_Hydrology.zip"))
  if (file.exists(out_file)) {
    out_file<-file.path(save_dir,paste0(basename(tempfile()), "Processed_Hydrology.zip"))
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
    warning(paste0("Target .zip file already exists. Saving as: ",out_file))
  }

  if (verbose) print("Generating Output")

  zip(out_file,
      unlist(dist_list_out),
      flags = '-r9Xjq'
  )

  output<-list(
    outfile=out_file
  )

  if (return_products){
    output<-c(
<<<<<<< HEAD
      lapply(setNames(dist_list_out,basename(unlist(dist_list_out))),function(x) wrap(rast(x))),
=======
      lapply(dist_list_out,function(x) wrap(rast(x))),
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
      output
    )
  }
  file.remove(list.files(temp_dir,full.names = T))

  return(output)
}
