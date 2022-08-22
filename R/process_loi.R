
#' processes layers of interest (loi) for attribution of stream network
#'
#' @param input \code{NULL} or output from `process_hydrology()` or `process_flowdir()`. If \code{NULL}, `dem` must be specified.
#' @param dem \code{NULL} or character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/dem.tif"), \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of GeoTiFF type. Digital elevation model raster. If \code{NULL}, input must be specified.
#' @param clip_region character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/lu.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of ESRI Shapefile type or GeoTiFF type only. Region over which loi are returned. defaults to `dem` extent.
#' @param num_inputs named list containing file paths (with extensions), and/or R GIS objects to be coerced into \code{SpatRaster}. Layers of interest to be summarized numerically (e.g., with mean, SD, etc.)
#' @param cat_inputs named list containing file paths (with extensions), and/or R GIS objects to be coerced into \code{SpatRaster}. Layers of interest to be summarized categorical (e.g., with with proportions or weighted proportions in the catchment)
#' @param output_filename \code{NULL} or character (full file path with extensions). If `input` is provided, loi are saved in the same .zip output folder regadless of `output_filename`. If `input` not provided, file path to write resulting .zip file.
#' @param variable_names named list containing variables from `num_inputs` and `cat_inputs` to include in output. If not specied, all variables are used.
#' @param return_products logical. If \code{TRUE}, a list containing all geospatial analysis products. If \code{FALSE}, folder path to resulting .zip file.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'
#' @examples
#'
#'
process_loi<-function(
    input=NULL,
    dem=NULL,
    clip_region=NULL,
    num_inputs=list(),
    cat_inputs=list(),
    output_filename=NULL,
    variable_names=NULL,
    return_products=F,
    temp_dir=NULL,
    verbose=F
) {

  # require(sf)
  # require(terra)
  # require(whitebox)
  # require(tidyverse)

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")
  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  if (!inherits(num_inputs,"list")) stop("'num_inputs' should be a named list")
  if (length(num_inputs)>0 & is.null(names(num_inputs))) stop("objects in 'num_inputs' must be named")
  if (!inherits(cat_inputs,"list")) stop("'cat_inputs' should be a named list")
  if (length(cat_inputs)>0 & is.null(names(cat_inputs))) stop("objects in 'cat_inputs' must be named")

  if (any(names(num_inputs) %in% names(cat_inputs)) |
      any(names(cat_inputs) %in% names(num_inputs))) stop("'num_inputs' and 'cat_inputs' cannot share names")

  if (!is.null(variable_names)) {
    # Add a check to make sure names match
  } else {
    message("No variables specified in 'variable_names', all variables in 'num_inputs' and 'cat_inputs' will be used")
  }

  if (is.null(input) & is.null(dem)) stop("Either 'input' or 'dem' must be specifed")
  if (is.null(input) & is.null(output_filename)) stop("Either 'input' or 'output_filename' must be specifed")

  if (!is.null(input) & is.null(output_filename)) {
    output_filename<-input$outfile
  }
  output_filename<-normalizePath(output_filename,mustWork =F)
  if (!grepl("\\.zip$",output_filename)) stop("output_filename must be a character ending in '.zip'")

  if (gsub(basename(output_filename),"",output_filename) == temp_dir) stop("'output_filename' should not be in the same directory as 'temp_dir'")


  # Prepare DEM -------------------------------------------------------------
  if (verbose) print("Preparing DEM")
  if (!is.null(input)){
    zip_loc<-input$outfile
    fl<-unzip(list=T,zip_loc)

    unzip(zip_loc,
          c("dem_final.tif"),
          exdir=temp_dir,
          overwrite=T,
          junkpaths=T)

    dem<-rast(file.path(temp_dir,"dem_final.tif"))
  } else {
    dem<-hydroweight::process_input(dem,input_name="dem",working_dir=temp_dir)
    if (!inherits(dem,"SpatRaster")) stop("dem must be a class 'SpatRaster'")
    writeRaster(dem,file.path(temp_dir,"dem_final.tif"),overwrite=T,gdal="COMPRESS=NONE")
  }

  if (is.null(clip_region)) {
    clip_region<-as.polygons(rast(ext(dem),crs=crs(dem)))
    writeVector(clip_region,file.path(temp_dir,"clip_region.shp"),overwrite=T)
  } else {
    clip_region<-hydroweight::process_input(clip_region,
                                            input_name="clip_region",
                                            working_dir=temp_dir,
                                            target=terra::vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",crs=terra::crs(dem)))

    writeVector(clip_region,file.path(temp_dir,"clip_region.shp"),overwrite=T)
  }

  target_crs<-crs(dem)

  # Prepare loi -------------------------------------------------------------
  ## Numeric loi

  if (verbose) print("Preparing Numeric Inputs")
  num_inputs_list<-list(lyr_nms=names(num_inputs),
                         lyr=num_inputs,
                         lyr_variables=variable_names[names(num_inputs)])

  num_inputs<-future_pmap(num_inputs_list,
                    function(lyr_nms,lyr,lyr_variables){
                      output<-hydroweight::process_input(
                        input=lyr,
                        input_name = lyr_nms,
                        variable_name=lyr_variables,
                        target=file.path(temp_dir,"dem_final.tif"),
                        clip_region = file.path(temp_dir,"clip_region.shp"),
                        resample_type = "bilinear"
                      )

                      ot<-writeRaster(output,file.path(temp_dir,paste0(lyr_nms,".tif")),overwrite=T,gdal="COMPRESS=NONE")
                      return(file.path(temp_dir,paste0(lyr_nms,".tif")))
                    })

  ## Categorical loi

  if (verbose) print("Preparing Categorical Inputs")
  cat_inputs_list<-list(lyr_nms=names(cat_inputs),
                        lyr=cat_inputs,
                        lyr_variables=variable_names[names(cat_inputs)])

  cat_inputs<-future_pmap(cat_inputs_list,
                   function(lyr_nms,lyr,lyr_variables){
                     output<-hydroweight::process_input(
                       input=lyr,
                       input_name = lyr_nms,
                       variable_name=lyr_variables,
                       target=file.path(temp_dir,"dem_final.tif"),
                       clip_region = file.path(temp_dir,"clip_region.shp"),
                       resample_type = "near"
                     )

                     ot<-writeRaster(output,file.path(temp_dir,paste0(lyr_nms,".tif")),overwrite=T,gdal="COMPRESS=NONE")
                     return(file.path(temp_dir,paste0(lyr_nms,".tif")))
                   })


  # Combine loi -------------------------------------------------------------
  if (verbose) print("Combining Numeric Inputs")
  num_inputs<-rast(map(num_inputs,rast))
  writeRaster(num_inputs,file.path(temp_dir,"num_rast.tif"),overwrite=T,gdal="COMPRESS=NONE")

  if (verbose) print("Combining Categorical Inputs")
  cat_inputs<-rast(map(cat_inputs,rast))
  writeRaster(cat_inputs,file.path(temp_dir,"cat_rast.tif"),overwrite=T,gdal="COMPRESS=NONE")


  # Generate Output ---------------------------------------------------------
  if (verbose) print("Generating Outputs")

  dist_list_out<-list(
    "num_rast.tif",
    "cat_rast.tif"
  )

  dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

  zip(output_filename,
      unlist(dist_list_out),
      flags = '-r9Xjq'
  )

  output<-input

  if (return_products){
    output<-c(
      list(
        num_inputs=wrap(num_inputs),
        cat_inputs=wrap(cat_inputs)
      ),
      output
    )
  }

  file.remove(list.files(temp_dir,full.names = T))

  return(output)

}
