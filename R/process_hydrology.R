#' Process DEM by generating various terrain products
#'
#' #' \code{IntegratedHydrology::process_hydrology()} processes a DEM into various geospatial analysis products.
#'
#' @param dem character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/dem.tif"), \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of GeoTiFF type. Digital elevation model raster.
#' @param output_filename character. File path to write resulting .zip file.
#' @param threshold numeric. Threshold in flow accumulation values for channelization.
#' @param return_products logical. If \code{TRUE}, a list containing all geospatial analysis products. If \code{FALSE}, folder path to resulting .zip file.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical. If \code{FALSE}, the function will not print output prints.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'
#' @examples
#'
#' library(whitebox)
#' library(terra)
#' library(sf)
#' library(hydroweight)
#'
#' ## Import toy_dem from whitebox package
#' toy_file<-sample_dem_data()
#' toy_file <- system.file("extdata", "DEM.tif", package = "whitebox")
#' toy_dem <- rast(raster::raster(x = toy_file)) # reading the file from terra directly sometimes crashes R for some reason
#' crs(toy_dem) <- "epsg:3161"
#'
#' ## Generate hydroweight_dir as a temporary directory
#' save_dir <- tempdir()
#'
#' hydro_out<-process_hydrology(
#'   dem=toy_dem,
#'   threshold=1000,
#'   return_products=T,
#'   save_dir=save_dir,
#'   temp_dir=NULL,
#'   verbose=F
#' )
#'
#
process_hydrology<-function(
    dem,
    output_filename,
    threshold,
    extra_attr=c( # Additional attributes to calculate and add to outputs
      "link_slope",
      "cont_slope",
      "USChnLn_To",
      "Elevation",
      "StOrd_Hack",
      "StOrd_Str",
      "StOrd_Hort",
      "StOrd_Shr"
    ),
    points=NULL,
    snap_distance=NULL,
    site_id_col=NULL,
    return_products=F,
    temp_dir=NULL,
    verbose=F
) {

  # require(sf)
  # require(terra)
  # require(whitebox)
  # require(tidyverse)

  if (!is.integer(threshold)) stop("'threshold' must be an integer value")
  if (!is.integer(snap_distance)) stop("'snap_distance' must be an integer value")

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

  if (!is.null(points)){
    points<-hydroweight::process_input(points,input_name="points")
    points<-st_as_sf(points)
    if (!is.character(site_id_col)) stop("'site_id_col' must be a single character")
    if (length(site_id_col)>1 | length(site_id_col)==0) stop("length 'site_id_col' must be 1")
    if (!site_id_col %in% names(points)) stop("'site_id_col' must be a variable name in 'points'")
  }

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)
  output_filename<-normalizePath(output_filename,mustWork =F)
  if (!grepl("\\.zip$",output_filename)) stop("output_filename must be a character ending in '.zip'")

  extra_attr<-match.arg(extra_attr,several.ok = T)

  if (verbose) print("Processing Flow Direction")
  hydro_out<-process_flowdir(
    dem=dem,
    threshold=threshold,
    return_products=return_products,
    output_filename=output_filename,
    temp_dir=temp_dir,
    verbose=verbose
  )

  if (verbose) print("Generate Subbasins")
  hydro_out<-generate_subbasins(
    input=hydro_out,
    return_products=return_products,
    temp_dir=temp_dir,
    verbose=verbose
  )

  if (verbose) print("Attribute Streamlines")
  hydro_out<-attrib_streamline(
    input=hydro_out,
    extra_attr=extra_attr,
    return_products=return_products,
    temp_dir=temp_dir,
    verbose=verbose
  )

  if (!is.null(points)) {
    if (verbose) print("Inserting Sampling Points")
    hydro_out<-insert_points(
      input=hydro_out,
      points=points,
      site_id_col=site_id_col,
      snap_distance=snap_distance,
      return_products=return_products,
      temp_dir=temp_dir,
      verbose=verbose
    )
  }

  if (verbose) print("Tracing Flowpaths")
  hydro_out<-trace_flowpaths(
    input=hydro_out,
    return_products=return_products,
    temp_dir=temp_dir,
    verbose=verbose
  )

  if (verbose) print("Generating Pairwise Distances")
  hydro_out<-generate_pwisedist(
    input=hydro_out,
    return_products=return_products,
    temp_dir=temp_dir,
    verbose=verbose
  )

  return(hydro_out)
}
