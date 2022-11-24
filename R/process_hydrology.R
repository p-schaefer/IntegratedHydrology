#' Generate geospatial analysis products
#'
#' `ihydro::process_hydrology` processes a DEM into various geospatial analysis products.
#'
#' @param dem character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/dem.tif"), \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of GeoTiFF type. Digital elevation model raster.
#' @param output_filename character. File path to write resulting .zip file.
#' @param threshold integer. Flow accumulation threshold for stream initiation.
#' @param burn_streams character. (full file path with extension, e.g., "C:/Users/Administrator/Desktop/input.shp"), sf, SpatVector, PackedSpatVector, RasterLayer, SpatRaster, or PackedSpatRaster. Stream vector to burn into DEM.
#' @param burn_depth numeric. Depth (in meters) to burn stream into the DEM.
#' @param min_length numeric. Minimum tributary length, shorter 1st order tributaries are removed.
#' @param depression_corr NULL or character. One of c("fill","breach"), specifying whether depressions should be filled or breached. NULL will perform neither, if DEM is already corrected.
#' @param extra_attr character. One or more of c("link_slope", "cont_slope", "USChnLn_To", "Elevation", "StOrd_Hack", "StOrd_Str", "StOrd_Hort", "StOrd_Shr"). Optional attributes to add to stream vector outputs.
#' @param points character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/points.shp"), or any GIS data object that will be converted to spatial points. Points representing sampling locations.
#' @param snap_distance integer. Maximum distance which points will be snapped to stream lines in map units
#' @param site_id_col character. Variable name in `points` that corresponds to unique site identifiers. This column will be included in all vector geospatial analysis products. Note, if multiple points have the same `site_id_col`, their centroid will be used and returned; if multiple points overlap after snapping, only the first is used.
#' @param pwise_dist logical. Calculate pairwise distances.
#' @param pwise_all_links logical. Should all pairwise distances be calculate, or only those originating from sampling points
#' @param return_products logical. If \code{TRUE}, a list containing all geospatial analysis products. If \code{FALSE}, folder path to resulting .zip file.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param compress logical. Should output rasters be compressed, slower but more space efficient.
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'

#' @importFrom hydroweight process_input

process_hydrology<-function(
    dem,
    output_filename,
    threshold,
    burn_streams=NULL,
    burn_depth=5,
    min_length=NULL,
    depression_corr=c(NULL,"fill","breach"),
    #calc_catch=c("all","none","sample_points"),
    extra_attr=c(
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
    break_on_noSnap=T,
    site_id_col=NULL,
    pwise_dist=F,
    pwise_all_links=F,
    return_products=F,
    temp_dir=NULL,
    compress=F,
    verbose=F
) {
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)

  # calc_catch<-calc_catch[1]
  # match.arg(calc_catch,several.ok = F)

  if (!is.integer(threshold)) stop("'threshold' must be an integer value")
  if (!is.null(snap_distance) && !is.integer(snap_distance)) stop("'snap_distance' must be an integer value")
  if (!is.numeric(burn_depth)) stop("'burn_depth' must be an numeric value")
  burn_depth<-as.integer(burn_depth)

  if (!is.logical(pwise_dist)) stop("'pwise_dist' must be logical")
  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(compress)) stop("'compress' must be logical")
  if (!is.logical(break_on_noSnap)) stop("'break_on_noSnap' must be logical")
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
  if (!grepl("\\.gpkg$",output_filename)) stop("output_filename must be a character ending in '.gpkg'")

  extra_attr<-match.arg(extra_attr,several.ok = T)

  if (verbose) message("Processing Flow Direction")
  hydro_out<-process_flowdir(
    dem=dem,
    threshold=threshold,
    burn_streams=burn_streams,
    burn_depth=burn_depth,
    min_length=min_length,
    depression_corr=depression_corr,
    return_products=return_products,
    output_filename=output_filename,
    temp_dir=temp_dir,
    compress=compress,
    verbose=verbose
  )

  if (verbose) message("Generate Vector")
  hydro_out<-generate_vectors(
    input=hydro_out,
    extra_attr=extra_attr,
    points=points,
    site_id_col=site_id_col,
    snap_distance=snap_distance,
    break_on_noSnap=break_on_noSnap,
    return_products=return_products,
    temp_dir=temp_dir,
    compress=compress,
    verbose=verbose
  )

  if (verbose) message("Tracing Flowpaths")
  hydro_out<-trace_flowpaths(
    input=hydro_out,
    return_products=return_products,
    temp_dir=temp_dir,
    pwise_dist=pwise_dist,
    pwise_all_links=pwise_all_links,
    verbose=verbose
  )

  return(hydro_out)
}
