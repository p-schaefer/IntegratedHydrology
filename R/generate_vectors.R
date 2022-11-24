#' Generate polygon subbasins, and attributed stream lines and points
#'
#' @param input resulting object from `generate_subbasins()`
#' @param extra_attr character. One or more of c("link_slope", "cont_slope", "USChnLn_To", "Elevation", "StOrd_Hack", "StOrd_Str", "StOrd_Hort", "StOrd_Shr"). Optional attributes to add to stream vector outputs.
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param points character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/points.shp"), or any GIS data object that will be converted to spatial points. Points representing sampling locations.
#' @param site_id_col character. Variable name in `points` that corresponds to unique site identifiers. This column will be included in all vector geospatial analysis products. Note, if multiple points have the same `site_id_col`, their centroid will be used and returned; if multiple points overlap after snapping, only the first is used.
#' @param snap_distance integer. Maximum distance which points will be snapped to stream lines in map units
#' @param break_on_noSnap logical. Should the function stop if any points are not snapped to a stream segment (i.e., are beyon snap_distance)
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param compress logical. Should output rasters be compressed, slower but more space efficient.
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'

generate_vectors<-function(
    input,
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
    site_id_col=NULL,
    snap_distance=NULL,
    break_on_noSnap=T,
    return_products=F,
    temp_dir=NULL,
    compress=F,
    verbose=F
) {
  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  options(future.rng.onMisuse="ignore")
  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")
  if (!is.logical(compress)) stop("'compress' must be logical")

  if (!is.null(points)) {
    if (!is.character(site_id_col)) stop("'site_id_col' must be a single character")
    if (length(site_id_col)>1 | length(site_id_col)==0) stop("length 'site_id_col' must be 1")
    if (site_id_col=="link_id") stop("'site_id_col' cannot be 'link_id'")
  }

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  extra_attr<-match.arg(extra_attr,several.ok = T)

  hydro_out<-input

  if (verbose) message("Attribute Streamlines")
  hydro_out<-attrib_streamline(
    input=hydro_out,
    extra_attr=extra_attr,
    return_products=return_products,
    points=points,
    site_id_col=site_id_col,
    snap_distance=snap_distance,
    break_on_noSnap=break_on_noSnap,
    temp_dir=temp_dir,
    compress=compress,
    verbose=verbose
  )

  if (verbose) message("Generate Subbasins")
  hydro_out<-generate_subbasins(
    input=hydro_out,
    points=points,
    #site_id_col=site_id_col,
    return_products=return_products,
    temp_dir=temp_dir,
    verbose=verbose
  )

  class(hydro_out)<-"ihydro"
  return(hydro_out)
}
