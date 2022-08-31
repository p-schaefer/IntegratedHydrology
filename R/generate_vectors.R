#' Generate polygon subbasins, and attributed stream lines and points
#'
#' @param input resulting object from `generate_subbasins()`
#' @param extra_attr character. One or more of c("link_slope", "cont_slope", "USChnLn_To", "Elevation", "StOrd_Hack", "StOrd_Str", "StOrd_Hort", "StOrd_Shr"). Optional attributes to add to stream vector outputs.
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
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
    return_products=F,
    temp_dir=NULL,
    verbose=F
) {

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")
  extra_attr<-match.arg(extra_attr,several.ok = T)

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  options(dplyr.summarise.inform = FALSE)

  wbt_options(exe_path=wbt_exe_path(),
              verbose=verbose,
              wd=temp_dir)

  terra::terraOptions(verbose = verbose,
                      tempdir = temp_dir
  )

  options(dplyr.summarise.inform = FALSE)

  extra_attr<-match.arg(extra_attr,several.ok = T)

  if (verbose) print("Generate Subbasins")
  hydro_out<-generate_subbasins(
    input=input,
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

  return(hydro_out)
}
