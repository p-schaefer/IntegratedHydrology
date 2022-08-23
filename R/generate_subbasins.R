
#' Extracts polygon subbasins from 'process_flowdir()'
#'
#' @param input resulting object from `process_flowdir()`
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param verbose logical.
#'
#' @seealso [whitebox::wbt_subbasins]
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'
#' @examples
#'

generate_subbasins<-function(
    input,
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

  wbt_options(exe_path=wbt_exe_path(),
              verbose=verbose,
              wd=temp_dir)

  terra::terraOptions(verbose = verbose,
                      tempdir = temp_dir
  )

  zip_loc<-input$outfile
  fl<-unzip(list=T,zip_loc)

  unzip(zip_loc,
        c("dem_d8.tif","dem_streams_d8.tif"),
        exdir=temp_dir,
        overwrite=T,
        junkpaths=T)

  # Generate subbasin polygons ----------------------------------------------
  if (verbose) print("Generating subbasins")
  wbt_subbasins(
    d8_pntr="dem_d8.tif",
    streams="dem_streams_d8.tif",
    output="Subbasins.tif"
  )

  if (verbose) print("Converting subbasins to polygons")
  subb<-rast(file.path(temp_dir,"Subbasins.tif"))
  subb<-as.polygons(subb,dissolve = TRUE)
  subb<-st_as_sf(subb) %>%
    mutate(sbbsn_area=st_area(.))
  names(subb)[1]<-"link_id"

  write_sf(subb,file.path(temp_dir,"Subbasins_poly.shp"))

  # Prepare Output ----------------------------------------------------------

  dist_list_out<-list(
    list.files(temp_dir,"Subbasins_poly")
  )

  dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

  out_file<-zip_loc

  if (verbose) print("Generating Output")

  zip(out_file,
      unlist(dist_list_out),
      flags = '-r9Xjq'
  )

  output<-input

  if (return_products){
    output<-c(
      list(subbasins=subb),
      output
    )
  }
  file.remove(list.files(temp_dir,full.names = T,recursive=T))

  return(output)
}
