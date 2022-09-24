
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
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")

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

  wbt_options(exe_path=wbt_exe_path(),
              verbose=verbose,
              wd=temp_dir)

  terra::terraOptions(verbose = verbose,
                      tempdir = temp_dir
  )

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

  num_inputs<-map(num_inputs,function(x){
    if (inherits(x,"character")) return(x)
    if (inherits(x,c("sf","sfc","sfg"))) return(wrap(vect(x)))
    if (inherits(x,c("SpatRaster","SpatVector"))) return(wrap(x))
    return(x)
  })

  cat_inputs<-map(cat_inputs,function(x){
    if (inherits(x,"character")) return(x)
    if (inherits(x,c("sf","sfc","sfg"))) return(wrap(vect(x)))
    if (inherits(x,c("SpatRaster","SpatVector"))) return(wrap(x))
    return(x)
  })


  inputs_list<-list(
    num_rast.tif=list(lyr_nms=as.list(names(num_inputs)),
                      lyr=num_inputs,
                      lyr_variables=variable_names[names(num_inputs)]),
    cat_rast.tif=list(lyr_nms=names(cat_inputs),
                      lyr=cat_inputs,
                      lyr_variables=variable_names[names(cat_inputs)])
  )

  inputs_list<-map2(inputs_list,names(inputs_list),~c(.x,
                                                      list(rln=as.list(rep(.y,length(.x[[1]])))),
                                                      list(temp_dir=as.list(rep(temp_dir,length(.x[[1]]))))
  ))

  loi_fn<-function(lyr_nms,lyr,lyr_variables,gdal_arg,p,rln,temp_dir){
    pmap(list(lyr_nms=lyr_nms,
              lyr=lyr,
              lyr_variables=lyr_variables,
              rln=rln,
              temp_dir=temp_dir),
         function(lyr_nms,
                  lyr,
                  lyr_variables,
                  rln,temp_dir){
           #print(lyr)
           #browser()

           temp_temp_dir<-file.path(temp_dir,basename(tempfile()))
           dir.create(temp_temp_dir)

           terra::terraOptions(verbose = verbose,
                               tempdir = temp_temp_dir
           )

           resaml<-ifelse(grepl("num_rast",rln),"bilinear","near")

           output<-hydroweight::process_input(
             input=unlist(lyr),
             input_name = unlist(lyr_nms),
             variable_name=unlist(lyr_variables),
             target=file.path(temp_dir,"dem_final.tif"),
             clip_region = file.path(temp_dir,"clip_region.shp"),
             resample_type = resaml,
             working_dir=temp_temp_dir
           )

           out_file<-file.path(temp_dir,rln)

           if (file.exists(out_file)) {
             output<-rast(list(rast(out_file),output))
             out_file<-file.path(temp_dir,paste0("new_",rln))
           }

           ot<-writeRaster(output,out_file,overwrite=T,gdal="COMPRESS=NONE")

           suppressWarnings(terra::tmpFiles(current = T,orphan=T,old=T,remove = T))
           file.remove(list.files(temp_temp_dir,full.names = T,recursive = T))
           unlink(temp_temp_dir,recursive = T, force = T)

           if (file.exists(file.path(temp_dir,paste0("new_",rln)))){
             file.remove(file.path(temp_dir,paste0(rln)))
             file.rename(file.path(temp_dir,paste0("new_",rln)),
                         file.path(temp_dir,paste0(rln)))
           }

           p()

           return(NULL)
         })
  }

  inputs_list<-inputs_list[sapply(inputs_list,function(x) length(x$lyr))>0]

  with_progress(enable=T,{
    print("Processing loi")
    p <- progressor(steps = (length(num_inputs) + length(cat_inputs)))

    ot<-future_map(inputs_list,~loi_fn(lyr_nms=.$lyr_nms,
                                       lyr=.$lyr,
                                       lyr_variables=.$lyr_variables,
                                       gdal_arg=.$gdal_arg,
                                       p=p,
                                       rln=.$rln,
                                       temp_dir=.$temp_dir))
  })

  # Generate Output ---------------------------------------------------------
  if (verbose) print("Generating Outputs")

  dist_list_out<-list(
    "num_rast.tif",
    "cat_rast.tif"
  )

  dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

  dist_list_out<-dist_list_out[sapply(dist_list_out,file.exists)]

  zip(output_filename,
      unlist(dist_list_out),
      flags = '-r9Xjq'
  )

  output<-input

  if (return_products){

    ot<-map(dist_list_out,~wrap(rast(.)))
    names(ot)[grepl("num_rast",unlist(dist_list_out))]<-"num_inputs"
    names(ot)[grepl("cat_rast",unlist(dist_list_out))]<-"cat_inputs"

    output<-c(
      ot,
      output
    )
  }

  file.remove(list.files(temp_dir,full.names = T,recursive = T))

  return(output)

}
