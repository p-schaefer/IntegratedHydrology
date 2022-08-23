
#' Title
#'
#' @param input ouptut of `trace_flowpaths()`
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'
#' @examples

generate_pwisedist<-function(
    input,
    return_products=F,
    temp_dir=NULL,
    verbose=F
){
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
        c("ds_flowpaths.rds","us_flowpaths.rds"),
        exdir=temp_dir,
        overwrite=T,
        junkpaths=T)

  ds_flowpaths<-readRDS(file.path(temp_dir,"ds_flowpaths.rds"))
  us_flowpaths<-readRDS(file.path(temp_dir,"us_flowpaths.rds"))

  pwise_dist<-pairwise_dist_fn(ds_flowpaths=ds_flowpaths,
                               us_flowpaths=us_flowpaths,
                               verbose=verbose)

  saveRDS(pwise_dist,file.path(temp_dir,"pwise_dist.rds"))

  dist_list_out<-list(
    "pwise_dist.rds"
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
      list(pwise_dist=pwise_dist),
      output
    )
  }
  file.remove(list.files(temp_dir,full.names = T))

  return(output)
}

#' @export
pairwise_dist_fn<-function(
    ds_flowpaths=NULL,
    us_flowpaths=NULL,
    stream_lines=NULL,
    verbose=F
) {
  # require(tidyverse)

  if (is.null(stream_lines) & is.null(ds_flowpaths)) stop("Either 'ds_flowpaths' or 'stream_lines' must be provided")
  if (is.null(stream_lines) & is.null(us_flowpaths)) stop("Either 'us_flowpaths' or 'stream_lines' must be provided")

  if (is.null(ds_flowpaths)) ds_flowpaths<-trace_ds_flowpath(stream_lines)
  if (is.null(us_flowpaths)) us_flowpaths<-trace_us_flowpath(stream_lines)
  us_catchment_areas<-lapply(us_flowpaths,function(x) sum(x$sbbsn_area,na.rm=T))

  if (any(duplicated(names(ds_flowpaths)))) stop("'ds_flowpaths' cannot contain duplicate names")
  if (any(duplicated(names(us_flowpaths)))) stop("'us_flowpaths' cannot contain duplicate names")
  if (any(!names(ds_flowpaths) %in% names(us_flowpaths))) stop("'ds_flowpaths' contains link_ds not present in 'us_flowpaths'")
  if (any(!names(us_flowpaths) %in% names(ds_flowpaths))) stop("'us_flowpaths' contains link_ds not present in 'ds_flowpaths'")

  if (verbose) print("Generating Pairwise Distances")

  pl_fn<-function(.x,.y,p){
    out<-mutate(.x,is_target=link_id==.y) %>%
      filter(!is.na(link_lngth)) %>%
      mutate(link_lngth=cumsum(link_lngth)) %>%
      filter(is_target) %>%
      pull(link_lngth) %>%
      ifelse(length(.)==0,0,.)

    p()

    return(out)
  }

  with_progress({
    print("Generating Pairwise Distances")
    p <- progressor(steps = length(ds_flowpaths)*length(ds_flowpaths))

    out_tbl<-tibble(
      origin=rep(names(ds_flowpaths),each=length(ds_flowpaths)),
      destination=rep(names(ds_flowpaths),length.out=length(ds_flowpaths)*length(ds_flowpaths))) %>%
      mutate(
        directed_path_length=ds_flowpaths[origin]
      ) %>%
      mutate(directed_path_length=future_map2_dbl(directed_path_length,destination,pl_fn,p)) %>%
      mutate(
        origin_catchment=unlist(us_catchment_areas[origin]),
        destination_catchment=unlist(us_catchment_areas[destination])
      ) %>%
      mutate(prop_shared_catchment=case_when(
        directed_path_length>0 ~ origin_catchment/destination_catchment,
        T ~ 0
      )) %>%
      select(-origin_catchment,-destination_catchment)
  })

  out_prop<-out_tbl %>%
    select(-directed_path_length) %>%
    mutate(origin=paste0("prop_link_id_",origin)) %>%
    rename(link_id=destination) %>%
    pivot_wider(names_from=origin,values_from=prop_shared_catchment)

  out_dist<-out_tbl %>%
    select(-prop_shared_catchment) %>%
    mutate(origin=paste0("dist_link_id_",origin)) %>%
    rename(link_id=destination) %>%
    pivot_wider(names_from=origin,values_from=directed_path_length)

  return(list(
    long_pwise=out_tbl,
    wide_prop_shared_catchment=out_prop,
    wide_directed_path_length=out_dist
  ))

}

