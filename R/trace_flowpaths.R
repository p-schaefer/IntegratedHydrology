
#' Title
#'
#' @param input
#' @param return_products
#' @param temp_dir
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples

trace_flowpaths<-function(
    input,
    return_products=F,
    temp_dir=NULL,
    verbose=F
){
  require(sf)
  require(terra)
  require(whitebox)
  require(tidyverse)

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

  final_lines<-read_sf(file.path("/vsizip",zip_loc,"stream_lines.shp"))

  ds_flowpaths<-trace_ds_flowpath(input=final_lines,verbose=verbose)
  us_flowpaths<-trace_us_flowpath(input=final_lines,verbose=verbose)

  saveRDS(ds_flowpaths,file.path(temp_dir,"ds_flowpaths.rds"))
  saveRDS(us_flowpaths,file.path(temp_dir,"us_flowpaths.rds"))

  dist_list_out<-list(
    "ds_flowpaths.rds",
    "us_flowpaths.rds"
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
      list(ds_flowpaths=ds_flowpaths,
           us_flowpaths=us_flowpaths),
      output
    )
  }
  file.remove(list.files(temp_dir,full.names = T))

  return(output)

}

#' @export
trace_ds_flowpath<-function(
    input,
    verbose=F
) {
  require(tidyverse)

  input<-as_tibble(input)
  input$link_id<-as.character(input$link_id)

  unique_link_id<-unique(input$link_id)
  unique_link_id<-split(unique_link_id,unique_link_id)
  if (verbose) print("Generating Downstream flowpaths")

  pb <- txtProgressBar(min = 0,
                       max = length(unique_link_id),
                       style = 2,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,
                       char = "=")

  pb <- txtProgressBar(min = 0,
                       max = length(unique_link_id),
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,
                       char = "=")

  unique_link_id<-lapply(unique_link_id,function(x){
    out<-x
    repeat {
      ds_id <- input %>%
        filter(link_id %in% out) %>%
        select(starts_with("dslink_id")) %>%
        unlist() %>%
        as.character() %>%
        unique()

      if (all(ds_id %in% out)) break

      out<-unique(c(out,ds_id))
    }

    out<-out[!is.na(out)]
    out<-tibble(link_id=out,
                link_lngth=input$link_lngth[match(out,input$link_id)])

    if (verbose) setTxtProgressBar(pb, which(names(unique_link_id)==out$link_id[[1]]))

    setTxtProgressBar(pb, which(names(unique_link_id)==out$link_id[[1]]))

    return(out)
  })

  if (verbose) print("")

  return(unique_link_id)

}

#' @export
trace_us_flowpath<-function(
    input,
    verbose=F
) {
  require(tidyverse)

  input<-as_tibble(input)
  input$link_id<-as.character(input$link_id)

  unique_link_id<-unique(input$link_id)
  unique_link_id<-split(unique_link_id,unique_link_id)
  if (verbose) print("Generating Upstream flowpaths")

  pb <- txtProgressBar(min = 0,
                       max = length(unique_link_id),
                       style = 2,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,
                       char = "=")

  pb <- txtProgressBar(min = 0,
                       max = length(unique_link_id),
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,
                       char = "=")

  unique_link_id<-lapply(unique_link_id,function(x){
    out<-x
    repeat {
      ds_id <- input %>%
        filter(link_id %in% out) %>%
        select(starts_with("uslink_id")) %>%
        unlist() %>%
        as.character() %>%
        unique()

      if (all(ds_id %in% out)) break

      out<-unique(c(out,ds_id))
    }

    out<-out[!is.na(out)]
    out<-tibble(link_id=out,
                link_lngth=input$link_lngth[match(out,input$link_id)],
                sbbsn_area=input$sbbsn_area[match(out,input$link_id)])

    if (verbose) setTxtProgressBar(pb, which(names(unique_link_id)==out$link_id[[1]]))

    setTxtProgressBar(pb, which(names(unique_link_id)==out$link_id[[1]]))

    return(out)
  })

  if (verbose) print("")

  return(unique_link_id)
}

