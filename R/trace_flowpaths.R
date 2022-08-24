
#' Generate upstream and downstream flow direction lists originating from each subbasin (and sampling point)
#'
#' @param input resulting object from `attrib_streamline()` or `insert_points()`
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'
#' @examples

trace_flowpaths<-function(
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

  final_links<-read_sf(file.path("/vsizip",zip_loc,"stream_links.shp"))

  ds_flowpaths<-trace_ds_flowpath(input=final_links,verbose=verbose)
  us_flowpaths<-trace_us_flowpath(input=final_links,verbose=verbose)

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

  input_tib<-input %>%
    as_tibble() %>%
    select(link_id,trib_id,link_lngth,sbbsn_area,starts_with("uslink_id"),starts_with("dslink_id")) %>%
    filter(!is.na(link_id)) %>%
    mutate(link_id=as.character(link_id))

  unique_link_id<-input_tib %>%
    filter(!is.na(link_id)) %>%
    filter(if_all(starts_with("uslink_id"),is.na)) %>%
    pull(link_id) %>%
    unique()

  unique_link_id<-split(unique_link_id,unique_link_id)
  print("Generating Downstream flowpaths")

  id_fn<-function(x,input_tib,p){
    out<-x
    #print(x)
    repeat {
      ds_id <- input_tib %>%
        filter(link_id %in% out) %>%
        dplyr::select(starts_with("dslink_id")) %>%
        unlist() %>%
        as.character() %>%
        unique() %>%
        .[!is.na(.)]

      if (all(ds_id %in% out)) break

      out<-unique(c(out,ds_id))
    }

    out<-out[!is.na(out)]
    out<-tibble(link_id=out) %>%
      left_join(input_tib %>% dplyr::select(link_id,trib_id,link_lngth,sbbsn_area),by="link_id") %>%
      distinct() %>%
      filter(!is.na(link_id))

    p()
    return(out)
  }

  input_tib_list<-as.list(rep(list(input_tib),length(unique_link_id)))

  with_progress({
    p <- progressor(steps = length(unique_link_id))

    unique_link_id<-future_map2(unique_link_id,input_tib_list,~id_fn(x=.x,input_tib=.y,p=p))
  })

  # Get remaining distances by unnesting the list
  final_out<-unique_link_id

  for (i in final_out){
    nrw<-nrow(i)
    new_entries<-lapply(2:nrw,function(x) i[seq(x,nrw),])
    names(new_entries)<-sapply(new_entries,function(x) head(x$link_id,1))

    keep_entries<-new_entries[!names(new_entries) %in% names(final_out)]

    final_out<-c(final_out,keep_entries)
  }

  final_out<-final_out[order(names(final_out))]
  final_out<-final_out[!is.na(names(final_out))]

  return(final_out)

}

#' @export
trace_us_flowpath<-function(
    input,
    verbose=F
) {

  input_tib<-input %>%
    as_tibble() %>%
    select(link_id,trib_id,link_lngth,sbbsn_area,
           starts_with("uslink_id"),starts_with("dslink_id"),starts_with("dstrib_id"),starts_with("ustrib_id")) %>%
    filter(!is.na(link_id)) %>%
    mutate(link_id=as.character(link_id))

  unique_link_id<-input_tib %>%
    # filter(if_all(starts_with("dslink_id"),is.na) | trib_id != dstrib_id1) %>% # heep any without a downstream link, or where the downstream trib has a different ID
    filter(!is.na(link_id)) %>%
    pull(link_id) %>%
    unique()

  unique_link_id<-split(unique_link_id,unique_link_id)
  print("Generating Upstream flowpaths")

  id_fn<-function(x,input_tib,p){
    out<-x
    repeat {
      us_id <- input_tib %>%
        filter(link_id %in% out) %>%
        dplyr::select(starts_with("uslink_id")) %>%
        unlist() %>%
        as.character() %>%
        unique() %>%
        .[!is.na(.)]

      if (all(us_id %in% out)) break

      out<-unique(c(out,us_id))
    }

    out<-out[!is.na(out)]
    out<-tibble(link_id=out)%>%
      left_join(input_tib %>% dplyr::select(link_id,trib_id,link_lngth,sbbsn_area),by="link_id") %>%
      distinct() %>%
      filter(!is.na(link_id))

    p()
    return(out)
  }

  #browser()

  input_tib_list<-as.list(rep(list(input_tib),length(unique_link_id)))

  with_progress({
    p <- progressor(steps = length(unique_link_id))

    unique_link_id<-future_map2(unique_link_id,input_tib_list,~id_fn(x=.x,input_tib=.y,p=p))
  })

  final_out<-unique_link_id

  # for (i in final_out){ # there is probably some way of doing this similar to ds_tracing but I can't figure it out right now
  #   nrw<-nrow(i)
  #   new_entries<-lapply(2:nrw,function(x) i[seq(x,nrw),])
  #   names(new_entries)<-sapply(new_entries,function(x) head(x$link_id,1))
  #
  #   keep_entries<-new_entries[!names(new_entries) %in% names(final_out)]
  #
  #   final_out<-c(final_out,keep_entries)
  # }

  final_out<-final_out[order(names(final_out))]
  final_out<-final_out[!is.na(names(final_out))]

  return(final_out)
}

