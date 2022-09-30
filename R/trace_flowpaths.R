
#' Generate upstream and downstream flow direction lists originating from each subbasin (and sampling point)
#'
#' @param input resulting object from `attrib_streamline()` or `insert_points()`
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export

trace_flowpaths<-function(
    input,
    return_products=F,
    temp_dir=NULL,
    verbose=F
){
  options(future.rng.onMisuse="ignore")
  options(scipen = 999)

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  options(dplyr.summarise.inform = FALSE)

  zip_loc<-input$outfile

  site_id_col<-paste0(data.table::fread(cmd=paste("unzip -p ",zip_loc,"site_id_col.csv")))

  final_links<-as_tibble(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv"))) %>%
    mutate(across(any_of(site_id_col),na_if,""))

  #browser()

  fp<-trace_flowpath_fn(input=final_links,verbose=verbose)
  ds_flowpaths<-fp$final_out_ds
  us_flowpaths<-fp$final_out_us

  saveRDS(ds_flowpaths,file.path(temp_dir,"ds_flowpaths.rds"),compress = F)
  saveRDS(us_flowpaths,file.path(temp_dir,"us_flowpaths.rds"),compress = F)

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

  output<-input[!names(input) %in% c("catchment_poly")]

  all_catch<-get_catchment(input = output,
                           target_points = final_links[["link_id"]]) %>%
    select(link_id)

  write_sf(all_catch %>% select(link_id),file.path(temp_dir,"Catchment_poly.shp"))

  zip(out_file,
      list.files(temp_dir,"Catchment_poly",full.names = T),
      flags = '-r9Xjq'
  )

  if (return_products){
    output<-c(
      list(ds_flowpaths=ds_flowpaths,
           us_flowpaths=us_flowpaths,
           catchments=all_catch %>% select(link_id)
      ),
      output
    )
  }
  suppressWarnings(file.remove(list.files(temp_dir,full.names = T)))

  return(output)

}


#' @export
trace_flowpath_fn<-function(
    input,
    verbose=F
) {
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")

  input_tib<-input %>%
    as_tibble() %>%
    select(link_id,trib_id,link_lngth,sbbsn_area,USChnLn_Fr,
           starts_with("uslink_id"),starts_with("dslink_id"),starts_with("dstrib_id"),starts_with("ustrib_id")) %>%
    filter(!is.na(link_id)) %>%
    mutate(link_id=as.character(link_id)) %>%
    pivot_longer(cols=c(starts_with("uslink_id"),starts_with("ustrib_id"))) %>%
    filter(case_when(
      name %in% c("uslink_id1","ustrib_id1") ~ T,
      is.na(value) ~ F,
      T ~ T
    )) %>%
    mutate(name=ifelse(grepl("trib",name),"ustrib_id1","uslink_id1")) %>%
    pivot_wider(names_from=name,values_from=value,values_fn=list) %>%
    unnest(cols=c(uslink_id1,ustrib_id1)) %>%
    distinct()

  # Downstream paths
  if (verbose) print("Identifying Exit Tributaries")

  int_tribs<-input_tib %>%
    mutate(trib_id2=trib_id) %>%
    group_by(trib_id) %>%
    nest() %>%
    ungroup() %>%
    #mutate(exit_trib=future_map_lgl(data,~any(is.na(.$dslink_id1)))) %>%
    mutate(join_id=future_map(data,
                              .options = furrr_options(globals = FALSE),
                              carrier::crate(function(x){
                                options(scipen = 999)
                                `%>%` <- magrittr::`%>%`

                                dplyr::arrange(x,USChnLn_Fr) %>%
                                  utils::tail(1) %>%
                                  dplyr::select(dslink_id=dslink_id1,dstrib_id=dstrib_id1)
                              } ))) %>%
    unnest(join_id) %>%
    mutate(join_USChnLn_Fr=input_tib$USChnLn_Fr[match(.$dslink_id,input_tib$link_id)])

  # Exiting tribs
  final_ds_paths<-int_tribs %>%
    filter(is.na(dslink_id)) %>%
    select(trib_id,data)

  int_tribs<-int_tribs %>%
    filter(!trib_id %in% final_ds_paths$trib_id)

  if (verbose) print("Tracing Downstream Flowpaths")
  while(nrow(int_tribs)>0){
    int_tribs_int<-int_tribs %>%
      filter(dstrib_id %in% final_ds_paths$trib_id) %>%
      mutate(final_paths=list(final_ds_paths)) %>%
      mutate(data=future_pmap(list(
        final_paths=final_paths,
        data=data,
        dstrib_id=dstrib_id,
        join_USChnLn_Fr=join_USChnLn_Fr),
        .options = furrr_options(globals = FALSE),
        carrier::crate(
          function(data,dstrib_id,join_USChnLn_Fr,final_paths) {
            options(scipen = 999)
            `%>%` <- magrittr::`%>%`

            dplyr::bind_rows(data,
                             final_paths$data[final_paths$trib_id==dstrib_id][[1]] %>% dplyr::arrange(USChnLn_Fr) %>% dplyr::filter(USChnLn_Fr>=join_USChnLn_Fr))
          }
        )
      ))

    int_tribs<-int_tribs %>% filter(!trib_id %in% int_tribs_int$trib_id)
    final_ds_paths<-bind_rows(final_ds_paths,int_tribs_int %>% select(trib_id,data))
  }

  #Finalize DS flowpaths
  if (verbose) print("Finalizing Downstream Flowpaths")
  final_ds_paths<-final_ds_paths %>%
    mutate(data2=future_map(data,
                            .options = furrr_options(globals = FALSE),
                            carrier::crate(function(x){
                              options(scipen = 999)
                              `%>%` <- magrittr::`%>%`

                              x<- x%>%
                                dplyr::select(link_id,trib_id=trib_id2,link_lngth,sbbsn_area,USChnLn_Fr) %>%
                                dplyr::distinct()
                              u<-x$link_id
                              u<-stats::setNames(u,u)
                              purrr:::map(u, function (y) x[which(x$link_id==y):nrow(x),] %>% dplyr::mutate(origin_id=y))
                            })))

  final_ds_paths<-final_ds_paths %>%
    select(data2) %>%
    unnest(cols=data2) %>%
    unnest(cols=data2) %>%
    distinct() %>%
    select(origin_id,everything())

  #US Paths

  if (verbose) print("Identifying Headwater Tributaries")
  final_us_paths<-final_ds_paths %>%
    group_by(link_id) %>%
    nest() %>%
    ungroup() %>%
    mutate(us_links=future_map(data,
                               .options = furrr_options(globals = FALSE),
                               carrier::crate(function(x){
                                 options(scipen = 999)
                                 `%>%` <- magrittr::`%>%`

                                 x %>% dplyr::pull(origin_id) %>% unique()
                               })
    )) %>%
    mutate(full_data=list(input[,c("link_id","trib_id","link_lngth","sbbsn_area","USChnLn_Fr")]))

  if (verbose) print("Finalizing Upstream Flowpaths")
  final_us_paths<-final_us_paths%>%
    mutate(data2=future_map2(us_links,full_data,
                             .options = furrr_options(globals = FALSE),
                             carrier::crate(function(x,y){
                               y[y$link_id %in% x,c("link_id","trib_id","link_lngth","sbbsn_area","USChnLn_Fr")]
                             })
    )) %>%
    select(source_ID=link_id, data2) %>%
    unnest(cols=data2)

  final_out_ds<-split(final_ds_paths[,-c(1)],final_ds_paths$origin_id)

  final_out_us<-split(final_us_paths[,-c(1)],final_us_paths$source_ID)

  return(list(final_out_us=final_out_us,final_out_ds=final_out_ds))

}

