
#' Generate pairwise distances between flow-connected and flow-unconnected site pairs, as well proportions of shared catchments
#'
#' @param input ouptut of `trace_flowpaths()`
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'


generate_pwisedist<-function(
    input,
    return_products=F,
    temp_dir=NULL,
    chunck_size=250,
    verbose=F
){
  stop("Function is depreciated, pairwise distances now calculated with 'trace_flowpaths()'")
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

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

  zip_loc<-input$outfile
  fl<-unzip(list=T,zip_loc)

  unzip(zip_loc,
        c("flowpaths_out.db"),
        exdir=temp_dir,
        overwrite=T,
        junkpaths=T)
  db_fp<-file.path(temp_dir,"flowpaths_out.db")

  # ds_flowpaths<-readRDS(file.path(temp_dir,"ds_flowpaths.rds"))
  # us_flowpaths<-readRDS(file.path(temp_dir,"us_flowpaths.rds"))

  pwise_dist<-pairwise_dist_fn(db_fp=db_fp,
                               chunck_size=chunck_size,
                               #ds_flowpaths_file=file.path(temp_dir,"ds_flowpaths.rds"),
                               #us_flowpaths_file=file.path(temp_dir,"us_flowpaths.rds"),
                               verbose=verbose,
                               temp_dir=temp_dir)

  # saveRDS(pwise_dist,file.path(temp_dir,"pwise_dist.rds"))

  dist_list_out<-list(
    db_fp
  )
  #browser()

  #dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

  out_file<-zip_loc

  if (verbose) print("Generating Output")

  zip(out_file,
      unlist(dist_list_out),
      flags = '-r9Xjq'
  )
  # dist_list_out<-list(
  #   "pwise_dist.rds"
  # )
  #
  # dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))
  #
  # out_file<-zip_loc
  #
  # if (verbose) print("Generating Output")
  #
  # zip(out_file,
  #     unlist(dist_list_out),
  #     flags = '-r9Xjq'
  # )

  output<-input

  if (return_products){
    output<-c(
      list(pwise_dist=pwise_dist),
      output
    )
  }
  suppressWarnings(file.remove(list.files(temp_dir,full.names = T)))

  return(output)
}

#' @export
pairwise_dist_fn<-function(
    db_fp,
    chunck_size,
    #ds_flowpaths_file=NULL,
    #us_flowpaths_file=NULL,
    #stream_links=NULL,
    verbose=F,
    temp_dir=NULL
) {
  library(dbplyr)
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")

  #browser()

  # if (is.null(stream_links) & is.null(ds_flowpaths_file)) stop("Either 'ds_flowpaths' or 'stream_links' must be provided")
  # if (is.null(stream_links) & is.null(us_flowpaths_file)) stop("Either 'us_flowpaths' or 'stream_links' must be provided")
  #
  # if (is.null(ds_flowpaths_file)) {
  #   ds_flowpaths<-trace_ds_flowpath(stream_links)
  # } else {
  #   ds_flowpaths<-readRDS(ds_flowpaths_file)
  # }
  # if (is.null(us_flowpaths_file)) {
  #   us_flowpaths<-trace_us_flowpath(stream_links)
  # } else {
  #   us_flowpaths<-readRDS(us_flowpaths_file)
  # }

  # if (any(duplicated(names(ds_flowpaths)))) stop("'ds_flowpaths' cannot contain duplicate names")
  # if (any(duplicated(names(us_flowpaths)))) stop("'us_flowpaths' cannot contain duplicate names")
  # if (any(!names(ds_flowpaths) %in% names(us_flowpaths))) stop("'ds_flowpaths' contains link_ds not present in 'us_flowpaths'")
  #if (any(!names(us_flowpaths) %in% names(ds_flowpaths))) stop("'us_flowpaths' contains link_ds not present in 'ds_flowpaths'") # this can happen at raster boundaries

  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)

  link_id_all<-unique(DBI::dbGetQuery(con, "SELECT DISTINCT link_id FROM ds_flowpaths")$link_id)
  us_catchment_areas<-DBI::dbGetQuery(con, paste0("SELECT source_id, SUM(sbbsn_area) AS sbbsn_area FROM us_flowpaths GROUP BY source_id"))

  chunks<-split(link_id_all,ceiling(seq_along(link_id_all)/chunck_size))

  # copy_to(con,
  #         tibble(origin =NA,
  #                destination =NA,
  #                undirected_path_length =NA,
  #                ),
  #         name="pairwise_pathlengths",overwrite = T)

  origin_catchment_tbl<-tbl(con,"us_flowpaths") %>%
    group_by(source_id) %>%
    summarize(origin_catchment=sum(sbbsn_area)) %>%
    rename(origin=source_id)

  destination_catchment_tbl<-tbl(con,"us_flowpaths") %>%
    group_by(source_id) %>%
    summarize(destination_catchment=sum(sbbsn_area)) %>%
    rename(destination=source_id)

  # DS directed path-lengths
  flowconn_out<-tbl(con,"ds_flowpaths") %>%
    filter(origin_id %in% link_id_all) %>%
    mutate(origin=origin_id) %>%
    mutate(destination=link_id) %>%
    group_by(origin_id) %>%
    window_order(USChnLn_Fr) %>%
    mutate(directed_path_length=cumsum(link_lngth)) %>%
    ungroup() %>%
    select(origin,destination,directed_path_length) %>%
    distinct() %>%
    left_join(
      origin_catchment_tbl,by="origin"
    ) %>%
    left_join(
      destination_catchment_tbl,by="destination"
    ) %>%
    mutate(prop_shared_catchment=case_when(
      directed_path_length>0 ~ as.numeric(origin_catchment) / as.numeric(destination_catchment),
      T ~ 0
    )) %>%
    mutate(undirected_path_length=directed_path_length)

  copy_to(con,flowconn_out,"pairwise_pathlengths")

  # Reverse are not directly connected, but undrirectly are connected
  flowconn_out_rev<-flowconn_out %>%
    rename(origin = destination,
           destination = origin) %>%
    mutate(directed_path_length =0,
           prop_shared_catchment = 0) %>%
    collect()

  DBI::dbAppendTable(con, "pairwise_pathlengths", flowconn_out_rev)


  out<-map(chunks,function(link_id_all){
    # US un-directed path-lengths
    ot<-map(setNames(link_id_all,link_id_all),~filter(flowconn_out,origin == .))
    ot<-map2(ot,
             rep(list(flowconn_out),length(ot)),
             #.options = furrr::furrr_options(globals = FALSE),
             #carrier::crate(
               function(x,flowconn_out) {
                 library(dbplyr)
                 #options(scipen = 999)
                 # `%>%` <- magrittr::`%>%`
                 # select<-dplyr::select
                 # mutate<-dplyr::mutate
                 # left_join<-dplyr::left_join
                 # case_when<-dplyr::case_when
                 # filter<-dplyr::filter

                 oo<-dplyr::full_join(x %>% select(origin,destination,directed_path_length),
                                      flowconn_out %>% select(origin,destination,directed_path_length),
                                      by=c("destination"="origin"),
                                      suffix=c("_p1","_p2")) %>%
                   filter(!is.na(origin)) %>%
                   left_join(flowconn_out %>% select(origin,destination,directed_path_length_dir1=directed_path_length),
                             by=c("destination_p2"="origin",
                                  "destination_p1"="destination")) %>%
                   mutate(undirected_path_length=case_when(
                     is.na(directed_path_length_dir1) ~ as.numeric(directed_path_length_p1)+as.numeric(directed_path_length_p2),
                     T ~ NA_real_
                   )) %>%
                   select(origin,destination=destination_p2,undirected_path_length) %>%
                   mutate(directed_path_length=0,
                          prop_shared_catchment=0) %>%
                   left_join(
                     origin_catchment_tbl,by="origin"
                   ) %>%
                   left_join(
                     destination_catchment_tbl,by="destination"
                   ) %>%
                   collect()

                 return(oo)

                 # bind_rows(
                 #   oo,
                 #   oo %>% rename(origin=destination,
                 #                 destination=origin)
                 # )
               })
    #)

    #browser()
    ot<-bind_rows(ot)

    DBI::dbAppendTable(con, "pairwise_pathlengths", ot)


  })

  # left_join(
  #   origin_catchment_tbl,by="origin"
  # ) %>%
  #   left_join(
  #     destination_catchment_tbl,by="destination"
  #   )
  #














  # o2<-tibble(
  #   origin=rep(ot,each=length(ot)),
  #   destination=rep(ot,length.out=length(ot)*length(ot))
  # )
  #
  #
  # ds<-ds %>%
  #   mutate(link_lngth=cumsum(link_lngth))
  # target<-ds$link_id[1]
  # target_trib<-ds$trib_id[1]
  # path<-tail(ds$link_id,1)
  #
  # ds_part<-ds %>% filter(link_id==path) %>% pull(link_lngth)
  # target_part<-ds %>% filter(link_id==target) %>% pull(link_lngth)
  #
  # diff_trib_part<-us[[path]] %>%
  #   filter(trib_id != target_trib) %>%
  #   mutate(link_lngth=cumsum(link_lngth)+ds_part)
  #
  # same_trib_part<-us[[path]] %>%
  #   filter(trib_id == target_trib) %>%
  #   mutate(link_lngth=cumsum(link_lngth)-ds_part+target_part) %>%
  #   filter(link_lngth>0)
  #
  # out<-bind_rows(
  #   diff_trib_part,
  #   same_trib_part
  # ) %>%
  #   mutate(origin=target,
  #          destination=link_id,
  #          undirected_path_length=link_lngth) %>%
  #   select(origin,destination,undirected_path_length)
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  #
  # us_fp_fun<-function(link_id,db_fp=db_fp){
  #   con <- DBI::dbConnect(RSQLite::SQLite(), db_fp[[1]])
  #   out<-DBI::dbGetQuery(con, paste0("SELECT * FROM us_flowpaths WHERE source_id IN (",paste0(link_id,collapse = ","),")")) %>%
  #     group_by(source_id) %>%
  #     nest() %>%
  #     ungroup()
  #
  #   out2<-out$data
  #   names(out2)<-out$source_id
  #
  #   out2<-out2[as.character(link_id)]
  #
  #   DBI::dbDisconnect(con)
  #   return(out2)
  # }
  #
  #
  #
  #
  # out<-DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths")) %>%
  #   mutate(origin=origin_id) %>%
  #   mutate(destination=link_id) %>%
  #   group_by(origin_id) %>%
  #   arrange(USChnLn_Fr) %>%
  #   mutate(directed_path_length=cumsum(link_lngth)) %>%
  #   select(origin_id,origin,destination,link_lngth,directed_path_length) %>%
  #   nest() %>%
  #   ungroup()
  #
  # out2<-out$data
  # names(out2)<-out$source_id
  #
  #
  # out2<-out$data
  # names(out2)<-out$source_id
  #
  # out2<-out2[link_id]
  #
  # DBI::dbDisconnect(con)
  # return(out2)
  #
  # ds_flowpaths
  #
  # #us_catchment_areas<-map(us_flowpaths,function(x) sum(x$sbbsn_area,na.rm=T))
  # us_catchment_areas<-DBI::dbGetQuery(con, paste0("SELECT source_id, SUM(sbbsn_area) AS sbbsn_area FROM us_flowpaths GROUP BY source_id"))
  # us_list<-rep(us_flowpaths_file,length(ds_flowpaths))
  # names(us_list)<-names(ds_flowpaths)
  #
  # if (verbose) print("Generating Pairwise Distances")
  #
  # out_tbl<-tibble(
  #   origin=rep(names(ds_flowpaths),each=as.numeric(length(ds_flowpaths))),
  #   destination=rep(names(ds_flowpaths),length.out=as.numeric(length(ds_flowpaths))*as.numeric(length(ds_flowpaths)))
  # )
  #
  # # Function - Flow Connected Distances
  # ds_pwise<-function(x,p) {
  #   out<-mutate(x,origin=link_id[1]) %>%
  #     rename(destination=link_id) %>%
  #     mutate(directed_path_length=cumsum(link_lngth)) %>%
  #     select(origin,destination,directed_path_length)
  #
  #   p()
  #   return(out)
  # }
  #
  # # Function - non-Flow Connected Distances
  # us_pwise<-function(ds,us_path,p) {
  #
  #   us<-readRDS(us_path)
  #
  #   ds<-ds %>%
  #     mutate(link_lngth=cumsum(link_lngth))
  #   target<-ds$link_id[1]
  #   target_trib<-ds$trib_id[1]
  #   path<-tail(ds$link_id,1)
  #
  #   ds_part<-ds %>% filter(link_id==path) %>% pull(link_lngth)
  #   target_part<-ds %>% filter(link_id==target) %>% pull(link_lngth)
  #
  #   diff_trib_part<-us[[path]] %>%
  #     filter(trib_id != target_trib) %>%
  #     mutate(link_lngth=cumsum(link_lngth)+ds_part)
  #
  #   same_trib_part<-us[[path]] %>%
  #     filter(trib_id == target_trib) %>%
  #     mutate(link_lngth=cumsum(link_lngth)-ds_part+target_part) %>%
  #     filter(link_lngth>0)
  #
  #   out<-bind_rows(
  #     diff_trib_part,
  #     same_trib_part
  #   ) %>%
  #     mutate(origin=target,
  #            destination=link_id,
  #            undirected_path_length=link_lngth) %>%
  #     select(origin,destination,undirected_path_length)
  #
  #   p()
  #
  #   return(out)
  # }
  #
  #
  # print("Generating Flow Connected Distances")
  # with_progress(enable=T,{
  #   p <- progressor(steps = length(ds_flowpaths))
  #
  #
  #
  #   #ds_out<-future_map_dfr(ds_flowpaths,~ds_pwise(.,p=p))
  # })
  #
  # print("Generating Flow Unconnected Distances")
  # with_progress(enable=T,{
  #   p <- progressor(steps = length(ds_flowpaths))
  #   us_out<-future_map2_dfr(ds_flowpaths,
  #                           us_list,
  #                           ~us_pwise(ds=.x,us=.y,p=p))
  # })
  #
  # out_tbl<-out_tbl %>%
  #   left_join(ds_out,by=c("origin","destination")) %>%
  #   left_join(us_out,by=c("origin","destination")) %>%
  #   # mutate(directed_path_length=if_else(is.na(directed_path_length),0,directed_path_length)) %>%
  #   # mutate(directed_path_length=if_else(is.na(directed_path_length),0,directed_path_length)) %>%
  #   mutate(
  #     origin_catchment=unlist(us_catchment_areas[origin]),
  #     destination_catchment=unlist(us_catchment_areas[destination])
  #   ) %>%
  #   mutate(prop_shared_catchment=case_when(
  #     directed_path_length>0 ~ origin_catchment/destination_catchment,
  #     T ~ 0
  #   ))

  return("pairwise_pathlengths")
}

