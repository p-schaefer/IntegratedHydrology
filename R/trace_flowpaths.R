
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
  # out_file<-zip_loc

  db_loc<-file.path(gsub(basename(zip_loc),"",zip_loc),gsub(".zip",".db",basename(zip_loc)))

  site_id_col<-paste0(data.table::fread(cmd=paste("unzip -p ",zip_loc,"site_id_col.csv")))

  final_links<-as_tibble(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv"))) %>%
    mutate(across(any_of(site_id_col),na_if,""))


  fp<-trace_flowpath_fn(input=final_links,
                        verbose=verbose,
                        db_loc=db_loc,
                        temp_dir=temp_dir)

  # dist_list_out<-list(
  #   fp
  # )
  #
  #
  # if (verbose) print("Generating Output")
  #
  # zip(out_file,
  #     unlist(dist_list_out),
  #     flags = '-r9Xjq'
  # )

  output<-input[!names(input) %in% c("catchment_poly")]

  output$db_loc<-db_loc

  all_catch<-get_catchment(input = output,
                           target_points = final_links[["link_id"]]) %>%
    select(link_id)

  write_sf(all_catch %>% select(link_id),file.path(temp_dir,"Catchment_poly.shp"))

  zip(zip_loc,
      list.files(temp_dir,"Catchment_poly",full.names = T),
      flags = '-r9Xjq'
  )

  if (return_products){
    output<-c(
      list(catchments=all_catch %>% select(link_id)
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
    db_loc,
    verbose=F,
    temp_dir=NULL
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
    select(-uslink_id1,-ustrib_id1) %>%
    distinct() %>%
    group_by(trib_id) %>%
    arrange(USChnLn_Fr) %>%
    ungroup()

  # Downstream paths
  if (verbose) print("Calculating Flowpaths")
  ds_fp<-db_loc
  # if (file.exists(ds_fp)) stop(paste0("sqlite database: ",ds_fp," Already Exists, please delete the file before proceeding."))
  if (file.exists(ds_fp)) {
    warning(paste0("sqlite database: ",ds_fp," Already Exists, it was deleted and replaced."))
    file.remove(ds_fp)
    }

  with_progress(enable=T,{

    int_tribs<-input_tib %>%
      mutate(trib_id2=trib_id) %>%
      group_by(trib_id) %>%
      nest() %>%
      ungroup() %>%
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

    p <- progressor(steps = nrow(int_tribs))

    # Exiting tribs
    final_ds_paths<-int_tribs %>%
      filter(is.na(dslink_id)) %>%
      select(trib_id,data) %>%
      mutate(p=list(p)) %>%
      mutate(data2=future_map2(data,p,
                               .options = furrr_options(globals = FALSE),
                               carrier::crate(function(x,p){
                                 options(scipen = 999)
                                 `%>%` <- magrittr::`%>%`

                                 x<- x%>%
                                   dplyr::select(link_id,trib_id=trib_id2,link_lngth,sbbsn_area,USChnLn_Fr)
                                 u<-x$link_id
                                 u<-stats::setNames(u,u)
                                 out<-purrr:::map(u, function (y) x[which(x$link_id==y):nrow(x),] %>% dplyr::mutate(origin_id=y)) %>%
                                   dplyr::bind_rows()
                                 p()
                                 return(out)
                               })))



    final_ds_paths_out<-final_ds_paths %>%
      select(data2) %>%
      unnest(data2)

    con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
    #DBI::dbExecute(con, "PRAGMA busy_timeout = 10000")
    ot<-DBI::dbCreateTable(con, "ds_flowpaths", final_ds_paths_out[F,])
    build_view<-DBI::dbExecute(con,"CREATE INDEX inx_ds_flowpaths ON ds_flowpaths (link_id, origin_id)")

    ot<-DBI::dbAppendTable(con, "ds_flowpaths", final_ds_paths_out)
    DBI::dbDisconnect(con)

    final_ds_paths_out<-unique(final_ds_paths_out$trib_id)

    int_tribs<-int_tribs %>%
      filter(!trib_id %in% final_ds_paths_out)

    ds_fp_fun<-function(dslink_id,db_fp=ds_fp){
      con <- DBI::dbConnect(RSQLite::SQLite(), db_fp[[1]])
      out<-DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths WHERE origin_id IN (",paste0(dslink_id,collapse = ","),")")) %>%
        group_by(origin_id) %>%
        nest() %>%
        ungroup()

      out2<-out$data
      names(out2)<-out$origin_id

      out2<-out2[as.character(dslink_id)]

      DBI::dbDisconnect(con)
      return(out2)
    }

    while(nrow(int_tribs)>0){
      int_tribs_int<-int_tribs %>%
        filter(dstrib_id %in% final_ds_paths_out) %>%
        mutate(p=list(p)) %>%
        mutate(ds_path=ds_fp_fun(dslink_id,db_fp=ds_fp))

      int_tribs_int<-int_tribs_int %>%
        mutate(ds_path=future_pmap(list(
          data=data,
          ds_path=ds_path,
          p=p
        ),
        .options = furrr_options(globals = FALSE),
        carrier::crate(
          function(data,ds_path,p){
            options(scipen = 999)
            `%>%` <- magrittr::`%>%`

            out<-dplyr::bind_rows(data %>%
                                    dplyr::rename(trib_id=trib_id2),
                                  ds_path %>%
                                    dplyr::mutate(link_id=as.character(link_id))
            ) %>%
              dplyr::distinct()

            x<- out%>%
              dplyr::select(link_id,trib_id,link_lngth,sbbsn_area,USChnLn_Fr) %>%
              dplyr::distinct()
            u<-data$link_id
            u<-stats::setNames(u,u)
            x<-rep(list(x),length(u))
            out<-purrr::map2(u,x,
                             carrier::crate(function (y,x) dplyr::mutate(x[which(x$link_id==y):nrow(x),],origin_id=y)))
            out<-dplyr::bind_rows(out)

            p()
            return(out)
          })
        ))

      int_tribs_int<-int_tribs_int %>%
        select(ds_path) %>%
        unnest(ds_path)

      con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
      ot<-DBI::dbAppendTable(con, "ds_flowpaths", int_tribs_int)
      int_tribs<-int_tribs %>% filter(!trib_id %in% int_tribs_int$trib_id)
      final_ds_paths_out<-unique(DBI::dbGetQuery(con, "SELECT DISTINCT trib_id FROM ds_flowpaths")$trib_id)
      DBI::dbDisconnect(con)

    }
  })


  con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)

  # Upstream portion is saved as a view
  sql_code<-dplyr::tbl(con,"ds_flowpaths") %>%
    dplyr::select(link_id,origin_id) %>%
    dplyr::distinct() %>%
    #window_order(link_id) %>% # up to here, is a list of every origin ID draining into each link_id
    left_join(dplyr::tbl(con,"ds_flowpaths") %>% #This adds the extra info from each orgin link to the table
                select(-origin_id),
              by=c("origin_id"="link_id")) %>%
    #window_order(link_id) %>%
    #dplyr::distinct() %>%
    dplyr::rename(source_id=link_id,
                  link_id=origin_id) %>%
    sql_render()

  sql_code2<-paste0("CREATE TABLE us_flowpaths AS ",gsub("\n","",sql_code))

  build_view<-DBI::dbExecute(con,sql_code2)
  build_view<-DBI::dbExecute(con,"CREATE INDEX inx_us_flowpaths ON us_flowpaths (link_id, source_id)")

  # Pairwise distance View --------------------------

  # DS directed path-lengths
  flowconn_out<-tbl(con,"ds_flowpaths") %>%
    mutate(origin=origin_id) %>%
    mutate(destination=link_id) %>%
    group_by(origin_id) %>%
    window_order(USChnLn_Fr) %>%
    mutate(directed_path_length=cumsum(link_lngth)) %>%
    ungroup() %>%
    select(origin,destination,directed_path_length) %>%
    distinct() %>%
    left_join(
      tbl(con,"us_flowpaths") %>%
        group_by(source_id) %>%
        summarize(origin_catchment=sum(sbbsn_area)) %>%
        rename(origin=source_id),
      by="origin"
    ) %>%
    left_join(
      tbl(con,"us_flowpaths") %>%
        group_by(source_id) %>%
        summarize(destination_catchment=sum(sbbsn_area)) %>%
        rename(destination=source_id),
      by="destination"
    ) %>%
    mutate(prop_shared_catchment=case_when(
      directed_path_length>0 ~ as.numeric(origin_catchment) / as.numeric(destination_catchment),
      T ~ 0
    )) %>%
    mutate(prop_shared_logcatchment=case_when(
      directed_path_length>0 ~ log(as.numeric(origin_catchment)) / log(as.numeric(destination_catchment)),
      T ~ 0
    )) %>%
    mutate(undirected_path_length=directed_path_length) %>%
    select(-origin_catchment,-destination_catchment)

  # Reverse are not directly connected, but indirectly are connected
  # flowconn_out_rev<-tbl(con,"ds_flowpaths") %>%
  #   mutate(origin=origin_id) %>%
  #   mutate(destination=link_id) %>%
  #   group_by(origin_id) %>%
  #   window_order(USChnLn_Fr) %>%
  #   mutate(directed_path_length=cumsum(link_lngth)) %>%
  #   ungroup() %>%
  #   select(origin,destination,directed_path_length) %>%
  #   distinct() %>%
  #   mutate(undirected_path_length=directed_path_length) %>%
  #   rename(origin = destination,
  #          destination = origin) %>%
  #   mutate(directed_path_length = 0,
  #          prop_shared_catchment = 0,
  #          prop_shared_logcatchment = 0)

  #un-directed path-lengths

  flowUNconn_out<-dplyr::full_join(flowconn_out %>% select(origin,destination,directed_path_length),
                                   flowconn_out %>% select(origin,destination,directed_path_length),
                                   by=c("destination"="origin"),
                                   suffix=c("_p1","_p2")) %>%
    filter(!is.na(origin)) %>%
    left_join(flowconn_out %>% select(origin,destination,directed_path_length_dir1=directed_path_length),
              by=c("destination_p2"="origin",
                   "destination_p1"="destination")) %>%
    mutate(undirected_path_length=case_when(
      is.na(directed_path_length_dir1) ~ as.numeric(directed_path_length_p1)+as.numeric(directed_path_length_p2),
      T ~ 0
    )) %>%
    select(origin,destination=destination_p2,undirected_path_length) %>%
    mutate(directed_path_length=0,
           prop_shared_catchment=0,
           prop_shared_logcatchment=0)

  final_sql<-flowconn_out %>%
    #rows_append(flowconn_out_rev) %>%
    rows_append(flowUNconn_out) %>%
    distinct() %>%
    sql_render()

  sql_code2<-paste0("CREATE TABLE pairwise_dist AS ",gsub("\n"," ",final_sql))

  build_view<-DBI::dbExecute(con,sql_code2)
  build_view<-DBI::dbExecute(con,"CREATE INDEX inx_pairwise_dist ON pairwise_dist (origin, destination)")


  DBI::dbDisconnect(con)

  return(ds_fp)

}

