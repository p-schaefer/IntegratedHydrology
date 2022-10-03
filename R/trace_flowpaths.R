
#' Generate upstream and downstream flow direction lists originating from each subbasin (and sampling point)
#'
#' @param input resulting object from `attrib_streamline()` or `insert_points()`
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param calc_catch character. One of "none", "sample_points", or "all" indicating which if any catchments should be calculated and included in the zip output
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export

trace_flowpaths<-function(
    input,
    return_products=F,
    temp_dir=NULL,
    calc_catch=c("all","none","sample_points"),
    verbose=F
){
  options(future.rng.onMisuse="ignore")
  options(scipen = 999)

  calc_catch<-calc_catch[1]
  match.arg(calc_catch,several.ok = F)

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

  # final_links<-as_tibble(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv"))) %>%
  #   mutate(across(any_of(site_id_col),na_if,""))

  con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
  final_links<-collect(tbl(con,"stream_links")) %>%
    mutate(across(c(link_id,any_of(site_id_col)),as.character)) %>%
    mutate(across(any_of(site_id_col),na_if,""))
  DBI::dbDisconnect(con)

  fp<-trace_flowpath_fn(input=final_links,
                        verbose=verbose,
                        db_loc=db_loc,
                        site_id_col=site_id_col,
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


  if (calc_catch=="none"){
    all_catch<-NULL
  } else {
    if (calc_catch=="all") {
      all_catch<-get_catchment(input = output,
                               target_points = final_links[["link_id"]]) %>%
        select(link_id)
    }
    if (calc_catch=="sample_points"){
      all_catch<-get_catchment(input = output,
                               site_id_col=site_id_col,
                               target_points = final_links[[site_id_col]]) %>%
        select(link_id)

    }

    write_sf(all_catch,file.path(temp_dir,"Catchment_poly.shp"))

    zip(zip_loc,
        list.files(temp_dir,"Catchment_poly",full.names = T),
        flags = '-r9Xjq'
    )
  }





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
    site_id_col,
    verbose=F,
    temp_dir=NULL
) {
  con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
  DBI::dbExecute(con,"DROP TABLE IF EXISTS ds_flowpaths")

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

    # con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
    # DBI::dbExecute(con, "PRAGMA busy_timeout = 10000")

    ot<-final_ds_paths_out %>%
      select(destination_link_id=link_id,origin_link_id=origin_id) %>%
      copy_to(df=.,
              con,
              "ds_flowpaths",
              overwrite =T,
              temporary =F,
              indexes=c("destination_link_id","origin_link_id"),
              analyze=T,
              in_transaction=T)

    # ot<-DBI::dbCreateTable(con, "ds_flowpaths", final_ds_paths_out[F,])
    # build_view<-DBI::dbExecute(con,"CREATE INDEX inx_ds_flowpaths ON ds_flowpaths (link_id, origin_id)")

    # ot<-DBI::dbAppendTable(con, "ds_flowpaths", final_ds_paths_out)
    # DBI::dbDisconnect(con)

    final_ds_paths_out<-unique(final_ds_paths_out$trib_id)

    int_tribs<-int_tribs %>%
      filter(!trib_id %in% final_ds_paths_out)

    ds_fp_fun<-function(dslink_id,db_fp=ds_fp){
      #con <- DBI::dbConnect(RSQLite::SQLite(), db_fp[[1]])
      out<-DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths WHERE origin_link_id IN (",paste0(dslink_id,collapse = ","),")")) %>%
        group_by(origin_link_id) %>%
        nest() %>%
        ungroup()

      out2<-out$data
      names(out2)<-out$origin_link_id

      out2<-out2[as.character(dslink_id)]

      #DBI::dbDisconnect(con)
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
            #browser()
            options(scipen = 999)
            `%>%` <- magrittr::`%>%`

            out<-dplyr::bind_rows(data %>%
                                    dplyr::select(link_id) %>%
                                    dplyr::rename(destination_link_id=link_id),
                                  ds_path
            ) #%>%
            #dplyr::distinct()

            x<- out%>%
              dplyr::distinct()
            u<-data$link_id
            u<-stats::setNames(u,u)
            x<-rep(list(x),length(u))
            out<-purrr::map2(u,x,
                             carrier::crate(function (y,x) dplyr::mutate(x[which(x$destination_link_id==y):nrow(x),],origin_link_id=y)))
            out<-dplyr::bind_rows(out)

            p()
            return(out)
          })
        ))

      int_tribs_int<-int_tribs_int %>%
        select(ds_path) %>%
        unnest(ds_path)

      #con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
      ot<-DBI::dbAppendTable(con, "ds_flowpaths", int_tribs_int)

      ot<-tbl(con,"ds_flowpaths") %>%
        left_join(
          tbl(con,"stream_links") %>%
            select(link_id,trib_id),
          by=c("destination_link_id"="link_id")
        ) %>%
        collect()

      int_tribs<-int_tribs %>% filter(!trib_id %in% ot$trib_id)

      final_ds_paths_out<-tbl(con,"ds_flowpaths") %>%
        left_join(
          tbl(con,"stream_links") %>%
            select(link_id,trib_id),
          by=c("destination_link_id"="link_id")
        ) %>%
        select(trib_id) %>%
        distinct() %>%
        collect() %>%
        pull(trib_id)
      #final_ds_paths_out<-unique(DBI::dbGetQuery(con, "SELECT DISTINCT trib_id FROM ds_flowpaths")$trib_id)
      #DBI::dbDisconnect(con)

    }
  })

  #browser()

  #con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)

  # Upstream portion is saved as a view

  DBI::dbExecute(con,"DROP TABLE IF EXISTS us_flowpaths")

  sql_code<-tbl(con,"ds_flowpaths") %>%
    select(pour_point_id=origin_link_id) %>%
    distinct() %>%
    full_join(
      tbl(con,"ds_flowpaths") ,
      by=c("pour_point_id"="destination_link_id")
    ) %>%
    copy_to(df=.,
            con,
            "us_flowpaths",
            overwrite =T,
            temporary =F,
            indexes=c("pour_point_id","origin_link_id"),
            analyze=T,
            in_transaction=T)
  #%>%
  #dbplyr::sql_render()

  #sql_code2<-paste0("CREATE TABLE us_flowpaths AS ",gsub("\n","",sql_code))

  #build_view<-DBI::dbExecute(con,sql_code2)

  # sql_code<-dplyr::tbl(con,"ds_flowpaths") %>%
  #   dplyr::select(link_id,origin_id) %>%
  #   dplyr::distinct() %>%
  #   #window_order(link_id) %>% # up to here, is a list of every origin ID draining into each link_id
  #   left_join(dplyr::tbl(con,"ds_flowpaths") %>% #This adds the extra info from each orgin link to the table
  #               select(-origin_id),
  #             by=c("origin_id"="link_id")) %>%
  #   #window_order(link_id) %>%
  #   #dplyr::distinct() %>%
  #   dplyr::rename(source_id=link_id,
  #                 link_id=origin_id) %>%
  #   sql_render()
  #
  # sql_code2<-paste0("CREATE TABLE us_flowpaths AS ",gsub("\n","",sql_code))
  #
  # build_view<-DBI::dbExecute(con,sql_code2)
  # build_view<-DBI::dbExecute(con,"CREATE INDEX inx_us_flowpaths ON us_flowpaths (link_id, source_id)")

  # Pairwise distance View --------------------------

  if (site_id_col != "link_id"){

    #browser()
    # UID<-tbl(con,"stream_links") %>%
    #   filter(!is.na(!!sym(site_id_col))) %>%
    #   select(link_id) %>%
    #   collect() %>%
    #   pull(link_id)

    # DS directed path-lengths
    flowconn_out<-tbl(con,"stream_links") %>%
      filter(!is.na(!!sym(site_id_col))) %>%
      select(link_id) %>%
      rename(origin_link_id=link_id) %>%
      left_join(
        tbl(con,"ds_flowpaths"),
        by=c("origin_link_id")
      ) %>%
      left_join(
        tbl(con,"stream_links") %>%
          select(link_id,link_lngth,USChnLn_Fr,sbbsn_area),
        by=c("destination_link_id"="link_id")
      ) %>%
      mutate(origin=origin_link_id) %>%
      mutate(destination=destination_link_id) %>%
      group_by(origin) %>%
      dbplyr::window_order(USChnLn_Fr) %>%
      mutate(directed_path_length=cumsum(link_lngth)) %>%
      ungroup() %>%
      select(origin,destination,directed_path_length) %>%
      distinct() %>%
      left_join(
        tbl(con,"us_flowpaths") %>%
          left_join(
            tbl(con,"stream_links") %>%
              select(link_id,sbbsn_area),
            by=c("pour_point_id"="link_id")
          ) %>%
          group_by(pour_point_id) %>%
          summarize(origin_catchment=sum(sbbsn_area)) %>%
          rename(origin=pour_point_id),
        by="origin"
      ) %>%
      left_join(
        tbl(con,"us_flowpaths") %>%
          left_join(
            tbl(con,"stream_links") %>%
              select(link_id,sbbsn_area),
            by=c("pour_point_id"="link_id")
          ) %>%
          group_by(pour_point_id) %>%
          summarize(destination_catchment=sum(sbbsn_area)) %>%
          rename(destination=pour_point_id),
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

    DBI::dbExecute(con,"DROP TABLE IF EXISTS fcon_pwise_dist")
    t1<-flowconn_out %>%
      copy_to(df=.,
              con,
              "fcon_pwise_dist",
              overwrite =T,
              temporary =F,
              indexes=c("origin","destination"),
              analyze=T,
              in_transaction=T)


    #un-directed path-lengths

    flowUNconn_out<-dplyr::full_join(tbl(con,"fcon_pwise_dist") %>% select(origin,destination,directed_path_length),
                                     tbl(con,"fcon_pwise_dist") %>% select(origin,destination,directed_path_length),
                                     by=c("destination"="origin"),
                                     suffix=c("_p1","_p2")) %>%
      filter(!is.na(origin)) %>%
      left_join(tbl(con,"fcon_pwise_dist") %>% select(origin,destination,directed_path_length_dir1=directed_path_length),
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

    DBI::dbExecute(con,"DROP TABLE IF EXISTS funcon_pwise_dist")



    t1<-flowconn_out %>%
      copy_to(df=.,
              con,
              "funcon_pwise_dist",
              overwrite =T,
              temporary =F,
              indexes=c("origin","destination"),
              analyze=T,
              in_transaction=T)

    # final_sql<-flowconn_out %>%
    #   rows_append(flowUNconn_out) %>%
    #   distinct() %>%
    #   copy_to(df=.,
    #           con,
    #           "pairwise_dist",
    #           overwrite =T,
    #           temporary =F,
    #           indexes=c("origin","destination"),
    #           analyze=T,
    #           in_transaction=T)


    # final_sql<-flowconn_out %>%
    #   rows_append(flowUNconn_out) %>%
    #   distinct() %>%
    #   dbplyr::sql_render()
    #
    # sql_code2<-paste0("CREATE VIEW pairwise_dist AS ",gsub("\n"," ",final_sql))
    #
    # build_view<-DBI::dbExecute(con,sql_code2)
  }



  DBI::dbDisconnect(con)

  return(db_loc)

}

