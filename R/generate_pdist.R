#' Generate pairwise distances
#'
#' @param input resulting object from `trace_flowpaths()`
#' @param pwise_all_links logical. Should all pairwise distances be calculate, or only those originating from sampling points
#' @param verbose logical.
#'
#' @return input with table 'fcon_pwise_dist' and 'funcon_pwise_dist' added to the database
#' @export

generate_pdist<-function(
    input,
    verbose=F,
    pwise_all_links=F
) {
  if (!is.logical(verbose)) stop("'verbose' must be logical")
  if (!is.logical(pwise_all_links)) stop("'pwise_all_links' must be logical")
  options(future.rng.onMisuse="ignore")
  options(scipen = 999)

  if (pwise_all_links) warning("Using pwise_all_links=T can be very slow for large datasets")

  zip_loc<-input$outfile
  db_loc<-input$db_loc

  site_id_col<-paste0(data.table::fread(cmd=paste("unzip -p ",zip_loc,"site_id_col.csv")))
  if (pwise_all_links) site_id_col<-"link_id"

  con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)

  #browser()
  # DS directed path-lengths
  if (verbose) print("Calculating Flow Connected Distances")
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
    rename(origin=origin_link_id) %>%
    rename(destination=destination_link_id) %>%
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
          by=c("origin_link_id"="link_id")
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
          by=c("origin_link_id"="link_id")
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
  DBI::dbExecute(con,"DROP TABLE IF EXISTS funcon_pwise_dist")

  if (pwise_all_links) {
    if (verbose) print("Calculating Flow Unconnected Distances")

    flowUNconn_out<-dplyr::full_join(tbl(con,"fcon_pwise_dist") %>%
                                       select(midpoint=destination,origin_p1=origin,directed_path_length_p1=directed_path_length),
                                     tbl(con,"fcon_pwise_dist") %>%
                                       select(midpoint=destination,origin_p2=origin,directed_path_length_p2=directed_path_length),
                                     by=c("midpoint")) %>% #this is everything  not directly flow connected to each midpoint
      filter(origin_p1 != origin_p2) %>%
      anti_join(
        tbl(con,"ds_flowpaths") %>% select(origin_p2=destination_link_id,origin_p1=origin_link_id),
        by = c("origin_p1", "origin_p2")
      )%>%
      anti_join(
        tbl(con,"ds_flowpaths") %>% select(origin_p1=destination_link_id,origin_p2=origin_link_id),
        by = c("origin_p1", "origin_p2")
      ) %>%
    # left_join(tbl(con,"fcon_pwise_dist") %>% select(origin,destination,directed_path_length_dir1=directed_path_length),
    #           by=c("origin_p2"="origin",
    #                "origin_p1"="destination")) %>%
    # mutate(undirected_path_length=case_when(
    #   is.na(directed_path_length_dir1) ~ as.numeric(directed_path_length_p1)+as.numeric(directed_path_length_p2),
    #   T ~ 0
    # )) %>%
    mutate(undirected_path_length=as.numeric(directed_path_length_p1)+as.numeric(directed_path_length_p2)) %>%
      select(origin=origin_p1,destination=origin_p2,undirected_path_length) %>%
      mutate(directed_path_length=0,
             prop_shared_catchment=0,
             prop_shared_logcatchment=0)


    t1<-flowUNconn_out %>%
      copy_to(df=.,
              con,
              "funcon_pwise_dist",
              overwrite =T,
              temporary =F,
              indexes=c("origin","destination"),
              analyze=T,
              in_transaction=T)

  }

  DBI::dbDisconnect(con)

  return(input)

}
