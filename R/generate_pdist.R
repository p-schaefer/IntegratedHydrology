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
  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  if (!is.logical(verbose)) stop("'verbose' must be logical")
  if (!is.logical(pwise_all_links)) stop("'pwise_all_links' must be logical")
  options(future.rng.onMisuse="ignore")
  options(scipen = 999)

  if (pwise_all_links) message("Using pwise_all_links=T can be very slow for large datasets")

  db_loc<-db_fp<-zip_loc<-input$outfile
  #db_loc<-input$db_loc

  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)

  site_id_col<-dplyr::collect(dplyr::tbl(con,"site_id_col"))$site_id_col

  if (pwise_all_links) site_id_col<-"link_id"

  #con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)

  #browser()
  # DS directed path-lengths
  DBI::dbExecute(con,"DROP TABLE IF EXISTS fcon_pwise_dist")

  if (verbose) message("Calculating Flow Connected Distances")
  flowconn_out<-dplyr::tbl(con,"stream_links_attr") %>%
    dplyr::filter(!is.na(dbplyr::sql(site_id_col))) %>%
    dplyr::select(link_id) %>%
    dplyr::rename(origin_link_id=link_id) %>%
    dplyr::left_join(
      dplyr::tbl(con,"ds_flowpaths"),
      by=c("origin_link_id")
    ) %>%
    dplyr::left_join(
      dplyr::tbl(con,"stream_links_attr") %>%
        dplyr::select(link_id,link_lngth,USChnLn_Fr,sbbsn_area),
      by=c("destination_link_id"="link_id")
    ) %>%
    dplyr::rename(origin=origin_link_id) %>%
    dplyr::rename(destination=destination_link_id) %>%
    dplyr::group_by(origin) %>%
    dbplyr::window_order(USChnLn_Fr) %>% # dbplyr::sql("USChnLn_Fr") # used to need this
    dplyr::mutate(directed_path_length=cumsum(link_lngth)) %>%
    dplyr::ungroup() %>%
    dplyr::select(origin,destination,directed_path_length) %>%
    distinct() %>%
    dplyr::left_join(
      dplyr::tbl(con,"us_flowpaths") %>%
        dplyr::left_join(
          dplyr::tbl(con,"stream_links_attr") %>%
            dplyr::select(link_id,sbbsn_area),
          by=c("origin_link_id"="link_id")
        ) %>%
        dplyr::group_by(pour_point_id) %>%
        dplyr::summarize(origin_catchment=sum(sbbsn_area,na.rm=T)) %>%
        dplyr::rename(origin=pour_point_id),
      by="origin"
    ) %>%
    dplyr::left_join(
      dplyr::tbl(con,"us_flowpaths") %>%
        dplyr::left_join(
          dplyr::tbl(con,"stream_links_attr") %>%
            dplyr::select(link_id,sbbsn_area),
          by=c("origin_link_id"="link_id")
        ) %>%
        dplyr::group_by(pour_point_id) %>%
        dplyr::summarize(destination_catchment=sum(sbbsn_area,na.rm=T)) %>%
        dplyr::rename(destination=pour_point_id),
      by="destination"
    ) %>%
    dplyr::mutate(prop_shared_catchment=dplyr::case_when(
      directed_path_length>0 ~ as.numeric(origin_catchment) / as.numeric(destination_catchment),
      T ~ 0
    )) %>%
    dplyr::mutate(prop_shared_logcatchment=dplyr::case_when(
      directed_path_length>0 ~ log(as.numeric(origin_catchment)) / log(as.numeric(destination_catchment)),
      T ~ 0
    )) %>%
    dplyr::mutate(undirected_path_length=directed_path_length) %>%
    dplyr::select(-origin_catchment,-destination_catchment) %>%
    dplyr::compute(name = "fcon_pwise_dist",
                   temporary = F,
                   overwrite=T,
                   indexes = c("origin","destination"))


  #un-directed path-lengths
  DBI::dbExecute(con,"DROP TABLE IF EXISTS funcon_pwise_dist")

  #if (pwise_all_links) {
    if (verbose) message("Calculating Flow Unconnected Distances")

    flowUNconn_out<-dplyr::full_join(dplyr::tbl(con,"fcon_pwise_dist") %>%
                                       dplyr::select(midpoint=destination,origin_p1=origin,directed_path_length_p1=directed_path_length),
                                     dplyr::tbl(con,"fcon_pwise_dist") %>%
                                       dplyr::select(midpoint=destination,origin_p2=origin,directed_path_length_p2=directed_path_length),
                                     by=c("midpoint")) %>% #this is everything  not directly flow connected to each midpoint
      dplyr::filter(origin_p1 != origin_p2) %>%
      dplyr::anti_join(
        dplyr::tbl(con,"ds_flowpaths") %>% dplyr::select(origin_p2=destination_link_id,origin_p1=origin_link_id),
        by = c("origin_p1", "origin_p2")
      )%>%
      dplyr::anti_join(
        dplyr::tbl(con,"ds_flowpaths") %>% dplyr::select(origin_p1=destination_link_id,origin_p2=origin_link_id),
        by = c("origin_p1", "origin_p2")
      ) %>%
      dplyr::mutate(undirected_path_length=as.numeric(directed_path_length_p1)+as.numeric(directed_path_length_p2)) %>%
      dplyr::select(origin=origin_p1,destination=origin_p2,undirected_path_length) %>%
      dplyr::mutate(directed_path_length=0,
                    prop_shared_catchment=0,
                    prop_shared_logcatchment=0)%>%
      dplyr::compute(name = "funcon_pwise_dist",
                     temporary = F,
                     overwrite=T,
                     indexes = c("origin","destination"))



  #}

  DBI::dbDisconnect(con)

  class(input)<-"ihydro"
  return(input)

}
