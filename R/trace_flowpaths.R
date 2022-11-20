
#' Generate upstream and downstream flow direction lists originating from each subbasin (and sampling point)
#'
#' @param input resulting object from `generate_vectors()`
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param pwise_dist logical. Calculate pairwise distances.
#' @param pwise_all_links logical. Should all pairwise distances be calculate, or only those originating from sampling points
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export

#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom dplyr collect tbl mutate across na_if
#' @importFrom RSQLite SQLite
#' @importFrom tidyselect any_of

trace_flowpaths<-function(
    input,
    return_products=F,
    temp_dir=NULL,
    #calc_catch=c("all","none","sample_points"),
    pwise_dist=F,
    pwise_all_links=F,
    verbose=F
){
  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  options(future.rng.onMisuse="ignore")
  options(scipen = 999)

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")
  if (!is.logical(pwise_dist)) stop("'pwise_dist' must be logical")
  if (!is.logical(pwise_all_links)) stop("'pwise_all_links' must be logical")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  options(dplyr.summarise.inform = FALSE)

  db_loc<-db_fp<-zip_loc<-input$outfile
  #db_loc<-input$db_loc

  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)

  site_id_col<-dplyr::collect(dplyr::tbl(con,"site_id_col"))$site_id_col

  # final_links<-as_tibble(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv"))) %>%
  #   mutate(across(any_of(site_id_col),na_if,""))

  #con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
  final_links<-dplyr::collect(dplyr::tbl(con,"stream_links_attr")) %>%
    dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
    dplyr::mutate(dplyr::across(tidyselect::any_of(site_id_col),~dplyr::na_if(.,"")))
  DBI::dbDisconnect(con)

  fp<-trace_flowpath_fn(input=final_links,
                        verbose=verbose,
                        db_loc=db_loc,
                        site_id_col=site_id_col,
                        temp_dir=temp_dir)

  if (pwise_dist){
    input<-generate_pdist(
      input=input,
      verbose=verbose,
      pwise_all_links=pwise_all_links
    )
  }


  suppressWarnings(file.remove(list.files(temp_dir,full.names = T)))

  class(input)<-"ihydro"
  return(input)

}

#' @export
#' @importFrom carrier crate
#' @importFrom DBI dbConnect dbExecute dbGetQuery dbAppendTable dbDisconnect
#' @importFrom dplyr select filter mutate case_when distinct group_by arrange ungroup bind_rows copy_to rename tbl left_join collect pull full_join
#' @importFrom furrr future_map furrr_options future_map2 future_pmap
#' @importFrom progressr with_progress progressor
#' @importFrom purrr map2
#' @importFrom RSQLite SQLite
#' @importFrom stats setNames
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer pivot_wider unnest nest
#' @importFrom tidyselect starts_with
#' @importFrom utils tail
#' @keywords internal
#' @noRd
#' @importFrom carrier crate
#' @importFrom DBI dbConnect dbExecute dbGetQuery dbAppendTable dbDisconnect
#' @importFrom dplyr select filter mutate case_when distinct group_by arrange ungroup bind_rows copy_to rename tbl left_join collect pull full_join
#' @importFrom furrr future_map furrr_options future_map2 future_pmap
#' @importFrom progressr with_progress progressor
#' @importFrom purrr map2
#' @importFrom RSQLite SQLite
#' @importFrom stats setNames
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer pivot_wider unnest nest
#' @importFrom tidyselect starts_with
#' @importFrom utils tail
#' @importFrom carrier crate
#' @importFrom DBI dbConnect dbExecute dbGetQuery dbAppendTable dbDisconnect
#' @importFrom dplyr select filter mutate case_when distinct group_by arrange ungroup bind_rows copy_to rename tbl left_join collect pull full_join
#' @importFrom furrr future_map furrr_options future_map2 future_pmap
#' @importFrom progressr with_progress progressor
#' @importFrom purrr map2
#' @importFrom RSQLite SQLite
#' @importFrom stats setNames
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer pivot_wider unnest nest
#' @importFrom tidyselect starts_with
#' @importFrom utils tail
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
    tibble::as_tibble() %>%
    dplyr::select(link_id,trib_id,link_lngth,sbbsn_area,USChnLn_Fr,
           tidyselect::starts_with("uslink_id"),tidyselect::starts_with("dslink_id"),tidyselect::starts_with("dstrib_id"),tidyselect::starts_with("ustrib_id")) %>%
    dplyr::filter(!is.na(link_id)) %>%
    dplyr::mutate(link_id=as.character(link_id)) %>%
    tidyr::pivot_longer(cols=c(tidyselect::starts_with("uslink_id"),tidyselect::starts_with("ustrib_id"))) %>%
    dplyr::filter(dplyr::case_when(
      name %in% c("uslink_id1","ustrib_id1") ~ T,
      is.na(value) ~ F,
      T ~ T
    )) %>%
    dplyr::mutate(name=ifelse(grepl("trib",name),"ustrib_id1","uslink_id1")) %>%
    tidyr::pivot_wider(names_from=name,values_from=value,values_fn=list) %>%
    tidyr::unnest(cols=c(uslink_id1,ustrib_id1)) %>%
    dplyr::select(-uslink_id1,-ustrib_id1) %>%
    dplyr::distinct() %>%
    dplyr::group_by(trib_id) %>%
    dplyr::arrange(USChnLn_Fr) %>%
    dplyr::ungroup()

  # Downstream paths
  if (verbose) message("Calculating Downstream Flowpaths")

  progressr::with_progress(enable=T,{

    int_tribs<-input_tib %>%
      dplyr::mutate(trib_id2=trib_id) %>%
      dplyr::group_by(trib_id) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(join_id=furrr::future_map(data,
                                .options = furrr::furrr_options(globals = FALSE),
                                carrier::crate(function(x){
                                  options(scipen = 999)
                                  `%>%` <- magrittr::`%>%`

                                  dplyr::arrange(x,USChnLn_Fr) %>%
                                    utils::tail(1) %>%
                                    dplyr::select(dslink_id=dslink_id1,dstrib_id=dstrib_id1)
                                } ))) %>%
      tidyr::unnest(join_id) %>%
      dplyr::mutate(join_USChnLn_Fr=input_tib$USChnLn_Fr[match(.$dslink_id,input_tib$link_id)])

    p <- progressr::progressor(steps = nrow(int_tribs))

    # Exiting tribs
    final_ds_paths<-int_tribs %>%
      dplyr::filter(is.na(dslink_id)) %>%
      dplyr::select(trib_id,data) %>%
      dplyr::mutate(p=list(p)) %>%
      dplyr::mutate(data2=furrr::future_map2(data,p,
                               .options = furrr::furrr_options(globals = FALSE),
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
      dplyr::select(data2) %>%
      tidyr::unnest(data2)

    # con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
    # DBI::dbExecute(con, "PRAGMA busy_timeout = 10000")

    ot<-final_ds_paths_out %>%
      dplyr::select(destination_link_id=link_id,origin_link_id=origin_id) %>%
      dplyr::copy_to(df=.,
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
      dplyr::filter(!trib_id %in% final_ds_paths_out)

    ds_fp_fun<-function(dslink_id,db_fp=ds_fp){
      #con <- DBI::dbConnect(RSQLite::SQLite(), db_fp[[1]])
      out<-DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths WHERE origin_link_id IN (",paste0(dslink_id,collapse = ","),")")) %>%
        dplyr::group_by(origin_link_id) %>%
        tidyr::nest() %>%
        dplyr::ungroup()

      out2<-out$data
      names(out2)<-out$origin_link_id

      out2<-out2[as.character(dslink_id)]

      #DBI::dbDisconnect(con)
      return(out2)
    }


    while(nrow(int_tribs)>0){
      int_tribs_int<-int_tribs %>%
        dplyr::filter(dstrib_id %in% final_ds_paths_out) %>%
        dplyr::mutate(p=list(p)) %>%
        dplyr::mutate(ds_path=ds_fp_fun(dslink_id,db_fp=ds_fp))

      int_tribs_int<-int_tribs_int %>%
        dplyr::mutate(ds_path=furrr::future_pmap(list(
          data=data,
          ds_path=ds_path,
          p=p
        ),
        .options = furrr::furrr_options(globals = FALSE),
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
        dplyr::select(ds_path) %>%
        tidyr::unnest(ds_path)

      #con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
      ot<-DBI::dbAppendTable(con, "ds_flowpaths", int_tribs_int)


      ot<-dplyr::tbl(con,"ds_flowpaths") %>%
        dplyr::left_join(
          dplyr::tbl(con,"stream_links_attr") %>%
            dplyr::select(link_id,trib_id),
          by=c("destination_link_id"="link_id")
        ) %>%
        dplyr::collect()

      int_tribs<-int_tribs %>% dplyr::filter(!trib_id %in% ot$trib_id)

      final_ds_paths_out<-dplyr::tbl(con,"ds_flowpaths") %>%
        dplyr::left_join(
          dplyr::tbl(con,"stream_links_attr") %>%
            dplyr::select(link_id,trib_id),
          by=c("destination_link_id"="link_id")
        ) %>%
        dplyr::select(trib_id) %>%
        dplyr::distinct() %>%
        dplyr::collect() %>%
        dplyr::pull(trib_id)
      #final_ds_paths_out<-unique(DBI::dbGetQuery(con, "SELECT DISTINCT trib_id FROM ds_flowpaths")$trib_id)
      #DBI::dbDisconnect(con)

    }
  })

  #browser()

  #con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
  if (verbose) message("Calculating Upstream Flowpaths")
  DBI::dbExecute(con,"DROP TABLE IF EXISTS us_flowpaths")

  sql_code<-dplyr::tbl(con,"ds_flowpaths") %>%
    dplyr::select(pour_point_id=origin_link_id) %>%
    dplyr::distinct() %>%
    dplyr::full_join(
      dplyr::tbl(con,"ds_flowpaths") ,
      by=c("pour_point_id"="destination_link_id")
    ) %>%
    dplyr::copy_to(df=.,
            con,
            "us_flowpaths",
            overwrite =T,
            temporary =F,
            indexes=c("pour_point_id","origin_link_id"),
            analyze=T,
            in_transaction=T)

  DBI::dbDisconnect(con)


  return(db_loc)

}

