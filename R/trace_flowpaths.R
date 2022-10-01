
#' Generate upstream and downstream flow direction lists originating from each subbasin (and sampling point)
#'
#' @param input resulting object from `attrib_streamline()` or `insert_points()`
#' @param chunck_size numeric. Size of chunks to process in flow directions in the database, when too large, function may crash due to lack of memory.
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
    chunck_size=2000,
    verbose=F
){
  options(future.rng.onMisuse="ignore")
  options(scipen = 999)

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

  if (!is.numeric(chunck_size)) stop("'chunck_size' must be numeric")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  options(dplyr.summarise.inform = FALSE)

  zip_loc<-input$outfile

  site_id_col<-paste0(data.table::fread(cmd=paste("unzip -p ",zip_loc,"site_id_col.csv")))

  final_links<-as_tibble(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv"))) %>%
    mutate(across(any_of(site_id_col),na_if,""))

  #browser()

  fp<-trace_flowpath_fn(input=final_links,
                        verbose=verbose,
                        temp_dir=temp_dir,
                        chunck_size=chunck_size)
  # ds_flowpaths<-fp$final_out_ds
  # us_flowpaths<-fp$final_out_us
  #
  # saveRDS(ds_flowpaths,file.path(temp_dir,"ds_flowpaths.rds"),compress = F)
  # saveRDS(us_flowpaths,file.path(temp_dir,"us_flowpaths.rds"),compress = F)

  dist_list_out<-list(
    fp
  )
  #browser()

  #dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

  out_file<-zip_loc

  if (verbose) print("Generating Output")

  zip(out_file,
      unlist(dist_list_out),
      flags = '-r9Xjq'
  )

  output<-input[!names(input) %in% c("catchment_poly")]

  #browser()

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
      list(flowpaths_db=basename(fp),
           #ds_flowpaths=ds_flowpaths,
           #us_flowpaths=us_flowpaths,
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
    verbose=F,
    chunck_size, #100 IDs at once shouldn't explode memory, right?
    temp_dir=NULL
) {

  chunck_size<-as.integer(chunck_size)

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

  #browser()
  # Downstream paths
  if (verbose) print("Calculating Downstream Flowpaths")
  ds_fp<-file.path(temp_dir,"flowpaths_out.db")
  if (file.exists(ds_fp)) stop(paste0("sqlite database: ",ds_fp,"Already Exists, please delete the file before proceeding."))
  with_progress(enable=T,{

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

    #p <- progressor(steps = length(unique(input_tib$link_id)))
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
                                   dplyr::select(link_id,trib_id=trib_id2,link_lngth,sbbsn_area,USChnLn_Fr) #%>%
                                 #dplyr::distinct()
                                 u<-x$link_id
                                 u<-stats::setNames(u,u)
                                 out<-purrr:::map(u, function (y) x[which(x$link_id==y):nrow(x),] %>% dplyr::mutate(origin_id=y)) %>%
                                   dplyr::bind_rows() #%>%
                                 #dplyr::distinct()
                                 p()
                                 return(out)
                               })))



    final_ds_paths_out<-final_ds_paths %>%
      select(data2) %>%
      unnest(data2)

    #data.table::fwrite(final_ds_paths_out,ds_fp)
    con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
    DBI::dbExecute(con, "PRAGMA busy_timeout = 10000")
    ot<-DBI::dbCreateTable(con, "ds_flowpaths", final_ds_paths_out[F,])
    ot<-DBI::dbAppendTable(con, "ds_flowpaths", final_ds_paths_out)
    DBI::dbDisconnect(con)

    final_ds_paths_out<-unique(final_ds_paths_out$trib_id)

    int_tribs<-int_tribs %>%
      filter(!trib_id %in% final_ds_paths_out)

    ds_fp_fun<-function(dslink_id,db_fp=ds_fp){
      con <- DBI::dbConnect(RSQLite::SQLite(), db_fp[[1]])
      out<-DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths WHERE origin_id IN (",paste0(dslink_id,collapse = ","),")")) %>%
        #distinct() %>%
        #mutate(link_id2=link_id) %>%
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
                             #.options = furrr::furrr_options(globals = FALSE),
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


      # int_tribs_int<-int_tribs %>%
      #   filter(dstrib_id %in% final_ds_paths_out) %>%
      #   mutate(ds_fp=ds_fp,
      #          p=list(p)) %>%
      #   mutate(data=pmap(list( #future_
      #     data=data,
      #     ds_fp=ds_fp,
      #     dstrib_id=dstrib_id,
      #     dslink_id=dslink_id,
      #     join_USChnLn_Fr=join_USChnLn_Fr,
      #     p=p),
      #     #.options = furrr_options(globals = FALSE),
      #     carrier::crate(
      #       function(data,dstrib_id,dslink_id,join_USChnLn_Fr,ds_fp,p) {
      #         #browser()
      #         options(scipen = 999)
      #         `%>%` <- magrittr::`%>%`
      #
      #         con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
      #
      #         #final_paths<-data.table::fread(cmd=paste0('awk -F, \'$1 == "',dstrib_id,'"\' ',ds_fp))
      #         #final_paths<-data.table::fread(ds_fp) %>% dplyr::filter(trib_id==dstrib_id)
      #         final_paths<-try(DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths WHERE origin_id='",dslink_id,"'")),silent=T)
      #
      #         while(inherits(final_paths,"try-error")){
      #           Sys.sleep(0.2)
      #           final_paths<-try(DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths WHERE origin_id='",dslink_id,"'")),silent=T)
      #         }
      #
      #         out<-dplyr::bind_rows(data %>%
      #                                 dplyr::rename(trib_id=trib_id2),
      #                               final_paths %>%
      #                                 # dplyr::arrange(USChnLn_Fr) %>%
      #                                 # dplyr::filter(USChnLn_Fr>=join_USChnLn_Fr) %>%
      #                                 dplyr::mutate(link_id=as.character(link_id))
      #         ) %>%
      #           dplyr::distinct()
      #
      #         x<- out%>%
      #           dplyr::select(link_id,trib_id,link_lngth,sbbsn_area,USChnLn_Fr) %>%
      #           dplyr::distinct()
      #         u<-data$link_id
      #         u<-stats::setNames(u,u)
      #         x<-rep(list(x),length(u))
      #         out<-furrr::future_map2(u,x,
      #                                 .options = furrr::furrr_options(globals = FALSE),
      #                                 carrier::crate(function (y,x) dplyr::mutate(x[which(x$link_id==y):nrow(x),],origin_id=y)))
      #         out<-dplyr::bind_rows(out)
      #
      #         ot<-try(DBI::dbAppendTable(con, "ds_flowpaths", out),silent=T)
      #
      #         while(inherits(ot,"try-error")){
      #           Sys.sleep(0.2)
      #           ot<-try(DBI::dbAppendTable(con, "ds_flowpaths", out),silent=T)
      #         }
      #
      #         DBI::dbDisconnect(con)
      #
      #         p()
      #
      #         return(NULL)
      #
      #       }
      #     )
      #   ))
      #
      # con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
      # int_tribs<-int_tribs %>% filter(!trib_id %in% int_tribs_int$trib_id)
      # final_ds_paths_out<-unique(DBI::dbGetQuery(con, "SELECT trib_id FROM ds_flowpaths")$trib_id)
      #
      # DBI::dbDisconnect(con)

    }

    # con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
    # final_ds_paths<-DBI::dbGetQuery(con, "SELECT * FROM ds_flowpaths") %>%
    #   group_by(trib_id) %>%
    #   nest() %>%
    #   ungroup()
    # DBI::dbDisconnect(con)

    #Finalize DS flowpaths
    # if (verbose) print("Finalizing Downstream Flowpaths")
    # final_ds_paths<-final_ds_paths %>%
    #   mutate(data2=future_map(data,
    #                           .options = furrr_options(globals = FALSE),
    #                           carrier::crate(function(x){
    #                             options(scipen = 999)
    #                             `%>%` <- magrittr::`%>%`
    #
    #                             x<- x%>%
    #                               dplyr::select(link_id,trib_id=trib_id2,link_lngth,sbbsn_area,USChnLn_Fr) %>%
    #                               dplyr::distinct()
    #                             u<-x$link_id
    #                             u<-stats::setNames(u,u)
    #                             purrr:::map(u, function (y) x[which(x$link_id==y):nrow(x),] %>% dplyr::mutate(origin_id=y))
    #                           })))
    #
    # final_ds_paths<-final_ds_paths %>%
    #   select(data2) %>%
    #   unnest(cols=data2) %>%
    #   unnest(cols=data2) %>%
    #   distinct() %>%
    #   select(origin_id,everything())
  })

  #browser()
  #US Paths
  if (verbose) print("Calculating Upstream Flowpaths")
  with_progress(enable=T,{

    con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
    final_us_paths<-unique(DBI::dbGetQuery(con, "SELECT DISTINCT link_id FROM ds_flowpaths")$link_id)
    colnms<-colnames(DBI::dbFetch(DBI::dbSendQuery(con, "SELECT * FROM ds_flowpaths"),n=1) %>% rename(source_id=origin_id))
    final_us<-data.frame(matrix(ncol=length(colnms)))[F,]
    colnames(final_us)<- colnms
    ot<-suppressWarnings(DBI::dbCreateTable(con, "us_flowpaths", final_us))

    DBI::dbDisconnect(con)

    chunks<-split(final_us_paths,ceiling(seq_along(final_us_paths)/chunck_size))

    p <- progressor(steps = length(chunks))



    final_us_paths<-pmap(list(
      final_us_paths=chunks,
      ds_fp=rep(list(ds_fp),length(chunks)),
      p=rep(list(p),length(chunks))
    ),
    #.options = furrr::furrr_options(globals = FALSE),
    carrier::crate(
      function(final_us_paths,ds_fp,p){
        #browser()
        options(scipen = 999)
        `%>%` <- magrittr::`%>%`
        con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
        DBI::dbExecute(con, "PRAGMA busy_timeout = 10000")

        base<-dplyr::tbl(con,"us_flowpaths")

        out <- dplyr::tbl(con,"ds_flowpaths") %>%
          dplyr::filter(link_id %in% final_us_paths) %>%
          dplyr::select(link_id,origin_id) %>%
          dplyr::distinct()

        oid<-suppressWarnings(dplyr::pull(out,origin_id))

        out<-out%>%
          dplyr::arrange(link_id)

        out<-out%>%
          dplyr::full_join(
            dplyr::tbl(con,"ds_flowpaths") %>%
              dplyr::filter(link_id %in% oid) %>%
              dplyr::select(-origin_id) %>%
              dplyr::rename(origin_id=link_id) %>%
              dplyr::distinct(),
            by="origin_id") %>%
          dplyr::distinct() %>%
          dbplyr::window_order(link_id) %>%
          dplyr::rename(source_id=link_id,
                        link_id=origin_id) %>%
          dplyr::collect()

        ot<-try(DBI::dbAppendTable(con, "us_flowpaths", out),silent=T)

        # out<-try(suppressWarnings(dplyr::rows_insert(base,out,conflict = "ignore",in_place=T)),silent=T)
        #
        while(inherits(ot,"try-error")){
          Sys.sleep(stats::runif(1,0,5))
          ot<-try(DBI::dbAppendTable(con, "us_flowpaths", out),silent=T)
        }

        DBI::dbDisconnect(con)

        # us_fp_fun<-function(link_id,db_fp=ds_fp){
        #   con <- DBI::dbConnect(RSQLite::SQLite(), db_fp[[1]])
        #   out<-DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths WHERE link_id IN (",paste0(link_id,collapse = ","),")")) %>%
        #     dplyr::select(link_id,origin_id) %>%
        #     dplyr::distinct() %>%
        #     dplyr::arrange(link_id)
        #
        #   out<-DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths WHERE link_id IN (",paste0(out$origin_id,collapse = ","),")")) %>%
        #     dplyr::select(-origin_id) %>%
        #     dplyr::rename(origin_id=link_id) %>%
        #     dplyr::distinct() %>%
        #     dplyr::full_join(out,by="origin_id") %>%
        #     dplyr::distinct() %>%
        #     dplyr::arrange(link_id) %>%
        #     dplyr::rename(source_id=link_id,
        #                   link_id=origin_id)
        #
        #   ot<-try(DBI::dbAppendTable(con, "us_flowpaths", out),silent=T)
        #
        #   while(inherits(ot,"try-error")){
        #     Sys.sleep(0.2)
        #     ot<-try(DBI::dbAppendTable(con, "us_flowpaths", out),silent=T)
        #   }
        #
        #   DBI::dbDisconnect(con)
        #   return(NULL)
        # }
        #
        # out<-us_fp_fun(final_us_paths,ds_fp)

        # con <- DBI::dbConnect(RSQLite::SQLite(), ds_fp)
        #
        # final_paths<-try(unique(DBI::dbGetQuery(con, paste0("SELECT origin_id FROM ds_flowpaths WHERE link_id='",final_us_paths,"'"))$origin_id),
        #                  silent=T)
        #
        # # while(inherits(final_paths,"try-error")){
        # #   Sys.sleep(0.2)
        # #   final_paths<-try(unique(DBI::dbGetQuery(con, paste0("SELECT origin_id FROM ds_flowpaths WHERE link_id='",final_us_paths,"'"))$origin_id),
        # #                    silent=T)
        # # }
        #
        # out<-try(DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths WHERE link_id IN (",paste0(final_paths,collapse = ","),")")),
        #          silent=T)
        #
        # # while(inherits(out,"try-error")){
        # #   Sys.sleep(0.2)
        # #   out<-try(DBI::dbGetQuery(con, paste0("SELECT * FROM ds_flowpaths WHERE link_id IN (",paste0(final_paths,collapse = ","),")")),
        # #            silent=T)
        # # }
        #
        # out<-out%>%
        #   dplyr::select(-origin_id) %>%
        #   dplyr::distinct() %>%
        #   dplyr::mutate(source_id=final_us_paths)
        #
        # ot<-try(DBI::dbAppendTable(con, "us_flowpaths", out),silent=T)
        #
        # # while(inherits(ot,"try-error")){
        # #   Sys.sleep(0.2)
        # #   ot<-try(DBI::dbAppendTable(con, "ds_flowpaths", out),silent=T)
        # # }
        #
        # DBI::dbDisconnect(con)

        p()
        return(NULL)

      }))
  })

  #browser()

  # final_us_paths<-final_ds_paths %>%
  #   group_by(link_id) %>%
  #   nest() %>%
  #   ungroup() %>%
  #   mutate(us_links=future_map(data,
  #                              .options = furrr_options(globals = FALSE),
  #                              carrier::crate(function(x){
  #                                options(scipen = 999)
  #                                `%>%` <- magrittr::`%>%`
  #
  #                                x %>% dplyr::pull(origin_id) %>% unique()
  #                              })
  #   )) %>%
  #   mutate(full_data=list(input[,c("link_id","trib_id","link_lngth","sbbsn_area","USChnLn_Fr")]))
  #
  # if (verbose) print("Finalizing Upstream Flowpaths")
  # final_us_paths<-final_us_paths%>%
  #   mutate(data2=future_map2(us_links,full_data,
  #                            .options = furrr_options(globals = FALSE),
  #                            carrier::crate(function(x,y){
  #                              y[y$link_id %in% x,c("link_id","trib_id","link_lngth","sbbsn_area","USChnLn_Fr")]
  #                            })
  #   )) %>%
  #   select(source_id=link_id, data2) %>%
  #   unnest(cols=data2)
  #
  # final_out_ds<-split(final_ds_paths[,-c(1)],final_ds_paths$origin_id)
  #
  # final_out_us<-split(final_us_paths[,-c(1)],final_us_paths$source_id)

  return(ds_fp)

}

