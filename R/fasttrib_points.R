#' Quickly attribute stream segments/sampling points with layers of interest (loi)
#'
#' @param input output from `process_hydrology()` (if `process_loi()` was not run on `process_hydrology()`, `loi_file` must be specified)
#' @param loi_file filepath of `process_loi()` output (optional, will overwrite data in `process_hydrology()` output if present).
#' @param loi_cols character or NULL. Names of loi layers to include in summary. If NULL, all layers used.
#' @param sample_points character or NULL. IDs of unique station identifiers priveded in 'site_id_col' of `generate_vectors()`
#' @param link_id character or NULL. 'link_id's of reaches to calculate attributes for.
#' @param target_o_type character. One of: c("point","segment_point","segment_whole"). Target for iEucO" "iFLO", and "HAiFLO" weighting schemes. "Point" represents the sampling point on the stream, "segment_point" represents the upstream segment of the sampling points, and "segment_whole" will target the entire reach, regardless of where sampling occurred.
#' @param weighting_scheme character. One or more weighting schemes: c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS")
#' @param loi_numeric_stats character. One or more of c("mean", "sd", "median", "min", "max", "sum"). Only distance-weighted versions of mean and SD are returned for all weighting schemes except lumped.
#' @param approx_distwtdsd coarsly approximates weighted
#' @param inv_function function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.
#' @param out_filename Output file name.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#'
#' @return A data.frame of weighted attributes for the requested areas
#' @export
#'

fasttrib_points<-function(
    input,
    loi_file=NULL,
    loi_cols=NULL,
    sample_points=NULL,
    link_id=NULL,
    target_o_type=c("point","segment_point","segment_whole"),
    weighting_scheme =  c("lumped", "iEucS", "iFLS", "HAiFLS","iEucO","iFLO",  "HAiFLO"),
    loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
    approx_distwtdsd=F,
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    out_filename=NULL,
    # return_products=F, # This is not possible using this faster method
    temp_dir=NULL,
    verbose=F
){

  n_cores<-nbrOfWorkers()-1
  if (is.infinite(n_cores)) n_cores<-availableCores(logical = F)
  if (n_cores==0) n_cores<-1

  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)

  weighting_scheme_s<-weighting_scheme[grepl("FLS|iEucS",weighting_scheme)]
  weighting_scheme_o<-weighting_scheme[!grepl("lumped|FLS|iEucS",weighting_scheme)]
  lumped_scheme<-"lumped" %in% weighting_scheme
  if (length(weighting_scheme_o)>0) warning("Calculation for iEucO, iFLO, and HAiFLO are slow")

  if (is.null(target_o_type)) target_o_type<-"point"
  if (length(target_o_type)>1) target_o_type<-target_o_type[[1]]
  match.arg(target_o_type,several.ok = F)
  match.arg(weighting_scheme,several.ok = T)
  match.arg(loi_numeric_stats,several.ok = T)

  if (is.null(out_filename)) out_filename<-paste0("attrib_",floor(as.numeric(Sys.time())),".csv")
  if (length(out_filename)>1) out_filename<-out_filename[[1]]
  if (!grepl("\\.csv$",out_filename)) stop("'out_filename' must be a .csv file")

  loi_numeric_stats<-setNames(loi_numeric_stats,loi_numeric_stats)

  zip_loc<-input$outfile
  db_loc<-input$db_loc
  loi_loc<-loi_file
  if (is.null(loi_loc)) loi_loc<-zip_loc

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  wbt_options(exe_path=wbt_exe_path(),
              verbose=verbose,
              wd=temp_dir)

  terra::terraOptions(verbose = verbose,
                      tempdir = temp_dir
  )

  fl<-unzip(list=T,zip_loc)
  fl_loi<-unzip(list=T,loi_loc)

  if (verbose) print("Reading in data")

  site_id_col<-paste0(data.table::fread(cmd=paste("unzip -p ",zip_loc,"site_id_col.csv")))

  db_loc<-input$db_loc
  con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
  stream_links<-collect(tbl(con,"stream_links")) %>%
    mutate(across(c(link_id,any_of(site_id_col)),as.character)) %>%
    mutate(across(any_of(site_id_col),na_if,""))
  DBI::dbDisconnect(con)

  all_points<-read_sf(file.path("/vsizip",zip_loc,"stream_links.shp"))%>%
    mutate(across(c(link_id,any_of(site_id_col)),as.character)) %>% #stream_links.shp
    left_join(stream_links, by = c("link_id"))

  all_subb<-read_sf(file.path("/vsizip",zip_loc,"Subbasins_poly.shp"))
  # all_points<-read_sf(file.path("/vsizip",zip_loc,"stream_links.shp")) %>%
  #   left_join(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv")) %>%
  #               mutate(across(any_of(site_id_col),na_if,"")),
  #             by="link_id")
  all_catch<-read_sf(file.path("/vsizip",zip_loc,"Catchment_poly.shp"))

  # Get target link_id ------------------------------------------------------
  sample_points<-as.character(sample_points)
  link_id<-as.character(link_id)
  if (length(sample_points)==0 & length(link_id)==0) {
    warning("`sample_points` and `link_id` are NULL, all `link_id`s will evaluated")
    target_IDs<-all_points %>%
      as_tibble() %>%
      select(link_id,any_of(site_id_col))
  } else {
    if (site_id_col!="link_id" & length(sample_points)>0){
      target_IDs<-all_points %>%
        as_tibble() %>%
        select(link_id,any_of(site_id_col)) %>%
        filter(!!sym(site_id_col) %in% sample_points)
    } else {
      target_IDs<-NULL
    }

    if (length(link_id)>0){
      target_IDs<-bind_rows(
        target_IDs,
        all_points %>%
          as_tibble() %>%
          select(link_id,any_of(site_id_col)) %>%
          filter(link_id %in% link_id)
      )
    }
  }

  if (target_o_type=="segment_whole") {
    target_IDs<-target_IDs %>%
      select(link_id) %>%
      mutate(link_id=as.character(floor(as.numeric(link_id))))
  }

  target_IDs<-distinct(target_IDs)


  # Setup loi  --------------------------------------------------------------
  loi_rasts_exists<-c("num_rast.tif","cat_rast.tif")
  if (!any(loi_rasts_exists %in% fl_loi$Name)) stop("No 'loi' present in input, please run 'process_loi()' first, or specify location of process_loi() ouput")
  loi_rasts_exists<-fl_loi$Name[grepl("num_rast|cat_rast",fl_loi$Name)]
  loi_rasts_exists<-map(loi_rasts_exists,~file.path("/vsizip",loi_loc,.))
  names(loi_rasts_exists)<-gsub("\\.tif","",sapply(loi_rasts_exists,basename))

  loi_rasts<-map(loi_rasts_exists,rast)

  loi_rasts_comb<-rast(loi_rasts)

  names(loi_rasts_comb)<-unlist(sapply(loi_rasts,names))
  if (is.null(loi_cols)) loi_cols<-names(loi_rasts_comb)

  if (any(!loi_cols %in% names(loi_rasts_comb))) stop(paste0("The following `loi_cols` are not present in the `loi`:",
                                                             paste0(loi_cols[!loi_cols %in% names(loi_rasts_comb)],collapse = ", ")
  ))

  loi_rasts_comb<-terra::subset(loi_rasts_comb,loi_cols)
  loi_rasts_names<-map(loi_rasts,names)
  loi_rasts_names<-loi_rasts_names[loi_rasts_names %in% loi_cols]
  loi_rasts_names<-map(loi_rasts_names,~map(.,~setNames(as.list(.),.)) %>% unlist(recursive=T))
  loi_rasts_names<-map(loi_rasts_names,~map(.,~loi_numeric_stats))
  loi_rasts_names$cat_rast<-map(loi_rasts_names$cat_rast,~NULL)

  target_crs<-crs(vect(all_subb[1,]))

  # Get Upstream flowpaths --------------------------------------------------

  # unzip(zip_loc,files =c("flowpaths_out.db"),exdir=temp_dir)
  # db_loc<-file.path(temp_dir,"flowpaths_out.db")

  us_fp_fun<-function(link_id_in,db_loc=db_loc){
    con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
    out<-tbl(con,"us_flowpaths") %>%
      filter(pour_point_id %in% link_id_in) %>%
      rename(link_id=origin_link_id) %>%
      collect() %>%
      group_by(pour_point_id) %>%
      nest() %>%
      ungroup()

    out2<-out$data
    names(out2)<-out$pour_point_id

    out2<-out2[link_id_in]

    DBI::dbDisconnect(con)
    return(out2)


    # con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
    # out<-DBI::dbGetQuery(con, paste0("SELECT * FROM us_flowpaths WHERE source_id IN (",paste0(link_id,collapse = ","),")")) %>%
    #   group_by(source_id) %>%
    #   nest() %>%
    #   ungroup()
    #
    # out2<-out$data
    # names(out2)<-out$source_id
    #
    # out2<-out2[link_id]
    #
    # DBI::dbDisconnect(con)
    # return(out2)
  }

  # browser()
  us_flowpaths_out<-target_IDs %>%
    select(link_id) %>%
    mutate(link_id=as.character(link_id)) %>%
    mutate(us_flowpaths=us_fp_fun(link_id,db_loc=db_loc))

  # Select correct target for O -------------------------------------
  if (target_o_type=="point"){
    target_O<-all_points
  } else {
    if (target_o_type=="segment_point"){
      target_O<-read_sf(file.path("/vsizip",zip_loc,"stream_lines.shp"))
    } else {
      target_O<-read_sf(file.path("/vsizip",zip_loc,"stream_lines.shp")) %>%
        select(link_id) %>%
        mutate(link_id=as.character(floor(as.numeric(link_id)))) %>%
        group_by(link_id) %>%
        summarize(geometry=st_union(geometry)) %>%
        ungroup()


    }
  }

  #browser()

  # Sort everything by target_IDs
  target_O<-target_O[match(target_IDs[["link_id"]],target_O[["link_id"]],nomatch = 0),]
  all_points<-all_points[match(target_IDs[["link_id"]],all_points[["link_id"]],nomatch = 0),]
  all_catch<-all_catch[match(target_IDs[["link_id"]],all_catch[["link_id"]],nomatch = 0),]

  target_IDs<-target_IDs[match(target_O[["link_id"]],target_IDs[["link_id"]],nomatch = 0),]
  target_IDs<-target_IDs[match(all_points[["link_id"]],target_IDs[["link_id"]],nomatch = 0),]
  target_IDs<-target_IDs[match(all_catch[["link_id"]],target_IDs[["link_id"]],nomatch = 0),]

  target_O<-target_O[match(target_IDs[["link_id"]],target_O[["link_id"]],nomatch = 0),]
  all_points<-all_points[match(target_IDs[["link_id"]],all_points[["link_id"]],nomatch = 0),]
  all_catch<-all_catch[match(target_IDs[["link_id"]],all_catch[["link_id"]],nomatch = 0),]

  us_flowpaths_out<-us_flowpaths_out[match(target_IDs[["link_id"]],us_flowpaths_out[["link_id"]],nomatch = 0),]

  all_subb<-all_subb %>%
    filter(link_id %in% unlist(map(us_flowpaths_out$us_flowpaths,~.$link_id)))

  # Calculate Lumped Stats --------------------------------------------------

  if (lumped_scheme){
    #browser()

    custfun <- function(x, type) {
      switch(type,
             mean = mean(x,na.rm=T),
             median = median(x,na.rm=T),
             sd = mean(x,na.rm=T),
             min = min(x,na.rm=T),
             max = max(x,na.rm=T),
             sum = sum(x,na.rm=T)
      )
    }

    if (verbose) print("Calculating Lumped Attributes")

    lumped_numeric<-map(loi_numeric_stats,function(x){
      out<-terra::extract(loi_rasts$num_rast,all_catch,fun=x,na.rm=T,method="simple",touches=F,ID=F)
      names(out)<-paste0(names(out),"_lumped_",x)
      mutate(out,link_id=all_catch$link_id)
    })

    lumped_cat_part1<-left_join(us_flowpaths_out %>%
                                  transmute(link_id_upper=as.numeric(link_id),us_flowpaths=us_flowpaths) %>%
                                  unnest(us_flowpaths) %>%
                                  mutate(link_id=as.numeric(link_id)),
                                terra::extract(loi_rasts$cat_rast,all_subb,
                                               fun=function(x,na.rm=T) sum(x,na.rm=na.rm),
                                               na.rm=T,method="simple",touches=F,ID=F) %>%
                                  mutate(link_id=all_subb$link_id) %>%
                                  as_tibble(),
                                by="link_id"
    ) %>%
      group_by(link_id_upper) %>%
      summarize(
        across(c(everything(),-contains("link_id")),sum)
      ) %>%
      rename(link_id=link_id_upper)


    lumped_cat_part2<-terra::extract(loi_rasts$cat_rast[[1]],all_subb,
                                     fun=function(x,na.rm=T) length(x),
                                     na.rm=T,method="simple",touches=F,ID=F) %>%
      mutate(link_id=all_subb$link_id)%>%
      as_tibble() %>%
      left_join(us_flowpaths_out %>%
                  transmute(link_id_upper=as.numeric(link_id),us_flowpaths=us_flowpaths) %>%
                  unnest(us_flowpaths) %>%
                  mutate(link_id=as.numeric(link_id)),
                by="link_id"
      ) %>%
      group_by(link_id_upper) %>%
      summarize(
        across(c(everything(),-contains("link_id")),sum)
      ) %>%
      setNames(c("link_id","cell_sum"))

    lumped_cat<-lumped_cat_part1 %>%
      mutate(across(c(everything(),-contains("link_id")),~./lumped_cat_part2$cell_sum)) %>%
      dplyr::rename_with(.cols =c(everything(),-contains("link_id")),~paste0(.x,"_lumped_prop"))


  } else {
    lumped_numeric<-NULL
    lumped_cat<-NULL
  }


  # Calculate weighted distances -------------------------------------
  if (verbose) print("Generating Stream Targeted Weights")
  hw_streams<-hydroweight::hydroweight(hydroweight_dir=temp_dir,
                                       target_O = target_O,
                                       target_S = file.path("/vsizip",zip_loc,"dem_streams_d8.tif"),
                                       target_uid = 'ALL',
                                       OS_combine = FALSE,
                                       dem=file.path("/vsizip",zip_loc,"dem_final.tif"),
                                       flow_accum = file.path("/vsizip",zip_loc,"dem_accum_d8.tif"),
                                       weighting_scheme = weighting_scheme_s,
                                       inv_function = inv_function,
                                       clean_tempfiles=T,
                                       return_products = F,
                                       wrap_return_products=F,
                                       save_output=T)

  hw_streams_nm<-unzip(list=T,hw_streams)$Name
  # unzip(hw_streams,
  #       exdir = temp_dir)

  # dw_s_sum<-terra::extract(rast(hw_streams),all_catch,fun="sum",na.rm=T,method="simple",touches=F,ID=F) %>%
  #   mutate(link_id=all_catch$link_id)

  hw_streams_lo<-map(hw_streams_nm,function(x){
    file.path("/vsizip",hw_streams,x)
    #writeRaster(x,filename = file.path(temp_dir,paste0(names(x),".tif")),overwrite=T)
    #return(file.path(temp_dir,paste0(names(x),".tif")))
  })

  # Separate target_o into non-overlapping groups ---------------------------
  if (length(weighting_scheme_o)>0){
    if (verbose) print("Generating Site Targeted Weights")
    if (verbose) print("Unnesting Basins")

    temp_dir_sub<-file.path(temp_dir,basename(tempfile()))
    dir.create(temp_dir_sub)

    unzip(zip_loc,
          c("dem_d8.tif"),
          exdir=temp_dir_sub,
          overwrite=T,
          junkpaths=T)
    write_sf(all_points %>% select(link_id),
             file.path(temp_dir_sub,"pour_points.shp"),
             overwrite=T)

    #browser()

    future_unnest<-future::future({
      wbt_unnest_basins(
        wd=temp_dir_sub,
        d8_pntr="dem_d8.tif",
        pour_pts="pour_points.shp",
        output="unnest.tif"
      )
    })

    future_unnest_status <- future::futureOf(future_unnest)

    rast_out<-list()
    while(!future::resolved(future_unnest_status)){
      Sys.sleep(0.2)
      fl_un<-list.files(temp_dir_sub,"unnest_",full.names = T)

      if (length(fl_un)==0) next

      rast_all<-map(fl_un,function(x) try(rast(x),silent=T))
      rast_all<-rast_all[!sapply(rast_all,function(x) inherits(x,"try-error"))]

      if (length(rast_all)>0){
        rast_out<-c(rast_out,map(rast_all,terra::unique))
        suppressWarnings(file.remove(unlist(map(rast_all,terra::sources))))
      }
    }

    fl_un<-list.files(temp_dir_sub,"unnest_",full.names = T)
    rast_all<-map(fl_un,function(x) try(rast(x),silent=T))
    rast_all<-rast_all[!sapply(rast_all,function(x) inherits(x,"try-error"))]
    if (length(rast_all)>0){
      rast_out<-c(rast_out,map(rast_all,terra::unique))
      suppressWarnings(file.remove(unlist(map(rast_all,terra::sources))))
    }

    #browser()
    target_O_sub<-map(rast_out,~target_O[unlist(.x),] %>% select(link_id))

  } else {
    target_O_sub<-NULL
  }

  # Calculate weighted loi --------------------------------------------------
  if (verbose) print("Calculating Weighted Attributes")

  #browser()

  # if (n_cores==1|
  #     length(target_O_sub)<n_cores |
  #     length(hw_streams_lo)<n_cores){
  #
  #   #spt1<-which.max(sapply(target_O_sub,nrow))
  #
  # }
  #
  max_ittr<-length(target_O_sub)+length(hw_streams_lo)

  if (n_cores==1) {
    splt<-length(hw_streams_lo)
  } else {
    splt<-n_cores
  }

  #browser()

  with_progress(enable=T,{
    p <- progressor(steps = max_ittr)

    #loi_dw_out<-pmap(
    loi_dw_out<-future_pmap(
      .options = furrr_options(globals = FALSE),
      list(
        target_O_subs=suppressWarnings(split(target_O_sub,1:splt)),
        all_catch=rep(list(all_catch),splt),
        hw=suppressWarnings(rev(split(unlist(hw_streams_lo),1:splt))),
        a_subb=rep(list(all_subb),splt),
        us_flowpaths=rep(list(us_flowpaths_out),splt),
        loi_rasts_exists=rep(list(loi_rasts_exists),splt),
        loi_cols_sub=rep(list(loi_cols),splt),
        weighting_scheme_o=rep(list(weighting_scheme_o),splt),
        loi_numeric_stats=rep(list(loi_numeric_stats),splt),
        approx_distwtdsd=rep(list(approx_distwtdsd),splt),
        inv_function=rep(list(inv_function),splt),
        temp_dir=rep(list(temp_dir),splt),
        zip_loc=rep(list(zip_loc),splt),
        p=rep(list(p),splt)
      ),
      carrier::crate(function(target_O_subs,
                              all_catch,
                              hw,
                              a_subb,
                              us_flowpaths,
                              loi_rasts_exists,
                              loi_cols_sub,
                              weighting_scheme_o,
                              loi_numeric_stats,
                              approx_distwtdsd,
                              inv_function,
                              temp_dir,
                              zip_loc,
                              p
      ) {
        if (is.null(hw) & is.null(target_O_subs)) return(NULL)
        options(scipen = 999)
        #browser()
        `%>%` <- magrittr::`%>%`

        temp_dir_sub<-file.path(temp_dir,basename(tempfile()))

        # Read in LOI -------------------------------------------------------------
        loi_rasts<-purrr::map(loi_rasts_exists,terra::rast)
        loi_rasts_nms<-purrr::map(loi_rasts,names)

        loi_rasts_nms<-purrr::map(loi_rasts_nms,~.[.%in%loi_cols_sub])

        loi_rasts_comb<-terra::rast(loi_rasts)
        names(loi_rasts_comb)<-unlist(sapply(loi_rasts,names))
        loi_rasts_comb<-terra::subset(loi_rasts_comb,loi_cols_sub)

        # Perform Stream Weighting ------------------------------------------------

        if (length(hw)>0 & any(c("mean","sd") %in% loi_numeric_stats)) {
          sumfun<-function(x,us_flowpaths=us_flowpaths){
            out<-purrr::map_dfr(us_flowpaths$us_flowpaths,
                                ~x[x$link_id %in% .$link_id,] %>%
                                  dplyr::select(-tidyselect::any_of("link_id")) %>%
                                  dplyr::summarise(dplyr::across(tidyselect::everything(),sum))
                                ,.id = "link_id")

            out<-out[match(us_flowpaths$link_id, out$link_id,nomatch = 0),]
            out<-dplyr::select(out,-link_id)
            return(out)
          }

          #browser()
          hw2<-terra::rast(hw[[1]])

          loi_dist <- loi_rasts_comb * hw2

          if (any(c("mean") %in% loi_numeric_stats) | length(loi_rasts_nms$cat_rast)>0){
            # distance_weighted mean
            dw_s_sum<-terra::extract(hw2,a_subb,fun="sum",na.rm=T,method="simple",touches=F,ID=F) %>%
              dplyr::mutate(link_id=a_subb$link_id)
            dw_s_sum<-sumfun(dw_s_sum,us_flowpaths)

            loi_dw_s_sum<-terra::extract(loi_dist,a_subb,fun="sum",na.rm=T,method="simple",touches=F,ID=F)%>%
              dplyr::mutate(link_id=a_subb$link_id)
            loi_dw_s_sum<-sumfun(loi_dw_s_sum,us_flowpaths)

            s_dwmean<-loi_dw_s_sum/unlist(dw_s_sum)

            colnames(s_dwmean)[colnames(s_dwmean) %in% loi_rasts_nms$num_rast]<-paste0(
              colnames(s_dwmean)[colnames(s_dwmean) %in% loi_rasts_nms$num_rast],
              "_",names(hw2),"_mean"
            )

            colnames(s_dwmean)[colnames(s_dwmean) %in% loi_rasts_nms$cat_rast]<-paste0(
              colnames(s_dwmean)[colnames(s_dwmean) %in% loi_rasts_nms$cat_rast],
              "_",names(hw2),"_prop"
            )
          } else {
            s_dwmean<-NULL
          }


          # distance_weighted SD
          if (length(loi_rasts_nms$num_rast)>0 & any(c("sd") %in% loi_numeric_stats)){

            distwtd_sum<-terra::extract(hw2,
                                        a_subb,
                                        fun="sum",na.rm=T,method="simple",touches=F,ID=F)%>%
              dplyr::mutate(link_id=a_subb$link_id)
            distwtd_sum<-sumfun(distwtd_sum,us_flowpaths)

            distwtd_M<-terra::extract(hw2,
                                      a_subb ,
                                      fun=function(x) sum(x!=0,na.rm=T),
                                      method="simple",touches=F,ID=F)%>%
              dplyr::mutate(link_id=a_subb$link_id)
            distwtd_M<-sumfun(distwtd_M,us_flowpaths)

            term2<-((distwtd_M - 1) / distwtd_M) * distwtd_sum

            if (approx_distwtdsd) {
              # this is an incorrect approximation, but is much faster
              term1<-terra::extract(terra::subset(loi_dist,loi_rasts_nms$num_rast),
                                    a_subb ,
                                    fun=function(x) sum((x-mean(x,na.rm=T))^2,na.rm=T) ,
                                    method="simple",touches=F,ID=F)%>%
                dplyr::mutate(link_id=a_subb$link_id)

              term1<-sumfun(term1,us_flowpaths)

            } else {
              # this is the correct way, but slower
              #browser()
              t0<-terra::extract(hw2,
                                 a_subb,
                                 fun=NULL ,
                                 method="simple",touches=F,ID=T,cells=T)

              t1<-terra::extract(terra::subset(loi_dist,loi_rasts_nms$num_rast),
                                 a_subb,
                                 fun=NULL,
                                 method="simple",touches=F,ID=T,cells=T)

              # t2<-terra::extract(terra::subset(loi_dist,loi_rasts_nms$num_rast)*hw2,
              #                    a_subb,
              #                    fun=NULL ,
              #                    method="simple",touches=F,ID=T,cells=T)

              t2<-t0 %>%
                dplyr::left_join(t1,by=c("ID","cell")) %>%
                dplyr::group_by(ID) %>%
                dplyr::mutate(dplyr::across(tidyselect::any_of(loi_rasts_nms$num_rast), ~sum(.x,na.rm=T),.names = "sum_{.col}")) %>%
                dplyr::mutate(dplyr::across(tidyselect::any_of(names(hw2)), ~sum(.x,na.rm=T),.names = "dw_tbl")) %>%
                dplyr::mutate(dplyr::across(tidyselect::any_of(paste0("sum_",loi_rasts_nms$num_rast)), ~.x/dw_tbl,.names = "mean_{.col}")) %>%
                dplyr::ungroup()

              t3<-purrr:::map_dfc(loi_rasts_nms$num_rast,function(x){
                t2[,names(dw)]*((t2[,x]-t2[,paste0("mean_sum_",x)])^2)
              })
              colnames(t3)<-loi_rasts_nms$num_rast
              t3$ID<-t2$ID

              # t2<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast)*dw,
              #                    all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
              #                    fun=NULL ,
              #                    method="simple",touches=F,ID=T,cells=T)
              #
              # t2<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast)*dw,
              #                    all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
              #                    fun=NULL ,
              #                    method="simple",touches=F,ID=T,cells=T) %>%
              #   dplyr::group_by(ID) %>%
              #   dplyr::mutate(dplyr::across(tidyselect::any_of(loi_rasts_nms$num_rast),mean)) %>%
              #   dplyr::ungroup()

              # t3<-(t0[[names(dw)]]*(t1-t2)^2)
              # t3$ID<-t2$ID
              term1<-t3 %>%
                dplyr::group_by(ID) %>%
                dplyr::summarise(dplyr::across(tidyselect::everything(),sum,na.rm=T)) %>%
                dplyr::ungroup() %>%
                dplyr::select(-tidyselect::any_of("ID"),-tidyselect::any_of("cell"))%>%
                dplyr::mutate(link_id=a_subb$link_id)

              # t3<-(t0[[names(hw2)]]*(t1-t2)^2)
              # t3$ID<-t2$ID
              # term1<-t3 %>%
              #   dplyr::group_by(ID) %>%
              #   dplyr::summarise(dplyr::across(tidyselect::everything(),sum,na.rm=T)) %>%
              #   dplyr::ungroup() %>%
              #   dplyr::select(-tidyselect::any_of("ID"),-tidyselect::any_of("cell"))%>%
              #   dplyr::mutate(link_id=a_subb$link_id)

              term1<-sumfun(term1,us_flowpaths)

            }

            s_dwSD <- suppressMessages(purrr::map_dfc(colnames(term1), ~tibble::tibble(a=sqrt(term1[[.]] / unlist(term2)))))
            colnames(s_dwSD)<-paste0(colnames(term1),"_",names(hw2),"_sd")
          } else {
            s_dwSD<-NULL
          }

          s_out<-dplyr::bind_cols(
            s_dwmean,
            s_dwSD
          ) %>%
            dplyr::mutate(link_id=us_flowpaths$link_id)

          p()
        } else {
          s_out<-NULL
        }

        # Perform Target Weighting ------------------------------------------------

        if (length(target_O_subs)>0){
          o_out<-purrr::map(target_O_subs,function(x){
            hw<-hydroweight::hydroweight(hydroweight_dir=temp_dir_sub,
                                         target_O = x,
                                         target_S = file.path("/vsizip",zip_loc,"dem_streams_d8.tif"),
                                         target_uid = basename(tempfile()),
                                         OS_combine = FALSE,
                                         dem=file.path("/vsizip",zip_loc,"dem_final.tif"),
                                         flow_accum = file.path("/vsizip",zip_loc,"dem_accum_d8.tif"),
                                         weighting_scheme = weighting_scheme_o,
                                         inv_function = inv_function,
                                         clean_tempfiles=T,
                                         return_products = T,
                                         wrap_return_products=F,
                                         save_output=F)

            if (any(c("mean") %in% loi_numeric_stats) | length(loi_rasts_nms$cat_rast)>0){
              # distance_weighted mean
              dw_o_sum<-purrr::map(hw,~terra::extract(.,all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
                                                      fun="sum",na.rm=T,method="simple",touches=F,ID=F))

              loi_dw_o_sum<-purrr::map2(hw,names(hw),function(hhw,nm){
                terra::extract(loi_rasts_comb*hhw,all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
                               fun="sum",na.rm=T,method="simple",touches=F,ID=F)  #%>%
                #dplyr::rename_with(~paste0(.,"_",nm))
              })

              o_dwmean<-purrr::map2(loi_dw_o_sum,dw_o_sum,~(.x/unlist(.y)))

              o_dwmean<-purrr::map2_dfc(o_dwmean,names(o_dwmean),function(xx,yy){
                colnames(xx)[colnames(xx) %in% loi_rasts_nms$num_rast]<-paste0(
                  colnames(xx)[colnames(xx) %in% loi_rasts_nms$num_rast],
                  "_",yy,"_mean"
                )

                colnames(xx)[colnames(xx) %in% loi_rasts_nms$cat_rast]<-paste0(
                  colnames(xx)[colnames(xx) %in% loi_rasts_nms$cat_rast],
                  "_",yy,"_prop"
                )

                return(xx)
              })

            } else {
              o_dwmean<-NULL
            }

            if (length(loi_rasts_nms$num_rast)>0 & any(c("sd") %in% loi_numeric_stats)){
              #browser()
              o_dwSD<-purrr::map(hw,function(dw){
                #browser()
                distwtd_sum<-terra::extract(dw,
                                            all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
                                            fun="sum",na.rm=T,method="simple",touches=F,ID=F)
                distwtd_M<-terra::extract(dw,
                                          all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
                                          fun=function(x) sum(x!=0,na.rm=T),
                                          method="simple",touches=F,ID=F)

                term2<-((distwtd_M - 1) / distwtd_M) * distwtd_sum

                if (approx_distwtdsd){

                  # this is an incorrect approximation, but is much faster
                  term1<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast)*dw,
                                        all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
                                        fun=function(x) sum((x-mean(x,na.rm=T))^2,na.rm=T) ,
                                        method="simple",touches=F,ID=F)
                } else {
                  # this is the correct way, but slower
                  t0<-terra::extract(dw,
                                     all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
                                     fun=NULL ,
                                     method="simple",touches=F,ID=T,cells=T)

                  t1<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast),
                                     all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
                                     fun=NULL ,
                                     method="simple",touches=F,ID=T,cells=T)

                  t2<-t0 %>%
                    dplyr::left_join(t1,by=c("ID","cell")) %>%
                    dplyr::group_by(ID) %>%
                    dplyr::mutate(dplyr::across(tidyselect::any_of(loi_rasts_nms$num_rast), ~sum(.x,na.rm=T),.names = "sum_{.col}")) %>%
                    dplyr::mutate(dplyr::across(tidyselect::any_of(names(dw)), ~sum(.x,na.rm=T),.names = "dw_tbl")) %>%
                    dplyr::mutate(dplyr::across(tidyselect::any_of(paste0("sum_",loi_rasts_nms$num_rast)), ~.x/dw_tbl,.names = "mean_{.col}")) %>%
                    dplyr::ungroup()

                  t3<-purrr:::map_dfc(loi_rasts_nms$num_rast,function(x){
                    t2[,names(dw)]*((t2[,x]-t2[,paste0("mean_sum_",x)])^2)
                  })
                  colnames(t3)<-loi_rasts_nms$num_rast
                  t3$ID<-t2$ID

                  # t2<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast)*dw,
                  #                    all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
                  #                    fun=NULL ,
                  #                    method="simple",touches=F,ID=T,cells=T)
                  #
                  # t2<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast)*dw,
                  #                    all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
                  #                    fun=NULL ,
                  #                    method="simple",touches=F,ID=T,cells=T) %>%
                  #   dplyr::group_by(ID) %>%
                  #   dplyr::mutate(dplyr::across(tidyselect::any_of(loi_rasts_nms$num_rast),mean)) %>%
                  #   dplyr::ungroup()

                  # t3<-(t0[[names(dw)]]*(t1-t2)^2)
                  # t3$ID<-t2$ID
                  term1<-t3 %>%
                    dplyr::group_by(ID) %>%
                    dplyr::summarise(dplyr::across(tidyselect::everything(),sum,na.rm=T)) %>%
                    dplyr::ungroup() %>%
                    dplyr::select(-tidyselect::any_of("ID"),-tidyselect::any_of("cell"))
                }

                loi_distwtd_sd <- suppressMessages(purrr::map_dfc(colnames(term1), ~tibble::tibble(a=sqrt(term1[[.]] / unlist(term2)))))
                colnames(loi_distwtd_sd)<-paste0(colnames(term1),"_",names(dw),"_sd")

                return(loi_distwtd_sd)

              })
            } else {
              o_dwSD<-NULL
            }
            # distance weighted SD is costly to calculate like this

            o_out<-dplyr::bind_cols(
              o_dwmean,
              o_dwSD
            ) %>%
              dplyr::mutate(link_id=all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]) %>% dplyr::pull(link_id))

            p()
            return(o_out)

          })

          o_out<-dplyr::bind_rows(o_out)
        } else {
          o_out<-NULL
        }


        return(
          list(
            s_out=s_out,
            o_out=o_out
          )
        )

      }
      ))
  })

  if (inherits(lumped_numeric,"list")){
    lumped_numeric<-reduce(lumped_numeric,left_join,by="link_id")
  }
  if (inherits(lumped_cat,"list")){
    lumped_cat<-reduce(lumped_cat,left_join,by="link_id")
  }

  s_out<-map(loi_dw_out,~.$s_out)
  s_out<-s_out[!sapply(s_out,is.null)]
  s_out<-reduce(s_out,left_join)

  o_out<-map_dfr(loi_dw_out,~.$o_out)

  final_out<-target_IDs %>%
    mutate(link_id=as.character(link_id)) %>%
    left_join(lumped_numeric %>% mutate(link_id=as.character(link_id)),by="link_id") %>%
    left_join(lumped_cat %>% mutate(link_id=as.character(link_id)),by="link_id") %>%
    left_join(s_out %>% mutate(link_id=as.character(link_id)),by="link_id") %>%
    left_join(o_out %>% mutate(link_id=as.character(link_id)),by="link_id")

  data.table::fwrite(final_out,file.path(temp_dir,out_filename))

  zip(zip_loc,
      file.path(temp_dir,out_filename),
      flags = '-r9Xjq'
  )

  # output<-input
  #
  # output[[gsub(".csv","",out_filename)]]<-final_out

  suppressWarnings(file.remove(list.files(temp_dir,full.names = T,recursive = T)))

  return(final_out)
}
