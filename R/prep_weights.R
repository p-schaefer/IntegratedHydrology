
#' Prepares hydroweights and adds them to existing zip file
#'
#' @param input output from `process_hydrology()` (if `process_loi()` was not run on `process_hydrology()`, `loi_file` must be specified)
#' @param sample_points character or NULL. IDs of unique station identifiers priveded in 'site_id_col' of `generate_vectors()`
#' @param link_id character or NULL. 'link_id's of reaches to calculate attributes for.
#' @param target_o_type character. One of: c("point","segment_point","segment_whole"). Target for iEucO" "iFLO", and "HAiFLO" weighting schemes. "Point" represents the sampling point on the stream, "segment_point" represents the upstream segment of the sampling points, and "segment_whole" will target the entire reach, regardless of where sampling occurred.
#' @param weighting_scheme character. One or more weighting schemes: c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS")
#' @param inv_function function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#'
#' @return input
#' @export
#'

prep_weights<-function(
    input,
    sample_points=NULL,
    link_id=NULL,
    target_o_type=c("point","segment_point","segment_whole"),
    weighting_scheme =  c("iEucS", "iFLS", "HAiFLS","iEucO","iFLO",  "HAiFLO"),
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    temp_dir=NULL,
    verbose=F
){

  n_cores<-nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-availableCores(logical = F)
  if (n_cores==0) n_cores<-1

  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)

  weighting_scheme_s<-weighting_scheme[grepl("FLS|iEucS",weighting_scheme)]
  weighting_scheme_o<-weighting_scheme[!grepl("lumped|FLS|iEucS",weighting_scheme)]
  if (length(weighting_scheme_o)>0) warning("Calculation for iEucO, iFLO, and HAiFLO are slow")

  if (is.null(target_o_type)) target_o_type<-"point"
  if (length(target_o_type)>1) target_o_type<-target_o_type[[1]]
  match.arg(target_o_type,several.ok = F)
  match.arg(weighting_scheme,several.ok = T)

  zip_loc<-input$outfile

  out_zip_loc<-zip_loc
  out_zip_loc<-file.path(gsub(basename(out_zip_loc),"",out_zip_loc),paste0(gsub(".zip","_DW.zip",basename(out_zip_loc))))
  if (file.exists(out_zip_loc)) file.remove(out_zip_loc)

  db_loc<-input$db_loc

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

  # Get Upstream flowpaths --------------------------------------------------

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


  # Calculate weighted S-target distances -------------------------------------

  temp_dir_sub<-file.path(temp_dir,basename(tempfile()))
  dir.create(temp_dir_sub)

  if (verbose) print("Generating Stream Targeted Weights")
  hw_streams<-hydroweight::hydroweight(hydroweight_dir=temp_dir_sub,
                                       target_O = NULL,
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

  uz_fls<-unzip(list=T,hw_streams)$Name
  unzip(hw_streams)

  rout<-sapply(uz_fls,function(x) file.path(temp_dir_sub,x))

  zip(out_zip_loc,
      unlist(rout),
      flags = '-r9Xjq'
  )

  t1<-try(file.remove(unlist(rout)),silent=T)
  t1<-try(file.remove(file.path(temp_dir_sub,paste0("All_inv_distances.zip"))),silent=T)

  # Calculate weighted O-target distances -------------------------------------

  unzip(zip_loc,
        c("dem_d8.tif"),
        exdir=temp_dir_sub,
        overwrite=T,
        junkpaths=T)
  write_sf(all_points %>% select(link_id),
           file.path(temp_dir_sub,"pour_points.shp"),
           overwrite=T)

  # Use unnest basins to find catchments that don't overlap
  # Use asynchronous evaluation (if parallel backend registered)
  # to clear up rasters as they are made
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

  target_O_sub<-map2(rast_out,seq_along(rast_out),~target_O[unlist(.x),] %>% select(link_id) %>% mutate(unn_group=.y))
  splt_lst<-suppressWarnings(split(target_O_sub,1:n_cores))
  splt_lst<-splt_lst[!sapply(splt_lst,is.null)]

  #browser()

  write_sf(bind_rows(target_O_sub),file.path(temp_dir_sub,"unnest_group_target_O.shp"))

  zip(out_zip_loc,
      list.files(temp_dir_sub,"unnest_group_target_O",full.names = T),
      flags = '-r9Xjq'
  )

  hw_o_targ<-future_map(splt_lst,function(x){
    target_S <- terra::rast(file.path("/vsizip",zip_loc,"dem_streams_d8.tif"))
    dem <- terra::rast(file.path("/vsizip",zip_loc,"dem_final.tif"))
    flow_accum <- terra::rast(file.path("/vsizip",zip_loc,"dem_accum_d8.tif"))

    o_out<-purrr::map(x,function(y){
      hw_o<-hydroweight::hydroweight(hydroweight_dir=temp_dir_sub,
                                     target_O = y,
                                     target_S = target_S,
                                     target_uid = paste0("unnest_group_",y$unn_group[[1]]),
                                     OS_combine = FALSE,
                                     dem=dem,
                                     flow_accum = flow_accum,
                                     weighting_scheme = weighting_scheme_o,
                                     inv_function = inv_function,
                                     clean_tempfiles=T,
                                     return_products = F,
                                     wrap_return_products=F,
                                     save_output=T)

      uz_fls<-unzip(list=T,hw_streams)$Name
      unzip(hw_streams)

      rout<-sapply(uz_fls,function(x) file.path(temp_dir_sub,x))

      zip(out_zip_loc,
          unlist(rout),
          flags = '-r9Xjq'
      )

      t1<-try(file.remove(unlist(rout)),silent=T)
      t1<-try(file.remove(file.path(temp_dir_sub,paste0("unnest_group_",y$unn_group[[1]],"_inv_distances.zip"))),silent=T)

      return(NULL)
    })

  })

  tt<-file.remove(list.files(temp_dir_sub,full.names = T))

  input$dw_dir<-out_zip_loc

  return(input)

}

