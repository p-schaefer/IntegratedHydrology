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
#' @param inv_function function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.
#' @param use_exising_hw logical. Should the function look for existing hydroweight layers in the zip file?
#' @param use_existing_attr logical. Should the function look for existing attribute layers in the database file?
#' @param store_hw logical. Should hydroweight layer be stored and added to the zip file?
#' @param out_filename Output file name.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#'
#' @return A data.frame of weighted attributes for the requested areas
#' @export
#'


# t1<-final_attributes_slow %>%
#   arrange(site_id) %>%
#   select(-any_of("distance_weights"),-any_of("weighted_attr"),-any_of("link_id"),-any_of("site_id")) %>%
#   dplyr::rename_with(~gsub("distwtd_","",.x))
# t2<-final_attributes%>%
#   arrange(as.numeric(site_id)) %>%
#   select(-any_of("distance_weights"),-any_of("weighted_attr"),-any_of("link_id"),-any_of("site_id"))
#
# t1<-t1[,colnames(t2)]
#
# (t1)[,grepl("slope",colnames(t2))]
# (t2)[,grepl("slope",colnames(t2))]
#
# all_equal(
#   t1
#   ,
#   t2
# )

fasttrib_points<-function(
    input,
    loi_file=NULL,
    loi_cols=NULL,
    sample_points=NULL,
    link_id=NULL,
    target_o_type=c("point","segment_point","segment_whole"),
    weighting_scheme =  c("lumped", "iEucS", "iFLS", "HAiFLS","iEucO","iFLO",  "HAiFLO"),
    loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    use_exising_hw=F,
    use_existing_attr=F,
    store_hw=F,
    out_filename=NULL,
    # return_products=F, # This is not possible using this faster method
    temp_dir=NULL,
    verbose=F
){
  if (!is.logical(use_exising_hw)) stop("'use_exising_hw' must be logical")
  if (!is.logical(store_hw)) stop("'store_hw' must be logical")
  if (!is.logical(use_existing_attr)) stop("'use_existing_attr' must be logical")

  if (use_exising_hw & is.null(input$dw_dir)) stop("If 'use_exising_hw'=TRUE, input$dw_dir ust not be NULL. Has prep_weights() been run?")

  if (store_hw){
    use_exising_hw<-T

    input<-prep_weights(
      input=input,
      sample_points=sample_points,
      link_id=link_id,
      target_o_type=target_o_type,
      weighting_scheme =  weighting_scheme[!grepl("lumped",weighting_scheme)],
      inv_function = inv_function,
      temp_dir=temp_dir,
      verbose=verbose
    )
  }


  n_cores<-nbrOfWorkers()
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
  dw_dir<-input$dw_dir
  db_loc<-input$db_loc
  attr_db_loc<-gsub(basename(db_loc),gsub(".db","_Attr.db",basename(db_loc)),db_loc)

  if (!use_existing_attr){
    if (file.exists(attr_db_loc)) {
      file.remove(attr_db_loc)
    }

    file.copy(
      db_loc,
      attr_db_loc
    )
  } else {
    if (!file.exists(attr_db_loc)) stop(paste0(attr_db_loc," must exist if use_existing_attr=TRUE"))
  }



  con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc)

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

  if (!is.null(dw_dir)){
    fl_dw<-unzip(list=T,dw_dir)
  } else {
    fl_dw<-NULL
  }

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


  # Setup loi  --------------------------------------------------------------
  if (verbose) print("Reading in LOI")

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
  loi_rasts_names$cat_rast<-as.list(setNames(rep(NA,length(names(loi_rasts$cat_rast))),names(loi_rasts$cat_rast)))

  target_crs<-crs(vect(all_subb[1,]))


  # Upload loi rasters to attributes database -------------------------------
  if (!use_existing_attr){
    if (verbose) print("Writing LOI to attributes database")

    attrib_tbl<-exactextractr::exact_extract(
      loi_rasts_comb,
      all_subb,
      weights=NULL,
      include_cell=T,
      fun=NULL,
      include_cols="link_id",
      progress=F
    )%>%
      dplyr::bind_rows() %>%
      dplyr::select(-coverage_fraction) %>%
      stats::setNames(c("subb_link_id",names(loi_rasts_comb),"cell_number")) %>%
      select(cell_number,subb_link_id,everything()) %>%
      copy_to(df=.,
              con_attr,
              "attrib_tbl",
              overwrite =T,
              temporary =F,
              indexes=c("subb_link_id","cell_number"),
              analyze=T,
              in_transaction=T)
  }

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


  # Calculate weighted distances -------------------------------------
  if (!use_existing_attr){

    if (!use_exising_hw){
      if (verbose) print("Generating Stream Targeted Weights")
      hw_streams<-hydroweight::hydroweight(hydroweight_dir=temp_dir,
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

      hw_streams_nm<-unzip(list=T,hw_streams)$Name

      hw_streams_lo<-map(hw_streams_nm,function(x){
        file.path("/vsizip",hw_streams,x)
      })
    } else {
      #browser()
      trg_fl<-paste0("ALL_",weighting_scheme_s,"_inv_distances.tif")
      if (all(sapply(trg_fl,function(x) any(grepl(x,fl_dw$Name))))) {
        hw_streams_lo<-map(trg_fl,~file.path("/vsizip",dw_dir,.))
      } else {
        stop(paste0("Not all 'weighting_scheme' found in zip file"))
      }
    }

    if (verbose) print("Writing S-targeted weights to attributes database")

    hw2<-map(hw_streams_lo,terra::rast)
    names(hw2)<-sapply(hw2,names)

    s_trg_weights<-exactextractr::exact_extract(
      terra::rast(hw2),
      all_subb,
      weights=NULL,
      include_cell=T,
      fun=NULL,
      include_cols="link_id",
      progress=F
    )%>%
      dplyr::bind_rows() %>%
      dplyr::select(-coverage_fraction) %>%
      stats::setNames(c("subb_link_id",names(hw2),"cell_number")) %>%
      select(cell_number,subb_link_id,everything()) %>%
      copy_to(df=.,
              con_attr,
              "s_target_weights",
              overwrite =T,
              temporary =F,
              indexes=c("subb_link_id","cell_number"),
              analyze=T,
              in_transaction=T)
  }

  # Separate target_o into non-overlapping groups ---------------------------
  if (length(weighting_scheme_o)>0){
    if (!use_existing_attr){
      if (!use_exising_hw){

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
        target_O_sub<-map2(rast_out,seq_along(rast_out),~target_O[unlist(.x),] %>% select(link_id) %>% mutate(unn_group=.y))
      } else {

        #browser()
        target_O_sub<-read_sf(file.path("/vsizip",dw_dir,"unnest_group_target_O.shp")) %>%
          split(.,paste0("unnest_group_",.$unn_group))


        trg_fls<-sapply(weighting_scheme_o,function(x) paste0(names(target_O_sub),"_",x,"_inv_distances.tif"))
        if (any(!trg_fls %in% fl_dw$Name)){
          stop("Some unnested target_O groups are missing hydroweights")
        }

      }

      # extract O target weights
      if (verbose) print("Writing O-targeted weights to attributes database")

      with_progress(enable=T,{

        p <- progressor(steps = length(target_O_sub))

        splt<-1

        o_trg_weights<-tibble(catch_link_id="1.1",
                              cell_number=1L,
        ) %>%
          bind_cols(
            data.frame(matrix(ncol=length(weighting_scheme_o),nrow=1)) %>%
              setNames(weighting_scheme_o)
          ) %>%
          mutate(across(any_of(weighting_scheme_o),~1.1)) %>%
          .[F,]%>%
          copy_to(df=.,
                  con_attr,
                  "o_target_weights",
                  overwrite =T,
                  temporary =F,
                  indexes=c("catch_link_id","cell_number"),
                  analyze=T,
                  in_transaction=T)

        loi_dw_out<-pmap( # I don't think this can be parallel
          #loi_dw_out<-future_pmap_dfr(
          #  .options = furrr_options(globals = FALSE),
          list(
            target_O_subs=suppressWarnings(split(target_O_sub,1:splt)),
            weighting_scheme_o=rep(list(weighting_scheme_o),splt),
            all_catch=rep(list(all_catch),splt),
            inv_function=rep(list(inv_function),splt),
            use_exising_hw=rep(list(use_exising_hw),splt),
            temp_dir=rep(list(temp_dir),splt),
            con_attr_l=rep(list(con_attr),splt),
            new_tbl=rep(list(o_trg_weights),splt),
            zip_loc=rep(list(zip_loc),splt),
            dw_dir=rep(list(dw_dir),splt),
            p=rep(list(p),splt)
          ),
          carrier::crate(function(target_O_subs,
                                  weighting_scheme_o,
                                  all_catch,
                                  inv_function,
                                  use_exising_hw,
                                  temp_dir,
                                  con_attr_l,
                                  new_tbl,
                                  zip_loc,
                                  dw_dir,
                                  p
          ) {
            options(scipen = 999)
            `%>%` <- magrittr::`%>%`

            temp_dir_sub<-file.path(temp_dir,basename(tempfile()))

            target_S <- file.path("/vsizip",zip_loc,"dem_streams_d8.tif")
            dem <- file.path("/vsizip",zip_loc,"dem_final.tif")
            flow_accum <- file.path("/vsizip",zip_loc,"dem_accum_d8.tif")

            o_out<-purrr::map(target_O_subs, # I don't think this can be parallel
                              function(x){
                                #print(x$unn_group[[1]])
                                sub_catch<-all_catch %>%
                                  dplyr::filter(link_id %in% x$link_id)

                                if (!use_exising_hw){
                                  hw<-hydroweight::hydroweight(hydroweight_dir=temp_dir_sub,
                                                               target_O = x,
                                                               target_S = target_S,
                                                               target_uid = basename(tempfile()),
                                                               OS_combine = FALSE,
                                                               dem=dem,
                                                               flow_accum = flow_accum,
                                                               weighting_scheme = weighting_scheme_o,
                                                               inv_function = inv_function,
                                                               clean_tempfiles=T,
                                                               return_products = T,
                                                               wrap_return_products=F,
                                                               save_output=F)
                                } else {
                                  trg_fl<-paste0("unnest_group_",x$unn_group[[1]],"_",weighting_scheme_o,"_inv_distances.tif")
                                  hw<-purrr::map(trg_fl,~terra::rast(file.path("/vsizip",dw_dir,.)))
                                  names(hw)<-sapply(hw,names)
                                }

                                out<-exactextractr::exact_extract(
                                  terra::rast(hw),
                                  sub_catch,
                                  weights=NULL,
                                  #fun="sum",
                                  include_cell=T,
                                  fun=NULL,
                                  include_cols="link_id",
                                  progress=F
                                )%>%
                                  dplyr::bind_rows() %>%
                                  dplyr::select(-coverage_fraction) %>%
                                  stats::setNames(c("catch_link_id",names(hw),"cell_number")) %>%
                                  dplyr::rows_insert(y=.,
                                                     x=new_tbl,
                                                     by="catch_link_id",
                                                     copy=T,
                                                     in_place=T,
                                                     conflict = "ignore")

                                p()
                                return(NULL)
                              })

          }))
      })

    }
  } else {
    target_O_sub<-NULL
  }

  DBI::dbDisconnect(con_attr)

  # Calculate Attributes ----------------------------------------------------
  # o_trg_weights
  # s_trg_weights
  # attrib_tbl

  lumped_out<-NULL
  s_targ_out<-NULL
  o_targ_out<-NULL

  #browser()
  # Lumped Attributes -------------------------------------------------------

  if (lumped_scheme) {
    if (verbose) print("Calculating Lumped Attributes")

    with_progress(enable=T,{

      p <- progressor(steps = nrow(us_flowpaths_out))
      lumped_out<-us_flowpaths_out %>%
        dplyr::mutate(attr=furrr::future_pmap_dfr(
          list(
            link_id_in=link_id,
            attr_db_loc=list(attr_db_loc),
            loi_rasts_names=list(loi_rasts_names),
            p=list(p)
          ),
          .options = furrr_options(globals = FALSE),
          carrier::crate(
            function(link_id_in,
                     attr_db_loc,
                     loi_rasts_names,
                     p
            ){
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`

              con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc)

              out<-dplyr::semi_join(
                dplyr::tbl(con_attr,"attrib_tbl") %>%
                  dplyr::rename(link_id=subb_link_id),
                dplyr::tbl(con_attr,"us_flowpaths") %>%
                  dplyr::filter(pour_point_id %in% link_id_in) %>%
                  dplyr::rename(link_id=origin_link_id),
                by="link_id"
              )

              attrs<-sapply(sapply(loi_rasts_names$num_rast,unique),unique)

              mean_out<-NULL
              sd_out<-NULL
              min_out<-NULL
              max_out<-NULL
              count_out<-NULL
              median_out<-NULL
              sum_out<-NULL

              if (any("mean"==attrs)|length(loi_rasts_names$cat_rast) >0){
                mean_out<-out %>%
                  dplyr::select(-link_id,-cell_number) %>%
                  dplyr::summarise(dplyr::across(tidyselect::everything(),~sum(.,na.rm=T)/dplyr::n())) %>%
                  dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_mean")) %>%
                  dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$cat_rast)),~paste0(.x,"_lumped_prop")) %>%
                  dplyr::collect()
              }

              if (any(attrs %in% c("sd","stdev"))){
                sd_out<-out %>%
                  dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~sd(.,na.rm=T))) %>%
                  dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_sd"))%>%
                  dplyr::collect()

              }
              if (any(attrs=="min")){
                min_out<-out %>%
                  dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~min(.,na.rm=T))) %>%
                  dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_min"))%>%
                  dplyr::collect()

              }
              if (any(attrs=="max")){
                max_out<-out %>%
                  dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~max(.,na.rm=T)))%>%
                  dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_max")) %>%
                  dplyr::collect()

              }
              if (any(attrs=="count")){
                count_out<-out %>%
                  dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~sum(.[!is.na(.)],na.rm=T)))%>%
                  dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_count")) %>%
                  dplyr::collect()

              }
              if (any(attrs=="median")){
                median_out<-out %>%
                  dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~median(.,na.rm=T)))%>%
                  dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_median")) %>%
                  dplyr::collect()

              }
              if (any(attrs=="sum")){
                sum_out<-out %>%
                  dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~sum(.,na.rm=T)))%>%
                  dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_sum")) %>%
                  dplyr::collect()
              }

              final_out<-dplyr::bind_cols(
                list(
                  mean_out,
                  sd_out,
                  min_out,
                  max_out,
                  count_out,
                  median_out,
                  sum_out
                )
              )


              DBI::dbDisconnect(con_attr)

              p()

              return(final_out)

            })
        )) %>%
        select(everything(),-us_flowpaths) %>%
        unnest(cols=attr)

    })
  }

  # s-targeted Attributes -------------------------------------------------------

  if (length(weighting_scheme_s)>0) {
    if (verbose) print("Calculating s-targeted Attributes")

    with_progress(enable=T,{

      p <- progressor(steps = nrow(us_flowpaths_out))
      s_targ_out<-us_flowpaths_out %>%
        dplyr::mutate(attr=future_pmap(
          list(
            link_id_in=link_id,
            attr_db_loc=list(attr_db_loc),
            loi_rasts_names=list(loi_rasts_names),
            weighting_scheme_s=list(weighting_scheme_s),
            p=list(p)
          ),
          .options = furrr_options(globals = FALSE),
          carrier::crate(
            function(link_id_in,
                     attr_db_loc,
                     loi_rasts_names,
                     weighting_scheme_s,
                     p
            ){
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`
              #browser()

              con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc)

              attr_nms<-names(c(loi_rasts_names$num_rast,loi_rasts_names$cat_rast))
              names(attr_nms)<-attr_nms

              names(weighting_scheme_s)<-weighting_scheme_s

              out<-dplyr::semi_join(
                dplyr::tbl(con_attr,"attrib_tbl") %>%
                  dplyr::rename(link_id=subb_link_id),
                dplyr::tbl(con_attr,"us_flowpaths") %>%
                  dplyr::filter(pour_point_id %in% link_id_in) %>%
                  dplyr::rename(link_id=origin_link_id),
                by="link_id"
              ) %>%
                dplyr::left_join(
                  dplyr::tbl(con_attr,"s_target_weights") %>%
                    dplyr::select(cell_number,tidyselect::any_of(weighting_scheme_s)),
                  by="cell_number"
                )

              if ("iFLS" %in% weighting_scheme_s){ # This is the only way I could get around an error by iterating over weighting_scheme_s
                out<-out %>%
                  dplyr::left_join(
                    out %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("iFLS")),.names="{.col}_iFLS" )) %>%
                      dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_","iFLS"))),

                    by = c("cell_number", "link_id")
                  )
              }

              if ("HAiFLS" %in% weighting_scheme_s){
                out<-out %>%
                  dplyr::left_join(
                    out %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("HAiFLS")),.names="{.col}_HAiFLS" )) %>%
                      dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_","HAiFLS"))),
                    by = c("cell_number", "link_id")
                  )
              }


              attrs<-sapply(sapply(loi_rasts_names$num_rast,unique),unique)
              mean_out<-NULL
              sd_out<-NULL

              if (any("mean"==attrs)|length(loi_rasts_names$cat_rast) >0){


                mean_out<-purrr::map_dfc(weighting_scheme_s,function(x){
                  out_iFLS<-NULL
                  out_HAiFLS<-NULL

                  if ("iFLS" %in% x){
                    out_iFLS<-out %>%
                      dplyr::summarize(dplyr::across(tidyselect::ends_with(paste0("_","iFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T) ))  %>%
                      dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$num_rast),"_")),~paste0(.x,"_mean")) %>%
                      dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$cat_rast),"_")),~paste0(.x,"_prop")) %>%
                      dplyr::collect()
                  }

                  if ("HAiFLS" %in% x){
                    out_iFLS<-out %>%
                      dplyr::summarize(dplyr::across(tidyselect::ends_with(paste0("_","HAiFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T) ))  %>%
                      dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$num_rast),"_")),~paste0(.x,"_mean")) %>%
                      dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$cat_rast),"_")),~paste0(.x,"_prop")) %>%
                      dplyr::collect()
                  }

                  dplyr::bind_cols(out_iFLS,out_HAiFLS)

                }) %>%
                  dplyr::mutate(link_id=link_id_in) %>%
                  dplyr::select(link_id,tidyselect::everything())
              }

              if (any(attrs %in% c("sd","stdev"))) {
                sd_out<-purrr::map_dfc(weighting_scheme_s,function(x){
                  out_iFLS<-NULL
                  out_HAiFLS<-NULL

                  if ("iFLS" %in% x){
                    out_iFLS<-out %>%
                      dplyr::select(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                    tidyselect::any_of("iFLS")
                      ) %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                                  ~(!!rlang::sym("iFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)))^2)),
                                                  .names="{.col}_iFLS_term1"),
                                    dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                                  ~ ((sum(!!rlang::sym("iFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLS")!=0,na.rm=T)) * sum(!!rlang::sym("iFLS"),na.rm=T),
                                                  .names="{.col}_iFLS_term2"
                                    ))%>%
                      dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),sum,na.rm=T),
                                       dplyr::across(tidyselect::ends_with("_term2"),~.[1])
                      ) %>%
                      dplyr::collect()
                  }

                  if ("HAiFLS" %in% x){
                    out_HAiFLS<-out %>%
                      dplyr::select(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                    tidyselect::any_of("HAiFLS")
                      ) %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                                  ~(!!rlang::sym("HAiFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)))^2)),
                                                  .names="{.col}_HAiFLS_term1"),
                                    dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                                  ~ ((sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLS"),na.rm=T),
                                                  .names="{.col}_HAiFLS_term2"
                                    ))%>%
                      dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),sum,na.rm=T),
                                       dplyr::across(tidyselect::ends_with("_term2"),~.[1])
                      ) %>%
                      dplyr::collect()
                  }


                  out_final<-dplyr::bind_cols(out_iFLS,out_HAiFLS)

                  out_final %>%
                    tidyr::pivot_longer(tidyselect::everything()) %>%
                    dplyr::mutate(attr=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_",2)[,1],
                                  term=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_",2)[,2]) %>%
                    dplyr::rowwise() %>%
                    dplyr::mutate(hw=gsub(paste0(attr,"_","|","_",term,""),"",name)) %>%
                    dplyr::ungroup() %>%
                    dplyr::mutate(name=paste0(attr,"_",hw,"_sd")) %>%
                    dplyr::group_by(name) %>%
                    dplyr::summarize(sd=sqrt(value[term=="term1"]/value[term=="term2"])) %>%
                    dplyr::ungroup() %>%
                    tidyr::pivot_wider(names_from = name,values_from=sd)

                  #browser()

                  # purrr::map_dfc(names(loi_rasts_names$num_rast),function(y){
                  #
                  #   out_iFLS<-NULL
                  #   out_HAiFLS<-NULL
                  #
                  #   if ("iFLS" %in% x){
                  #     out_iFLS<-out_final %>%
                  #       dplyr::select(tidyselect::contains(y)) %>%
                  #       dplyr::mutate(a= sqrt( !!rlang::sym(paste0(y,"_","iFLS","_term1")) / !!rlang::sym(paste0(y,"_","iFLS","_term2")) )) %>%
                  #       dplyr::select(a) %>%
                  #       dplyr::rename_with(~paste0(y,"_","iFLS","_sd"))
                  #   }
                  #
                  #   if ("HAiFLS" %in% x){
                  #     out_HAiFLS<-out_final %>%
                  #       dplyr::select(tidyselect::contains(y)) %>%
                  #       dplyr::mutate(a= sqrt( !!rlang::sym(paste0(y,"_","HAiFLS","_term1")) / !!rlang::sym(paste0(y,"_","HAiFLS","_term2")) )) %>%
                  #       dplyr::select(a) %>%
                  #       dplyr::rename_with(~paste0(y,"_","HAiFLS","_sd"))
                  #   }
                  #
                  #   dplyr::bind_cols(out_iFLS,out_HAiFLS)
                  #
                  # })


                })

              }

              final_out<-dplyr::bind_cols(
                list(
                  mean_out,
                  sd_out
                )
              ) %>%
                dplyr::select(-tidyselect::any_of("link_id"))

              p()

              DBI::dbDisconnect(con_attr)
              return(final_out)

            })
        )) %>%
        select(everything(),-us_flowpaths) %>%
        unnest(cols=attr)
    })
  }

  # o-targeted Attributes -------------------------------------------------------

  if (length(weighting_scheme_o)>0) {
    if (verbose) print("Calculating o-targeted Attributes")

    with_progress(enable=T,{

      p <- progressor(steps = nrow(us_flowpaths_out))
      o_targ_out<-us_flowpaths_out %>%
        dplyr::mutate(attr=furrr::future_pmap_dfr(
          list(
            link_id_in=link_id,
            attr_db_loc=list(attr_db_loc),
            loi_rasts_names=list(loi_rasts_names),
            weighting_scheme_o=list(weighting_scheme_o),
            p=list(p)
          ),
          .options = furrr_options(globals = FALSE),
          carrier::crate(
            function(link_id_in,
                     attr_db_loc,
                     loi_rasts_names,
                     weighting_scheme_o,
                     p
            ){
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`

              con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc)

              attr_nms<-names(c(loi_rasts_names$num_rast,loi_rasts_names$cat_rast))
              names(attr_nms)<-attr_nms

              names(weighting_scheme_o)<-weighting_scheme_o

              attrs<-sapply(sapply(loi_rasts_names$num_rast,unique),unique)
              mean_out<-NULL
              sd_out<-NULL

              out<-dplyr::left_join(
                dplyr::tbl(con_attr,"o_target_weights") %>%
                  dplyr::select(catch_link_id,cell_number,tidyselect::any_of(weighting_scheme_o)) %>%
                  dplyr::filter(catch_link_id %in% link_id_in),
                dplyr::tbl(con_attr,"attrib_tbl") %>%
                  dplyr::rename(link_id=subb_link_id),
                by="cell_number"
              )

              if ("iFLO" %in% weighting_scheme_o){ # This is the only way I could get around an error by iterating over weighting_scheme_s
                out<-out %>%
                  dplyr::left_join(
                    out %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("iFLO")),.names="{.col}_iFLO" )) %>%
                      dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_","iFLO"))),

                    by = c("cell_number", "link_id")
                  )
              }

              if ("HAiFLO" %in% weighting_scheme_o){
                out<-out %>%
                  dplyr::left_join(
                    out %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("HAiFLO")),.names="{.col}_HAiFLO" )) %>%
                      dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_","HAiFLO"))),
                    by = c("cell_number", "link_id")
                  )
              }

              # out<-out %>%
              #   dplyr::left_join(
              #     purrr::map(weighting_scheme_o,function(x){
              #       out %>%
              #         dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms),~.*!!rlang::sym(x),.names="{.col}_{x}" )) %>%
              #         dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_",x)))
              #     }) %>%
              #       purrr::reduce(dplyr::left_join,by=c("cell_number","link_id")),
              #     by = c("cell_number", "link_id")
              #   )

              mean_out<-purrr::map_dfc(weighting_scheme_o,function(x){
                out_iFLO<-NULL
                out_HAiFLO<-NULL

                if ("iFLO" %in% x){
                  out_iFLO<-out %>%
                    dplyr::summarize(dplyr::across(tidyselect::ends_with(paste0("_","iFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T) ))  %>%
                    dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$num_rast),"_")),~paste0(.x,"_mean")) %>%
                    dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$cat_rast),"_")),~paste0(.x,"_prop")) %>%
                    dplyr::collect()
                }

                if ("HAiFLO" %in% x){
                  out_HAiFLO<-out %>%
                    dplyr::summarize(dplyr::across(tidyselect::ends_with(paste0("_","HAiFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T) ))  %>%
                    dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$num_rast),"_")),~paste0(.x,"_mean")) %>%
                    dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$cat_rast),"_")),~paste0(.x,"_prop")) %>%
                    dplyr::collect()
                }

                dplyr::bind_cols(out_iFLO,out_HAiFLO)

              }) %>%
                dplyr::mutate(link_id=link_id_in) %>%
                dplyr::select(link_id,tidyselect::everything())

              # if (any("mean"==attrs)|length(loi_rasts_names$cat_rast) >0){
              #
              #   mean_out<-purrr::map_dfc(weighting_scheme_o,function(x){
              #     out %>%
              #       dplyr::summarize(dplyr::across(tidyselect::ends_with(paste0("_",x)),~sum(.,na.rm=T)/sum(!!rlang::sym(x),na.rm=T) ))  %>%
              #       dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$num_rast),"_")),~paste0(.x,"_mean")) %>%
              #       dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$cat_rast),"_")),~paste0(.x,"_prop")) %>%
              #       dplyr::collect()
              #   }) %>%
              #     dplyr::mutate(link_id=link_id_in) %>%
              #     dplyr::select(link_id,tidyselect::everything())
              # }

              if (any(attrs %in% c("sd","stdev"))) {
                sd_out<-purrr::map_dfc(weighting_scheme_o,function(x){
                  out_iFLO<-NULL
                  out_HAiFLO<-NULL

                  if ("iFLO" %in% x){
                    out_iFLO<-out %>%
                      dplyr::select(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                    tidyselect::any_of("iFLO")
                      ) %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                                  ~(!!rlang::sym("iFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T)))^2)),
                                                  .names="{.col}_iFLO_term1"),
                                    dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                                  ~ ((sum(!!rlang::sym("iFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLO")!=0,na.rm=T)) * sum(!!rlang::sym("iFLO"),na.rm=T),
                                                  .names="{.col}_iFLO_term2"
                                    ))%>%
                      dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),sum,na.rm=T),
                                       dplyr::across(tidyselect::ends_with("_term2"),~.[1])
                      ) %>%
                      dplyr::collect()
                  }

                  if ("HAiFLO" %in% x){
                    out_HAiFLO<-out %>%
                      dplyr::select(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                    tidyselect::any_of("HAiFLO")
                      ) %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                                  ~(!!rlang::sym("HAiFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T)))^2)),
                                                  .names="{.col}_HAiFLO_term1"),
                                    dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                                                  ~ ((sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLO"),na.rm=T),
                                                  .names="{.col}_HAiFLO_term2"
                                    ))%>%
                      dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),sum,na.rm=T),
                                       dplyr::across(tidyselect::ends_with("_term2"),~.[1])
                      ) %>%
                      dplyr::collect()
                  }


                  out_final<-dplyr::bind_cols(out_iFLO,out_HAiFLO)

                  out_final %>%
                    tidyr::pivot_longer(tidyselect::everything()) %>%
                    dplyr::mutate(attr=stringr::str_split_fixed(name,"_iFLO_|_HAiFLO_",2)[,1],
                                  term=stringr::str_split_fixed(name,"_iFLO_|_HAiFLO_",2)[,2]) %>%
                    dplyr::rowwise() %>%
                    dplyr::mutate(hw=gsub(paste0(attr,"_","|","_",term,""),"",name)) %>%
                    dplyr::ungroup() %>%
                    dplyr::mutate(name=paste0(attr,"_",hw,"_sd")) %>%
                    dplyr::group_by(name) %>%
                    dplyr::summarize(sd=sqrt(value[term=="term1"]/value[term=="term2"])) %>%
                    dplyr::ungroup() %>%
                    tidyr::pivot_wider(names_from = name,values_from=sd)

                })

              }
              # if (any(attrs %in% c("sd","stdev"))) {
              #   sd_out<-purrr::map_dfc(weighting_scheme_o,function(x){
              #     out_t<-out %>%
              #       dplyr::select(tidyselect::any_of(names(loi_rasts_names$num_rast)),
              #                     tidyselect::any_of(x)
              #       ) %>%
              #       dplyr::mutate(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
              #                                   ~(!!rlang::sym(x) * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym(x),na.rm=T)))^2)),
              #                                   .names="{.col}_{x}_term1"),
              #                     dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
              #                                   ~ ((sum(!!rlang::sym(x)!=0,na.rm=T)-1)/sum(!!rlang::sym(x)!=0,na.rm=T)) * sum(!!rlang::sym(x),na.rm=T),
              #                                   .names="{.col}_{x}_term2"
              #                     ))%>%
              #       dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),sum,na.rm=T),
              #                        dplyr::across(tidyselect::ends_with("_term2"),~.[1])
              #       ) %>%
              #       dplyr::collect()
              #
              #     purrr::map_dfc(names(loi_rasts_names$num_rast),function(y){
              #       out_t %>%
              #         dplyr::select(tidyselect::contains(y)) %>%
              #         dplyr::mutate(a= sqrt( !!sym(paste0(y,"_",x,"_term1")) / !!sym(paste0(y,"_",x,"_term2")) )) %>%
              #         dplyr::select(a) %>%
              #         dplyr::rename_with(~paste0(y,"_",x,"_sd"))
              #     })
              #   })
              #
              # }

              final_out<-dplyr::bind_cols(
                list(
                  mean_out,
                  sd_out
                )
              ) %>%
                dplyr::select(-tidyselect::any_of("link_id"))

              p()

              DBI::dbDisconnect(con_attr)
              return(final_out)

            })
        )) %>%
        select(everything(),-us_flowpaths) %>%
        unnest(cols=attr)
    })
  }

  #browser()

  final_out<-target_IDs %>%
    mutate(link_id=as.character(link_id)) %>%
    left_join(lumped_out ,by="link_id") %>%
    left_join(s_targ_out ,by="link_id") %>%
    left_join(o_targ_out ,by="link_id") %>%
    mutate(across(ends_with("_prop"),~ifelse(is.na(.),0,.)))

  data.table::fwrite(final_out,file.path(temp_dir,out_filename))

  zip(zip_loc,
      file.path(temp_dir,out_filename),
      flags = '-r9Xjq'
  )


  suppressWarnings(file.remove(list.files(temp_dir,full.names = T,recursive = T)))


  return(final_out)



































  # # Calculate Lumped Stats --------------------------------------------------
  # custfun <- function(x, type) {
  #   #browser()
  #   sapply(type,function(type) switch(type,
  #                                     mean = mean(x,na.rm=T),
  #                                     median = median(x,na.rm=T),
  #                                     sd = sd(x,na.rm=T),
  #                                     stdev = sd(x,na.rm=T),
  #                                     min = min(x,na.rm=T),
  #                                     max = max(x,na.rm=T),
  #                                     sum = sum(x,na.rm=T),
  #                                     prop=sum(x,na.rm=T)/length(x)
  #   )
  #   )
  # }
  #
  # hw2<-purrr::map(hw_streams_lo,terra::rast)
  # names(hw2)<-sapply(hw2,names)
  #
  #
  # # data.table::data.table() %>%
  # # data.table::setkey("subb_link_id")
  #
  #
  # s_weighted<-purrr::map(hw2,function(hww) exactextractr::exact_extract(
  #   loi_rasts_comb*hww,
  #   all_subb,
  #   weights=NULL,
  #   fun="sum",
  #   append_cols="link_id",
  #   progress=F
  # ) %>%
  #   #dplyr::bind_rows() %>%
  #   #dplyr::select(-coverage_fraction) %>%
  #   stats::setNames(c("link_id",names(loi_rasts_comb))) %>%
  #   data.table::data.table() %>%
  #   data.table::setkey("link_id")#%>%
  # #dplyr::rename_with(.cols=c(-link_id),~paste0(.x,"_",names(hww)))
  # )
  #
  # s_dwmean<-us_flowpaths_out%>%
  #   dplyr::transmute(link_id_upper=as.numeric(link_id),us_flowpaths=us_flowpaths) %>%
  #   dplyr::group_by(link_id_upper) %>%
  #   dplyr::summarise(
  #     result=furrr::future_pmap_dfc(list(hww=hw_weighted,sww=s_weighted,lid=us_flowpaths,loi_rasts_nms=list(loi_rasts_names)),
  #                                   carrier::crate(function(hww,sww,lid,loi_rasts_nms){
  #                                     require(data.table)
  #                                     out<-sww[link_id %in% lid$link_id,lapply(.SD, sum, na.rm=TRUE),.SDcols=-1]/
  #                                       unlist(hww[link_id %in% lid$link_id,lapply(.SD, sum, na.rm=TRUE),.SDcols=-1])
  #
  #                                     stats::setNames(out,paste0(names(sww)[-c(1)],
  #                                                                "_",
  #                                                                names(hww)[-c(1)],
  #                                                                dplyr::case_when(
  #                                                                  names(sww)[-c(1)] %in% loi_rasts_nms$num_rast ~ "_mean",
  #                                                                  T ~ "_prop"
  #                                                                )))
  #
  #                                   }))
  #   ) %>%
  #   tidyr::unnest(cols=result) %>%
  #   dplyr::rename(link_id=link_id_upper)
  #
  #
  # # if (lumped_scheme){
  # #   if (verbose) print("Calculating Lumped Attributes")
  # #   lumped<-exactextractr::exact_extract(
  # #     loi_rasts_comb,
  # #     all_subb,
  # #     weights=NULL,
  # #     fun=NULL,
  # #     include_cols="link_id",
  # #     progress=F
  # #   ) %>%
  # #     dplyr::bind_rows() %>%
  # #     select(-coverage_fraction) %>%
  # #     stats::setNames(c("link_id",names(loi_rasts_comb)))
  # #
  # #   lumped_numeric<-left_join(us_flowpaths_out %>%
  # #                               transmute(link_id_upper=as.numeric(link_id),us_flowpaths=us_flowpaths) %>%
  # #                               unnest(us_flowpaths) %>%
  # #                               mutate(link_id=as.numeric(link_id)),
  # #                             lumped %>% select(link_id,any_of(names(loi_rasts$num_rast))),
  # #                             by="link_id"
  # #   ) %>%
  # #     group_by(link_id_upper) %>%
  # #     summarize(
  # #       across(c(any_of(names(loi_rasts$num_rast)),-contains("link_id")),~list(custfun(type=loi_numeric_stats,x=.x)))
  # #     ) %>%
  # #     mutate(
  # #       across(any_of(names(loi_rasts$num_rast)),~map_dfr(.,~.))
  # #     ) %>%
  # #     unnest(cols=c(everything(),-link_id_upper),names_sep="_lumped_") %>%
  # #     rename(link_id=link_id_upper)
  # #
  # #
  # #   lumped_cat<-left_join(us_flowpaths_out %>%
  # #                           transmute(link_id_upper=as.numeric(link_id),us_flowpaths=us_flowpaths) %>%
  # #                           unnest(us_flowpaths) %>%
  # #                           mutate(link_id=as.numeric(link_id)),
  # #                         lumped %>% select(link_id,any_of(names(loi_rasts$cat_rast))),
  # #                         by="link_id"
  # #   ) %>%
  # #     group_by(link_id_upper) %>%
  # #     summarize(
  # #       across(c(any_of(names(loi_rasts$cat_rast)),-contains("link_id")),~list(custfun(type="prop",x=.x)))
  # #     ) %>%
  # #     mutate(
  # #       across(any_of(names(loi_rasts$cat_rast)),~map_dfr(.,~.))
  # #     ) %>%
  # #     unnest(cols=c(everything(),-link_id_upper),names_sep="_lumped_") %>%
  # #     rename(link_id=link_id_upper)
  # #
  # # } else {
  # #   lumped_numeric<-NULL
  # #   lumped_cat<-NULL
  # # }
  #
  #
  #
  #
  #
  # # Calculate weighted loi --------------------------------------------------
  # if (verbose) print("Calculating Weighted Attributes")
  #
  # #browser()
  #
  # # if (n_cores==1|
  # #     length(target_O_sub)<n_cores |
  # #     length(hw_streams_lo)<n_cores){
  # #
  # #   #spt1<-which.max(sapply(target_O_sub,nrow))
  # #
  # # }
  # #
  # max_ittr<-length(target_O_sub)+length(hw_streams_lo)
  #
  # if (n_cores==1) {
  #   splt<-length(hw_streams_lo)
  # } else {
  #   splt<-n_cores
  # }
  #
  # browser()
  #
  # with_progress(enable=T,{
  #   p <- progressor(steps = max_ittr)
  #
  #   loi_dw_out<-pmap(
  #     #loi_dw_out<-future_pmap(
  #     #  .options = furrr_options(globals = FALSE),
  #     list(
  #       target_O_subs=suppressWarnings(split(target_O_sub,1:splt)),
  #       all_catch=rep(list(all_catch),splt),
  #       hw=suppressWarnings(rev(split(unlist(hw_streams_lo),1:splt))),
  #       a_subb=rep(list(all_subb),splt),
  #       us_flowpaths=rep(list(us_flowpaths_out),splt),
  #       loi_rasts_exists=rep(list(loi_rasts_exists),splt),
  #       loi_cols_sub=rep(list(loi_cols),splt),
  #       weighting_scheme_o=rep(list(weighting_scheme_o),splt),
  #       loi_numeric_stats=rep(list(loi_numeric_stats),splt),
  #       inv_function=rep(list(inv_function),splt),
  #       use_exising_hw=rep(list(use_exising_hw),splt),
  #       temp_dir=rep(list(temp_dir),splt),
  #       zip_loc=rep(list(zip_loc),splt),
  #       dw_dir=rep(list(dw_dir),splt),
  #       p=rep(list(p),splt)
  #     ),
  #     carrier::crate(function(target_O_subs,
  #                             all_catch,
  #                             hw,
  #                             a_subb,
  #                             us_flowpaths,
  #                             loi_rasts_exists,
  #                             loi_cols_sub,
  #                             weighting_scheme_o,
  #                             loi_numeric_stats,
  #                             inv_function,
  #                             use_exising_hw,
  #                             temp_dir,
  #                             zip_loc,
  #                             dw_dir,
  #                             p
  #     ) {
  #       if (is.null(hw) & is.null(target_O_subs)) return(NULL)
  #       options(scipen = 999)
  #       browser()
  #       `%>%` <- magrittr::`%>%`
  #
  #       temp_dir_sub<-file.path(temp_dir,basename(tempfile()))
  #
  #       # Read in LOI -------------------------------------------------------------
  #       loi_rasts<-purrr::map(loi_rasts_exists,terra::rast)
  #       loi_rasts_nms<-purrr::map(loi_rasts,names)
  #
  #       loi_rasts_nms<-purrr::map(loi_rasts_nms,~.[.%in%loi_cols_sub])
  #
  #       loi_rasts_comb<-terra::rast(loi_rasts)
  #       names(loi_rasts_comb)<-unlist(sapply(loi_rasts,names))
  #       loi_rasts_comb<-terra::subset(loi_rasts_comb,loi_cols_sub)
  #
  #       # Perform Stream Weighting ------------------------------------------------
  #
  #       if (length(hw)>0 & any(c("mean","sd") %in% loi_numeric_stats)) {
  #         custfun <- function(x,y, type) {
  #           #browser()
  #           sapply(type,function(type) switch(type,
  #                                             mean = mean(x,na.rm=T),
  #                                             weighted.mean = sum(x,na.rm=T)/sum(y,na.rm=T),
  #                                             median = median(x,na.rm=T),
  #                                             sd = sd(x,na.rm=T),
  #                                             stdev = sd(x,na.rm=T),
  #                                             min = min(x,na.rm=T),
  #                                             max = max(x,na.rm=T),
  #                                             sum = sum(x,na.rm=T),
  #                                             prop=sum(x,na.rm=T)/length(x)
  #           )
  #           )
  #         }
  #
  #         # sumfun<-function(x,us_flowpaths=us_flowpaths){
  #         #   out<-purrr::map_dfr(us_flowpaths$us_flowpaths,
  #         #                       ~x[x$link_id %in% .$link_id,] %>%
  #         #                         dplyr::select(-tidyselect::any_of("link_id")) %>%
  #         #                         dplyr::summarise(dplyr::across(tidyselect::everything(),sum))
  #         #                       ,.id = "link_id")
  #         #
  #         #   out<-out[match(us_flowpaths$link_id, out$link_id,nomatch = 0),]
  #         #   out<-dplyr::select(out,-link_id)
  #         #   return(out)
  #         # }
  #
  #         #browser()
  #         #hw2<-terra::rast(hw[[1]])
  #         hw2<-purrr::map(hw,terra::rast)
  #         names(hw2)<-sapply(hw2,names)
  #
  #         hw_weighted<-purrr::map(hw2,function(hww) exactextractr::exact_extract(
  #           hww,
  #           a_subb,
  #           weights=NULL,
  #           fun="sum",
  #           append_cols="link_id",
  #           progress=F
  #         )%>%
  #           #dplyr::bind_rows() %>%
  #           #dplyr::select(-coverage_fraction) %>%
  #           stats::setNames(c("link_id",names(hww))) %>%
  #           data.table::data.table() %>%
  #           data.table::setkey("link_id")
  #         )
  #
  #         #loi_dist <- loi_rasts_comb * hw2
  #
  #         if (any(c("mean","sd") %in% loi_numeric_stats) | length(loi_rasts_nms$cat_rast)>0){
  #
  #           s_weighted<-purrr::map(hw2,function(hww) exactextractr::exact_extract(
  #             loi_rasts_comb*hww,
  #             a_subb,
  #             weights=NULL,
  #             fun="sum",
  #             append_cols="link_id",
  #             progress=F
  #           ) %>%
  #             #dplyr::bind_rows() %>%
  #             #dplyr::select(-coverage_fraction) %>%
  #             stats::setNames(c("link_id",names(loi_rasts_comb))) %>%
  #             data.table::data.table() %>%
  #             data.table::setkey("link_id")#%>%
  #           #dplyr::rename_with(.cols=c(-link_id),~paste0(.x,"_",names(hww)))
  #           )
  #
  #           s_dwmean<-us_flowpaths%>%
  #             dplyr::transmute(link_id_upper=as.numeric(link_id),us_flowpaths=us_flowpaths) %>%
  #             dplyr::group_by(link_id_upper) %>%
  #             dplyr::summarise(
  #               result=purrr::pmap_dfc(list(hww=hw_weighted,sww=s_weighted,lid=us_flowpaths),function(hww,sww,lid){
  #                 out<-sww[link_id %in% lid$link_id,lapply(.SD, sum, na.rm=TRUE),.SDcols=-1]/
  #                   unlist(hww[link_id %in% lid$link_id,lapply(.SD, sum, na.rm=TRUE),.SDcols=-1])
  #
  #                 stats::setNames(out,paste0(names(sww)[-c(1)],
  #                                            "_",
  #                                            names(hww)[-c(1)],
  #                                            dplyr::case_when(
  #                                              names(sww)[-c(1)] %in% loi_rasts_nms$num_rast ~ "_mean",
  #                                              T ~ "_prop"
  #                                            )))
  #               })
  #             ) %>%
  #             tidyr::unnest(cols=result) %>%
  #             dplyr::rename(link_id=link_id_upper)
  #
  #
  #
  #           # lumped_cat<-left_join(us_flowpaths_out %>%
  #           #                         transmute(link_id_upper=as.numeric(link_id),us_flowpaths=us_flowpaths) %>%
  #           #                         unnest(us_flowpaths) %>%
  #           #                         mutate(link_id=as.numeric(link_id)),
  #           #                       lumped %>% select(link_id,any_of(names(loi_rasts$cat_rast))),
  #           #                       by="link_id"
  #           # ) %>%
  #           #   group_by(link_id_upper) %>%
  #           #   summarize(
  #           #     across(c(any_of(names(loi_rasts$cat_rast)),-contains("link_id")),~list(custfun(type="prop",x=.x)))
  #           #   ) %>%
  #           #   mutate(
  #           #     across(any_of(names(loi_rasts$cat_rast)),~map_dfr(.,~.))
  #           #   ) %>%
  #           #   unnest(cols=c(everything(),-link_id_upper),names_sep="_lumped_") %>%
  #           #   rename(link_id=link_id_upper)
  #           #
  #           #
  #           # distance_weighted mean
  #
  #           # dw_o_sum<-purrr::map(hw,~exactextractr::exact_extract(.,a_subb,
  #           #                                                       fun="sum",progress=F))
  #           #
  #           # loi_dw_o_sum<-purrr::map2(hw,names(hw),function(hhw,nm){
  #           #   exactextractr::exact_extract(loi_rasts_comb*hhw,a_subb,
  #           #                                fun="sum",progress=F)  #%>%
  #           #   #dplyr::rename_with(~paste0(.,"_",nm))
  #           # })
  #           #
  #           # o_dwmean<-purrr::map2(loi_dw_o_sum,dw_o_sum,~(.x/unlist(.y)))
  #           #
  #           # o_dwmean<-purrr::map2_dfc(o_dwmean,names(o_dwmean),function(xx,yy){
  #           #   colnames(xx)<-gsub("sum.","",colnames(xx))
  #           #   colnames(xx)[colnames(xx) %in% loi_rasts_nms$num_rast]<-paste0(
  #           #     colnames(xx)[colnames(xx) %in% loi_rasts_nms$num_rast],
  #           #     "_",yy,"_mean"
  #           #   )
  #           #
  #           #   colnames(xx)[colnames(xx) %in% loi_rasts_nms$cat_rast]<-paste0(
  #           #     colnames(xx)[colnames(xx) %in% loi_rasts_nms$cat_rast],
  #           #     "_",yy,"_prop"
  #           #   )
  #           #
  #           #   return(xx)
  #           # })
  #
  #
  #           # dw_s_sum<-exactextractr::exact_extract(hw2,a_subb,fun="sum",na.rm=T,progress=F) %>%
  #           #   dplyr::mutate(link_id=a_subb$link_id)
  #           # dw_s_sum<-sumfun(dw_s_sum,us_flowpaths)
  #           #
  #           # loi_dw_s_sum<-exactextractr::exact_extract(loi_dist,a_subb,fun="sum",na.rm=T,progress=F)%>%
  #           #   dplyr::mutate(link_id=a_subb$link_id)
  #           # loi_dw_s_sum<-sumfun(loi_dw_s_sum,us_flowpaths)
  #           #
  #           # s_dwmean<-loi_dw_s_sum/unlist(dw_s_sum)
  #           #
  #           # colnames(s_dwmean)[colnames(s_dwmean) %in% loi_rasts_nms$num_rast]<-paste0(
  #           #   colnames(s_dwmean)[colnames(s_dwmean) %in% loi_rasts_nms$num_rast],
  #           #   "_",names(hw2),"_mean"
  #           # )
  #           #
  #           # colnames(s_dwmean)[colnames(s_dwmean) %in% loi_rasts_nms$cat_rast]<-paste0(
  #           #   colnames(s_dwmean)[colnames(s_dwmean) %in% loi_rasts_nms$cat_rast],
  #           #   "_",names(hw2),"_prop"
  #           # )
  #         } else {
  #           s_dwmean<-NULL
  #         }
  #
  #
  #         # distance_weighted SD
  #         if (length(loi_rasts_nms$num_rast)>0 & any(c("sd") %in% loi_numeric_stats)){
  #
  #
  #           distance_weights_sum<-hw_weighted
  #
  #           M<-purrr::map(hw2,function(hww) exactextractr::exact_extract(
  #             hww!=0,
  #             a_subb,
  #             weights=NULL,
  #             fun="sum",
  #             append_cols="link_id",
  #             progress=F
  #           )%>%
  #             stats::setNames(c("link_id",names(hww))) %>%
  #             data.table::data.table() %>%
  #             data.table::setkey("link_id")
  #           )
  #
  #           term2 <- purrr::map(distance_weights_sum,~ M[[names(.)[2]]][,sapply(.SD, function(x) (x-1)/x),.SDcols=-1])  #((M[,sapply(.SD, function(x) (x-1)/x),.SDcols=-1]) / M) * distance_weights_sum
  #           term2 <- purrr::map(term2,~./distance_weights_sum[[colnames(.)]][,colnames(.),with=F])
  #           term2 <- purrr::map(term2,~cbind(hw_weighted[[1]][,1],.) %>% data.table::setkey("link_id"))
  #
  #           all_loi_vals<-exactextractr::exact_extract(
  #             loi_rasts$num_rast,
  #             a_subb,
  #             weights=NULL,
  #             fun=NULL,
  #             include_cols="link_id",
  #             progress=F
  #           ) %>%
  #             dplyr::bind_rows() %>%
  #             dplyr::select(-coverage_fraction) %>%
  #             stats::setNames(c("link_id",names(loi_rasts$num_rast))) %>%
  #             data.table::data.table() %>%
  #             data.table::setkey("link_id")
  #
  #           all_hw_vals<-exactextractr::exact_extract(
  #             terra::rast(hw2),
  #             a_subb,
  #             weights=NULL,
  #             fun=NULL,
  #             include_cols="link_id",
  #             progress=F
  #           ) %>%
  #             dplyr::bind_rows() %>%
  #             dplyr::select(-coverage_fraction) %>%
  #             stats::setNames(c("link_id",names(hw2))) %>%
  #             data.table::data.table() %>%
  #             data.table::setkey("link_id")
  #
  #           term1<-us_flowpaths%>%
  #             dplyr::transmute(link_id_upper=as.numeric(link_id),us_flowpaths=us_flowpaths) %>%
  #             dplyr::group_by(link_id_upper) %>%
  #             dplyr::summarise(
  #               result=purrr::pmap(
  #                 list(hww=list(all_hw_vals),
  #                      sww=list(all_loi_vals),
  #                      lid=list(us_flowpaths),
  #                      mean_s=list(s_dwmean),
  #                      term2=list(term2)),
  #                 function(hww,sww,lid,mean_s,term2){
  #                   browser()
  #                   lid<-lid[[1]]
  #
  #                   term2_fin<-purrr::map(term2,~.[link_id %in% lid$link_id,lapply(.SD, sum, na.rm=TRUE),.SDcols=-1])
  #
  #                   loi_distwtd_mean<-mean_s %>%
  #                     dplyr::filter(link_id %in% lid) %>%
  #                     dplyr::select(tidyselect::contains(colnames(sww)[-c(1)])) %>%
  #                     dplyr::select(tidyselect::contains(colnames(hww)[-c(1)])) %>%
  #                     as.list()
  #
  #                   sww_l<-as.list(sww[link_id %in% lid$link_id,colnames(sww)[-c(1)],with=F])
  #                   hww_l<-as.list(hww[link_id %in% lid$link_id,colnames(hww)[-c(1)],with=F])
  #
  #                   term1<- purrr::map(stats::setNames(names(sww_l),names(sww_l)),
  #                                      function(y) purrr::map_dfc(stats::setNames(names(hww_l),names(hww_l)),
  #                                                                 function(x) {
  #                                                                   term1<-sum(na.rm=T,hww_l[[which(names(hww_l)==x)]]*
  #                                                                                ((sww_l[[which(names(sww_l)==y)]]-
  #                                                                                    loi_distwtd_mean[[which(grepl(paste0("_",x,"_mean"),names(loi_distwtd_mean)) & grepl(paste0(y,"_"),names(loi_distwtd_mean)))]])
  #                                                                                 ^2))
  #
  #                                                                   sqrt(term1/unlist(term2_fin[[which(names(term2_fin)==x)]]))
  #
  #                                                                 }))
  #
  #                   tibble::enframe(term1) %>%
  #                     tidyr::unnest(value) %>%
  #                     tidyr::pivot_wider(names_from=name,values_from=names(hww_l),names_glue="{name}_{.value}_sd")
  #
  #                 })
  #             )
  #
  #
  #           s_dwsd<-us_flowpaths%>%
  #             dplyr::transmute(link_id_upper=as.numeric(link_id),us_flowpaths=us_flowpaths) %>%
  #             dplyr::group_by(link_id_upper) %>%
  #             dplyr::summarise(
  #               result=purrr::pmap_dfr(list(hww=hw_weighted,sww=s_weighted,lid=us_flowpaths),function(hww,sww,lid){
  #                 out<-sww[link_id %in% lid$link_id,lapply(.SD, sum, na.rm=TRUE),.SDcols=-1]/
  #                   unlist(hww[link_id %in% lid$link_id,lapply(.SD, sum, na.rm=TRUE),.SDcols=-1])
  #
  #                 stats::setNames(out,paste0(names(sww)[-c(1)],
  #                                            "_",
  #                                            names(hww)[-c(1)],
  #                                            dplyr::case_when(
  #                                              names(hww)[-c(1)] %in% loi_rasts_nms$num_rast ~ "_mean",
  #                                              T ~ "_prop"
  #                                            )))
  #               })
  #             ) %>%
  #             tidyr::unnest(cols=result)
  #
  #
  #
  #
  #           distwtd_sum<-terra::extract(hw2,
  #                                       a_subb,
  #                                       fun="sum",na.rm=T,method="simple",touches=F,ID=F)%>%
  #             dplyr::mutate(link_id=a_subb$link_id)
  #           distwtd_sum<-sumfun(distwtd_sum,us_flowpaths)
  #
  #           distwtd_M<-terra::extract(hw2,
  #                                     a_subb ,
  #                                     fun=function(x) sum(x!=0,na.rm=T),
  #                                     method="simple",touches=F,ID=F)%>%
  #             dplyr::mutate(link_id=a_subb$link_id)
  #           distwtd_M<-sumfun(distwtd_M,us_flowpaths)
  #
  #           term2<-((distwtd_M - 1) / distwtd_M) * distwtd_sum
  #
  #           # this is the correct way, but slower
  #           #browser()
  #           t0<-terra::extract(hw2,
  #                              a_subb,
  #                              fun=NULL ,
  #                              method="simple",touches=F,ID=T,cells=T)
  #
  #           t1<-terra::extract(terra::subset(loi_dist,loi_rasts_nms$num_rast),
  #                              a_subb,
  #                              fun=NULL,
  #                              method="simple",touches=F,ID=T,cells=T)
  #
  #           # t2<-terra::extract(terra::subset(loi_dist,loi_rasts_nms$num_rast)*hw2,
  #           #                    a_subb,
  #           #                    fun=NULL ,
  #           #                    method="simple",touches=F,ID=T,cells=T)
  #
  #           t2<-t0 %>%
  #             dplyr::left_join(t1,by=c("ID","cell")) %>%
  #             dplyr::group_by(ID) %>%
  #             dplyr::mutate(dplyr::across(tidyselect::any_of(loi_rasts_nms$num_rast), ~sum(.x,na.rm=T),.names = "sum_{.col}")) %>%
  #             dplyr::mutate(dplyr::across(tidyselect::any_of(names(hw2)), ~sum(.x,na.rm=T),.names = "dw_tbl")) %>%
  #             dplyr::mutate(dplyr::across(tidyselect::any_of(paste0("sum_",loi_rasts_nms$num_rast)), ~.x/dw_tbl,.names = "mean_{.col}")) %>%
  #             dplyr::ungroup()
  #
  #           t3<-purrr:::map_dfc(loi_rasts_nms$num_rast,function(x){
  #             t2[,names(hw2)]*((t2[,x]-t2[,paste0("mean_sum_",x)])^2)
  #           })
  #           colnames(t3)<-loi_rasts_nms$num_rast
  #           t3$ID<-t2$ID
  #
  #           # t2<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast)*dw,
  #           #                    all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
  #           #                    fun=NULL ,
  #           #                    method="simple",touches=F,ID=T,cells=T)
  #           #
  #           # t2<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast)*dw,
  #           #                    all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
  #           #                    fun=NULL ,
  #           #                    method="simple",touches=F,ID=T,cells=T) %>%
  #           #   dplyr::group_by(ID) %>%
  #           #   dplyr::mutate(dplyr::across(tidyselect::any_of(loi_rasts_nms$num_rast),mean)) %>%
  #           #   dplyr::ungroup()
  #
  #           # t3<-(t0[[names(dw)]]*(t1-t2)^2)
  #           # t3$ID<-t2$ID
  #           term1<-t3 %>%
  #             dplyr::group_by(ID) %>%
  #             dplyr::summarise(dplyr::across(tidyselect::everything(),sum,na.rm=T)) %>%
  #             dplyr::ungroup() %>%
  #             dplyr::select(-tidyselect::any_of("ID"),-tidyselect::any_of("cell"))%>%
  #             dplyr::mutate(link_id=a_subb$link_id)
  #
  #           # t3<-(t0[[names(hw2)]]*(t1-t2)^2)
  #           # t3$ID<-t2$ID
  #           # term1<-t3 %>%
  #           #   dplyr::group_by(ID) %>%
  #           #   dplyr::summarise(dplyr::across(tidyselect::everything(),sum,na.rm=T)) %>%
  #           #   dplyr::ungroup() %>%
  #           #   dplyr::select(-tidyselect::any_of("ID"),-tidyselect::any_of("cell"))%>%
  #           #   dplyr::mutate(link_id=a_subb$link_id)
  #
  #           term1<-sumfun(term1,us_flowpaths)
  #
  #
  #           s_dwSD <- suppressMessages(purrr::map_dfc(colnames(term1), ~tibble::tibble(a=sqrt(term1[[.]] / unlist(term2)))))
  #           colnames(s_dwSD)<-paste0(colnames(term1),"_",names(hw2),"_sd")
  #         } else {
  #           s_dwSD<-NULL
  #         }
  #
  #         s_out<-dplyr::bind_cols(
  #           s_dwmean,
  #           s_dwSD
  #         ) %>%
  #           dplyr::mutate(link_id=us_flowpaths$link_id)
  #
  #         p()
  #       } else {
  #         s_out<-NULL
  #       }
  #
  #       # Perform Target Weighting ------------------------------------------------
  #
  #       if (length(target_O_subs)>0){
  #         target_S <- file.path("/vsizip",zip_loc,"dem_streams_d8.tif")
  #         dem <- file.path("/vsizip",zip_loc,"dem_final.tif")
  #         flow_accum <- file.path("/vsizip",zip_loc,"dem_accum_d8.tif")
  #
  #         o_out<-purrr::map(target_O_subs,function(x){
  #           if (!use_exising_hw){
  #             hw<-hydroweight::hydroweight(hydroweight_dir=temp_dir_sub,
  #                                          target_O = x,
  #                                          target_S = target_S,
  #                                          target_uid = basename(tempfile()),
  #                                          OS_combine = FALSE,
  #                                          dem=dem,
  #                                          flow_accum = flow_accum,
  #                                          weighting_scheme = weighting_scheme_o,
  #                                          inv_function = inv_function,
  #                                          clean_tempfiles=T,
  #                                          return_products = T,
  #                                          wrap_return_products=F,
  #                                          save_output=F)
  #           } else {
  #             trg_fl<-paste0("unnest_group_",x$unn_group[[1]],"_",weighting_scheme_o,"_inv_distances.tif")
  #             hw<-purrr::map(trg_fl,~terra::rast(file.path("/vsizip",dw_dir,.)))
  #             names(hw)<-sapply(hw,names)
  #           }
  #
  #           if (any(c("mean") %in% loi_numeric_stats) | length(loi_rasts_nms$cat_rast)>0){
  #             # distance_weighted mean
  #
  #             dw_o_sum<-purrr::map(hw,~exactextractr::exact_extract(.,all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
  #                                                                   fun="sum",progress=F))
  #
  #             loi_dw_o_sum<-purrr::map2(hw,names(hw),function(hhw,nm){
  #               exactextractr::exact_extract(loi_rasts_comb*hhw,all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
  #                                            fun="sum",progress=F)  #%>%
  #               #dplyr::rename_with(~paste0(.,"_",nm))
  #             })
  #
  #             o_dwmean<-purrr::map2(loi_dw_o_sum,dw_o_sum,~(.x/unlist(.y)))
  #
  #             o_dwmean<-purrr::map2_dfc(o_dwmean,names(o_dwmean),function(xx,yy){
  #               colnames(xx)<-gsub("sum.","",colnames(xx))
  #               colnames(xx)[colnames(xx) %in% loi_rasts_nms$num_rast]<-paste0(
  #                 colnames(xx)[colnames(xx) %in% loi_rasts_nms$num_rast],
  #                 "_",yy,"_mean"
  #               )
  #
  #               colnames(xx)[colnames(xx) %in% loi_rasts_nms$cat_rast]<-paste0(
  #                 colnames(xx)[colnames(xx) %in% loi_rasts_nms$cat_rast],
  #                 "_",yy,"_prop"
  #               )
  #
  #               return(xx)
  #             })
  #
  #           } else {
  #             o_dwmean<-NULL
  #           }
  #
  #           if (length(loi_rasts_nms$num_rast)>0 & any(c("sd") %in% loi_numeric_stats)){
  #             #browser()
  #             o_dwSD<-purrr::map(hw,function(dw){
  #               #browser()
  #               distwtd_sum<-terra::extract(dw,
  #                                           all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
  #                                           fun="sum",na.rm=T,method="simple",touches=F,ID=F)
  #               distwtd_M<-terra::extract(dw,
  #                                         all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
  #                                         fun=function(x) sum(x!=0,na.rm=T),
  #                                         method="simple",touches=F,ID=F)
  #
  #               term2<-((distwtd_M - 1) / distwtd_M) * distwtd_sum
  #
  #               # this is the correct way, but slower
  #               t0<-terra::extract(dw,
  #                                  all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
  #                                  fun=NULL ,
  #                                  method="simple",touches=F,ID=T,cells=T)
  #
  #               t1<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast),
  #                                  all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
  #                                  fun=NULL ,
  #                                  method="simple",touches=F,ID=T,cells=T)
  #
  #               t2<-t0 %>%
  #                 dplyr::left_join(t1,by=c("ID","cell")) %>%
  #                 dplyr::group_by(ID) %>%
  #                 dplyr::mutate(dplyr::across(tidyselect::any_of(loi_rasts_nms$num_rast), ~sum(.x,na.rm=T),.names = "sum_{.col}")) %>%
  #                 dplyr::mutate(dplyr::across(tidyselect::any_of(names(dw)), ~sum(.x,na.rm=T),.names = "dw_tbl")) %>%
  #                 dplyr::mutate(dplyr::across(tidyselect::any_of(paste0("sum_",loi_rasts_nms$num_rast)), ~.x/dw_tbl,.names = "mean_{.col}")) %>%
  #                 dplyr::ungroup()
  #
  #               t3<-purrr:::map_dfc(loi_rasts_nms$num_rast,function(x){
  #                 t2[,names(dw)]*((t2[,x]-t2[,paste0("mean_sum_",x)])^2)
  #               })
  #               colnames(t3)<-loi_rasts_nms$num_rast
  #               t3$ID<-t2$ID
  #
  #               # t2<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast)*dw,
  #               #                    all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
  #               #                    fun=NULL ,
  #               #                    method="simple",touches=F,ID=T,cells=T)
  #               #
  #               # t2<-terra::extract(terra::subset(loi_rasts_comb,loi_rasts_nms$num_rast)*dw,
  #               #                    all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]),
  #               #                    fun=NULL ,
  #               #                    method="simple",touches=F,ID=T,cells=T) %>%
  #               #   dplyr::group_by(ID) %>%
  #               #   dplyr::mutate(dplyr::across(tidyselect::any_of(loi_rasts_nms$num_rast),mean)) %>%
  #               #   dplyr::ungroup()
  #
  #               # t3<-(t0[[names(dw)]]*(t1-t2)^2)
  #               # t3$ID<-t2$ID
  #               term1<-t3 %>%
  #                 dplyr::group_by(ID) %>%
  #                 dplyr::summarise(dplyr::across(tidyselect::everything(),sum,na.rm=T)) %>%
  #                 dplyr::ungroup() %>%
  #                 dplyr::select(-tidyselect::any_of("ID"),-tidyselect::any_of("cell"))
  #
  #               loi_distwtd_sd <- suppressMessages(purrr::map_dfc(colnames(term1), ~tibble::tibble(a=sqrt(term1[[.]] / unlist(term2)))))
  #               colnames(loi_distwtd_sd)<-paste0(colnames(term1),"_",names(dw),"_sd")
  #
  #               return(loi_distwtd_sd)
  #
  #             })
  #           } else {
  #             o_dwSD<-NULL
  #           }
  #           # distance weighted SD is costly to calculate like this
  #
  #           o_out<-dplyr::bind_cols(
  #             o_dwmean,
  #             o_dwSD
  #           ) %>%
  #             dplyr::mutate(link_id=all_catch %>% dplyr::filter(link_id %in% x[["link_id"]]) %>% dplyr::pull(link_id))
  #
  #           p()
  #           return(o_out)
  #
  #         })
  #
  #         o_out<-dplyr::bind_rows(o_out)
  #       } else {
  #         o_out<-NULL
  #       }
  #
  #       return(
  #         list(
  #           s_out=s_out,
  #           o_out=o_out
  #         )
  #       )
  #
  #     }
  #     ))
  # })
  #
  # if (inherits(lumped_numeric,"list")){
  #   lumped_numeric<-reduce(lumped_numeric,left_join,by="link_id")
  # }
  # if (inherits(lumped_cat,"list")){
  #   lumped_cat<-reduce(lumped_cat,left_join,by="link_id")
  # }
  #
  # s_out<-map(loi_dw_out,~.$s_out)
  # s_out<-s_out[!sapply(s_out,is.null)]
  # s_out<-reduce(s_out,left_join)
  #
  # o_out<-map_dfr(loi_dw_out,~.$o_out)
  #
  # final_out<-target_IDs %>%
  #   mutate(link_id=as.character(link_id)) %>%
  #   left_join(lumped_numeric %>% mutate(link_id=as.character(link_id)),by="link_id") %>%
  #   left_join(lumped_cat %>% mutate(link_id=as.character(link_id)),by="link_id") %>%
  #   left_join(s_out %>% mutate(link_id=as.character(link_id)),by="link_id") %>%
  #   left_join(o_out %>% mutate(link_id=as.character(link_id)),by="link_id")
  #
  # data.table::fwrite(final_out,file.path(temp_dir,out_filename))
  #
  # zip(zip_loc,
  #     file.path(temp_dir,out_filename),
  #     flags = '-r9Xjq'
  # )
  #
  # # output<-input
  # #
  # # output[[gsub(".csv","",out_filename)]]<-final_out
  #
  # DBI::dbDisconnect(con_attr)
  #
  # suppressWarnings(file.remove(list.files(temp_dir,full.names = T,recursive = T)))
  #
  #
  # return(final_out)
}


