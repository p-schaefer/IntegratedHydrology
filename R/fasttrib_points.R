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



  con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)
  # DBI::dbSendStatement(con_attr, "PRAGMA busy_timeout = 10000")
  # DBI::dbSendStatement(con_attr,"PRAGMA journal_mode = OFF")
  # DBI::dbSendStatement(con_attr,"PRAGMA synchronous = 0")
  # DBI::dbSendStatement(con_attr,"PRAGMA cache_size = 1000000")
  # # DBI::dbSendStatement(con_attr,"PRAGMA locking_mode = EXCLUSIVE")
  # DBI::dbSendStatement(con_attr,"PRAGMA temp_store = MEMORY")
  # DBI::dbSendStatement(con_attr,"PRAGMA mmap_size = 30000000000")
  # DBI::dbSendStatement(con_attr,"PRAGMA page_size = 32768")

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

  all_subb_v<-terra::vect(all_subb %>% select(link_id))
  all_catch_v<-terra::vect(all_catch %>% select(link_id))

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

    attrib_tbl<-copy_to(df=as_tibble(matrix(ncol = length(names(loi_rasts_comb))+2,nrow=1),.name_repair="minimal") %>%
                          setNames(c("subb_link_id","cell_number",names(loi_rasts_comb))) %>%
                          mutate(across(everything(),~1.1)) %>%
                          .[F,],
                        con_attr,
                        "attrib_tbl",
                        overwrite =T,
                        temporary =F,
                        #indexes=c("subb_link_id","cell_number"),
                        analyze=T,
                        in_transaction=T)


    n_cores_2<-n_cores
    if (n_cores_2>1) n_cores_2<-n_cores_2-1

    all_subb_v$core<-rep(1:(n_cores_2),length.out=nrow(all_subb_v))
    splt<-terra::split(all_subb_v,"core")
    splt<-lapply(splt,function(x){
      x$split<-rep(1:10,length.out=nrow(x))
      terra::split(x,"split")
    })
    splt<-map(splt,~map(.,terra::wrap))

    total_outs<-sum(unlist(map(splt,~map(.,length))))

    out<-furrr::future_pmap(list(x=splt,
                                 loi_rasts_exists=list(loi_rasts_exists),
                                 loi_cols=list(loi_cols),
                                 temp_dir=list(temp_dir)),
                            .options = furrr_options(globals = FALSE),
                            carrier::crate(
                              function(x,
                                       loi_rasts_exists,
                                       loi_cols,
                                       temp_dir
                              ){
                                #browser()

                                options(scipen = 999)
                                `%>%` <- magrittr::`%>%`

                                loi_rasts<-purrr::map(loi_rasts_exists,terra::rast)
                                loi_rasts_comb<-terra::rast(loi_rasts)
                                names(loi_rasts_comb)<-unlist(sapply(loi_rasts,names))
                                loi_rasts_comb<-terra::subset(loi_rasts_comb,loi_cols)


                                splt<-x

                                attrib_tbl<-future::future(
                                  packages = c("future","furrr","purrr","terra","sf","dplyr","data.table","carrier","magrittr","stats","base","utils"),
                                  globals = c("splt","loi_rasts_comb","temp_dir"),
                                  {
                                    options(future.rng.onMisuse = "ignore")

                                    out<-purrr::pmap(list(xx=splt,
                                                          loi_rasts_comb=list(loi_rasts_comb),
                                                          temp_dir=list(temp_dir)
                                    ),
                                    carrier::crate(
                                      function(xx,
                                               loi_rasts_comb,
                                               temp_dir
                                      ){
                                        #browser()
                                        options(scipen = 999)
                                        `%>%` <- magrittr::`%>%`

                                        xx<-terra::unwrap(xx)

                                        out<-exactextractr::exact_extract(
                                          loi_rasts_comb,
                                          xx,
                                          weights=NULL,
                                          include_cell=T,
                                          fun=NULL,
                                          include_cols="link_id",
                                          progress=F
                                        ) %>%
                                          dplyr::bind_rows() %>%
                                          dplyr::select(-coverage_fraction) %>%
                                          stats::setNames(c("subb_link_id",names(loi_rasts_comb),"cell_number")) %>%
                                          dplyr::select(cell_number,subb_link_id,everything()) %>%
                                          data.table::fwrite(file=file.path(temp_dir,paste0("attr_sub_s_",xx$core[[1]],"_",xx$split[[1]],".csv")))

                                        # out<-terra::extract(
                                        #   loi_rasts_comb,
                                        #   xx,
                                        #   cells=T,
                                        #   ID=T,
                                        #   fun=NULL
                                        # ) %>%
                                        #   tibble::as_tibble() %>%
                                        #   dplyr::left_join(sf::st_as_sf(xx) %>%
                                        #                      tibble::as_tibble() %>%
                                        #                      dplyr::select(link_id) %>%
                                        #                      dplyr::mutate(ID=dplyr::row_number()),
                                        #                    by="ID") %>%
                                        #   dplyr::rename(cell_number=cell,
                                        #                 subb_link_id=link_id) %>%
                                        #   dplyr::select(-ID) %>%
                                        #   data.table::fwrite(file=file.path(temp_dir,paste0("attr_sub_s_",xx$core[[1]],"_",xx$split[[1]],".csv")))

                                        file.rename(
                                          file.path(temp_dir,paste0("attr_sub_s_",xx$core[[1]],"_",xx$split[[1]],".csv")),
                                          file.path(temp_dir,paste0("attr_sub_",xx$core[[1]],"_",xx$split[[1]],".csv"))
                                        )

                                        return(NA)

                                      }
                                    ))
                                  })

                                return(attrib_tbl)
                              })

    )

    #browser()
    #future_attrib_tbl <- lapply(out,future::futureOf)

    total_procs<-0

    with_progress(enable=T,{
      p <- progressor(steps = total_outs)

      while(total_procs<total_outs){
        Sys.sleep(0.2)

        fl_attr<-list.files(temp_dir,"attr_sub_",full.names = T)
        fl_attr<-fl_attr[grepl(".csv",fl_attr)]

        if (length(fl_attr)>0) {
          df<-purrr::map(fl_attr,data.table::fread) %>%
            dplyr::bind_rows()

          out<-DBI::dbAppendTable(conn=con_attr,
                                  name="attrib_tbl",
                                  value=df)

          fr<-file.remove(fl_attr)
          for (i in seq_along(fr)){
            p()
          }
          #eval(parse(text=paste0(rep("p();",sum(fr)))))
          #rep(p,sum(fr))
          total_procs<-total_procs+sum(fr)
        }
      }

    })

    Sys.sleep(0.2)

    fl_attr<-list.files(temp_dir,"attr_sub_",full.names = T)
    fl_attr<-fl_attr[grepl(".csv",fl_attr)]

    if (length(fl_attr)>0) {
      df<-purrr::map(fl_attr,data.table::fread) %>%
        dplyr::bind_rows()

      out<-DBI::dbAppendTable(conn=con_attr,
                              name="attrib_tbl",
                              value=df)

      fr<-file.remove(fl_attr)
    }


    #})

    attrib_tbl<-tbl(con_attr,"attrib_tbl")


    DBI::dbSendStatement(con_attr,"CREATE INDEX inx_attrib_tbl ON attrib_tbl (subb_link_id, cell_number)")
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
      if (verbose) print("Merging stream segments")

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
      names(hw_streams_nm)<-gsub(".tif","",hw_streams_nm)

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

    # hw2<-map(hw_streams_lo,terra::rast)
    # names(hw2)<-sapply(hw2,names)

    s_trg_weights<-copy_to(df=as_tibble(matrix(ncol = length(names(hw_streams_lo))+2),nrow=1) %>%
                             setNames(c("subb_link_id","cell_number",names(hw_streams_lo))) %>%
                             mutate(across(everything(),~1.1)) %>%
                             .[F,],
                           con_attr,
                           "s_target_weights",
                           overwrite =T,
                           temporary =F,
                           #indexes=c("subb_link_id","cell_number"),
                           analyze=T,
                           in_transaction=T)

    # n_cores_2<-n_cores # already run above
    # if (n_cores_2>1) n_cores_2<-n_cores_2-1
    #
    # all_subb_v$core<-rep(1:(n_cores_2),length.out=nrow(all_subb_v))
    # splt<-terra::split(all_subb_v,"core")
    # splt<-lapply(splt,function(x){
    #   x$split<-rep(1:10,length.out=nrow(x))
    #   terra::split(x,"split")
    # })
    # splt<-map(splt,~map(.,terra::wrap))

    #total_outs<-sum(unlist(map(splt,~map(.,length))))

    out<-furrr::future_pmap(list(x=splt,
                                 loi_rasts_exists=list(hw_streams_lo),
                                 temp_dir=list(temp_dir)),
                            .options = furrr_options(globals = FALSE),
                            carrier::crate(
                              function(x,
                                       loi_rasts_exists,
                                       temp_dir
                              ){
                                #browser()

                                options(scipen = 999)
                                `%>%` <- magrittr::`%>%`

                                loi_rasts<-purrr::map(loi_rasts_exists,terra::rast)
                                loi_rasts_comb<-terra::rast(loi_rasts)
                                names(loi_rasts_comb)<-unlist(sapply(loi_rasts,names))


                                splt<-x

                                attrib_tbl<-future::future(
                                  packages = c("future","furrr","purrr","terra","sf","dplyr","data.table","carrier","magrittr","stats","base","utils"),
                                  globals = c("splt","loi_rasts_comb","temp_dir"),
                                  {
                                    options(future.rng.onMisuse = "ignore")

                                    out<-purrr::pmap(list(xx=splt,
                                                          loi_rasts_comb=list(loi_rasts_comb),
                                                          temp_dir=list(temp_dir)
                                    ),
                                    carrier::crate(
                                      function(xx,
                                               loi_rasts_comb,
                                               temp_dir
                                      ){
                                        #browser()
                                        options(scipen = 999)
                                        `%>%` <- magrittr::`%>%`

                                        xx<-terra::unwrap(xx)

                                        out<-exactextractr::exact_extract(
                                          loi_rasts_comb,
                                          xx,
                                          weights=NULL,
                                          include_cell=T,
                                          fun=NULL,
                                          include_cols="link_id",
                                          progress=F
                                        ) %>%
                                          dplyr::bind_rows() %>%
                                          dplyr::select(-coverage_fraction) %>%
                                          stats::setNames(c("subb_link_id",names(loi_rasts_comb),"cell_number")) %>%
                                          dplyr::select(cell_number,subb_link_id,tidyselect::everything()) %>%
                                          data.table::fwrite(file=file.path(temp_dir,paste0("s_target_weights_sub_s_",xx$core[[1]],"_",xx$split[[1]],".csv")))


                                        # out<-terra::extract(
                                        #   loi_rasts_comb,
                                        #   xx,
                                        #   cells=T,
                                        #   ID=T,
                                        #   fun=NULL
                                        # ) %>%
                                        #   tibble::as_tibble() %>%
                                        #   dplyr::left_join(sf::st_as_sf(xx) %>%
                                        #                      tibble::as_tibble() %>%
                                        #                      dplyr::select(link_id) %>%
                                        #                      dplyr::mutate(ID=dplyr::row_number()),
                                        #                    by="ID") %>%
                                        #   dplyr::rename(cell_number=cell,
                                        #                 subb_link_id=link_id) %>%
                                        #   dplyr::select(-ID) %>%
                                        #   data.table::fwrite(file=file.path(temp_dir,paste0("s_target_weights_sub_s_",xx$core[[1]],"_",xx$split[[1]],".csv")))

                                        file.rename(
                                          file.path(temp_dir,paste0("s_target_weights_sub_s_",xx$core[[1]],"_",xx$split[[1]],".csv")),
                                          file.path(temp_dir,paste0("s_target_weights_sub_",xx$core[[1]],"_",xx$split[[1]],".csv"))
                                        )

                                        return(NA)

                                      }
                                    ))
                                  })

                                return(attrib_tbl)
                              })

    )

    #browser()
    #future_attrib_tbl <- lapply(out,future::futureOf)

    total_procs<-0

    with_progress(enable=T,{
      p <- progressor(steps = total_outs)

      while(total_procs<total_outs){
        Sys.sleep(0.2)

        fl_attr<-list.files(temp_dir,"s_target_weights_sub_",full.names = T)
        fl_attr<-fl_attr[grepl(".csv",fl_attr)]

        if (length(fl_attr)>0) {
          df<-purrr::map(fl_attr,data.table::fread) %>%
            dplyr::bind_rows()

          out<-DBI::dbAppendTable(conn=con_attr,
                                  name="s_target_weights",
                                  value=df)

          fr<-file.remove(fl_attr)
          for (i in seq_along(fr)){
            p()
          }
          total_procs<-total_procs+sum(fr)
        }
      }

    })

    Sys.sleep(0.2)

    fl_attr<-list.files(temp_dir,"s_target_weights_sub_",full.names = T)
    fl_attr<-fl_attr[grepl(".csv",fl_attr)]

    if (length(fl_attr)>0) {
      df<-purrr::map(fl_attr,data.table::fread) %>%
        dplyr::bind_rows()

      out<-DBI::dbAppendTable(conn=con_attr,
                              name="s_target_weights",
                              value=df)

      fr<-file.remove(fl_attr)
    }

    s_trg_weights<-tbl(con_attr,"s_target_weights")

    DBI::dbSendStatement(con_attr,"CREATE INDEX inx_s_target_weights ON s_target_weights (subb_link_id, cell_number)")

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
                  #indexes=c("catch_link_id","cell_number"),
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
            n_cores=rep(list(n_cores),splt),
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
                                  n_cores,
                                  con_attr_l,
                                  new_tbl,
                                  zip_loc,
                                  dw_dir,
                                  p
          ) {
            options(scipen = 999)
            `%>%` <- magrittr::`%>%`

            temp_dir_sub<-file.path(temp_dir,basename(tempfile()))
            dir.create(temp_dir_sub)

            target_S <- file.path("/vsizip",zip_loc,"dem_streams_d8.tif")
            dem <- file.path("/vsizip",zip_loc,"dem_final.tif")
            flow_accum <- file.path("/vsizip",zip_loc,"dem_accum_d8.tif")

            o_out<-purrr::map(target_O_subs, # I don't think this can be parallel
                              function(x){
                                #print(x$unn_group[[1]])
                                sub_catch<-all_catch %>%
                                  dplyr::filter(link_id %in% x$link_id)

                                sub_catch_v<-terra::vect(sub_catch %>% dplyr::select(link_id))

                                if (!use_exising_hw){
                                  hw<-hydroweight::hydroweight(hydroweight_dir=temp_dir_sub,
                                                               target_O = x,
                                                               target_S = target_S,
                                                               target_uid = paste0("unnest_group_",x$unn_group[[1]]),
                                                               OS_combine = FALSE,
                                                               dem=dem,
                                                               flow_accum = flow_accum,
                                                               weighting_scheme = weighting_scheme_o,
                                                               inv_function = inv_function,
                                                               clean_tempfiles=T,
                                                               return_products = F,
                                                               wrap_return_products=F,
                                                               save_output=T)

                                  hw_o_nm<-utils::unzip(list=T,hw)$Name
                                  names(hw_o_nm)<-gsub(".tif","",hw_o_nm)

                                  hw_o_lo<-purrr::map(hw_o_nm,function(x){
                                    file.path("/vsizip",hw,x)
                                  })

                                } else {
                                  trg_fl<-paste0("unnest_group_",x$unn_group[[1]],"_",weighting_scheme_o,"_inv_distances.tif")
                                  hw<-purrr::map(trg_fl,~terra::rast(file.path("/vsizip",dw_dir,.)))
                                  names(hw)<-sapply(hw,names)
                                }

                                #browser()

                                n_cores_2<-n_cores # already run above
                                if (n_cores_2>1) n_cores_2<-n_cores_2-1

                                sub_catch_v$core<-rep(1:(n_cores_2),length.out=nrow(sub_catch_v))
                                splt<-terra::split(sub_catch_v,"core")
                                splt<-lapply(splt,function(x){
                                  x$split<-rep(1:10,length.out=nrow(x))
                                  terra::split(x,"split")
                                })
                                splt<-purrr::map(splt,~purrr::map(.,terra::wrap))

                                total_outs<-sum(unlist(purrr::map(splt,~purrr::map(.,length))))

                                out<-furrr::future_pmap(list(x=splt,
                                                             loi_rasts_exists=list(hw_o_lo),
                                                             temp_dir=list(temp_dir_sub)),
                                                        .options = furrr::furrr_options(globals = FALSE),
                                                        carrier::crate(
                                                          function(x,
                                                                   loi_rasts_exists,
                                                                   temp_dir
                                                          ){
                                                            #browser()

                                                            options(scipen = 999)
                                                            `%>%` <- magrittr::`%>%`

                                                            loi_rasts<-purrr::map(loi_rasts_exists,terra::rast)
                                                            loi_rasts_comb<-terra::rast(loi_rasts)
                                                            names(loi_rasts_comb)<-unlist(sapply(loi_rasts,names))


                                                            splt<-x

                                                            attrib_tbl<-future::future(
                                                              packages = c("future","furrr","purrr","terra","sf","dplyr","data.table","carrier","magrittr","stats","base","utils"),
                                                              globals = c("splt","loi_rasts_comb","temp_dir"),
                                                              {
                                                                options(future.rng.onMisuse = "ignore")

                                                                out<-purrr::pmap(list(xx=splt,
                                                                                      loi_rasts_comb=list(loi_rasts_comb),
                                                                                      temp_dir=list(temp_dir)
                                                                ),
                                                                carrier::crate(
                                                                  function(xx,
                                                                           loi_rasts_comb,
                                                                           temp_dir
                                                                  ){
                                                                    #browser()
                                                                    options(scipen = 999)
                                                                    `%>%` <- magrittr::`%>%`

                                                                    xx<-terra::unwrap(xx)

                                                                    out<-exactextractr::exact_extract(
                                                                      loi_rasts_comb,
                                                                      xx,
                                                                      weights=NULL,
                                                                      include_cell=T,
                                                                      fun=NULL,
                                                                      include_cols="link_id",
                                                                      progress=F
                                                                    ) %>%
                                                                      dplyr::bind_rows() %>%
                                                                      dplyr::select(-coverage_fraction) %>%
                                                                      stats::setNames(c("catch_link_id",names(loi_rasts_comb),"cell_number")) %>%
                                                                      dplyr::select(cell_number,subb_link_id,tidyselect::everything()) %>%
                                                                      data.table::fwrite(file=file.path(temp_dir,paste0("o_target_weights_sub_s_",xx$core[[1]],"_",xx$split[[1]],".csv")))


                                                                    # out<-terra::extract(
                                                                    #   loi_rasts_comb,
                                                                    #   xx,
                                                                    #   cells=T,
                                                                    #   ID=T,
                                                                    #   fun=NULL
                                                                    # ) %>%
                                                                    #   tibble::as_tibble() %>%
                                                                    #   dplyr::left_join(sf::st_as_sf(xx) %>%
                                                                    #                      tibble::as_tibble() %>%
                                                                    #                      dplyr::select(link_id) %>%
                                                                    #                      dplyr::mutate(ID=dplyr::row_number()),
                                                                    #                    by="ID") %>%
                                                                    #   dplyr::rename(cell_number=cell,
                                                                    #                 catch_link_id=link_id) %>%
                                                                    #   dplyr::select(-ID) %>%
                                                                    #   data.table::fwrite(file=file.path(temp_dir,paste0("o_target_weights_sub_s_",xx$core[[1]],"_",xx$split[[1]],".csv")))

                                                                    file.rename(
                                                                      file.path(temp_dir,paste0("o_target_weights_sub_s_",xx$core[[1]],"_",xx$split[[1]],".csv")),
                                                                      file.path(temp_dir,paste0("o_target_weights_sub_",xx$core[[1]],"_",xx$split[[1]],".csv"))
                                                                    )

                                                                    return(NA)

                                                                  }
                                                                ))
                                                              })

                                                            return(attrib_tbl)
                                                          })

                                )

                                #browser()
                                #future_attrib_tbl <- lapply(out,future::futureOf)

                                total_procs<-0

                                while(total_procs<total_outs){
                                  Sys.sleep(0.2)

                                  fl_attr<-list.files(temp_dir_sub,"o_target_weights_sub_",full.names = T)
                                  fl_attr<-fl_attr[grepl(".csv",fl_attr)]

                                  if (length(fl_attr)>0) {
                                    df<-purrr::map(fl_attr,data.table::fread) %>%
                                      dplyr::bind_rows()

                                    out<-DBI::dbAppendTable(conn=con_attr_l,
                                                            name="o_target_weights",
                                                            value=df)

                                    fr<-file.remove(fl_attr)
                                    total_procs<-total_procs+sum(fr)
                                  }
                                }


                                Sys.sleep(0.2)

                                fl_attr<-list.files(temp_dir_sub,"o_target_weights_sub_",full.names = T)
                                fl_attr<-fl_attr[grepl(".csv",fl_attr)]

                                if (length(fl_attr)>0) {
                                  df<-purrr::map(fl_attr,data.table::fread) %>%
                                    dplyr::bind_rows()

                                  out<-DBI::dbAppendTable(conn=con_attr_l,
                                                          name="o_target_weights",
                                                          value=df)

                                  fr<-file.remove(fl_attr)
                                }

                                file.remove(hw)

                                # out<-purrr::pmap(list(x=terra::split(sub_catch_v,"link_id"),
                                #                       con_attr_l=list(con_attr_l),
                                #                       hw2_l=list(terra::rast(hw))
                                # ),
                                # carrier::crate(function(x,
                                #                         con_attr_l,
                                #                         hw2_l){
                                #   options(scipen = 999)
                                #   `%>%` <- magrittr::`%>%`
                                #
                                #   out<-terra::extract(
                                #     hw2_l,
                                #     x,
                                #     cells=T,
                                #     ID=F,
                                #     fun=NULL
                                #   ) %>%
                                #     tibble::as_tibble() %>%
                                #     dplyr::mutate(catch_link_id=x$link_id) %>%
                                #     dplyr::rename(cell_number=cell) %>%
                                #     DBI::dbAppendTable(conn=con_attr_l,
                                #                        name="o_target_weights",
                                #                        value=.)
                                #
                                #   return(NULL)
                                #
                                # }))

                                p()

                                return(NULL)
                              })

          }))

        DBI::dbSendStatement(con_attr,"CREATE INDEX idx_o_target_weights ON o_target_weights (catch_link_id, cell_number)")

      })

    }
  } else {
    target_O_sub<-NULL
  }

  DBI::dbSendStatement(con_attr,"PRAGMA analysis_limit=1000")
  DBI::dbSendStatement(con_attr,"PRAGMA vacuum")
  DBI::dbSendStatement(con_attr,"PRAGMA optimize")

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
        dplyr::mutate(attr=future_pmap_dfr(#furrr::future_
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
              #browser()
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`

              con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)
              # DBI::dbSendStatement(con_attr,"PRAGMA journal_mode = OFF")
              # DBI::dbSendStatement(con_attr,"PRAGMA synchronous = 0")
              # DBI::dbSendStatement(con_attr,"PRAGMA cache_size = 1000000")
              # # DBI::dbSendStatement(con_attr,"PRAGMA locking_mode = EXCLUSIVE")
              # DBI::dbSendStatement(con_attr,"PRAGMA temp_store = MEMORY")
              # DBI::dbSendStatement(con_attr,"PRAGMA mmap_size = 30000000000")
              # DBI::dbSendStatement(con_attr,"PRAGMA page_size = 32768")

              out<-dplyr::semi_join(
                dplyr::tbl(con_attr,"attrib_tbl") %>%
                  dplyr::rename(link_id=subb_link_id),
                dplyr::tbl(con_attr,"us_flowpaths") %>%
                  dplyr::filter(pour_point_id %in% link_id_in) %>%
                  dplyr::rename(link_id=origin_link_id),
                by="link_id"
              ) %>%
                dplyr::compute()

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

              con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)
              # DBI::dbSendStatement(con_attr,"PRAGMA journal_mode = OFF")
              # DBI::dbSendStatement(con_attr,"PRAGMA synchronous = 0")
              # DBI::dbSendStatement(con_attr,"PRAGMA cache_size = 1000000")
              # # DBI::dbSendStatement(con_attr,"PRAGMA locking_mode = EXCLUSIVE")
              # DBI::dbSendStatement(con_attr,"PRAGMA temp_store = MEMORY")
              # DBI::dbSendStatement(con_attr,"PRAGMA mmap_size = 30000000000")
              # DBI::dbSendStatement(con_attr,"PRAGMA page_size = 32768")

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
                ) %>%
                dplyr::compute()

              if ("iFLS" %in% weighting_scheme_s){ # This is the only way I could get around an error by iterating over weighting_scheme_s
                out<-out %>%
                  dplyr::left_join(
                    out %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("iFLS")),.names="{.col}_iFLS" )) %>%
                      dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_","iFLS"))),

                    by = c("cell_number", "link_id")
                  )%>%
                  dplyr::compute()
              }

              if ("HAiFLS" %in% weighting_scheme_s){
                out<-out %>%
                  dplyr::left_join(
                    out %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("HAiFLS")),.names="{.col}_HAiFLS" )) %>%
                      dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_","HAiFLS"))),
                    by = c("cell_number", "link_id")
                  )%>%
                  dplyr::compute()
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

              con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)
              # DBI::dbSendStatement(con_attr,"PRAGMA journal_mode = OFF")
              # DBI::dbSendStatement(con_attr,"PRAGMA synchronous = 0")
              # DBI::dbSendStatement(con_attr,"PRAGMA cache_size = 1000000")
              # # DBI::dbSendStatement(con_attr,"PRAGMA locking_mode = EXCLUSIVE")
              # DBI::dbSendStatement(con_attr,"PRAGMA temp_store = MEMORY")
              # DBI::dbSendStatement(con_attr,"PRAGMA mmap_size = 30000000000")
              # DBI::dbSendStatement(con_attr,"PRAGMA page_size = 32768")

              attr_nms<-names(c(loi_rasts_names$num_rast,loi_rasts_names$cat_rast))
              names(attr_nms)<-attr_nms

              names(weighting_scheme_o)<-weighting_scheme_o

              attrs<-sapply(sapply(loi_rasts_names$num_rast,unique),unique)
              mean_out<-NULL
              sd_out<-NULL

              out<-dplyr::semi_join(
                dplyr::tbl(con_attr,"attrib_tbl") %>%
                  dplyr::rename(link_id=subb_link_id),
                dplyr::tbl(con_attr,"us_flowpaths") %>%
                  dplyr::filter(pour_point_id %in% link_id_in) %>%
                  dplyr::rename(link_id=origin_link_id),
                by="link_id"
              ) %>%
                dplyr::left_join(
                  dplyr::tbl(con_attr,"o_target_weights") %>%
                    dplyr::select(cell_number,catch_link_id ,tidyselect::any_of(weighting_scheme_o)) %>%
                    dplyr::filter(catch_link_id %in% link_id_in) %>%
                    dplyr::select(-catch_link_id),
                  by=c("cell_number")
                )%>%
                dplyr::compute()


              # out<-dplyr::left_join(
              #   dplyr::tbl(con_attr,"o_target_weights") %>%
              #     dplyr::select(catch_link_id,cell_number,tidyselect::any_of(weighting_scheme_o)) %>%
              #     dplyr::rename(link_id=catch_link_id) %>%
              #     dplyr::filter(link_id %in% link_id_in),
              #   dplyr::tbl(con_attr,"attrib_tbl") %>%
              #     dplyr::select(-tidyselect::any_of("link_id")),
              #   by="cell_number"
              # ) %>%
              #   dplyr::compute()

              if ("iFLO" %in% weighting_scheme_o){ # This is the only way I could get around an error by iterating over weighting_scheme_s
                out<-out %>%
                  dplyr::left_join(
                    out %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("iFLO")),.names="{.col}_iFLO" )) %>%
                      dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_","iFLO"))),

                    by = c("cell_number", "link_id")
                  )%>%
                  dplyr::compute()
              }

              if ("HAiFLO" %in% weighting_scheme_o){
                out<-out %>%
                  dplyr::left_join(
                    out %>%
                      dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("HAiFLO")),.names="{.col}_HAiFLO" )) %>%
                      dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_","HAiFLO"))),
                    by = c("cell_number", "link_id")
                  )%>%
                  dplyr::compute()
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

}

#' @export
dbWithWriteTransaction <- function(conn, code) {
  dbExecute(conn, "BEGIN IMMEDIATE")
  rollback <- function(e) {
    call <- DBI::dbExecute(conn, "ROLLBACK")
    if (identical(call, FALSE)) {
      stop(paste(
        "Failed to rollback transaction.",
        "Tried to roll back because an error occurred:",
        conditionMessage(e)
      ), call. = FALSE)
    }
    if (inherits(e, "error")) stop(e)
  }
  tryCatch(
    {
      res <- force(code)
      DBI::dbExecute(conn, "COMMIT")
      res
    },
    db_abort = rollback,
    error = rollback,
    interrupt = rollback
  )
}
