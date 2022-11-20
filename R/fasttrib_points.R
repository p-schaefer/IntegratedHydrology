#' Quickly attribute stream segments/sampling points with layers of interest (loi)
#'
#' @param input output from `process_hydrology()` (if `process_loi()` was not run on `process_hydrology()`, `loi_file` must be specified)
#' @param out_filename Output file name, ending in '.csv'.
#' @param loi_file filepath of `process_loi()` output (optional, will search for loi layers in 'input' if not specified).
#' @param loi_cols character or NULL. Names of loi layers to include in summary. If NULL, all layers used.
#' @param iDW_file filepath of `prep_weights()` output (optional, will search for weight layers in 'input' if not specified).
#' @param store_iDW logical, if 'iDW_file' is NULL, weights will be added to 'input'
#' @param raw_file Optional, file path to store intermediate SQL table, if file doesn't exist, it will be created at the specified path, if it does exist it will be queried first, and only extra data that is needed will be added.
#' @param sample_points character or NULL. IDs of unique station identifiers provided in 'site_id_col' of `generate_vectors()` to calculate attributes for.
#' @param link_id character or NULL. 'link_id's of reaches to calculate attributes for.
#' @param target_o_type character. One of: c("point","segment_point","segment_whole"). Target for"iFLO", and "HAiFLO" weighting schemes. "Point" represents the sampling point on the stream, "segment_point" represents the upstream segment of the sampling points, and "segment_whole" will target the entire reach, regardless of where sampling occurred.
#' @param weighting_scheme character. One or more weighting schemes: c("lumped","iFLO", "iFLS", "HAiFLO", "HAiFLS")
#' @param loi_numeric_stats character. One or more of c("mean", "sd", "median", "min", "max", "sum"). Only distance-weighted versions of mean and SD are returned for all weighting schemes except lumped.
#' @param inv_function function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#'
#' @return A data.frame of weighted attributes for the requested areas
#' @export
#'
#' @importFrom carrier crate
#' @importFrom data.table fwrite
#' @importFrom DBI dbConnect dbDisconnect dbExecute dbRemoveTable
#' @importFrom dbplyr tbl_memdb sql
#' @importFrom dplyr filter tbl select distinct collect left_join mutate arrange copy_to group_by n desc ungroup bind_rows across collapse compute rows_append rename_with summarise summarize rowwise bind_cols
#' @importFrom furrr future_pmap furrr_options
#' @importFrom future nbrOfWorkers availableCores
#' @importFrom progressr with_progress progressor
#' @importFrom purrr map map2_dfr
#' @importFrom RSQLite SQLite
#' @importFrom sf read_sf
#' @importFrom stats setNames
#' @importFrom stringr str_split_fixed
#' @importFrom terra terraOptions
#' @importFrom tibble enframe tibble
#' @importFrom tidyr nest pivot_wider pivot_longer unnest
#' @importFrom tidyselect everything contains any_of ends_with starts_with
#' @importFrom whitebox wbt_options wbt_exe_path

fasttrib_points<-function(
    input,
    out_filename,
    loi_file=NULL,
    loi_cols=NULL,
    iDW_file=NULL,
    store_iDW=F,
    raw_file=NULL,
    sample_points=NULL,
    link_id=NULL,
    target_o_type=c("point","segment_point","segment_whole"),
    weighting_scheme =  c("lumped",  "iFLS", "HAiFLS","iFLO",  "HAiFLO"),
    loi_numeric_stats = c("mean", "sd", "median", "min", "max", "sum"),
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    temp_dir=NULL,
    verbose=F
){
  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  if (inherits(loi_file,"ihydro")) loi_file<-loi_file$outfile
  if (inherits(iDW_file,"ihydro")) iDW_file<-iDW_file$outfile

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  whitebox::wbt_options(exe_path=whitebox::wbt_exe_path(),
                        verbose=verbose>2,
                        wd=temp_dir)

  terra::terraOptions(verbose = verbose>3,
                      tempdir = temp_dir
  )

  target_o_type<-match.arg(target_o_type,several.ok = F)
  weighting_scheme<-match.arg(weighting_scheme,several.ok = T)
  loi_numeric_stats<-match.arg(loi_numeric_stats,several.ok = T)

  weighting_scheme_s<-weighting_scheme[grepl("FLS|iEucS",weighting_scheme)]
  weighting_scheme_o<-weighting_scheme[!grepl("lumped|FLS|iEucS",weighting_scheme)]
  lumped_scheme<-"lumped" %in% weighting_scheme
  #if (length(weighting_scheme_o)>0) message("Calculation for iFLO, and HAiFLO are slow")

  loi_numeric_stats<-stats::setNames(loi_numeric_stats,loi_numeric_stats)

  attr_db_loc<-db_loc<-dw_dir<-zip_loc<-input$outfile

  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)
  if (n_cores==0) n_cores<-1

  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)

  if (!is.logical(store_iDW)) stop("'store_iDW' must be logical")

  lyrs<-ihydro_layers(input)

  if (is.null(loi_file)) {
    if (!"loi_meta" %in% lyrs$layer_name) stop("no 'loi' data found in 'input'")
    loi_file<-input$outfile
    loi_file<-as.ihydro(loi_file)
  } else {
    loi_file<-as.ihydro(loi_file)
    loi_lyrs<-ihydro_layers(loi_file)
    if (!"loi_meta" %in% loi_lyrs$layer_name) stop("no 'loi' data found in 'loi_lyrs'")
  }
  # TODO: Add checks to make sure requested loi_cols are in input
  loi_meta<-sf::read_sf(loi_file,"loi_meta")
  if (is.null(loi_cols)) loi_cols<-loi_meta$loi_var_nms
  if (any(!loi_cols %in% loi_meta$loi_var_nms)) stop("some `loi_cols` not present in `loi_file`")
  loi_meta<-loi_meta %>%
    dplyr::filter(loi_var_nms %in% loi_cols)


  # Not 100% sure this is going to work
  if (is.null(iDW_file)) {
    if (store_iDW){
      iDW_file<-input$outfile

    } else {
      iDW_file<-tempfile()
      dir.create(iDW_file)
      iDW_file<-file.path(iDW_file,"temp_iDW.gpkg")
    }
  }

  if (verbose) message("Preparing inverse distance weights")
  iDW_out<-prep_weights(
    input=input,
    output_filename=iDW_file,
    sample_points=sample_points,
    link_id=link_id,
    target_o_type=target_o_type,
    weighting_scheme =   weighting_scheme[!grepl("lumped",weighting_scheme)],
    inv_function = inv_function,
    temp_dir=temp_dir,
    verbose=verbose
  )

  idw_lyrs<-ihydro_layers(iDW_out)

  # TODO: Add checks to make sure requested weighting_scheme are in input
  if (any(!sapply(weighting_scheme[!grepl("lumped",weighting_scheme)],function(x) any(grepl(x,idw_lyrs$layer_name))))) stop("no 'iDW' data found in 'input'")

  target_o_meta<-sf::read_sf(iDW_file,"target_o_meta")



  # Get Target IDs ----------------------------------------------------------
  target_IDs<-target_id_fun(
    db_fp=db_loc,
    sample_points=sample_points,
    link_id=link_id,
    segment_whole=target_o_type=="segment_whole",
    target_o_type=target_o_type
  )


  # Get subbasin IDs --------------------------------------------------------
  con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)

  subb_IDs<-dplyr::tbl(con,"us_flowpaths") %>%
    dplyr::filter(pour_point_id %in% local(target_IDs$link_id)) %>%
    dplyr::select(pour_point_id,origin_link_id) %>%
    dplyr::distinct() %>%
    dplyr::collect() %>%
    dplyr::left_join(
      target_o_meta %>%
        dplyr::mutate(link_id=as.character(link_id)),
      by=c("pour_point_id"="link_id")
    ) %>%
    dplyr::arrange(origin_link_id)

  DBI::dbDisconnect(con)

  # Setup output ------------------------------------------------------------
  temp_dir_sub<-file.path(temp_dir,basename(tempfile()))
  dir.create(temp_dir_sub)

  # if (is.null(raw_file)){
  #   raw_file<-file.path(temp_dir_sub,"temp_db.sqlite")
  # }
  #
  # if (!grepl("\\.sql$|\\.sqlite$|\\.db$",raw_file)) stop("`raw_file` must be a file path ending in '.sqlite', '.sql' or '.db'")
  #
  # con <- DBI::dbConnect(RSQLite::SQLite(), dbname = raw_file)
  # t1<-dplyr::copy_to(con,
  #                    subb_IDs,
  #                    "subb_IDs",
  #                    temporary=F,
  #                    overwrite=T,
  #                    indexes=c("pour_point_id","origin_link_id"))
  #
  # DBI::dbDisconnect(con)

  #browser()
  o1<-extract_raster_attributes(
    input=input,
    iDW_file=iDW_file,
    loi_file=loi_file,
    weighting_scheme=weighting_scheme,
    loi_numeric_stats=loi_numeric_stats,
    loi_cols=loi_cols,
    subb_IDs=subb_IDs,
    verbose=verbose
  )

  target_IDs_out<-target_id_fun(
    db_fp=db_loc,
    sample_points=sample_points,
    link_id=link_id,
    segment_whole=F,
    target_o_type=target_o_type
  )

  #browser()
  final_out<-target_IDs_out %>%
    dplyr::left_join(
      o1 %>% dplyr::mutate(pour_point_id=as.character(pour_point_id)),
      by=c("link_id"="pour_point_id"),
      multiple = "all"
    )

  data.table::fwrite(x=final_out,
                     file=out_filename,
                     buffMB = 128L,
                     nThread = 1,
                     showProgress = F)

  return(final_out)


  # Save extracted values to SQLite -----------------------------------------
  ev<-extract_raster_vals(
    input=input,
    raw_file=raw_file,
    iDW_file=iDW_file,
    loi_file=loi_file,
    weighting_scheme=weighting_scheme,
    loi_cols=loi_cols,
    subb_IDs=subb_IDs,
    temp_dir_sub=temp_dir_sub,
    verbose=verbose
  )

  # Calculate Endpoints -----------------------------------------------------
  #browser()
  if (verbose) message("Calculating Attributes")
  progressr::with_progress(enable=T,{
    p <- progressr::progressor(steps = (length(unique(subb_IDs$pour_point_id))))

    n_cores_sub<-n_cores

    #if (length(unique(subb_IDs$pour_point_id))<(n_cores_sub/2)) n_cores_sub<-1

    attr_out<-subb_IDs %>%
      dplyr::group_by(pour_point_id) %>%
      dplyr::mutate(gp_size=dplyr::n()) %>%
      dplyr::arrange(dplyr::desc(gp_size)) %>%
      dplyr::select(-gp_size) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(core=rep(1:n_cores_sub,length.out=dplyr::n())) %>%
      dplyr::group_by(core) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        out=furrr::future_pmap(
          #out=purrr::pmap(
          list(x=data,
               core=core,
               input=list(input),
               raw_file=list(raw_file),
               iDW_file=list(iDW_file),
               loi_file=list(loi_file),
               weighting_scheme=list(weighting_scheme),
               loi_cols=list(loi_cols),
               loi_numeric_stats=list(loi_numeric_stats),
               loi_meta=list(loi_meta),
               temp_dir=list(temp_dir),
               p=list(p)
          ),
          .options = furrr::furrr_options(globals = F),
          carrier::crate(
            function(x,
                     core,
                     input,
                     raw_file,
                     iDW_file,
                     loi_file,
                     weighting_scheme,
                     loi_cols,
                     loi_numeric_stats,
                     loi_meta,
                     temp_dir,
                     p){
              #browser()
              options(dplyr.summarise.inform = FALSE)
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`

              con<-DBI::dbConnect(RSQLite::SQLite(),raw_file)
              t1<-DBI::dbExecute(con,paste0("PRAGMA temp_store_directory = '",temp_dir,"'"))


              add_d<-dplyr::bind_rows(x$data)

              cols<-c(
                "link_id",
                ".rowid",
                loi_cols,
                weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")],
                unlist(purrr::map(paste0(weighting_scheme[weighting_scheme %in% c("iFLO","HAiFLO")],"_unn_group"),~paste0(.,unique(add_d$unn_group))))
              )

              cols<-c(cols,
                      unlist(purrr::map(cols[!cols %in% c("link_id",
                                                          ".rowid",
                                                          loi_cols)],
                                        ~paste0(loi_cols,"_",.)))
              )

              temp1<-tibble::enframe(cols) %>%
                dplyr::mutate(name=1) %>%
                tidyr::pivot_wider(names_from=value) %>%
                .[F,] %>%
                dplyr::mutate(dplyr::across(tidyselect::everything(),as.numeric)) %>%
                dplyr::select(-name)


              temp1<-temp1 %>%
                dplyr::copy_to(dest=con,
                               df=.,
                               name=paste0("temp_outer_",core),
                               overwrite = T,
                               temporary=T)

              add_d<-split(add_d,add_d$origin_link_id)

              for (i in add_d){

                sel_col<-unlist(purrr::map(paste0(weighting_scheme[weighting_scheme %in% c("iFLO","HAiFLO")],"_unn_group"),~paste0(.,unique(i$unn_group))))

                temp_sub<-dplyr::tbl(con,paste0(i$origin_link_id[[1]],"_raw")) %>%
                  dplyr::select(link_id,
                                tidyselect::contains(loi_cols),
                                tidyselect::any_of(weighting_scheme),#
                                tidyselect::any_of(sel_col)
                  ) %>%
                  dplyr::collapse()

                for (ii in c(weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")],sel_col)) {
                  temp_sub<-temp_sub %>%
                    dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols),
                                         ~.*(dbplyr::sql(ii)),.names=paste0("{.col}_",ii))) %>%
                    dplyr::collapse()
                  # dplyr::compute(analyze=F,
                  #                name=paste0("temp_temp_",core),
                  #                temporary=T,
                  #                cte=T,
                  #                overwrite=T)
                }

                # temp_sub<-dplyr::copy_to(dest=con,
                #                          df=temp_sub.,
                #                          name=paste0("temp_temp_",core),
                #                          overwrite = T,
                #                          temporary=T)
                # temp_sub<-dplyr::compute(temp_sub,
                #                          name=paste0("temp_temp_",core),
                #                          analyze=F,
                #                          temporary=T,
                #                          overwrite=T)

                temp1<-dplyr::rows_append(dplyr::tbl(con,paste0("temp_outer_",core)) ,
                                          temp_sub,
                                          copy=T,
                                          in_place =T
                )

                #DBI::dbRemoveTable(con,paste0("temp_temp_",core))
              }

              t1<-DBI::dbExecute(con,paste0("CREATE INDEX idx_temp ON ",paste0("temp_outer_",core),"(link_id)"))

              # raw_tbl_all2<-purrr::map(add_d,
              #                          function(z) {
              #                            dplyr::tbl(con,paste0(z$origin_link_id[[1]],"_raw")) %>%
              #                              dplyr::select(link_id,
              #                                            tidyselect::contains(loi_cols),
              #                                            tidyselect::any_of(weighting_scheme),
              #                                            tidyselect::any_of(paste0(paste0(weighting_scheme,"_unn_group"),z$unn_group))
              #                              ) %>%
              #                              dplyr::collapse()
              #                          })
              #
              # t1<-purrr::reduce(raw_tbl_all2,dplyr::full_join,by="link_id")
              #
              # raw_tbl_all<-raw_tbl_all2[[1]]
              #
              # if (length(raw_tbl_all2)>1){
              #   for (i in 2:length(raw_tbl_all2)){
              #     raw_tbl_all<-dplyr::union_all(raw_tbl_all,raw_tbl_all2[[i]]) %>%
              #       dplyr::collapse()
              #   }
              # }

              # for (i in cols[!cols %in% c("link_id",
              #                             ".rowid",
              #                             loi_cols)]) {
              #   temp1<-temp1 %>%
              #     dplyr::mutate(across(any_of(loi_cols), ~.*(dbplyr::sql(i)),.names=paste0("{.col}_",i))) %>%
              #     dplyr::collapse()
              # }
              #
              # temp2<- temp1 %>%
              #   dplyr::compute(temporary=T,
              #                  name=paste0("temp_fin_",core),
              #                  indexes=c("link_id"))

              out<-purrr::map2_dfr(x$data,x$pour_point_id,
                                   function(y,yy) {

                                     suffix<-paste0("_unn_group",y$unn_group[[1]])
                                     sel_cols<-c(
                                       loi_cols,
                                       weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")],
                                       paste0(weighting_scheme[weighting_scheme %in% c("iFLO","HAiFLO")],suffix),
                                       unlist(purrr::map(loi_cols,~paste0(.,paste0("_",weighting_scheme[weighting_scheme %in% c("iFLO","HAiFLO")],suffix))))
                                     )

                                     raw_tbl<-dplyr::tbl(con,paste0("temp_outer_",core)) %>%
                                       dplyr::filter(link_id %in% local(y$origin_link_id)) %>%
                                       dplyr::select(tidyselect::any_of(sel_cols)) %>%
                                       dplyr::rename_with(.cols=tidyselect::ends_with(y$unn_group[[1]]),
                                                          ~gsub(paste0("_unn_group",y$unn_group[[1]]),"",.x)) %>%
                                       dplyr::compute(
                                         name=paste0("temp_",yy),
                                         analyze=T,
                                         temporary=T
                                       )

                                     # Combine Tables ----------------------------------------------------------


                                     # raw_tbl_all<-purrr::map2(y$origin_link_id,y$unn_group,
                                     #                          function(z,zz) {
                                     #                            dplyr::tbl(con,paste0(z,"_raw")) %>%
                                     #                              dplyr::select(tidyselect::contains(loi_cols),
                                     #                                            tidyselect::any_of(weighting_scheme),
                                     #                                            tidyselect::any_of(paste0(paste0(weighting_scheme,"_unn_group"),zz))
                                     #                              ) %>%
                                     #                              dplyr::rename_with(.cols=tidyselect::contains(paste0("_unn_group",zz)), ~gsub(paste0("_unn_group",zz),"",.))
                                     #                          }) #%>%
                                     # #purrr::reduce(dplyr::union_all)
                                     #
                                     # raw_tbl<-raw_tbl_all[[1]]
                                     #
                                     # if (length(raw_tbl_all)>1){
                                     #   for (i in 2:length(raw_tbl_all)){
                                     #     raw_tbl<-dplyr::union_all(raw_tbl,raw_tbl_all[[i]]) %>%
                                     #       dplyr::collapse()
                                     #   }
                                     # }



                                     # for (i in weighting_scheme[!weighting_scheme %in% "lumped"]) {
                                     #   raw_tbl<-raw_tbl %>%
                                     #     dplyr::mutate(across(any_of(loi_cols), ~.*(dbplyr::sql(i)),.names=paste0("{.col}_",i))) %>%
                                     #     dplyr::collapse()
                                     # }

                                     # raw_tbl<- raw_tbl %>%
                                     #   dplyr::compute(temporary=T)

                                     # Calculate Summaries ------------------------------------------------

                                     weighted_mean_out<-NULL
                                     lumped_mean_out<-NULL
                                     weighted_sd_out<-NULL
                                     lumped_sd_out<-NULL
                                     min_out<-NULL
                                     max_out<-NULL
                                     count_out<-NULL
                                     median_out<-NULL
                                     sum_out<-NULL

                                     # Lumped Summaries --------------------------------------------------------

                                     if ("lumped" %in% weighting_scheme | any(loi_meta$loi_type=="cat_rast")) {
                                       lumped_mean_out<-raw_tbl %>%
                                         dplyr::select(tidyselect::any_of(loi_cols)) %>%
                                         dplyr::mutate(lumped=1) %>%
                                         dplyr::summarise(
                                           dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                                         ~sum(.,na.rm=T)/sum(!is.na(.),na.rm=T)
                                           ),
                                           dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="cat_rast"]),
                                                         ~sum(.,na.rm=T)/sum(dbplyr::sql("lumped"),na.rm=T)
                                           )
                                         ) %>%
                                         dplyr::rename_with(.cols=tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),~paste0(.,"_lumped_mean")) %>%
                                         dplyr::rename_with(.cols=tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="cat_rast"]),~paste0(.,"_lumped_prop")) %>%
                                         dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.))) %>%
                                         dplyr::collect()

                                     }

                                     if ("lumped" %in% weighting_scheme & any(loi_numeric_stats=="min")){
                                       min_out<-raw_tbl %>%
                                         dplyr::select(tidyselect::any_of(loi_cols)) %>%
                                         dplyr::summarise(dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),~min(.,na.rm=T))) %>%
                                         dplyr::rename_with(~paste0(.,"_lumped_min"))%>%
                                         dplyr::collect()

                                     }
                                     if ("lumped" %in% weighting_scheme & any(loi_numeric_stats=="max")){
                                       max_out<-raw_tbl %>%
                                         dplyr::select(tidyselect::any_of(loi_cols)) %>%
                                         dplyr::summarise(dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),~max(.,na.rm=T))) %>%
                                         dplyr::rename_with(~paste0(.,"_lumped_max"))%>%
                                         dplyr::collect()

                                     }
                                     if ("lumped" %in% weighting_scheme & any(loi_numeric_stats=="count")){
                                       count_out<-raw_tbl %>%
                                         dplyr::select(tidyselect::any_of(loi_cols)) %>%
                                         dplyr::summarise(dplyr::across(tidyselect::everything(),~sum(!is.na(.),na.rm=T))) %>%
                                         dplyr::rename_with(~paste0(.,"_lumped_count"))%>%
                                         dplyr::collect()

                                     }

                                     if ("lumped" %in% weighting_scheme & any(loi_numeric_stats=="sum")){
                                       sum_out<-raw_tbl %>%
                                         dplyr::select(tidyselect::any_of(loi_cols)) %>%
                                         dplyr::summarise(dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),~sum(.,na.rm=T))) %>%
                                         dplyr::rename_with(~paste0(.,"_lumped_sum"))%>%
                                         dplyr::collect()
                                     }

                                     if ("lumped" %in% weighting_scheme & any(loi_numeric_stats=="median")){
                                       median_out<-raw_tbl %>%
                                         dplyr::select(tidyselect::any_of(loi_cols)) %>%
                                         dplyr::summarise(dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),~median(.,na.rm=T))) %>%
                                         dplyr::rename_with(~paste0(.,"_lumped_median"))%>%
                                         dplyr::collect()
                                     }
                                     if ("lumped" %in% weighting_scheme & any(loi_numeric_stats %in% c("sd","stdev"))){
                                       lumped_sd_out<-raw_tbl %>%
                                         dplyr::select(tidyselect::any_of(loi_cols)) %>%
                                         dplyr::summarise(dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),~stats::sd(.,na.rm=T))) %>%
                                         dplyr::rename_with(~paste0(.,"_lumped_sd"))%>%
                                         dplyr::collect()
                                     }

                                     # Weighted summaries -----------------------------------------------------------

                                     if (length(weighting_scheme[!weighting_scheme %in% "lumped"])>0 & (any(loi_numeric_stats %in% c("mean"))|any(loi_meta$loi_type=="cat_rast"))) {
                                       weighted_mean_out<-raw_tbl %>%
                                         dplyr::summarize(
                                           dplyr::across(tidyselect::ends_with(paste0("_iFLS")),~sum(.,na.rm=T)/sum(dbplyr::sql("iFLS"),na.rm=T)),
                                           dplyr::across(tidyselect::ends_with(paste0("_HAiFLS")),~sum(.,na.rm=T)/sum(dbplyr::sql("HAiFLS"),na.rm=T)),
                                           dplyr::across(tidyselect::ends_with(paste0("_iFLO")),~sum(.,na.rm=T)/sum(dbplyr::sql("iFLO"),na.rm=T)),
                                           dplyr::across(tidyselect::ends_with(paste0("_HAiFLO")),~sum(.,na.rm=T)/sum(dbplyr::sql("HAiFLO"),na.rm=T))
                                         ) %>%
                                         dplyr::rename_with(.cols=tidyselect::starts_with(paste0(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"],"_")),~paste0(.,"_mean")) %>%
                                         dplyr::rename_with(.cols=tidyselect::starts_with(paste0(loi_meta$loi_var_nms[loi_meta$loi_type=="cat_rast"],"_")),~paste0(.,"_prop")) %>%
                                         dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.))) %>%
                                         dplyr::collect()
                                     }


                                     if (length(weighting_scheme[!weighting_scheme %in% "lumped"])>0 &
                                         any(loi_numeric_stats %in% c("sd","stdev")) &
                                         any(loi_meta$loi_type=="num_rast")
                                     ) {

                                       weighted_sd_out<-raw_tbl %>%
                                         dplyr::select(
                                           tidyselect::starts_with(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                           tidyselect::any_of(weighting_scheme)
                                         ) %>%
                                         dplyr::mutate(dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                                                     ~(dbplyr::sql("iFLS") * ((.-(sum(.,na.rm=T)/sum(dbplyr::sql("iFLS"),na.rm=T)))^2)),
                                                                     .names="{.col}_iFLS_term1"),
                                                       dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                                                     ~ ((sum(dbplyr::sql("iFLS")!=0,na.rm=T)-1)/sum(dbplyr::sql("iFLS")!=0,na.rm=T)) * sum(dbplyr::sql("iFLS"),na.rm=T),
                                                                     .names="{.col}_iFLS_term2"
                                                       )) %>%
                                         dplyr::mutate(dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                                                     ~(dbplyr::sql("HAiFLS") * ((.-(sum(.,na.rm=T)/sum(dbplyr::sql("HAiFLS"),na.rm=T)))^2)),
                                                                     .names="{.col}_HAiFLS_term1"),
                                                       dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                                                     ~ ((sum(dbplyr::sql("HAiFLS")!=0,na.rm=T)-1)/sum(dbplyr::sql("HAiFLS")!=0,na.rm=T)) * sum(dbplyr::sql("HAiFLS"),na.rm=T),
                                                                     .names="{.col}_HAiFLS_term2"
                                                       )) %>%
                                         dplyr::mutate(dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                                                     ~(dbplyr::sql("iFLO") * ((.-(sum(.,na.rm=T)/sum(dbplyr::sql("iFLO"),na.rm=T)))^2)),
                                                                     .names="{.col}_iFLO_term1"),
                                                       dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                                                     ~ ((sum(dbplyr::sql("iFLO")!=0,na.rm=T)-1)/sum(dbplyr::sql("iFLO")!=0,na.rm=T)) * sum(dbplyr::sql("iFLO"),na.rm=T),
                                                                     .names="{.col}_iFLO_term2"
                                                       )) %>%
                                         dplyr::mutate(dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                                                     ~(dbplyr::sql("HAiFLO") * ((.-(sum(.,na.rm=T)/sum(dbplyr::sql("HAiFLO"),na.rm=T)))^2)),
                                                                     .names="{.col}_HAiFLO_term1"),
                                                       dplyr::across(tidyselect::any_of(loi_meta$loi_var_nms[loi_meta$loi_type=="num_rast"]),
                                                                     ~ ((sum(dbplyr::sql("HAiFLO")!=0,na.rm=T)-1)/sum(dbplyr::sql("HAiFLO")!=0,na.rm=T)) * sum(dbplyr::sql("HAiFLO"),na.rm=T),
                                                                     .names="{.col}_HAiFLO_term2"
                                                       )) %>%
                                         dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),~sum(.,na.rm=T)),
                                                          dplyr::across(tidyselect::ends_with("_term2"),~.[1])
                                         ) %>%
                                         dplyr::collect() %>%
                                         # The below is some rearranging
                                         tidyr::pivot_longer(cols=c(tidyselect::everything())) %>%
                                         dplyr::mutate(attr=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_|_iFLO_|_HAiFLO_",2)[,1],
                                                       term=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_|_iFLO_|_HAiFLO_",2)[,2]) %>%
                                         dplyr::rowwise() %>%
                                         dplyr::mutate(hw=gsub(paste0(attr,"_","|","_",term,""),"",name)) %>%
                                         dplyr::ungroup() %>%
                                         dplyr::mutate(name=paste0(attr,"_",hw,"_sd")) %>%
                                         dplyr::group_by(name) %>%
                                         dplyr::summarize(sd=sqrt(value[term=="term1"]/value[term=="term2"])) %>%
                                         dplyr::ungroup() %>%
                                         tidyr::pivot_wider(names_from = name,values_from=sd)

                                     }


                                     p()

                                     final_list<-list(
                                       weighted_mean_out,
                                       lumped_mean_out,
                                       weighted_sd_out,
                                       lumped_sd_out,
                                       min_out,
                                       max_out,
                                       count_out,
                                       median_out,
                                       sum_out
                                     )
                                     final_list<-final_list[!sapply(final_list,is.null)]

                                     final_out<-dplyr::bind_cols(tibble::tibble(pour_point_id=yy),final_list)

                                     final_out<-final_out %>%
                                       dplyr::select(
                                         pour_point_id,
                                         tidyselect::contains(loi_meta$loi_var_nms)
                                       )

                                     DBI::dbRemoveTable(con,paste0("temp_",yy))

                                     return(final_out)

                                   })
              DBI::dbDisconnect(con)

              return(out)

            })
        ))

  })

  #browser()
  final_out<-attr_out %>%
    dplyr::select(out) %>%
    tidyr::unnest(cols = c(out))

  # Important to keep this
  if (!store_iDW) unlink(gsub(basename(iDW_file),"",iDW_file),force=T,recursive = T)

  unlink(temp_dir_sub,recursive = T,force = T)

  target_IDs_out<-target_id_fun(
    db_fp=db_loc,
    sample_points=sample_points,
    link_id=link_id,
    segment_whole=F,
    target_o_type=target_o_type
  )

  #browser()
  final_out<-target_IDs_out %>%
    dplyr::left_join(
      final_out %>% dplyr::mutate(pour_point_id=as.character(pour_point_id)),
      by=c("link_id"="pour_point_id")
    )

  data.table::fwrite(x=final_out,
                     file=out_filename,
                     buffMB = 128L,
                     nThread = 1,
                     showProgress = F)

  return(final_out)

}




#' @importFrom carrier crate
#' @importFrom data.table fwrite fread
#' @importFrom DBI dbConnect dbListTables dbDisconnect
#' @importFrom dplyr group_by ungroup mutate n filter select row_number tbl rows_update copy_to
#' @importFrom exactextractr exact_extract
#' @importFrom furrr future_pmap furrr_options
#' @importFrom future nbrOfWorkers availableCores plan tweak multisession future futureOf resolved
#' @importFrom progressr with_progress progressor
#' @importFrom purrr map2
#' @importFrom RSQLite SQLite
#' @importFrom sf read_sf
#' @importFrom terra rast
#' @importFrom tidyr nest
#' @importFrom tidyselect everything
#' @importFrom carrier crate
#' @importFrom data.table fwrite fread
#' @importFrom DBI dbConnect dbListTables dbExecute dbDisconnect
#' @importFrom dplyr group_by ungroup mutate n filter select row_number tbl rows_update copy_to
#' @importFrom exactextractr exact_extract
#' @importFrom furrr future_pmap furrr_options
#' @importFrom future nbrOfWorkers availableCores plan tweak multisession future futureOf resolved
#' @importFrom progressr with_progress progressor
#' @importFrom purrr map2
#' @importFrom RSQLite SQLite
#' @importFrom sf read_sf
#' @importFrom terra rast
#' @importFrom tidyr nest
#' @importFrom tidyselect everything
extract_raster_vals<-function(
    input,
    raw_file,
    iDW_file,
    loi_file,
    weighting_scheme,
    loi_cols,
    subb_IDs,
    temp_dir_sub,
    verbose
){
  if (inherits(loi_file,"ihydro")) loi_file<-loi_file$outfile
  if (inherits(iDW_file,"ihydro")) iDW_file<-iDW_file$outfile
  if (inherits(raw_file,"ihydro")) raw_file<-raw_file$outfile

  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)
  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)
  if (n_cores==0) n_cores<-1


  n_cores_2<-n_cores

  if (n_cores>1) {
    n_cores_2<-n_cores_2-1
    oplan <- future::plan(list(future::tweak(future::multisession, workers = 2), future::tweak(future::multisession, workers = n_cores_2)))
    on.exit(future::plan(oplan), add = TRUE)
  }


  # Check tables in database -------------------------------------------------

  # con<-DBI::dbConnect(RSQLite::SQLite(), raw_file,cache_size=200000)
  #
  # tl<-DBI::dbListTables(con)
  #
  # are_equal<-F
  # if (paste0(yy,"_raw") %in% tl) {
  #   dbt<-dplyr::collect(dplyr::tbl(con,paste0(yy,"_raw")))
  #
  #   are_equal<-all.equal(dbt,ot,tolerance=sqrt(.Machine$double.eps))
  # }
  #
  # DBI::dbDisconnect(con)
  #
  # if (!are_equal) data.table::fwrite(ot,file=file.path(temp_dir_sub,paste0(yy,"_raw.csv")))


  # Save extracted values to SQLite -----------------------------------------

  if (verbose) message("Extracting Raster Values")

  future_proc<-future::future({
    ot<-subb_IDs %>%
      dplyr::group_by(origin_link_id) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(core=rep(1:n_cores,length.out=dplyr::n())) %>%
      dplyr::group_by(core) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        out=furrr::future_pmap(list(x=data,
                                    input=list(as.ihydro(input$outfile)),
                                    raw_file=list(raw_file),
                                    iDW_file=list(iDW_file),
                                    loi_file=list(loi_file),
                                    weighting_scheme=list(weighting_scheme),
                                    loi_cols=list(loi_cols),
                                    temp_dir_sub=list(temp_dir_sub)
        ),
        .options = furrr::furrr_options(globals = F),
        carrier::crate(
          function(x,
                   input,
                   raw_file,
                   iDW_file,
                   loi_file,
                   weighting_scheme,
                   loi_cols,
                   temp_dir_sub){
            options(dplyr.summarise.inform = FALSE)
            options(scipen = 999)
            `%>%` <- magrittr::`%>%`

            iDWs_rasts<-NULL
            loi_rasts<-terra::rast(loi_file,loi_cols)
            if (length(weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")])>0){
              iDWs_rasts<-terra::rast(iDW_file,weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")])
            }

            purrr::map2(x$data,x$origin_link_id,function(y,yy){
              iDWo_rasts<-NULL

              input_poly<-sf::read_sf(input$outfile,
                                      query=paste0("SELECT `link_id`, `geom` FROM `Subbasins_poly` WHERE (`link_id` IN (",paste0("'",yy,collapse = "',"),"'))")
              )

              if (length(weighting_scheme[weighting_scheme %in% c("iFLO","HAiFLO")])>0){
                iDWo_rasts<-terra::rast(iDW_file,
                                        unlist(purrr:::map(unique(y$unn_group),
                                                           ~paste0(
                                                             paste0(
                                                               weighting_scheme[weighting_scheme %in% c("iFLO","HAiFLO")],
                                                               "_unn_group"
                                                             ),
                                                             .
                                                           )))
                )
              }


              input_rasts<-c(loi_rasts,iDWs_rasts,iDWo_rasts)

              ot<-exactextractr::exact_extract(input_rasts,input_poly)[[1]] %>%
                dplyr::filter(coverage_fraction>0.5) %>%
                dplyr::select(-coverage_fraction) %>%
                dplyr::mutate(.rowid=dplyr::row_number()) %>%
                dplyr::mutate(link_id=yy) %>%
                dplyr::select(.rowid,link_id,tidyselect::everything())

              data.table::fwrite(ot,file=file.path(temp_dir_sub,paste0(yy,"_raw.csv")))

              # con<-DBI::dbConnect(RSQLite::SQLite(), raw_file,cache_size=200000)
              #
              # tl<-DBI::dbListTables(con)
              #
              # are_equal<-F
              # if (paste0(yy,"_raw") %in% tl) {
              #   dbt<-dplyr::collect(dplyr::tbl(con,paste0(yy,"_raw")))
              #
              #   are_equal<-all.equal(dbt,ot,tolerance=sqrt(.Machine$double.eps))
              # }
              #
              # DBI::dbDisconnect(con)
              #
              # if (!are_equal) data.table::fwrite(ot,file=file.path(temp_dir_sub,paste0(yy,"_raw.csv")))

              return(NULL)

            })
          }
        )))
    return(ot)
  })

  # Write raw data to database ----------------------------------------------

  #browser()
  progressr::with_progress(enable=T,{
    p <- progressr::progressor(steps = (length(unique(subb_IDs$origin_link_id))))

    future_proc_status <- future::futureOf(future_proc)

    while(!future::resolved(future_proc_status)){
      Sys.sleep(0.2)

      fl<-list.files(temp_dir_sub,"_raw.csv",full.names = T)

      for (x in fl) {

        y<-try(data.table::fread(x),silent = T)
        if (inherits(y,"try-error") || nrow(y)==0 || is.null(y)) next()

        con<-DBI::dbConnect(RSQLite::SQLite(), raw_file,cache_size=200000)

        tl<-DBI::dbListTables(con)

        if (paste0(gsub(".csv","",basename(x))) %in% tl){
          ot<-dplyr::tbl(con,paste0(gsub(".csv","",basename(x))))

          if (any(!colnames(y) %in% colnames(ot))) {

            add_col_statement<-paste0(
              paste0("ALTER TABLE '",gsub(".csv","",basename(x)),"' ADD "),
              paste0(
                colnames(y)[!colnames(y) %in% colnames(ot)]
              ),
              " REAL "
            )

            for (xx in add_col_statement) DBI::dbExecute(con,xx)
          }

          ot<-dplyr::tbl(con,paste0(gsub(".csv","",basename(x))))

          ot<-dplyr::rows_update(
            x=ot,
            y=y,
            by = ".rowid",
            unmatched = "ignore",
            #...,
            copy = T,
            in_place = T
          )
        } else {
          ot<-dplyr::copy_to(
            dest=con,
            df=y,
            name=paste0(gsub(".csv","",basename(x))),
            temporary=F,
            overwrite=T,
            unique_indexes = ".rowid",
            analyze=F
          )

        }

        DBI::dbDisconnect(con)

        file.remove(x)

        p()
      }

    }

    Sys.sleep(2)

    fl<-list.files(temp_dir_sub,"_raw.csv",full.names = T)

    for (x in fl) {
      Sys.sleep(0.2)

      y<-try(data.table::fread(x),silent = T)
      if (inherits(y,"try-error")  || nrow(y)==0 || is.null(y)) stop(paste0("Error reading: ",x))

      con<-DBI::dbConnect(RSQLite::SQLite(), raw_file,cache_size=200000)

      tl<-DBI::dbListTables(con)

      if (paste0(gsub(".csv","",basename(x))) %in% tl){
        ot<-dplyr::tbl(con,paste0(gsub(".csv","",basename(x))))

        if (any(!colnames(y) %in% colnames(ot))) {

          add_col_statement<-paste0(
            paste0("ALTER TABLE '",gsub(".csv","",basename(x)),"' ADD "),
            paste0(
              colnames(y)[!colnames(y) %in% colnames(ot)]
            ),
            " REAL "
          )

          for (xx in add_col_statement) DBI::dbExecute(con,xx)
        }

        ot<-dplyr::tbl(con,paste0(gsub(".csv","",basename(x))))

        ot<-dplyr::rows_update(
          x=ot,
          y=y,
          by = ".rowid",
          unmatched = "ignore",
          #...,
          copy = T,
          in_place = T
        )
      } else {
        ot<-dplyr::copy_to(
          dest=con,
          df=y,
          name=paste0(gsub(".csv","",basename(x))),
          temporary=F,
          overwrite=T,
          unique_indexes = ".rowid",
          analyze=F
        )

      }

      DBI::dbDisconnect(con)

      file.remove(x)

      p()
    }

  })

  return(NULL)

}


extract_raster_attributes<-function(
    input,
    iDW_file,
    loi_file,
    weighting_scheme,
    loi_cols,
    subb_IDs,
    loi_numeric_stats,
    verbose
){
  if (inherits(loi_file,"ihydro")) loi_file<-loi_file$outfile
  if (inherits(iDW_file,"ihydro")) iDW_file<-iDW_file$outfile

  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)
  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)
  if (n_cores==0) n_cores<-1

  loi_meta<-sf::read_sf(loi_file,"loi_meta")
  loi_meta<-loi_meta %>% dplyr::filter(loi_var_nms %in% loi_cols)

  if (verbose) message("Calculating Catchments")

  input_poly<-ihydro::get_catchment(input,link_id=unique(subb_IDs$pour_point_id))

  if (verbose) message("Calculating Weighted Attributes")

  progressr::with_progress(enable=T,{
    p <- progressr::progressor(steps = (length(unique(subb_IDs$unn_group))))

    ot<-subb_IDs %>%
      dplyr::group_by(unn_group) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(core=rep(1:n_cores,length.out=dplyr::n())) %>%
      dplyr::group_by(core) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        out=furrr::future_pmap(
          #out=purrr::pmap(
          list(x=data,
               input=list(as.ihydro(input$outfile)),
               iDW_file=list(iDW_file),
               loi_file=list(loi_file),
               loi_meta=list(loi_meta),
               weighting_scheme=list(weighting_scheme),
               loi_numeric_stats=list(loi_numeric_stats),
               loi_cols=list(loi_cols),
               p=list(p)
          ),
          .options = furrr::furrr_options(globals = F),
          carrier::crate(
            function(x,
                     input,
                     iDW_file,
                     loi_file,
                     loi_meta,
                     weighting_scheme,
                     loi_numeric_stats,
                     loi_cols,
                     p,
                     temp_dir_sub){
              options(dplyr.summarise.inform = FALSE)
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`

              iDWs_rasts<-NULL
              loi_rasts<-terra::rast(loi_file,loi_cols)
              if (length(weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")])>0){
                iDWs_rasts<-terra::rast(iDW_file,weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")])
              }

              purrr::map2_dfr(x$data,x$unn_group,function(y,yy){
                iDWo_rasts<-NULL

                input_poly<-ihydro::get_catchment(input,link_id=unique(y$pour_point_id))

                if (length(weighting_scheme[weighting_scheme %in% c("iFLO","HAiFLO")])>0){
                  iDWo_rasts<-terra::rast(iDW_file,
                                          unlist(purrr:::map(unique(yy),
                                                             ~paste0(
                                                               paste0(
                                                                 weighting_scheme[weighting_scheme %in% c("iFLO","HAiFLO")],
                                                                 "_unn_group"
                                                               ),
                                                               .
                                                             )))
                  )

                  names(iDWo_rasts)<-gsub(paste0("_unn_group",yy),"",names(iDWo_rasts))

                }


                input_rasts<-c(loi_rasts,iDWs_rasts,iDWo_rasts)


                ot<-exactextractr::exact_extract(input_rasts,
                                                 input_poly,
                                                 include_cols="link_id",
                                                 summarize_df=T,
                                                 progress=F,
                                                 force_df=T,
                                                 fun=function(df,
                                                              weighting_scheme2=weighting_scheme,
                                                              loi_meta2=loi_meta,
                                                              loi_cols2=loi_cols,
                                                              loi_numeric_stats2=loi_numeric_stats
                                                 ){
                                                   #browser()
                                                   options(dplyr.summarise.inform = FALSE)
                                                   options(scipen = 999)
                                                   `%>%` <- magrittr::`%>%`

                                                   #con<-DBI::dbConnect(RSQLite::SQLite(),temp_dir_sub2)

                                                   pour_point_id<-df$link_id[[1]]

                                                   df<-df %>%
                                                     dplyr::filter(coverage_fraction>0.5) %>%
                                                     dplyr::select(-coverage_fraction)

                                                   df<-df %>%
                                                     dplyr::mutate(dplyr::across(tidyselect::where(is.numeric),~ifelse(is.nan(.),NA_real_,.)))

                                                   df<-df %>%
                                                     dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2),
                                                                          ~.*(!!rlang::sym("iFLS")),.names=paste0("{.col}_","iFLS")))
                                                   df<-df %>%
                                                     dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2),
                                                                          ~.*(!!rlang::sym("iFLO")),.names=paste0("{.col}_","iFLO")))
                                                   df<-df %>%
                                                     dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2),
                                                                          ~.*(!!rlang::sym("HAiFLS")),.names=paste0("{.col}_","HAiFLS")))
                                                   df<-df %>%
                                                     dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2),
                                                                          ~.*(!!rlang::sym("HAiFLO")),.names=paste0("{.col}_","HAiFLO")))

                                                   # for (ii in weighting_scheme2[!weighting_scheme2 %in% c("lumped")]) {
                                                   #   df<-df %>%
                                                   #     dplyr::mutate(across(tidyselect::any_of(loi_cols2),
                                                   #                          ~.*(!!rlang::sym(ii)),.names=paste0("{.col}_",ii)))
                                                   # }

                                                   # raw_tbl<-df2 %>%
                                                   #   dplyr::copy_to(dest=con,
                                                   #                  df=.,
                                                   #                  name=paste0("raw_",df$link_id[[1]]))


                                                   weighted_mean_out<-NULL
                                                   lumped_mean_out<-NULL
                                                   weighted_sd_out<-NULL
                                                   lumped_sd_out<-NULL
                                                   min_out<-NULL
                                                   max_out<-NULL
                                                   count_out<-NULL
                                                   median_out<-NULL
                                                   sum_out<-NULL


                                                   # Lumped Summaries --------------------------------------------------------

                                                   if ("lumped" %in% weighting_scheme2 | any(loi_meta2$loi_type=="cat_rast")) {
                                                     lumped_mean_out<-df %>%
                                                       dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                       #dplyr::mutate(lumped=1) %>%
                                                       dplyr::summarise(
                                                         dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),
                                                                       ~sum(.,na.rm=T)/sum(!is.na(.),na.rm=T)
                                                         ),
                                                         dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="cat_rast"]),
                                                                       ~sum(.,na.rm=T)/dplyr::n()#sum(!!rlang::sym("lumped"),na.rm=T)
                                                         )
                                                       ) %>%
                                                       dplyr::rename_with(.cols=tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta$loi_type=="num_rast"]),~paste0(.,"_lumped_mean")) %>%
                                                       dplyr::rename_with(.cols=tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta$loi_type=="cat_rast"]),~paste0(.,"_lumped_prop")) %>%
                                                       dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.))) #%>%
                                                     #dplyr::collect()

                                                   }

                                                   if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="min")){
                                                     min_out<-df %>%
                                                       dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                       dplyr::summarise(dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),~min(.,na.rm=T))) %>%
                                                       dplyr::rename_with(~paste0(.,"_lumped_min"))#%>%
                                                     #dplyr::collect()

                                                   }
                                                   if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="max")){
                                                     max_out<-df %>%
                                                       dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                       dplyr::summarise(dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),~max(.,na.rm=T))) %>%
                                                       dplyr::rename_with(~paste0(.,"_lumped_max"))#%>%
                                                     #dplyr::collect()

                                                   }
                                                   if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="count")){
                                                     count_out<-df %>%
                                                       dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                       dplyr::summarise(dplyr::across(tidyselect::everything(),~sum(!is.na(.),na.rm=T))) %>%
                                                       dplyr::rename_with(~paste0(.,"_lumped_count"))#%>%
                                                     #dplyr::collect()

                                                   }

                                                   if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="sum")){
                                                     sum_out<-df %>%
                                                       dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                       dplyr::summarise(dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),~sum(.,na.rm=T))) %>%
                                                       dplyr::rename_with(~paste0(.,"_lumped_sum"))#%>%
                                                     #dplyr::collect()
                                                   }

                                                   if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="median")){
                                                     median_out<-df %>%
                                                       dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                       dplyr::summarise(dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),~stats::median(.,na.rm=T))) %>%
                                                       dplyr::rename_with(~paste0(.,"_lumped_median"))#%>%
                                                     #dplyr::collect()
                                                   }
                                                   if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2 %in% c("sd","stdev"))){
                                                     lumped_sd_out<-df %>%
                                                       dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                       dplyr::summarise(dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),~stats::sd(.,na.rm=T))) %>%
                                                       dplyr::rename_with(~paste0(.,"_lumped_sd"))#%>%
                                                     #dplyr::collect()
                                                   }

                                                   # Weighted summaries -----------------------------------------------------------

                                                   if (length(weighting_scheme2[!weighting_scheme2 %in% "lumped"])>0 & (any(loi_numeric_stats2 %in% c("mean"))|any(loi_meta2$loi_type=="cat_rast"))) {
                                                     weighted_mean_out<-df %>%
                                                       dplyr::summarize(
                                                         dplyr::across(tidyselect::ends_with(paste0("_iFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)),
                                                         dplyr::across(tidyselect::ends_with(paste0("_HAiFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)),
                                                         dplyr::across(tidyselect::ends_with(paste0("_iFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T)),
                                                         dplyr::across(tidyselect::ends_with(paste0("_HAiFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T))
                                                       ) %>%
                                                       dplyr::rename_with(.cols=tidyselect::starts_with(paste0(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"],"_")),~paste0(.,"_mean")) %>%
                                                       dplyr::rename_with(.cols=tidyselect::starts_with(paste0(loi_meta2$loi_var_nms[loi_meta2$loi_type=="cat_rast"],"_")),~paste0(.,"_prop")) %>%
                                                       dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.))) #%>%
                                                     #dplyr::collect()
                                                   }


                                                   if (length(weighting_scheme2[!weighting_scheme2 %in% "lumped"])>0 &
                                                       any(loi_numeric_stats2 %in% c("sd","stdev")) &
                                                       any(loi_meta2$loi_type=="num_rast")
                                                   ) {

                                                     weighted_sd_out<-df %>%
                                                       dplyr::select(
                                                         tidyselect::starts_with(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),
                                                         tidyselect::any_of(weighting_scheme2)
                                                       ) %>%
                                                       dplyr::mutate(dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),
                                                                                   ~(!!rlang::sym("iFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)))^2)),
                                                                                   .names="{.col}_iFLS_term1"),
                                                                     dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),
                                                                                   ~ ((sum(!!rlang::sym("iFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLS")!=0,na.rm=T)) * sum(!!rlang::sym("iFLS"),na.rm=T),
                                                                                   .names="{.col}_iFLS_term2"
                                                                     )) %>%
                                                       dplyr::mutate(dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),
                                                                                   ~(!!rlang::sym("HAiFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)))^2)),
                                                                                   .names="{.col}_HAiFLS_term1"),
                                                                     dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),
                                                                                   ~ ((sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLS"),na.rm=T),
                                                                                   .names="{.col}_HAiFLS_term2"
                                                                     )) %>%
                                                       dplyr::mutate(dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),
                                                                                   ~(!!rlang::sym("iFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T)))^2)),
                                                                                   .names="{.col}_iFLO_term1"),
                                                                     dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),
                                                                                   ~ ((sum(!!rlang::sym("iFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLO")!=0,na.rm=T)) * sum(!!rlang::sym("iFLO"),na.rm=T),
                                                                                   .names="{.col}_iFLO_term2"
                                                                     )) %>%
                                                       dplyr::mutate(dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),
                                                                                   ~(!!rlang::sym("HAiFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T)))^2)),
                                                                                   .names="{.col}_HAiFLO_term1"),
                                                                     dplyr::across(tidyselect::any_of(loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]),
                                                                                   ~ ((sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLO"),na.rm=T),
                                                                                   .names="{.col}_HAiFLO_term2"
                                                                     )) %>%
                                                       dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),~sum(.,na.rm=T)),
                                                                        dplyr::across(tidyselect::ends_with("_term2"),~.[1])
                                                       ) %>%
                                                       #dplyr::collect() %>%
                                                       # The below is some rearranging
                                                       tidyr::pivot_longer(cols=c(tidyselect::everything())) %>%
                                                       dplyr::mutate(attr=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_|_iFLO_|_HAiFLO_",2)[,1],
                                                                     term=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_|_iFLO_|_HAiFLO_",2)[,2]) %>%
                                                       dplyr::rowwise() %>%
                                                       dplyr::mutate(hw=gsub(paste0(attr,"_","|","_",term,""),"",name)) %>%
                                                       dplyr::ungroup() %>%
                                                       dplyr::mutate(name=paste0(attr,"_",hw,"_sd")) %>%
                                                       dplyr::group_by(name) %>%
                                                       dplyr::summarize(sd=sqrt(value[term=="term1"]/value[term=="term2"])) %>%
                                                       dplyr::ungroup() %>%
                                                       tidyr::pivot_wider(names_from = name,values_from=sd)

                                                   }


                                                   final_list<-list(
                                                     weighted_mean_out,
                                                     lumped_mean_out,
                                                     weighted_sd_out,
                                                     lumped_sd_out,
                                                     min_out,
                                                     max_out,
                                                     count_out,
                                                     median_out,
                                                     sum_out
                                                   )
                                                   final_list<-final_list[!sapply(final_list,is.null)]

                                                   final_out<-dplyr::bind_cols(tibble::tibble(pour_point_id=pour_point_id),final_list)

                                                   final_out<-final_out %>%
                                                     dplyr::select(
                                                       pour_point_id,
                                                       tidyselect::contains(loi_meta2$loi_var_nms)
                                                     )

                                                   #DBI::dbDisconnect(con)

                                                   return(final_out)

                                                 }

                )

                p()
                return(ot)

              })
            }
          )))

  })

  out<-ot %>%
    dplyr::select(out) %>%
    tidyr::unnest(cols=out)

  return(out)

}
