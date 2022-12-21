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

  #browser()

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

  o1<-extract_raster_attributes(
    input=input,
    iDW_file=iDW_file,
    loi_file=loi_file,
    weighting_scheme=weighting_scheme,
    loi_numeric_stats=loi_numeric_stats,
    loi_cols=loi_cols,
    subb_IDs=subb_IDs,
    temp_dir_sub=temp_dir_sub,
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
extract_raster_attributes<-function(
    input,
    iDW_file,
    loi_file,
    weighting_scheme,
    loi_cols,
    subb_IDs,
    loi_numeric_stats,
    temp_dir_sub,
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

  #browser()

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
             p=list(p),
             temp_dir_sub=list(temp_dir_sub)
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
            Sys.setenv(GTIFF_DIRECT_IO='ON')
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

              #browser()


              # Summary Function --------------------------------------------------------


              #browser()

              ot<-try(exactextractr::exact_extract(input_rasts,
                                                   input_poly,
                                                   include_cols="link_id",
                                                   summarize_df=T,
                                                   progress=F,
                                                   force_df=T,
                                                   fun=function(df,
                                                                weighting_scheme2=weighting_scheme,
                                                                loi_meta2=loi_meta,
                                                                loi_cols2=loi_cols,
                                                                loi_numeric_stats2=loi_numeric_stats){
                                                     #browser()
                                                     options(dplyr.summarise.inform = FALSE)
                                                     options(scipen = 999)
                                                     `%>%` <- magrittr::`%>%`

                                                     #con<-DBI::dbConnect(RSQLite::SQLite(),temp_dir_sub2)

                                                     loi_meta2<-dplyr::filter(loi_meta2,loi_var_nms == loi_cols2)

                                                     pour_point_id<-df$link_id[[1]]

                                                     #browser()

                                                     df_class<-sapply(df,class)

                                                     df<-df %>%
                                                       dtplyr::lazy_dt() %>%
                                                       dplyr::filter(coverage_fraction>0.5) %>%
                                                       dplyr::select(-coverage_fraction)


                                                     df<-df %>%
                                                       dplyr::mutate(dplyr::across(tidyselect::any_of(names(df_class[df_class=="numeric"])),~ifelse(is.nan(.),NA_real_,.))) #tidyselect::where(is.numeric)

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

                                                     df<-df %>%
                                                       tibble::as_tibble() %>%
                                                       dtplyr::lazy_dt()


                                                     weighted_mean_out<-NULL
                                                     lumped_mean_out<-NULL
                                                     weighted_sd_out<-NULL
                                                     lumped_sd_out<-NULL
                                                     min_out<-NULL
                                                     max_out<-NULL
                                                     count_out<-NULL
                                                     median_out<-NULL
                                                     sum_out<-NULL

                                                     numb_rast<-loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]
                                                     cat_rast<-loi_meta2$loi_var_nms[loi_meta2$loi_type=="cat_rast"]


                                                     # Lumped Summaries --------------------------------------------------------

                                                     if ("lumped" %in% weighting_scheme2 | any(loi_meta2$loi_type=="cat_rast")) {
                                                       lumped_mean_out<-df %>%
                                                         dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                         #dplyr::mutate(lumped=1) %>%
                                                         dplyr::summarise(
                                                           dplyr::across(tidyselect::any_of(numb_rast),
                                                                         ~sum(.,na.rm=T)/sum(!is.na(.),na.rm=T)
                                                           ),
                                                           dplyr::across(tidyselect::any_of(cat_rast),
                                                                         ~sum(.,na.rm=T)/dplyr::n()#sum(!!rlang::sym("lumped"),na.rm=T)
                                                           )
                                                         ) %>%
                                                         tibble::as_tibble()

                                                       if (length(numb_rast)>0) {
                                                         lumped_mean_out<-lumped_mean_out %>%
                                                           dplyr::rename_with(.cols=tidyselect::any_of(numb_rast),~paste0(.,"_lumped_mean"))
                                                       }

                                                       if (length(cat_rast)>0) {
                                                         lumped_mean_out<-lumped_mean_out %>%
                                                           dplyr::rename_with(.cols=tidyselect::any_of(cat_rast),~paste0(.,"_lumped_prop")) %>%
                                                           dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.)))
                                                       }


                                                     }

                                                     if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="min") & length(numb_rast)>0){
                                                       min_out<-df %>%
                                                         dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                         dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~min(.,na.rm=T))) %>%
                                                         tibble::as_tibble()%>%
                                                         dplyr::rename_with(~paste0(.,"_lumped_min"))#%>%
                                                       #dplyr::collect()

                                                     }
                                                     if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="max" & length(numb_rast)>0)){
                                                       max_out<-df %>%
                                                         dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                         dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~max(.,na.rm=T))) %>%
                                                         tibble::as_tibble()%>%
                                                         dplyr::rename_with(~paste0(.,"_lumped_max"))#%>%
                                                       #dplyr::collect()

                                                     }
                                                     if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="count")){
                                                       count_out<-df %>%
                                                         dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                         dplyr::summarise(dplyr::across(tidyselect::everything(),~sum(!is.na(.),na.rm=T))) %>%
                                                         tibble::as_tibble()%>%
                                                         dplyr::rename_with(~paste0(.,"_lumped_count"))#%>%
                                                       #dplyr::collect()

                                                     }

                                                     if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="sum" & length(numb_rast)>0)){
                                                       sum_out<-df %>%
                                                         dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                         dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~sum(.,na.rm=T))) %>%
                                                         tibble::as_tibble()%>%
                                                         dplyr::rename_with(~paste0(.,"_lumped_sum"))#%>%
                                                       #dplyr::collect()
                                                     }

                                                     if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="median" & length(numb_rast)>0)){
                                                       median_out<-df %>%
                                                         dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                         dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~stats::median(.,na.rm=T))) %>%
                                                         tibble::as_tibble()%>%
                                                         dplyr::rename_with(~paste0(.,"_lumped_median"))#%>%
                                                       #dplyr::collect()
                                                     }
                                                     if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2 %in% c("sd","stdev") & length(numb_rast)>0)){
                                                       lumped_sd_out<-df %>%
                                                         dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                         dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~stats::sd(.,na.rm=T))) %>%
                                                         tibble::as_tibble()%>%
                                                         dplyr::rename_with(~paste0(.,"_lumped_sd"))#%>%
                                                       #dplyr::collect()
                                                     }

                                                     # Weighted summaries -----------------------------------------------------------

                                                     if (length(weighting_scheme2[!weighting_scheme2 %in% "lumped"])>0 & (any(loi_numeric_stats2 %in% c("mean")) | length(cat_rast)>0)) {
                                                       weighted_mean_out<-df %>%
                                                         dplyr::summarize(
                                                           dplyr::across(tidyselect::ends_with(paste0("_iFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)),
                                                           dplyr::across(tidyselect::ends_with(paste0("_HAiFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)),
                                                           dplyr::across(tidyselect::ends_with(paste0("_iFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T)),
                                                           dplyr::across(tidyselect::ends_with(paste0("_HAiFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T))
                                                         ) %>%
                                                         tibble::as_tibble()

                                                       if (length(numb_rast)>0) {
                                                         weighted_mean_out<-weighted_mean_out %>%
                                                           dplyr::rename_with(.cols=tidyselect::starts_with(paste0(numb_rast,"_")),~paste0(.,"_mean"))
                                                       }

                                                       if (length(cat_rast)>0) {
                                                         weighted_mean_out<-weighted_mean_out %>%
                                                           dplyr::rename_with(.cols=tidyselect::starts_with(paste0(cat_rast,"_")),~paste0(.,"_prop")) %>%
                                                           dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.)))
                                                       }

                                                     }


                                                     if (length(weighting_scheme2[!weighting_scheme2 %in% "lumped"])>0 &
                                                         any(loi_numeric_stats2 %in% c("sd","stdev")) &
                                                         length(numb_rast)>0
                                                     ) {

                                                       weighted_sd_out<-df %>%
                                                         dplyr::select(
                                                           tidyselect::starts_with(numb_rast),
                                                           tidyselect::any_of(weighting_scheme2)
                                                         ) %>%
                                                         dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                                                                     ~(!!rlang::sym("iFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)))^2)),
                                                                                     .names="{.col}_iFLS_term1"),
                                                                       dplyr::across(tidyselect::any_of(numb_rast),
                                                                                     ~ ((sum(!!rlang::sym("iFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLS")!=0,na.rm=T)) * sum(!!rlang::sym("iFLS"),na.rm=T),
                                                                                     .names="{.col}_iFLS_term2"
                                                                       )) %>%
                                                         dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                                                                     ~(!!rlang::sym("HAiFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)))^2)),
                                                                                     .names="{.col}_HAiFLS_term1"),
                                                                       dplyr::across(tidyselect::any_of(numb_rast),
                                                                                     ~ ((sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLS"),na.rm=T),
                                                                                     .names="{.col}_HAiFLS_term2"
                                                                       )) %>%
                                                         dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                                                                     ~(!!rlang::sym("iFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T)))^2)),
                                                                                     .names="{.col}_iFLO_term1"),
                                                                       dplyr::across(tidyselect::any_of(numb_rast),
                                                                                     ~ ((sum(!!rlang::sym("iFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLO")!=0,na.rm=T)) * sum(!!rlang::sym("iFLO"),na.rm=T),
                                                                                     .names="{.col}_iFLO_term2"
                                                                       )) %>%
                                                         dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                                                                     ~(!!rlang::sym("HAiFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T)))^2)),
                                                                                     .names="{.col}_HAiFLO_term1"),
                                                                       dplyr::across(tidyselect::any_of(numb_rast),
                                                                                     ~ ((sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLO"),na.rm=T),
                                                                                     .names="{.col}_HAiFLO_term2"
                                                                       )) %>%
                                                         dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),~sum(.,na.rm=T)),
                                                                          dplyr::across(tidyselect::ends_with("_term2"),~.[1])
                                                         ) %>%
                                                         tibble::as_tibble() %>%
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
                                                     final_list<-final_list[sapply(final_list,nrow)>0]

                                                     final_out<-dplyr::bind_cols(tibble::tibble(pour_point_id=pour_point_id),final_list)

                                                     final_out<-final_out %>%
                                                       dplyr::select(
                                                         pour_point_id,
                                                         tidyselect::contains(loi_meta2$loi_var_nms)
                                                       )

                                                     #DBI::dbDisconnect(con)
                                                     rm(df)
                                                     gg<-gc()

                                                     return(final_out)
                                                   }

              ),silent=F)

              #browser()
              if (inherits(ot,"try-error")) {

                message("Insufficient memory to hold loi and IDW layers, calculation will be slower")

                rm(ot)
                gg<-gc()

                if (!dir.exists(temp_dir_sub)) dir.create(temp_dir_sub,recursive = T)

                #browser()

                ot<-try(
                  purrr::map(split(input_poly,seq_along(input_poly$link_id)),
                             function(xx) {
                               purrr::map(loi_cols,function(x){
                                 ot<-try(exactextractr::exact_extract(terra::subset(input_rasts,c(x,weighting_scheme)[c(x,weighting_scheme) %in% names(input_rasts)]),
                                                                      xx,
                                                                      max_cells_in_memory=(3e+07)/4,
                                                                      include_cols="link_id",
                                                                      summarize_df=T,
                                                                      progress=F,
                                                                      force_df=T,
                                                                      fun=function(df,
                                                                                   weighting_scheme2=weighting_scheme,
                                                                                   loi_meta2=loi_meta,
                                                                                   loi_cols2=x,
                                                                                   loi_numeric_stats2=loi_numeric_stats
                                                                      ){
                                                                        #browser()
                                                                        options(dplyr.summarise.inform = FALSE)
                                                                        options(scipen = 999)
                                                                        `%>%` <- magrittr::`%>%`

                                                                        #con<-DBI::dbConnect(RSQLite::SQLite(),temp_dir_sub2)

                                                                        loi_meta2<-dplyr::filter(loi_meta2,loi_var_nms == loi_cols2)

                                                                        pour_point_id<-df$link_id[[1]]

                                                                        #browser()

                                                                        df_class<-sapply(df,class)

                                                                        df<-df %>%
                                                                          dtplyr::lazy_dt() %>%
                                                                          dplyr::filter(coverage_fraction>0.5) %>%
                                                                          dplyr::select(-coverage_fraction)


                                                                        df<-df %>%
                                                                          dplyr::mutate(dplyr::across(tidyselect::any_of(names(df_class[df_class=="numeric"])),~ifelse(is.nan(.),NA_real_,.))) #tidyselect::where(is.numeric)

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

                                                                        df<-df %>%
                                                                          tibble::as_tibble() %>%
                                                                          dtplyr::lazy_dt()


                                                                        weighted_mean_out<-NULL
                                                                        lumped_mean_out<-NULL
                                                                        weighted_sd_out<-NULL
                                                                        lumped_sd_out<-NULL
                                                                        min_out<-NULL
                                                                        max_out<-NULL
                                                                        count_out<-NULL
                                                                        median_out<-NULL
                                                                        sum_out<-NULL

                                                                        numb_rast<-loi_meta2$loi_var_nms[loi_meta2$loi_type=="num_rast"]
                                                                        cat_rast<-loi_meta2$loi_var_nms[loi_meta2$loi_type=="cat_rast"]


                                                                        # Lumped Summaries --------------------------------------------------------

                                                                        if ("lumped" %in% weighting_scheme2 | any(loi_meta2$loi_type=="cat_rast")) {
                                                                          lumped_mean_out<-df %>%
                                                                            dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                                            #dplyr::mutate(lumped=1) %>%
                                                                            dplyr::summarise(
                                                                              dplyr::across(tidyselect::any_of(numb_rast),
                                                                                            ~sum(.,na.rm=T)/sum(!is.na(.),na.rm=T)
                                                                              ),
                                                                              dplyr::across(tidyselect::any_of(cat_rast),
                                                                                            ~sum(.,na.rm=T)/dplyr::n()#sum(!!rlang::sym("lumped"),na.rm=T)
                                                                              )
                                                                            ) %>%
                                                                            tibble::as_tibble()

                                                                          if (length(numb_rast)>0) {
                                                                            lumped_mean_out<-lumped_mean_out %>%
                                                                              dplyr::rename_with(.cols=tidyselect::any_of(numb_rast),~paste0(.,"_lumped_mean"))
                                                                          }

                                                                          if (length(cat_rast)>0) {
                                                                            lumped_mean_out<-lumped_mean_out %>%
                                                                              dplyr::rename_with(.cols=tidyselect::any_of(cat_rast),~paste0(.,"_lumped_prop")) %>%
                                                                              dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.)))
                                                                          }


                                                                        }

                                                                        if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="min") & length(numb_rast)>0){
                                                                          min_out<-df %>%
                                                                            dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                                            dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~min(.,na.rm=T))) %>%
                                                                            tibble::as_tibble()%>%
                                                                            dplyr::rename_with(~paste0(.,"_lumped_min"))#%>%
                                                                          #dplyr::collect()

                                                                        }
                                                                        if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="max" & length(numb_rast)>0)){
                                                                          max_out<-df %>%
                                                                            dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                                            dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~max(.,na.rm=T))) %>%
                                                                            tibble::as_tibble()%>%
                                                                            dplyr::rename_with(~paste0(.,"_lumped_max"))#%>%
                                                                          #dplyr::collect()

                                                                        }
                                                                        if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="count")){
                                                                          count_out<-df %>%
                                                                            dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                                            dplyr::summarise(dplyr::across(tidyselect::everything(),~sum(!is.na(.),na.rm=T))) %>%
                                                                            tibble::as_tibble()%>%
                                                                            dplyr::rename_with(~paste0(.,"_lumped_count"))#%>%
                                                                          #dplyr::collect()

                                                                        }

                                                                        if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="sum" & length(numb_rast)>0)){
                                                                          sum_out<-df %>%
                                                                            dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                                            dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~sum(.,na.rm=T))) %>%
                                                                            tibble::as_tibble()%>%
                                                                            dplyr::rename_with(~paste0(.,"_lumped_sum"))#%>%
                                                                          #dplyr::collect()
                                                                        }

                                                                        if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="median" & length(numb_rast)>0)){
                                                                          median_out<-df %>%
                                                                            dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                                            dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~stats::median(.,na.rm=T))) %>%
                                                                            tibble::as_tibble()%>%
                                                                            dplyr::rename_with(~paste0(.,"_lumped_median"))#%>%
                                                                          #dplyr::collect()
                                                                        }
                                                                        if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2 %in% c("sd","stdev") & length(numb_rast)>0)){
                                                                          lumped_sd_out<-df %>%
                                                                            dplyr::select(tidyselect::any_of(loi_cols2)) %>%
                                                                            dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~stats::sd(.,na.rm=T))) %>%
                                                                            tibble::as_tibble()%>%
                                                                            dplyr::rename_with(~paste0(.,"_lumped_sd"))#%>%
                                                                          #dplyr::collect()
                                                                        }

                                                                        # Weighted summaries -----------------------------------------------------------

                                                                        if (length(weighting_scheme2[!weighting_scheme2 %in% "lumped"])>0 & (any(loi_numeric_stats2 %in% c("mean")) | length(cat_rast)>0)) {
                                                                          weighted_mean_out<-df %>%
                                                                            dplyr::summarize(
                                                                              dplyr::across(tidyselect::ends_with(paste0("_iFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)),
                                                                              dplyr::across(tidyselect::ends_with(paste0("_HAiFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)),
                                                                              dplyr::across(tidyselect::ends_with(paste0("_iFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T)),
                                                                              dplyr::across(tidyselect::ends_with(paste0("_HAiFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T))
                                                                            ) %>%
                                                                            tibble::as_tibble()

                                                                          if (length(numb_rast)>0) {
                                                                            weighted_mean_out<-weighted_mean_out %>%
                                                                              dplyr::rename_with(.cols=tidyselect::starts_with(paste0(numb_rast,"_")),~paste0(.,"_mean"))
                                                                          }

                                                                          if (length(cat_rast)>0) {
                                                                            weighted_mean_out<-weighted_mean_out %>%
                                                                              dplyr::rename_with(.cols=tidyselect::starts_with(paste0(cat_rast,"_")),~paste0(.,"_prop")) %>%
                                                                              dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.)))
                                                                          }

                                                                        }


                                                                        if (length(weighting_scheme2[!weighting_scheme2 %in% "lumped"])>0 &
                                                                            any(loi_numeric_stats2 %in% c("sd","stdev")) &
                                                                            length(numb_rast)>0
                                                                        ) {

                                                                          weighted_sd_out<-df %>%
                                                                            dplyr::select(
                                                                              tidyselect::starts_with(numb_rast),
                                                                              tidyselect::any_of(weighting_scheme2)
                                                                            ) %>%
                                                                            dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                                                                                        ~(!!rlang::sym("iFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)))^2)),
                                                                                                        .names="{.col}_iFLS_term1"),
                                                                                          dplyr::across(tidyselect::any_of(numb_rast),
                                                                                                        ~ ((sum(!!rlang::sym("iFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLS")!=0,na.rm=T)) * sum(!!rlang::sym("iFLS"),na.rm=T),
                                                                                                        .names="{.col}_iFLS_term2"
                                                                                          )) %>%
                                                                            dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                                                                                        ~(!!rlang::sym("HAiFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)))^2)),
                                                                                                        .names="{.col}_HAiFLS_term1"),
                                                                                          dplyr::across(tidyselect::any_of(numb_rast),
                                                                                                        ~ ((sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLS"),na.rm=T),
                                                                                                        .names="{.col}_HAiFLS_term2"
                                                                                          )) %>%
                                                                            dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                                                                                        ~(!!rlang::sym("iFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T)))^2)),
                                                                                                        .names="{.col}_iFLO_term1"),
                                                                                          dplyr::across(tidyselect::any_of(numb_rast),
                                                                                                        ~ ((sum(!!rlang::sym("iFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLO")!=0,na.rm=T)) * sum(!!rlang::sym("iFLO"),na.rm=T),
                                                                                                        .names="{.col}_iFLO_term2"
                                                                                          )) %>%
                                                                            dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                                                                                        ~(!!rlang::sym("HAiFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T)))^2)),
                                                                                                        .names="{.col}_HAiFLO_term1"),
                                                                                          dplyr::across(tidyselect::any_of(numb_rast),
                                                                                                        ~ ((sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLO"),na.rm=T),
                                                                                                        .names="{.col}_HAiFLO_term2"
                                                                                          )) %>%
                                                                            dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),~sum(.,na.rm=T)),
                                                                                             dplyr::across(tidyselect::ends_with("_term2"),~.[1])
                                                                            ) %>%
                                                                            tibble::as_tibble() %>%
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
                                                                        final_list<-final_list[sapply(final_list,nrow)>0]

                                                                        final_out<-dplyr::bind_cols(tibble::tibble(pour_point_id=pour_point_id),final_list)

                                                                        final_out<-final_out %>%
                                                                          dplyr::select(
                                                                            pour_point_id,
                                                                            tidyselect::contains(loi_meta2$loi_var_nms)
                                                                          )

                                                                        #DBI::dbDisconnect(con)
                                                                        rm(df)
                                                                        gg<-gc()

                                                                        return(final_out)

                                                                      }

                                 ),silent=T)

                                 if (inherits(ot,"try-error")) {
                                   warning(paste0("Could not process the following loi: ",x,
                                                  " at link_id: ",xx))

                                   ot<-NULL
                                 }

                                 return(ot)
                               }) %>%
                               .[!sapply(.,is.null)] %>%
                               purrr::reduce(dplyr::left_join,by="pour_point_id")
                             }) %>%
                  .[!sapply(.,is.null)] %>%
                  dplyr::bind_rows()

                  ,silent=F)

                #file.remove(terra::sources(input_rasts))
              }

              if (inherits(ot,"try-error")) {
                ot<-tibble::tibble(link_id=input_poly$link_id)
                warning(paste0("Could not process the following link_id: ",paste0(input_poly$link_id,collapse = ", ")))
              }

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
