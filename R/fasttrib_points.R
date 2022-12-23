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
  if (!dir.exists(temp_dir)) dir.create(temp_dir,recursive = T)
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
      dir.create(iDW_file,recursive = T)
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
  dir.create(temp_dir_sub,recursive = T)

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

  input_poly<-ihydro::get_catchment(input,link_id=unique(subb_IDs$pour_point_id))

  if (verbose) message("Calculating Weighted Attributes")


  # Calculate in parallel ---------------------------------------------------

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
               temp_dir_sub=list(temp_dir_sub),
               n_cores=list(n_cores)
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
                     temp_dir_sub,
                     n_cores) {
              #browser()
              ihydro::.summ_fn(
                x=x,
                input=input,
                iDW_file=iDW_file,
                loi_file=loi_file,
                loi_meta=loi_meta,
                weighting_scheme=weighting_scheme,
                loi_numeric_stats=loi_numeric_stats,
                loi_cols=loi_cols,
                p=p,
                temp_dir_sub=temp_dir_sub,
                n_cores=n_cores,
                backend=c("data.table")
              )

            }
          )))

  })

  out<-ot %>%
    dplyr::select(out) %>%
    tidyr::unnest(cols=out)

  #browser()


  # Any that couldn't be calculated in parallel -----------------------------
  incomplete_proc<-out$pour_point_id[out$status=="Incomplete"]

  if (length(incomplete_proc)>0){
    if (n_cores>1) {
      message(
        paste0(
          "The following link_ids could not be processed in parallel due to insufficient memory: ",
          paste(out$pour_point_id[out$status=="Incomplete"],collapse = ","),
          ". Attempting to recalculate sequentially."
        )
      )

      progressr::with_progress(enable=T,{
        p <- progressr::progressor(steps = length(incomplete_proc))

        ot2<-subb_IDs %>%
          dplyr::filter(pour_point_id %in% incomplete_proc) %>%
          dplyr::group_by(unn_group) %>%
          tidyr::nest() %>%
          dplyr::ungroup() %>%
          dplyr::mutate(core=rep(1,length.out=dplyr::n())) %>%
          dplyr::group_by(core) %>%
          tidyr::nest() %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            out=purrr::pmap(
              list(x=data,
                   input=list(as.ihydro(input$outfile)),
                   iDW_file=list(iDW_file),
                   loi_file=list(loi_file),
                   loi_meta=list(loi_meta),
                   weighting_scheme=list(weighting_scheme),
                   loi_numeric_stats=list(loi_numeric_stats),
                   loi_cols=list(loi_cols),
                   p=list(p),
                   temp_dir_sub=list(temp_dir_sub),
                   n_cores=list(1)
              ),
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
                         temp_dir_sub,
                         n_cores) {

                  ihydro::.summ_fn(
                    x=x,
                    input=input,
                    iDW_file=iDW_file,
                    loi_file=loi_file,
                    loi_meta=loi_meta,
                    weighting_scheme=weighting_scheme,
                    loi_numeric_stats=loi_numeric_stats,
                    loi_cols=loi_cols,
                    p=p,
                    temp_dir_sub=temp_dir_sub,
                    n_cores=1,
                    backend=c("data.table")
                  )
                }
              )))
      })

      out2<-ot2 %>%
        dplyr::select(out) %>%
        tidyr::unnest(cols=out)

      out<-bind_rows(out,out2)

      incomplete_proc<-out$pour_point_id[out$status=="Incomplete"]
    }
  }

  if (length(incomplete_proc)>0) {
    warning(
      paste0(
        "The following link_ids could not be processed due to insufficient memory: ",
        paste(out$pour_point_id[out$status=="Incomplete"],collapse = ",")
      )
    )

    # TODO: If all else fails, extract the data using subcatchments, write to
    #        sqlite database and do calculations there...
  }


  return(out)

}

#' Process gpkg rasters for attribute calculation
#'
#' @keywords internal
#' @return a data.frame of attributes
#' @export
#'
.summ_fn<-function(x,
                   input,
                   iDW_file,
                   loi_file,
                   loi_meta,
                   weighting_scheme,
                   loi_numeric_stats,
                   loi_cols,
                   p,
                   temp_dir_sub,
                   n_cores,
                   backend=c("data.table","SQLite")
){
  backend<-match.arg(backend,several.ok=T)

  sys.mem<-(memuse::Sys.meminfo()$totalram/n_cores)*0.9
  sys.disk<-NULL

  if ("SQLite" %in% backend){
    sys.disk <- system(paste0("df -BM ",temp_dir_sub),intern=T)
    sys.disk <- strsplit(sys.disk[length(sys.disk)], "[ ]+")[[1]]
    sys.disk <- as.numeric(gsub("\\D", "", sys.disk[4]))/1024
  }

  options(dplyr.summarise.inform = FALSE)
  options(scipen = 999)
  `%>%` <- magrittr::`%>%`

  iDWs_rasts<-NULL
  iDWo_rasts<-NULL

  loi_rasts<-terra::rast(loi_file,loi_cols)

  if (length(weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")])>0){
    iDWs_rasts<-terra::rast(iDW_file,weighting_scheme[weighting_scheme %in% c("iFLS","HAiFLS")])
  }

  purrr::map2_dfr(x$data,x$unn_group,function(y,yy){
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

    loi_types<-table(loi_meta$loi_type[loi_meta$loi_var_nms %in% loi_cols])
    num_rast_analysis_cols<-loi_types["num_rast"]*length(weighting_scheme[weighting_scheme!="lumped"])*3
    cat_rast_analysis_cols<-loi_types["cat_rast"]*length(weighting_scheme[weighting_scheme!="lumped"])*2
    analysis_cols<-sum(c(num_rast_analysis_cols,cat_rast_analysis_cols),na.rm=T)

    #Estimates the number of rows that could fit into memory without analysis
    max.obj.fulldata<-memuse::howmany(sys.mem,
                                      ncol=length(loi_cols)+length(weighting_scheme[weighting_scheme!="lumped"]))

    # Summary Function --------------------------------------------------------
    ot<-purrr::map_dfr(split(input_poly,seq_along(input_poly$link_id)),
                       function(sub_poly) {
                         #browser()
                         sub_id<-sub_poly$link_id
                         sub_poly_rast<-terra::rasterize(sub_poly,input_rasts)
                         sub_poly_rast<-terra::cells(sub_poly_rast)


                         #Estimates the number of rows that could fit into memory and do analysis
                         max.obj.fullanalysis_row<-memuse::howmany(sys.mem,
                                                                   ncol=analysis_cols)
                         max.obj.fullanalysis_col<-memuse::howmany(sys.mem,
                                                                   nrow=length(sub_poly_rast))

                         # Estimate table size if reading entire dataframe into memory
                         obj.size.fulldata<-memuse::howbig(nrow=length(sub_poly_rast),
                                                           ncol=length(loi_cols)+length(weighting_scheme[weighting_scheme!="lumped"]),
                                                           unit="GiB")

                         # Estimate table size if reading entire dataframe into memory
                         obj.size.fullanalysis<-memuse::howbig(nrow=length(sub_poly_rast),
                                                               ncol=analysis_cols,
                                                               unit="GiB")

                         # Estimate table size if single loi and IDW into memory
                         obj.size.single<-memuse::howbig(nrow=length(sub_poly_rast),
                                                         ncol=5,
                                                         unit="GiB")

                         if (obj.size.single>sys.mem) { # object won't fit into memory
                           # TODO: process rasters in chunks and write to sqlite, then do calculations there
                           out<-tibble::tibble(pour_point_id=sub_id,status="Incomplete")
                           # ot2<-terra::extract(
                           #   input_rasts,
                           #   sub_poly_rast
                           # )

                         } else {
                           if (obj.size.fullanalysis<sys.mem) { # object will fit entirely into memory
                             # ot<-exactextractr::exact_extract(
                             #   input_rasts,
                             #   sub_poly,
                             #   progress=F
                             # )[[1]]
                             ot<-terra::extract(
                               input_rasts,
                               sub_poly_rast
                             )

                             out<-ihydro::.attr_fn(ot,
                                                   point_id=sub_id,
                                                   weighting_scheme2=weighting_scheme,
                                                   loi_meta2=loi_meta,
                                                   loi_cols2=loi_cols,
                                                   loi_numeric_stats2=loi_numeric_stats)

                           } else {
                             if (max.obj.fullanalysis_col[2]>(length(weighting_scheme[weighting_scheme!="lumped"])*3)) { # all IDW plus at least some loi will fit into memory

                               weighting_scheme2<-weighting_scheme
                               if (all(weighting_scheme2=="lumped")) {
                                 weighting_scheme2<-"lumped"
                               } else {
                                 weighting_scheme2<-weighting_scheme[weighting_scheme!="lumped"]
                               }

                               if (all(weighting_scheme2=="lumped")){
                                 ot_idw<-data.frame(lumped=1)
                               } else {
                                 # ot_idw<-exactextractr::exact_extract(
                                 #   terra::subset(input_rasts,c(weighting_scheme2)[c(weighting_scheme2) %in% names(input_rasts)]),
                                 #   sub_poly,
                                 #   progress=F
                                 # )[[1]] %>%
                                 #   dplyr::filter(dplyr::if_any(tidyselect::any_of("coverage_fraction"),~.x>0.5)) %>%
                                 #   dplyr::select(-tidyselect::any_of("coverage_fraction"))

                                 ot_idw<-terra::extract(
                                   terra::subset(input_rasts,c(weighting_scheme2)[c(weighting_scheme2) %in% names(input_rasts)]),
                                   sub_poly_rast
                                 )


                                 if (ncol(ot_idw)==1) colnames(ot_idw)<-weighting_scheme2
                               }

                               remaining_cols<-max.obj.fullanalysis_col[2]-(length(weighting_scheme[weighting_scheme!="lumped"])*3)
                               remaining_cols_ratio<-remaining_cols/analysis_cols

                               if (remaining_cols_ratio>1){
                                 remaining_cols<-length(loi_cols)
                               } else {
                                 remaining_cols<-floor(length(loi_cols)*remaining_cols_ratio*0.9)
                               }

                               loi_cols_split<-split(loi_cols,
                                                     rep(1:length(loi_cols),
                                                         length.out=length(loi_cols),
                                                         each =remaining_cols)
                               )

                               out<-purrr::map(loi_cols_split,function(loi_sub){

                                 # ot<-exactextractr::exact_extract(
                                 #   terra::subset(input_rasts,c(loi_sub)[c(loi_sub) %in% names(input_rasts)]),
                                 #   sub_poly,
                                 #   progress=F
                                 # )[[1]] %>%
                                 #   dplyr::filter(dplyr::if_any(tidyselect::any_of("coverage_fraction"),~.x>0.5)) %>%
                                 #   dplyr::select(-tidyselect::any_of("coverage_fraction"))
                                 ot<-terra::extract(
                                   terra::subset(input_rasts,c(loi_sub)[c(loi_sub) %in% names(input_rasts)]),
                                   sub_poly_rast
                                 )

                                 if (ncol(ot)==1) colnames(ot)<-loi_sub

                                 out<-ihydro::.attr_fn(dplyr::bind_cols(ot,ot_idw),
                                                       point_id=sub_id,
                                                       weighting_scheme2=weighting_scheme,
                                                       loi_meta2=loi_meta,
                                                       loi_cols2=loi_sub,
                                                       loi_numeric_stats2=loi_numeric_stats)

                                 return(out)
                               }) %>%
                                 purrr::reduce(dplyr::left_join,by=c("pour_point_id","status"))

                             } else {
                               # one IDW plus at one loi will fit into memory

                               weighting_scheme2<-weighting_scheme
                               if (all(weighting_scheme2=="lumped")) {
                                 weighting_scheme2<-"lumped"
                               } else {
                                 weighting_scheme2<-weighting_scheme[weighting_scheme!="lumped"]
                               }

                               out<- purrr::map(weighting_scheme2,function(sub_weighting_scheme){

                                 if (sub_weighting_scheme=="lumped"){
                                   ot_idw<-data.frame(lumped=1)
                                 } else {
                                   # ot_idw<-exactextractr::exact_extract(
                                   #   terra::subset(input_rasts,c(sub_weighting_scheme)[c(sub_weighting_scheme) %in% names(input_rasts)]),
                                   #   sub_poly,
                                   #   progress=F
                                   # )[[1]] %>%
                                   #   dplyr::filter(dplyr::if_any(tidyselect::any_of("coverage_fraction"),~.x>0.5)) %>%
                                   #   dplyr::select(-tidyselect::any_of("coverage_fraction"))
                                   ot_idw<-terra::extract(
                                     terra::subset(input_rasts,c(sub_weighting_scheme)[c(sub_weighting_scheme) %in% names(input_rasts)]),
                                     sub_poly_rast
                                   )

                                   if (ncol(ot_idw)==1) colnames(ot_idw)<-sub_weighting_scheme
                                 }


                                 purrr::map(loi_cols,function(loi_sub){
                                   # ot<-exactextractr::exact_extract(
                                   #   terra::subset(input_rasts,c(loi_sub)[c(loi_sub) %in% names(input_rasts)]),
                                   #   sub_poly,
                                   #   progress=F
                                   # )[[1]] %>%
                                   #   dplyr::filter(dplyr::if_any(tidyselect::any_of("coverage_fraction"),~.x>0.5)) %>%
                                   #   dplyr::select(-tidyselect::any_of("coverage_fraction"))
                                   ot<-terra::extract(
                                     terra::subset(input_rasts,c(loi_sub)[c(loi_sub) %in% names(input_rasts)]),
                                     sub_poly_rast
                                   )

                                   if (ncol(ot)==1) colnames(ot)<-loi_sub

                                   with_lumped<-"lumped" %in% weighting_scheme

                                   if (with_lumped & sub_weighting_scheme==weighting_scheme2[1]){
                                     sub_weighting_scheme<-unique(c("lumped",sub_weighting_scheme))
                                   }

                                   out<-ihydro::.attr_fn(dplyr::bind_cols(ot,ot_idw),
                                                         point_id=sub_id,
                                                         weighting_scheme2=sub_weighting_scheme,
                                                         loi_meta2=loi_meta,
                                                         loi_cols2=loi_sub,
                                                         loi_numeric_stats2=loi_numeric_stats)
                                   return(out)
                                 })%>%
                                   purrr::reduce(dplyr::left_join,by=c("pour_point_id","status"))
                               }) %>%
                                 purrr::reduce(dplyr::left_join,by=c("pour_point_id","status"))

                             }
                           }
                         }

                         return(out)

                       })


    p()
    return(ot)
  })
}

#' Calculate loi attributes from input rasters
#'
#' @keywords internal
#' @return a data.frame of attributes
#' @export
#'
.attr_fn<-function(df=NULL,
                   db_path=NULL,
                   tbl_name=NULL,
                   point_id,
                   weighting_scheme2,
                   loi_meta2,
                   loi_cols2,
                   loi_numeric_stats2,
                   backend=c("data.table","SQLite")){
  options(dplyr.summarise.inform = FALSE)
  options(scipen = 999)
  `%>%` <- magrittr::`%>%`

  loi_meta2<-dplyr::filter(loi_meta2,loi_var_nms %in% loi_cols2)

  backend<-match.arg(backend)

  if (backend=="data.table") {
    stopifnot(!is.null(df))

    df_class<-sapply(df,class)

    df<-df %>%
      dtplyr::lazy_dt() %>%
      dplyr::filter(dplyr::if_any(tidyselect::any_of("coverage_fraction"),~.x>0.5)) %>%
      dplyr::select(-tidyselect::any_of("coverage_fraction"))
  } else {
    if (backend=="SQLite") {
      stopifnot(!is.null(db_path))
      stopifnot(!is.null(tbl_name))

      con<-DBI::dbConnect(RSQLite::SQLite(),db_path)

      df<-dplyr::tbl(con,tbl_name)

      df_class<-sapply(dplyr::collect(df,n=1),class)
    }
  }




  df<-df %>%
    dplyr::mutate(dplyr::across(tidyselect::any_of(names(df_class[df_class=="numeric"])),~ifelse(is.nan(.),NA_real_,.))) #tidyselect::where(is.numeric)

  if ("iFLS" %in% weighting_scheme2) {
    df<-df %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2),
                                  ~.*(!!rlang::sym("iFLS")),.names=paste0("{.col}_","iFLS")))
  }
  if ("iFLO" %in% weighting_scheme2) {
    df<-df %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2),
                                  ~.*(!!rlang::sym("iFLO")),.names=paste0("{.col}_","iFLO")))
  }
  if ("HAiFLS" %in% weighting_scheme2) {
    df<-df %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2),
                                  ~.*(!!rlang::sym("HAiFLS")),.names=paste0("{.col}_","HAiFLS")))
  }
  if ("HAiFLO" %in% weighting_scheme2) {
    df<-df %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(loi_cols2),
                                  ~.*(!!rlang::sym("HAiFLO")),.names=paste0("{.col}_","HAiFLO")))
  }

  df<-df %>%
    dplyr::compute()


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

  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="mean")) {
    lumped_mean_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(
        dplyr::across(tidyselect::any_of(numb_rast),
                      ~sum(.,na.rm=T)/sum(!is.na(.),na.rm=T)
        ),
        dplyr::across(tidyselect::any_of(cat_rast),
                      ~sum(.,na.rm=T)/dplyr::n()
        )
      ) %>%
      dplyr::collect()

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
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_min"))

  }
  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="max" & length(numb_rast)>0)){
    max_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~max(.,na.rm=T))) %>%
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_max"))

  }
  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="count")){
    count_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::everything(),~sum(!is.na(.),na.rm=T))) %>%
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_count"))

  }

  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="sum" & length(numb_rast)>0)){
    sum_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~sum(.,na.rm=T))) %>%
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_sum"))
  }

  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2=="median" & length(numb_rast)>0)){
    median_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~stats::median(.,na.rm=T))) %>%
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_median"))
  }
  if ("lumped" %in% weighting_scheme2 & any(loi_numeric_stats2 %in% c("sd","stdev") & length(numb_rast)>0)){
    lumped_sd_out<-df %>%
      dplyr::select(tidyselect::any_of(loi_cols2)) %>%
      dplyr::summarise(dplyr::across(tidyselect::any_of(numb_rast),~stats::sd(.,na.rm=T))) %>%
      dplyr::collect()%>%
      dplyr::rename_with(~paste0(.,"_lumped_sd"))
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
      dplyr::collect()

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


  if (
    length(weighting_scheme2[!weighting_scheme2 %in% "lumped"])>0 &
    any(loi_numeric_stats2 %in% c("sd","stdev")) &
    length(numb_rast)>0
  ) {

    weighted_sd_out<-df %>%
      dplyr::select(
        tidyselect::starts_with(numb_rast),
        tidyselect::any_of(weighting_scheme2)
      )

    if ("iFLS" %in% weighting_scheme2) {
      weighted_sd_out<-weighted_sd_out %>%
        dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                    ~(!!rlang::sym("iFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)))^2)),
                                    .names="{.col}_iFLS_term1"),
                      dplyr::across(tidyselect::any_of(numb_rast),
                                    ~ ((sum(!!rlang::sym("iFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLS")!=0,na.rm=T)) * sum(!!rlang::sym("iFLS"),na.rm=T),
                                    .names="{.col}_iFLS_term2"
                      ))
    }
    if ("HAiFLS" %in% weighting_scheme2) {
      weighted_sd_out<-weighted_sd_out %>%
        dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                    ~(!!rlang::sym("HAiFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)))^2)),
                                    .names="{.col}_HAiFLS_term1"),
                      dplyr::across(tidyselect::any_of(numb_rast),
                                    ~ ((sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLS"),na.rm=T),
                                    .names="{.col}_HAiFLS_term2"
                      ))
    }
    if ("iFLO" %in% weighting_scheme2) {
      weighted_sd_out<-weighted_sd_out %>%
        dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                    ~(!!rlang::sym("iFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T)))^2)),
                                    .names="{.col}_iFLO_term1"),
                      dplyr::across(tidyselect::any_of(numb_rast),
                                    ~ ((sum(!!rlang::sym("iFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLO")!=0,na.rm=T)) * sum(!!rlang::sym("iFLO"),na.rm=T),
                                    .names="{.col}_iFLO_term2"
                      ))
    }
    if ("HAiFLO" %in% weighting_scheme2) {
      weighted_sd_out<-weighted_sd_out %>%
        dplyr::mutate(dplyr::across(tidyselect::any_of(numb_rast),
                                    ~(!!rlang::sym("HAiFLO") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T)))^2)),
                                    .names="{.col}_HAiFLO_term1"),
                      dplyr::across(tidyselect::any_of(numb_rast),
                                    ~ ((sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLO")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLO"),na.rm=T),
                                    .names="{.col}_HAiFLO_term2"
                      ))
    }

    weighted_sd_out<- weighted_sd_out%>%
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


  final_list<-list(
    lumped_mean_out,
    lumped_sd_out,
    min_out,
    max_out,
    count_out,
    median_out,
    sum_out,
    weighted_mean_out,
    weighted_sd_out
  )
  final_list<-final_list[!sapply(final_list,is.null)]
  final_list<-final_list[sapply(final_list,nrow)>0]

  final_out<-dplyr::bind_cols(tibble::tibble(pour_point_id=point_id,status="Complete"),final_list)

  final_out<-final_out %>%
    dplyr::select(
      tidyselect::any_of("pour_point_id"),
      tidyselect::any_of("status"),
      tidyselect::contains(loi_meta2$loi_var_nms)
    )

  if (backend=="SQLite") {
    DBI::dbDisconnect(con)
  }

  return(final_out)
}



