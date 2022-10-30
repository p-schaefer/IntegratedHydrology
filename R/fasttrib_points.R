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
#' @param subb_per_core numeric, number of subbasins to evaluate per core cycle
#' @param catch_per_core numeric, number of catchments to evaluate per core cycle
#' @param store_hw logical. Should hydroweight layer be stored and added to the zip file?
#' @param out_filename Output file name.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param unzip_arg character. Arguments passed to unzip. In Unix systems, may need to set to unzip_arg = "/usr/bin/unzip" to avoid warning about corrupted zip file.
#' @param verbose logical.
#'
#' @return A data.frame of weighted attributes for the requested areas
#' @export
#'
#' @importFrom carrier crate
#' @importFrom data.table fread data.table fwrite
#' @importFrom DBI dbConnect dbListTables dbDisconnect dbExecute
#' @importFrom dplyr collect tbl mutate across na_if left_join copy_to select filter bind_rows distinct rename group_by ungroup summarize bind_cols arrange desc n collapse compute summarise rename_with rowwise
#' @importFrom furrr future_pmap furrr_options
#' @importFrom future nbrOfWorkers availableCores future futureOf resolved
#' @importFrom hydroweight hydroweight
#' @importFrom progressr with_progress progressor
#' @importFrom purrr map map2 pmap map_dbl pmap_dfr reduce
#' @importFrom rlang sym
#' @importFrom RSQLite SQLite
#' @importFrom sf read_sf st_union write_sf
#' @importFrom stats setNames sd median
#' @importFrom stringr str_split_fixed
#' @importFrom terra terraOptions rasterize rast values ncell rowFromCell colFromCell crs vect unique sources writeRaster
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr nest unnest pivot_longer pivot_wider
#' @importFrom tidyselect any_of everything ends_with contains
#' @importFrom utils unzip zip
#' @importFrom whitebox wbt_options wbt_exe_path wbt_unnest_basins

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
    subb_per_core=1000,
    catch_per_core=8,
    temp_dir=NULL,
    unzip_arg="internal",
    verbose=F
){

  if (!is.logical(use_exising_hw)) stop("'use_exising_hw' must be logical")
  if (!is.logical(store_hw)) stop("'store_hw' must be logical")
  if (!is.logical(use_existing_attr)) stop("'use_existing_attr' must be logical")

  if (use_exising_hw & is.null(input$dw_dir)) stop("If 'use_exising_hw'=TRUE, input$dw_dir ust not be NULL. Has prep_weights() been run?")

  if (store_hw){
    use_exising_hw<-T

    input<-ihydro::prep_weights(
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


  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)
  if (n_cores==0) n_cores<-1

  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)

  weighting_scheme_s<-weighting_scheme[grepl("FLS|iEucS",weighting_scheme)]
  weighting_scheme_o<-weighting_scheme[!grepl("lumped|FLS|iEucS",weighting_scheme)]
  lumped_scheme<-"lumped" %in% weighting_scheme
  if (length(weighting_scheme_o)>0) message("Calculation for iEucO, iFLO, and HAiFLO are slow")

  if (is.null(target_o_type)) target_o_type<-"point"
  if (length(target_o_type)>1) target_o_type<-target_o_type[[1]]
  match.arg(target_o_type,several.ok = F)
  match.arg(weighting_scheme,several.ok = T)
  match.arg(loi_numeric_stats,several.ok = T)

  if (is.null(out_filename)) out_filename<-paste0("attrib_",floor(as.numeric(Sys.time())),".csv")
  if (length(out_filename)>1) out_filename<-out_filename[[1]]
  if (!grepl("\\.csv$",out_filename)) stop("'out_filename' must be a .csv file")

  loi_numeric_stats<-stats::setNames(loi_numeric_stats,loi_numeric_stats)

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

    #if (!file.exists(attr_db_loc)) stop(paste0(attr_db_loc," must exist if use_existing_attr=TRUE"))
  }


  con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)
  # # DBI::dbSendStatement(con_attr, "PRAGMA busy_timeout = 10000")
  # # DBI::dbSendStatement(con_attr,"PRAGMA journal_mode = WAL")
  # # DBI::dbSendStatement(con_attr,"PRAGMA synchronous = 0")
  # DBI::dbSendStatement(con_attr,"PRAGMA cache_size = 1000000")
  # # # DBI::dbSendStatement(con_attr,"PRAGMA locking_mode = EXCLUSIVE")
  # DBI::dbSendStatement(con_attr,"PRAGMA temp_store = MEMORY")
  # DBI::dbSendStatement(con_attr,"PRAGMA mmap_size = 30000000000")
  # DBI::dbSendStatement(con_attr,"PRAGMA page_size = 32768")

  attr_tbl_list<-DBI::dbListTables(con_attr)
  DBI::dbDisconnect(con_attr)

  loi_loc<-loi_file
  if (is.null(loi_loc)) loi_loc<-zip_loc

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  whitebox::wbt_options(exe_path=whitebox::wbt_exe_path(),
                        verbose=verbose,
                        wd=temp_dir)

  terra::terraOptions(verbose = verbose,
                      tempdir = temp_dir
  )

  fl<-utils::unzip(list=T,zip_loc)
  fl_loi<-utils::unzip(list=T,loi_loc)

  if (!is.null(dw_dir)){
    fl_dw<-utils::unzip(list=T,dw_dir)
  } else {
    fl_dw<-NULL
  }

  if (verbose) print("Reading in data")

  site_id_col<-paste0(data.table::fread(cmd=paste("unzip -p ",zip_loc,"site_id_col.csv")))

  db_loc<-input$db_loc
  con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
  stream_links<-dplyr::collect(dplyr::tbl(con,"stream_links")) %>%
    dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
    dplyr::mutate(dplyr::across(tidyselect::any_of(site_id_col),dplyr::na_if,""))
  DBI::dbDisconnect(con)

  all_points<-sf::read_sf(file.path("/vsizip",zip_loc,"stream_links.shp"))%>%
    dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>% #stream_links.shp
    dplyr::left_join(stream_links, by = c("link_id"))

  all_subb<-sf::read_sf(file.path("/vsizip",zip_loc,"Subbasins_poly.shp"))
  all_catch<-sf::read_sf(file.path("/vsizip",zip_loc,"Catchment_poly.shp"))
  #browser()

  if (!"link_id_cellstats" %in% attr_tbl_list){
    if (verbose) print("Generating cell number tables")

    decimalplaces <- function(x) {
      if (abs(x - round(x)) > .Machine$double.eps^0.5) {
        nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
      } else {
        return(0)
      }
    }

    con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)
    n_dec<-max(sapply(all_subb$link_id,decimalplaces))
    all_subb_rast<-terra::rasterize(all_subb,
                                    terra::rast(file.path("/vsizip",zip_loc,"dem_final.tif")),
                                    field="link_id")

    all_subb_out<-data.table::data.table(
      subb_link_id=terra::values(all_subb_rast),
      cell_number=1:terra::ncell(all_subb_rast)
    ) %>%
      stats::setNames(c("subb_link_id","cell_number")) %>%
      dplyr::mutate(row=terra::rowFromCell(all_subb_rast,cell_number)) %>%
      dplyr::mutate(col=terra::colFromCell(all_subb_rast,cell_number)) %>%
      dplyr::mutate(subb_link_id=round(subb_link_id,n_dec)) %>%
      dplyr::mutate(subb_link_id=as.character(subb_link_id)) %>%
      dplyr::copy_to(df=.,
                     con_attr,
                     "link_id_cellstats",
                     overwrite =T,
                     temporary =F,
                     indexes=c("subb_link_id","cell_number"),
                     analyze=T,
                     in_transaction=T)

    all_subb_out<-NULL
    DBI::dbDisconnect(con_attr)

  }

  # Get target link_id ------------------------------------------------------
  sample_points<-as.character(sample_points)
  link_id<-as.character(link_id)
  if (length(sample_points)==0 & length(link_id)==0) {
    message("`sample_points` and `link_id` are NULL, all `link_id`s will evaluated")
    target_IDs<-all_points %>%
      tibble::as_tibble() %>%
      dplyr::select(link_id,tidyselect::any_of(site_id_col))
  } else {
    if (site_id_col!="link_id" & length(sample_points)>0){
      target_IDs<-all_points %>%
        tibble::as_tibble() %>%
        dplyr::select(link_id,tidyselect::any_of(site_id_col)) %>%
        dplyr::filter(!!rlang::sym(site_id_col) %in% sample_points)
    } else {
      target_IDs<-NULL
    }

    if (length(link_id)>0){
      link_id_sel<-link_id
      target_IDs<-dplyr::bind_rows(
        target_IDs,
        all_points %>%
          tibble::as_tibble() %>%
          dplyr::select(link_id,tidyselect::any_of(site_id_col)) %>%
          dplyr::filter(link_id %in% link_id_sel)
      )
    }
  }

  if (target_o_type=="segment_whole") {
    target_IDs<-target_IDs %>%
      dplyr::select(link_id) %>%
      dplyr::mutate(link_id=as.character(floor(as.numeric(link_id))))
  }

  target_IDs<-dplyr::distinct(target_IDs)


  target_crs<-terra::crs(terra::vect(all_subb[1,]))

  # Get Upstream flowpaths --------------------------------------------------

  us_fp_fun<-function(link_id_in,db_loc=db_loc){
    con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
    out<-dplyr::tbl(con,"us_flowpaths") %>%
      dplyr::filter(pour_point_id %in% link_id_in) %>%
      dplyr::rename(link_id=origin_link_id) %>%
      dplyr::collect() %>%
      dplyr::group_by(pour_point_id) %>%
      tidyr::nest() %>%
      dplyr::ungroup()

    out2<-out$data
    names(out2)<-out$pour_point_id

    out2<-out2[link_id_in]

    DBI::dbDisconnect(con)
    return(out2)
  }

  # browser()
  us_flowpaths_out<-target_IDs %>%
    dplyr::select(link_id) %>%
    dplyr::mutate(link_id=as.character(link_id)) %>%
    dplyr::mutate(us_flowpaths=us_fp_fun(link_id,db_loc=db_loc))

  # Select correct target for O -------------------------------------
  if (target_o_type=="point"){
    target_O<-all_points
  } else {
    if (target_o_type=="segment_point"){
      target_O<-sf::read_sf(file.path("/vsizip",zip_loc,"stream_lines.shp"))
    } else {
      if (verbose) print("Merging stream segments")

      target_O<-sf::read_sf(file.path("/vsizip",zip_loc,"stream_lines.shp")) %>%
        dplyr::select(link_id) %>%
        dplyr::mutate(link_id=as.character(floor(as.numeric(link_id)))) %>%
        dplyr::group_by(link_id) %>%
        dplyr::summarize(geometry=sf::st_union(geometry)) %>%
        dplyr::ungroup()

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
    dplyr::filter(link_id %in% unlist(purrr::map(us_flowpaths_out$us_flowpaths,~.$link_id)))

  # all_subb_v<-terra::vect(all_subb)
  # all_catch_v<-terra::vect(all_catch)


  # Calculate s-weighted distances -------------------------------------
  if (!use_existing_attr | !"s_target_weights" %in% attr_tbl_list){

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

      hw_streams_nm<-utils::unzip(list=T,hw_streams)$Name
      names(hw_streams_nm)<-gsub(".tif","",hw_streams_nm)

      hw_streams_lo<-purrr::map(hw_streams_nm,function(x){
        file.path("/vsizip",hw_streams,x)
      })
    } else {
      #browser()
      trg_fl<-paste0("ALL_",weighting_scheme_s,"_inv_distances.tif")
      if (all(sapply(trg_fl,function(x) any(grepl(x,fl_dw$Name))))) {
        hw_streams_lo<-purrr::map(trg_fl,~file.path("/vsizip",dw_dir,.))
        names(hw_streams_lo)<-weighting_scheme_s
      } else {
        stop(paste0("Not all 'weighting_scheme' found in zip file"))
      }
    }

    if (verbose) print("Writing S-targeted weights to attributes database")

    # hw2<-map(hw_streams_lo,terra::rast)
    # names(hw2)<-sapply(hw2,names)

    con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)

    s_trg_weights<-dplyr::copy_to(df=tibble::as_tibble(matrix(ncol = length(names(hw_streams_lo))+2),nrow=1,.name_repair="minimal") %>%
                                    stats::setNames(c("subb_link_id","cell_number",names(hw_streams_lo))) %>%
                                    dplyr::mutate(dplyr::across(tidyselect::everything(),~1.1)) %>%
                                    .[F,],
                                  con_attr,
                                  "s_target_weights",
                                  overwrite =T,
                                  temporary =F,
                                  #indexes=c("subb_link_id","cell_number"),
                                  analyze=T,
                                  in_transaction=T)


    ot<-parallel_layer_processing(n_cores=n_cores,
                                  attr_db_loc=attr_db_loc,
                                  polygons=all_subb,
                                  n_per_cycle=subb_per_core,
                                  rasts=hw_streams_lo,
                                  cols=names(hw_streams_lo),
                                  temp_dir=temp_dir,
                                  tbl_nm="s_target_weights",
                                  sub_nm="s_target_weights",
                                  link_id_nm="subb_link_id"
    )


    s_trg_weights<-NULL

    t1<-DBI::dbExecute(con_attr,"CREATE INDEX inx_s_target_weights ON s_target_weights (subb_link_id, cell_number)")

    DBI::dbDisconnect(con_attr)

  }


  # Separate target_o into non-overlapping groups ---------------------------
  if (length(weighting_scheme_o)>0){
    if (!use_existing_attr | !"o_target_weights" %in% attr_tbl_list){
      if (!use_exising_hw){

        if (verbose) print("Generating Site Targeted Weights")
        if (verbose) print("Unnesting Basins")

        temp_dir_sub<-file.path(temp_dir,basename(tempfile()))
        dir.create(temp_dir_sub)

        utils::unzip(zip_loc,
                     c("dem_d8.tif"),
                     exdir=temp_dir_sub,
                     overwrite=T,
                     junkpaths=T,
                     unzip=unzip_arg)
        sf::write_sf(all_points %>% dplyr::select(link_id),
                     file.path(temp_dir_sub,"pour_points.shp"),
                     overwrite=T)

        #browser()

        future_unnest<-future::future({
          whitebox::wbt_unnest_basins(
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

          rast_all<-purrr::map(fl_un,function(x) try(terra::rast(x),silent=T))
          rast_all<-rast_all[!sapply(rast_all,function(x) inherits(x,"try-error"))]

          if (length(rast_all)>0){
            rast_out<-c(rast_out,purrr::map(rast_all,terra::unique))
            suppressWarnings(file.remove(unlist(purrr::map(rast_all,terra::sources))))
          }
        }

        fl_un<-list.files(temp_dir_sub,"unnest_",full.names = T)
        rast_all<-purrr::map(fl_un,function(x) try(terra::rast(x),silent=T))
        rast_all<-rast_all[!sapply(rast_all,function(x) inherits(x,"try-error"))]
        if (length(rast_all)>0){
          rast_out<-c(rast_out,purrr::map(rast_all,terra::unique))
          suppressWarnings(file.remove(unlist(purrr::map(rast_all,terra::sources))))
        }

        #browser()
        target_O_sub<-purrr::map2(rast_out,seq_along(rast_out),~target_O[unlist(.x),] %>% dplyr::select(link_id) %>% dplyr::mutate(unn_group=.y))
      } else {

        #browser()
        target_O_sub<-sf::read_sf(file.path("/vsizip",dw_dir,"unnest_group_target_O.shp")) %>%
          split(.,paste0("unnest_group_",.$unn_group))


        trg_fls<-sapply(weighting_scheme_o,function(x) paste0(names(target_O_sub),"_",x,"_inv_distances.tif"))
        if (any(!trg_fls %in% fl_dw$Name)){
          stop("Some unnested target_O groups are missing hydroweights")
        }

      }

      # extract O target weights
      if (verbose) print("Writing O-targeted weights to attributes database")

      progressr::with_progress(enable=T,{

        p <- progressr::progressor(steps = length(target_O_sub))


        con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)

        o_trg_weights<-tibble::tibble(catch_link_id="1.1",
                                      cell_number=1L,
        ) %>%
          dplyr::bind_cols(
            data.frame(matrix(ncol=length(weighting_scheme_o),nrow=1)) %>%
              stats::setNames(weighting_scheme_o)
          ) %>%
          dplyr::mutate(dplyr::across(tidyselect::any_of(weighting_scheme_o),~1.1)) %>%
          .[F,]%>%
          dplyr::copy_to(df=.,
                         con_attr,
                         "o_target_weights",
                         overwrite =T,
                         temporary =F,
                         #indexes=c("catch_link_id","cell_number"),
                         analyze=T,
                         in_transaction=T)

        o_trg_weights<-NULL

        DBI::dbDisconnect(con_attr)

        spltl<-1

        loi_dw_out<-purrr::pmap( # I don't think this can be parallel
          #loi_dw_out<-future_pmap_dfr(
          #  .options = furrr_options(globals = FALSE),
          list(
            target_O_subs=suppressWarnings(split(target_O_sub,1:spltl)),
            weighting_scheme_o=rep(list(weighting_scheme_o),spltl),
            all_catch=rep(list(all_catch),spltl),
            inv_function=rep(list(inv_function),spltl),
            use_exising_hw=rep(list(use_exising_hw),spltl),
            temp_dir=rep(list(temp_dir),spltl),
            n_cores=rep(list(n_cores),spltl),
            new_tbl=rep(list(o_trg_weights),spltl),
            zip_loc=rep(list(zip_loc),spltl),
            dw_dir=rep(list(dw_dir),spltl),
            p=rep(list(p),spltl),
            attr_db_loc=rep(list(attr_db_loc),spltl),
            catch_per_core=rep(list(catch_per_core),spltl),
            subb_per_core=rep(list(subb_per_core),spltl)
          ),
          carrier::crate(function(target_O_subs,
                                  weighting_scheme_o,
                                  all_catch,
                                  inv_function,
                                  use_exising_hw,
                                  temp_dir,
                                  n_cores,
                                  new_tbl,
                                  zip_loc,
                                  dw_dir,
                                  p,
                                  attr_db_loc,
                                  catch_per_core,
                                  subb_per_core
          ) {
            #browser()
            options(scipen = 999)
            `%>%` <- magrittr::`%>%`

            temp_dir_sub<-file.path(temp_dir,basename(tempfile()))
            dir.create(temp_dir_sub)

            target_S <- file.path("/vsizip",zip_loc,"dem_streams_d8.tif")
            dem <- file.path("/vsizip",zip_loc,"dem_final.tif")
            flow_accum <- file.path("/vsizip",zip_loc,"dem_accum_d8.tif")

            o_out<-purrr::map(target_O_subs, # I don't think this can be parallel
                              function(x){
                                #browser()
                                #print(x$unn_group[[1]])
                                sub_catch<-all_catch %>%
                                  dplyr::filter(link_id %in% x$link_id)

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
                                  hw<-NULL
                                  trg_fl<-paste0("unnest_group_",x$unn_group[[1]],"_",weighting_scheme_o,"_inv_distances.tif")
                                  hw_o_lo<-purrr::map(trg_fl,~file.path("/vsizip",dw_dir,.))#terra::rast(
                                  names(hw_o_lo)<-weighting_scheme_o
                                }

                                temp_dir_sub_sub<-file.path(temp_dir_sub,basename(tempfile()))
                                dir.create(temp_dir_sub_sub)

                                sub_catch<-sub_catch %>%
                                  dplyr::mutate(split=rep(1:ceiling(nrow(x)/catch_per_core) ,length.out=nrow(sub_catch)))

                                sub_catch_split<-split(sub_catch,sub_catch$split)

                                for (i in sub_catch_split){
                                  ihydro::parallel_layer_processing(n_cores=n_cores,
                                                                    attr_db_loc=attr_db_loc,
                                                                    polygons=i,
                                                                    n_per_cycle=subb_per_core,
                                                                    rasts=hw_o_lo,
                                                                    cols=names(hw_o_lo),
                                                                    temp_dir=temp_dir_sub_sub,
                                                                    tbl_nm="o_target_weights",
                                                                    sub_nm="o_target_weights",
                                                                    link_id_nm="catch_link_id",
                                                                    progress=F
                                  )
                                }



                                if (!is.null(hw)) file.remove(hw)

                                p()

                                return(NULL)
                              })

          }))

        con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc)

        t1<-DBI::dbExecute(con_attr,"CREATE INDEX idx_o_target_weights ON o_target_weights (catch_link_id, cell_number)")

        DBI::dbDisconnect(con_attr)


      })

    }
  } else {
    target_O_sub<-NULL
  }

  utils::unzip(loi_loc,file="loi_meta.rds",exdir = temp_dir,unzip=unzip_arg)

  loi_meta<-readRDS(file.path(temp_dir,"loi_meta.rds"))

  loi_rasts_exists<-purrr::map(loi_meta,~unlist(purrr::map(.,~.$output_filename)))
  loi_rasts_names<-purrr::map(loi_meta,~unlist(purrr::map(.,~.$lyr_variables)))
  loi_rasts_names<-purrr::map(loi_rasts_names,~stats::setNames(.,.))

  if (is.null(loi_cols)) loi_cols<-unlist(loi_rasts_names,use.names=F)

  # Upload loi rasters to attributes database -------------------------------
  if (!use_existing_attr | !"attrib_tbl" %in% attr_tbl_list){
    #browser()


    if (F){
      loi_rasts<-purrr::map(loi_rasts_exists,terra::rast)
      loi_rasts_comb<-terra::rast(loi_rasts)
      names(loi_rasts_comb)<-unlist(sapply(loi_rasts,names))

      # try to see if this is faster with a big dataset
      terra::writeRaster(loi_rasts_comb,file.path(temp_dir,"all_preds.tif"))
      loi_rasts_comb<-list(file.path(temp_dir,"all_preds.tif"))
    }

    # loi_rasts_exists<-c("num_rast.tif","cat_rast.tif")
    # if (!any(loi_rasts_exists %in% fl_loi$Name)) stop("No 'loi' present in input, please run 'process_loi()' first, or specify location of process_loi() ouput")
    # loi_rasts_exists<-fl_loi$Name[grepl("num_rast|cat_rast",fl_loi$Name)]
    # loi_rasts_exists<-map(loi_rasts_exists,~file.path("/vsizip",loi_loc,.))
    # names(loi_rasts_exists)<-gsub("\\.tif","",sapply(loi_rasts_exists,basename))
    #
    # loi_rasts<-map(loi_rasts_exists,rast)
    #
    # loi_rasts_comb<-rast(loi_rasts)
    #
    # names(loi_rasts_comb)<-unlist(sapply(loi_rasts,names))
    # if (is.null(loi_cols)) loi_cols<-names(loi_rasts_comb)
    #
    # if (any(!loi_cols %in% names(loi_rasts_comb))) stop(paste0("The following `loi_cols` are not present in the `loi`:",
    #                                                            paste0(loi_cols[!loi_cols %in% names(loi_rasts_comb)],collapse = ", ")
    # ))
    #
    # loi_rasts_comb<-terra::subset(loi_rasts_comb,loi_cols)
    # loi_rasts_names<-map(loi_rasts,names)
    # loi_rasts_names<-loi_rasts_names[loi_rasts_names %in% loi_cols]
    # loi_rasts_names<-map(loi_rasts_names,~map(.,~setNames(as.list(.),.)) %>% unlist(recursive=T))
    # loi_rasts_names<-map(loi_rasts_names,~map(.,~loi_numeric_stats))
    # loi_rasts_names$cat_rast<-as.list(setNames(rep(NA,length(names(loi_rasts$cat_rast))),names(loi_rasts$cat_rast)))


    if (verbose) print("Writing LOI to attributes database")

    con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)

    attrib_tbl<-dplyr::copy_to(df=tibble::as_tibble(matrix(ncol = length(unlist(loi_rasts_names))+2,nrow=1),.name_repair="minimal") %>%
                                 stats::setNames(c("subb_link_id","cell_number",unlist(loi_rasts_names))) %>%
                                 dplyr::mutate(dplyr::across(tidyselect::everything(),~1.1)) %>%
                                 .[F,],
                               con_attr,
                               "attrib_tbl",
                               overwrite =T,
                               temporary =F,
                               analyze=T,
                               in_transaction=T)

    ot<-parallel_layer_processing(n_cores=n_cores,
                                  attr_db_loc=attr_db_loc,
                                  polygons=all_subb,
                                  n_per_cycle=subb_per_core,
                                  rasts=loi_rasts_exists,
                                  cols=unlist(loi_rasts_names),
                                  temp_dir=temp_dir,
                                  tbl_nm="attrib_tbl",
                                  sub_nm="attrib_tbl",
                                  link_id_nm="subb_link_id"
    )

    attrib_tbl<-NULL
    t1<-DBI::dbExecute(con_attr,"CREATE INDEX inx_attrib_tbl ON attrib_tbl (subb_link_id, cell_number)")

    DBI::dbDisconnect(con_attr)
  }

  con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc)

  #browser()
  #t1<-DBI::dbExecute(con_attr,"PRAGMA analysis_limit=1000")
  if (F)  t1<-DBI::dbExecute(con_attr,"VACUUM")
  #t1<-DBI::dbExecute(con_attr,"OPTIMIZE")

  #DBI::dbSendStatement(con_attr,"PRAGMA journal_mode = OFF")
  DBI::dbDisconnect(con_attr)

  # Calculate Attributes ----------------------------------------------------
  # o_trg_weights
  # s_trg_weights
  # attrib_tbl

  lumped_out<-NULL
  s_targ_out<-NULL
  o_targ_out<-NULL

  #
  # Lumped Attributes -------------------------------------------------------

  #browser() # here make calculations work with more than 1 link_id at a time

  us_flowpaths_out_o<-us_flowpaths_out
  us_flowpaths_out_o<-us_flowpaths_out_o %>%
    dplyr::arrange(dplyr::desc(purrr::map_dbl(us_flowpaths,nrow))) %>%
    dplyr::mutate(splt1=rep(1:n_cores,length.out=dplyr::n())) %>%
    dplyr::group_by(splt1) %>%
    # mutate(splt2=rep(1:catch_per_core,length.out=n())) %>%
    # ungroup() %>%
    # group_by(splt2) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(data=purrr::map(data,~dplyr::mutate(.,splt2=rep(1:ceiling(nrow(.)/catch_per_core),length.out=nrow(.))) %>% split(.,.$splt2))) %>%
    dplyr::mutate(
      attr_db_loc=list(attr_db_loc),
      loi_rasts_names=list(loi_rasts_names),
      loi_numeric_stats=list(loi_numeric_stats),
      loi_cols=list(loi_cols)
    ) %>%
    dplyr::group_by(splt1) %>%
    tidyr::nest() %>%
    dplyr::ungroup()

  nrep<-sum(purrr::map_dbl(us_flowpaths_out_o$data,~length(.$data[[1]])))

  #browser()

  if (F & lumped_scheme) {
    if (verbose) print("Calculating Lumped Attributes")

    progressr::with_progress(enable=T,{
      p <- progressr::progressor(steps = nrep)
      lumped_out<-us_flowpaths_out_o %>%
        dplyr::mutate(p=list(p)) %>%
        dplyr::mutate(attr=furrr::future_pmap(#
          list(
            data=data,
            p=p
          ),
          .options = furrr::furrr_options(globals = FALSE),
          carrier::crate(
            function(data,
                     p
            ){
              #browser()
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`

              purrr::pmap_dfr(
                list(
                  x=data$data[[1]],
                  attr_db_loc=data$attr_db_loc,
                  loi_rasts_names=data$loi_rasts_names,
                  loi_numeric_stats=data$loi_numeric_stats,
                  loi_cols=data$loi_cols,
                  p=list(p)
                ),
                carrier::crate(
                  function(x,
                           attr_db_loc,
                           loi_rasts_names,
                           loi_numeric_stats,
                           loi_cols,
                           p
                  ){
                    #browser()

                    link_id_in<-x$link_id

                    #browser()
                    options(scipen = 999)
                    `%>%` <- magrittr::`%>%`

                    con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)

                    out<-dplyr::tbl(con_attr,"us_flowpaths") %>%
                      dplyr::rename(link_id=origin_link_id) %>%
                      dplyr::filter(pour_point_id %in% local(link_id_in)) %>%
                      dplyr::collapse() %>%
                      dplyr::left_join(dplyr::tbl(con_attr,"attrib_tbl") %>%
                                         dplyr::select(subb_link_id,cell_number,tidyselect::any_of(loi_cols)) %>%
                                         dplyr::rename(link_id=subb_link_id )%>%
                                         dplyr::collapse(),
                                       by=c("link_id")
                      ) %>%
                      dplyr::group_by(pour_point_id) %>%
                      dplyr::compute()

                    # out<-dplyr::left_join(
                    #   dplyr::tbl(con_attr,"us_flowpaths") %>%
                    #     dplyr::filter(pour_point_id %in% link_id_in) %>%
                    #     dplyr::rename(link_id=origin_link_id),
                    #   dplyr::tbl(con_attr,"attrib_tbl") %>%
                    #     dplyr::rename(link_id=subb_link_id),
                    #   by="link_id"
                    # ) %>%
                    #   dplyr::compute()

                    #attrs<-sapply(sapply(loi_rasts_names$num_rast,unique),unique)
                    attrs<-loi_numeric_stats

                    mean_out<-NULL
                    sd_out<-NULL
                    min_out<-NULL
                    max_out<-NULL
                    count_out<-NULL
                    median_out<-NULL
                    sum_out<-NULL

                    if (any("mean"==attrs)|length(loi_rasts_names$cat_rast) > 0){
                      mean_out<-out %>%
                        dplyr::select(-link_id,-cell_number) %>%
                        dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~sum(.,na.rm=T)/sum(!is.na(.))),
                                         dplyr::across(tidyselect::any_of(names(loi_rasts_names$cat_rast)),~sum(.,na.rm=T)/dplyr::n())) %>%
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
                        dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~stats::median(.,na.rm=T)))%>%
                        dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_median")) %>%
                        dplyr::collect()

                    }
                    if (any(attrs=="sum")){
                      sum_out<-out %>%
                        dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~sum(.,na.rm=T)))%>%
                        dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_sum")) %>%
                        dplyr::collect()
                    }

                    final_list<-list(
                      mean_out,
                      sd_out,
                      min_out,
                      max_out,
                      count_out,
                      median_out,
                      sum_out
                    )
                    final_list<-final_list[!sapply(final_list,is.null)]

                    final_out<-purrr::reduce(
                      final_list,
                      dplyr::left_join,
                      by="pour_point_id"
                    )

                    final_out <- x %>%
                      dplyr::select(-us_flowpaths,-splt2) %>%
                      dplyr::left_join(final_out,by=c("link_id"="pour_point_id"))


                    DBI::dbDisconnect(con_attr)

                    p()

                    gg<-gc()

                    return(final_out)

                  })
              )

            })
        )) %>%
        dplyr::select(attr) %>%
        tidyr::unnest(cols = c(attr))
    })

  }

  # s-targeted Attributes -------------------------------------------------------

  if (lumped_scheme | length(weighting_scheme_s)>0) {
    if (verbose) print("Calculating s-targeted Attributes")

    progressr::with_progress(enable=T,{

      p <- progressr::progressor(steps = nrep)
      s_targ_out<-us_flowpaths_out_o %>%
        dplyr::mutate(p=list(p),
                      weighting_scheme_s=list(weighting_scheme_s),
                      lumped_scheme=list(lumped_scheme)) %>%
        dplyr::mutate(attr=furrr::future_pmap(
          list(
            data=data,
            weighting_scheme_s=weighting_scheme_s,
            lumped_scheme=lumped_scheme,
            p=p
          ),
          .options = furrr::furrr_options(globals = FALSE),
          carrier::crate(
            function(data,
                     weighting_scheme_s,
                     lumped_scheme,
                     p
            ){
              #browser()
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`

              purrr::pmap_dfr(
                list(
                  x=data$data[[1]],
                  attr_db_loc=data$attr_db_loc,
                  loi_rasts_names=data$loi_rasts_names,
                  weighting_scheme_s=list(weighting_scheme_s),
                  lumped_scheme=list(lumped_scheme),
                  loi_numeric_stats=data$loi_numeric_stats,
                  loi_cols=data$loi_cols,
                  p=list(p)
                ),
                carrier::crate(
                  function(x,
                           attr_db_loc,
                           loi_rasts_names,
                           weighting_scheme_s,
                           lumped_scheme,
                           loi_numeric_stats,
                           loi_cols,
                           p
                  ){
                    #browser()
                    attrs<-loi_numeric_stats

                    options(scipen = 999)
                    `%>%` <- magrittr::`%>%`
                    #browser()

                    con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)

                    attr_nms<-names(c(loi_rasts_names$num_rast,loi_rasts_names$cat_rast))
                    names(attr_nms)<-attr_nms

                    names(weighting_scheme_s)<-weighting_scheme_s

                    link_id_in<-x$link_id

                    sub_id_id<-dplyr::tbl(con_attr,"us_flowpaths") %>%
                      dplyr::filter(pour_point_id %in% local(link_id_in)) %>%
                      dplyr::select(origin_link_id) %>%
                      dplyr::distinct() %>%
                      dplyr::pull(origin_link_id)

                    out<-dplyr::tbl(con_attr,"attrib_tbl") %>%
                      dplyr::select(subb_link_id,cell_number,tidyselect::any_of(loi_cols)) %>%
                      dplyr::rename(link_id=subb_link_id) %>%
                      dplyr::filter(link_id %in% local(sub_id_id)) %>%
                      dplyr::collapse() %>%
                      dplyr::group_by(link_id)


                    if (length(weighting_scheme_s)>0) {
                      out_extra<-out%>%
                        dplyr::left_join(
                          dplyr::tbl(con_attr,"s_target_weights") %>%
                            dplyr::select(subb_link_id,cell_number,tidyselect::any_of(weighting_scheme_s)) %>%
                            dplyr::rename(link_id=subb_link_id ) %>%
                            dplyr::collapse(),
                          by=c("link_id","cell_number")
                        )

                      if ("iFLS" %in% weighting_scheme_s){ # This is the only way I could get around an error by iterating over weighting_scheme_s
                        out_extra<-out_extra %>%
                          dplyr::left_join(
                            out_extra %>%
                              dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("iFLS")),.names="{.col}_iFLS" )) %>%
                              dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_","iFLS"))),

                            by = c("cell_number", "link_id")
                          )
                      }

                      if ("HAiFLS" %in% weighting_scheme_s){
                        out_extra<-out_extra %>%
                          dplyr::left_join(
                            out_extra %>%
                              dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("HAiFLS")),.names="{.col}_HAiFLS" )) %>%
                              dplyr::select(cell_number,link_id,tidyselect::ends_with(paste0("_","HAiFLS"))),
                            by = c("cell_number", "link_id")
                          )
                      }
                    }

                    lumped_out<-NULL
                    if (lumped_scheme) {
                      if (length(weighting_scheme_s)>0 & any(attrs %in% c("sd","stdev"))) {
                        lumped_out<-dplyr::tbl(con_attr,"us_flowpaths") %>%
                          dplyr::rename(link_id=origin_link_id) %>%
                          dplyr::filter(pour_point_id %in% local(link_id_in)) %>%
                          dplyr::group_by(pour_point_id) %>%
                          dplyr::left_join(out_extra,by="link_id") %>%
                          dplyr::compute()
                      } else {
                        lumped_out<-dplyr::tbl(con_attr,"us_flowpaths") %>%
                          dplyr::rename(link_id=origin_link_id) %>%
                          dplyr::filter(pour_point_id %in% local(link_id_in)) %>%
                          dplyr::group_by(pour_point_id) %>%
                          dplyr::left_join(out,by="link_id") %>%
                          dplyr::compute()
                      }

                    }

                    attrs<-loi_numeric_stats

                    weighted_mean_out<-NULL
                    lumped_mean_out<-NULL
                    weighted_sd_out<-NULL
                    lumped_sd_out<-NULL
                    min_out<-NULL
                    max_out<-NULL
                    count_out<-NULL
                    median_out<-NULL
                    sum_out<-NULL

                    if (length(weighting_scheme_s)>0 & (any("mean"==attrs)|length(loi_rasts_names$cat_rast) >0)) {
                      dw_sum<-out_extra %>%
                        dplyr::summarize(dplyr::across(tidyselect::any_of(paste0(weighting_scheme_s)),~sum(.,na.rm=T)),
                                         dplyr::across(tidyselect::ends_with(paste0("_",weighting_scheme_s)),~sum(.,na.rm=T)))

                      dw_sum_comp<-dplyr::compute(dw_sum)

                      weighted_mean_out<-dplyr::tbl(con_attr,"us_flowpaths") %>%
                        dplyr::rename(link_id=origin_link_id) %>%
                        dplyr::filter(pour_point_id %in% local(link_id_in)) %>%
                        dplyr::group_by(pour_point_id) %>%
                        dplyr::left_join(dw_sum_comp,by="link_id") %>%
                        dplyr::summarize(
                          dplyr::across(tidyselect::ends_with(paste0("_","iFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)),
                          dplyr::across(tidyselect::ends_with(paste0("_","HAiFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T))
                        ) %>%
                        dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$num_rast),"_")),~paste0(.x,"_mean")) %>%
                        dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$cat_rast),"_")),~paste0(.x,"_prop")) %>%
                        dplyr::ungroup() %>%
                        dplyr::collect()
                    }

                    if (lumped_scheme & any(attrs=="min")){
                      min_out<-lumped_out %>%
                        dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~min(.,na.rm=T))) %>%
                        dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_min"))%>%
                        dplyr::ungroup() %>%
                        dplyr::collect()

                    }
                    if (lumped_scheme & any(attrs=="max")){
                      max_out<-lumped_out %>%
                        dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~max(.,na.rm=T)))%>%
                        dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_max")) %>%
                        dplyr::ungroup() %>%
                        dplyr::collect()

                    }
                    if (lumped_scheme & any(attrs=="count")){
                      count_out<-lumped_out %>%
                        dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~sum(.[!is.na(.)],na.rm=T)))%>%
                        dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_count")) %>%
                        dplyr::ungroup() %>%
                        dplyr::collect()

                    }
                    if (lumped_scheme & any(attrs=="median")){
                      median_out<-lumped_out %>%
                        dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~median(.,na.rm=T)))%>%
                        dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_median")) %>%
                        dplyr::ungroup() %>%
                        dplyr::collect()

                    }
                    if (lumped_scheme & any(attrs=="sum")){
                      sum_out<-lumped_out %>%
                        dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~sum(.,na.rm=T)))%>%
                        dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_sum")) %>%
                        dplyr::ungroup() %>%
                        dplyr::collect()
                    }
                    if (lumped_scheme & any(attrs %in% c("sd","stdev"))){
                      lumped_sd_out<-lumped_out %>%
                        dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~sd(.,na.rm=T))) %>%
                        dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_sd"))%>%
                        dplyr::ungroup() %>%
                        dplyr::collect()

                    }
                    if (lumped_scheme & (any(attrs %in% c("mean"))|length(loi_rasts_names$cat_rast)>0)){
                      lumped_sd_out<-lumped_out %>%
                        dplyr::summarise(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),~sum(.,na.rm=T)/sum(!is.na(.))),
                                         dplyr::across(tidyselect::any_of(names(loi_rasts_names$cat_rast)),~sum(.,na.rm=T)/dplyr::n())) %>%
                        dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$num_rast)),~paste0(.x,"_lumped_mean")) %>%
                        dplyr::rename_with(.cols=tidyselect::any_of(names(loi_rasts_names$cat_rast)),~paste0(.x,"_lumped_prop")) %>%
                        dplyr::ungroup() %>%
                        dplyr::collect()

                    }

                    if (length(weighting_scheme_s)>0 & any(attrs %in% c("sd","stdev"))) {

                      temp_sd<-lumped_out

                      weighted_sd_out<-purrr::map(weighting_scheme_s,function(x){

                        if ("iFLS" %in% x){
                          out2<-temp_sd %>%
                            dplyr::select(pour_point_id,
                                          tidyselect::any_of(names(loi_rasts_names$num_rast)),
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
                          out2<-temp_sd %>%
                            dplyr::select(pour_point_id,
                                          tidyselect::any_of(names(loi_rasts_names$num_rast)),
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

                        out2 %>%
                          tidyr::pivot_longer(cols=c(tidyselect::everything(),-pour_point_id)) %>%
                          dplyr::mutate(attr=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_",2)[,1],
                                        term=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_",2)[,2]) %>%
                          dplyr::rowwise() %>%
                          dplyr::mutate(hw=gsub(paste0(attr,"_","|","_",term,""),"",name)) %>%
                          dplyr::ungroup() %>%
                          dplyr::mutate(name=paste0(attr,"_",hw,"_sd")) %>%
                          dplyr::group_by(name,pour_point_id) %>%
                          dplyr::summarize(sd=sqrt(value[term=="term1"]/value[term=="term2"])) %>%
                          dplyr::ungroup() %>%
                          tidyr::pivot_wider(names_from = name,values_from=sd)

                      })

                      weighted_sd_out<-purrr::reduce(weighted_sd_out,dplyr::left_join,by="pour_point_id")


                    }










                    #
                    #                     #attrs<-sapply(sapply(loi_rasts_names$num_rast,unique),unique)
                    #                     attrs<-loi_numeric_stats
                    #
                    #                     mean_out<-NULL
                    #                     sd_out<-NULL
                    #
                    #                     if (any("mean"==attrs)|length(loi_rasts_names$cat_rast) >0){
                    #
                    #                       mean_out<-purrr::map(weighting_scheme_s,function(x){
                    #
                    #                         if ("iFLS" %in% x){
                    #                           out2<-out_comp %>%
                    #                             dplyr::summarize(dplyr::across(tidyselect::ends_with(paste0("_","iFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T) ))  %>%
                    #                             dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$num_rast),"_")),~paste0(.x,"_mean")) %>%
                    #                             dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$cat_rast),"_")),~paste0(.x,"_prop")) %>%
                    #                             dplyr::collect()
                    #                         }
                    #
                    #                         if ("HAiFLS" %in% x){
                    #                           out2<-out_comp %>%
                    #                             dplyr::summarize(dplyr::across(tidyselect::ends_with(paste0("_","HAiFLS")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T) ))  %>%
                    #                             dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$num_rast),"_")),~paste0(.x,"_mean")) %>%
                    #                             dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$cat_rast),"_")),~paste0(.x,"_prop")) %>%
                    #                             dplyr::collect()
                    #                         }
                    #
                    #                         return(out2)
                    #
                    #                       })
                    #                       mean_out<-purrr::reduce(mean_out,dplyr::left_join,by="pour_point_id")
                    #                     }
                    #
                    #                     if (any(attrs %in% c("sd","stdev"))) {
                    #                       sd_out<-purrr::map(weighting_scheme_s,function(x){
                    #
                    #                         if ("iFLS" %in% x){
                    #                           out2<-out_comp %>%
                    #                             dplyr::select(pour_point_id,
                    #                                           tidyselect::any_of(names(loi_rasts_names$num_rast)),
                    #                                           tidyselect::any_of("iFLS")
                    #                             ) %>%
                    #                             dplyr::mutate(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                    #                                                         ~(!!rlang::sym("iFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("iFLS"),na.rm=T)))^2)),
                    #                                                         .names="{.col}_iFLS_term1"),
                    #                                           dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                    #                                                         ~ ((sum(!!rlang::sym("iFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("iFLS")!=0,na.rm=T)) * sum(!!rlang::sym("iFLS"),na.rm=T),
                    #                                                         .names="{.col}_iFLS_term2"
                    #                                           ))%>%
                    #                             dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),sum,na.rm=T),
                    #                                              dplyr::across(tidyselect::ends_with("_term2"),~.[1])
                    #                             ) %>%
                    #                             dplyr::collect()
                    #                         }
                    #
                    #                         if ("HAiFLS" %in% x){
                    #                           out2<-out_comp %>%
                    #                             dplyr::select(pour_point_id,
                    #                                           tidyselect::any_of(names(loi_rasts_names$num_rast)),
                    #                                           tidyselect::any_of("HAiFLS")
                    #                             ) %>%
                    #                             dplyr::mutate(dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                    #                                                         ~(!!rlang::sym("HAiFLS") * ((.-(sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLS"),na.rm=T)))^2)),
                    #                                                         .names="{.col}_HAiFLS_term1"),
                    #                                           dplyr::across(tidyselect::any_of(names(loi_rasts_names$num_rast)),
                    #                                                         ~ ((sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)-1)/sum(!!rlang::sym("HAiFLS")!=0,na.rm=T)) * sum(!!rlang::sym("HAiFLS"),na.rm=T),
                    #                                                         .names="{.col}_HAiFLS_term2"
                    #                                           ))%>%
                    #                             dplyr::summarize(dplyr::across(tidyselect::ends_with("_term1"),sum,na.rm=T),
                    #                                              dplyr::across(tidyselect::ends_with("_term2"),~.[1])
                    #                             ) %>%
                    #                             dplyr::collect()
                    #                         }
                    #
                    #                         out2 %>%
                    #                           tidyr::pivot_longer(cols=c(tidyselect::everything(),-pour_point_id)) %>%
                    #                           dplyr::mutate(attr=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_",2)[,1],
                    #                                         term=stringr::str_split_fixed(name,"_iFLS_|_HAiFLS_",2)[,2]) %>%
                    #                           dplyr::rowwise() %>%
                    #                           dplyr::mutate(hw=gsub(paste0(attr,"_","|","_",term,""),"",name)) %>%
                    #                           dplyr::ungroup() %>%
                    #                           dplyr::mutate(name=paste0(attr,"_",hw,"_sd")) %>%
                    #                           dplyr::group_by(name,pour_point_id) %>%
                    #                           dplyr::summarize(sd=sqrt(value[term=="term1"]/value[term=="term2"])) %>%
                    #                           dplyr::ungroup() %>%
                    #                           tidyr::pivot_wider(names_from = name,values_from=sd)
                    #
                    #                       })
                    #
                    #                       sd_out<-purrr::reduce(sd_out,dplyr::left_join,by="pour_point_id")
                    #
                    #
                    #                     }

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

                    final_out<-purrr::reduce(
                      final_list,
                      dplyr::left_join,
                      by="pour_point_id"
                    )

                    final_out <- x %>%
                      dplyr::select(-us_flowpaths,-splt2) %>%
                      dplyr::left_join(final_out,by=c("link_id"="pour_point_id"))

                    p()

                    gg<-gc()

                    DBI::dbDisconnect(con_attr)
                    return(final_out)

                  })
              )

            })
        )) %>%
        dplyr::select(attr) %>%
        tidyr::unnest(cols = c(attr))
    })

  }

  # o-targeted Attributes -------------------------------------------------------

  if (length(weighting_scheme_o)>0) {
    if (verbose) print("Calculating o-targeted Attributes")

    progressr::with_progress(enable=T,{
      p <- progressr::progressor(steps = nrep)
      o_targ_out<-us_flowpaths_out_o %>%
        dplyr::mutate(p=list(p),
                      weighting_scheme_o=list(weighting_scheme_o)) %>%
        dplyr::mutate(attr=furrr::future_pmap(
          list(
            data=data,
            weighting_scheme_o=weighting_scheme_o,
            p=p
          ),
          .options = furrr::furrr_options(globals = FALSE),
          carrier::crate(
            function(data,
                     weighting_scheme_o,
                     p
            ){
              #browser()
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`

              purrr::pmap_dfr(
                list(
                  x=data$data[[1]],
                  attr_db_loc=data$attr_db_loc,
                  loi_rasts_names=data$loi_rasts_names,
                  weighting_scheme_o=list(weighting_scheme_o),
                  loi_numeric_stats=data$loi_numeric_stats,
                  loi_cols=data$loi_cols,
                  p=list(p)
                ),
                carrier::crate(
                  function(x,
                           attr_db_loc,
                           loi_rasts_names,
                           weighting_scheme_o,
                           loi_numeric_stats,
                           loi_cols,
                           p
                  ){
                    #browser()
                    link_id_in<-x$link_id

                    options(scipen = 999)
                    `%>%` <- magrittr::`%>%`
                    #browser()

                    con_attr<-DBI::dbConnect(RSQLite::SQLite(), attr_db_loc,cache_size=1000000)

                    attr_nms<-names(c(loi_rasts_names$num_rast,loi_rasts_names$cat_rast))
                    names(attr_nms)<-attr_nms

                    names(weighting_scheme_o)<-weighting_scheme_o

                    out<-dplyr::tbl(con_attr,"us_flowpaths") %>%
                      dplyr::rename(link_id=origin_link_id) %>%
                      dplyr::filter(pour_point_id %in% local(link_id_in)) %>%
                      dplyr::collapse() %>%
                      dplyr::left_join(dplyr::tbl(con_attr,"attrib_tbl") %>%
                                         dplyr::select(subb_link_id,cell_number,tidyselect::any_of(loi_cols)) %>%
                                         dplyr::rename(link_id=subb_link_id ) %>%
                                         dplyr::collapse(),
                                       by=c("link_id")
                      ) %>%
                      dplyr::group_by(pour_point_id) %>%
                      dplyr::left_join(
                        dplyr::tbl(con_attr,"o_target_weights") %>%
                          dplyr::filter(catch_link_id %in% local(link_id_in)) %>%
                          dplyr::select(catch_link_id,cell_number,tidyselect::any_of(weighting_scheme_o)) %>%
                          dplyr::rename(link_id=catch_link_id ) %>%
                          dplyr::collapse(),
                        by=c("link_id","cell_number")
                      ) %>%
                      dplyr::group_by(pour_point_id) #%>%
                    #dplyr::compute()

                    if ("iFLO" %in% weighting_scheme_o){ # This is the only way I could get around an error by iterating over weighting_scheme_s
                      out<-out %>%
                        dplyr::left_join(
                          out %>%
                            dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("iFLO")),.names="{.col}_iFLO" )) %>%
                            dplyr::select(pour_point_id,cell_number,link_id,tidyselect::ends_with(paste0("_","iFLO"))),

                          by = c("cell_number", "link_id","pour_point_id")
                        )#%>%
                      #dplyr::compute()
                    }

                    if ("HAiFLO" %in% weighting_scheme_o){
                      out<-out %>%
                        dplyr::left_join(
                          out %>%
                            dplyr::mutate(dplyr::across(tidyselect::any_of(attr_nms), ~.*(!!rlang::sym("HAiFLO")),.names="{.col}_HAiFLO" )) %>%
                            dplyr::select(pour_point_id,cell_number,link_id,tidyselect::ends_with(paste0("_","HAiFLO"))),
                          by = c("cell_number", "link_id","pour_point_id")
                        )#%>%
                      #dplyr::compute()
                    }


                    out<-dplyr::compute(out)
                    #attrs<-sapply(sapply(loi_rasts_names$num_rast,unique),unique)
                    attrs<-loi_numeric_stats

                    mean_out<-NULL
                    sd_out<-NULL

                    if (any("mean"==attrs)|length(loi_rasts_names$cat_rast) >0){

                      mean_out<-purrr::map(weighting_scheme_o,function(x){

                        if ("iFLO" %in% x){
                          out2<-out %>%
                            dplyr::summarize(dplyr::across(tidyselect::ends_with(paste0("_","iFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("iFLO"),na.rm=T) ))  %>%
                            dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$num_rast),"_")),~paste0(.x,"_mean")) %>%
                            dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$cat_rast),"_")),~paste0(.x,"_prop")) %>%
                            dplyr::collect()
                        }

                        if ("HAiFLO" %in% x){
                          out2<-out %>%
                            dplyr::summarize(dplyr::across(tidyselect::ends_with(paste0("_","HAiFLO")),~sum(.,na.rm=T)/sum(!!rlang::sym("HAiFLO"),na.rm=T) ))  %>%
                            dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$num_rast),"_")),~paste0(.x,"_mean")) %>%
                            dplyr::rename_with(.cols=tidyselect::contains(paste0(names(loi_rasts_names$cat_rast),"_")),~paste0(.x,"_prop")) %>%
                            dplyr::collect()
                        }

                        return(out2)

                      })
                      mean_out<-purrr::reduce(mean_out,dplyr::left_join,by="pour_point_id")
                    }

                    if (any(attrs %in% c("sd","stdev"))) {
                      sd_out<-purrr::map(weighting_scheme_o,function(x){

                        if ("iFLO" %in% x){
                          out2<-out %>%
                            dplyr::select(pour_point_id,
                                          tidyselect::any_of(names(loi_rasts_names$num_rast)),
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
                          out2<-out %>%
                            dplyr::select(pour_point_id,
                                          tidyselect::any_of(names(loi_rasts_names$num_rast)),
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

                        out2 %>%
                          tidyr::pivot_longer(cols=c(tidyselect::everything(),-pour_point_id)) %>%
                          dplyr::mutate(attr=stringr::str_split_fixed(name,"_iFLO_|_HAiFLO_",2)[,1],
                                        term=stringr::str_split_fixed(name,"_iFLO_|_HAiFLO_",2)[,2]) %>%
                          dplyr::rowwise() %>%
                          dplyr::mutate(hw=gsub(paste0(attr,"_","|","_",term,""),"",name)) %>%
                          dplyr::ungroup() %>%
                          dplyr::mutate(name=paste0(attr,"_",hw,"_sd")) %>%
                          dplyr::group_by(name,pour_point_id) %>%
                          dplyr::summarize(sd=sqrt(value[term=="term1"]/value[term=="term2"])) %>%
                          dplyr::ungroup() %>%
                          tidyr::pivot_wider(names_from = name,values_from=sd)

                      })

                      sd_out<-purrr::reduce(sd_out,dplyr::left_join,by="pour_point_id")


                    }

                    final_list<-list(
                      mean_out,
                      sd_out
                    )
                    final_list<-final_list[!sapply(final_list,is.null)]

                    final_out<-purrr::reduce(
                      final_list,
                      dplyr::left_join,
                      by="pour_point_id"
                    )

                    final_out <- x %>%
                      dplyr::select(-us_flowpaths,-splt2) %>%
                      dplyr::left_join(final_out,by=c("link_id"="pour_point_id"))

                    p()

                    gg<-gc()

                    DBI::dbDisconnect(con_attr)
                    return(final_out)
                  })
              )

            })
        )) %>%
        dplyr::select(attr) %>%
        tidyr::unnest(cols = c(attr))
    })

  }

  #browser()

  final_out<-target_IDs %>%
    dplyr::mutate(link_id=as.character(link_id))

  if (!is.null(lumped_out)) final_out<-dplyr::left_join(final_out,lumped_out ,by="link_id")
  if (!is.null(s_targ_out)) final_out<-dplyr::left_join(final_out,s_targ_out ,by="link_id")
  if (!is.null(o_targ_out)) final_out<-dplyr::left_join(final_out,o_targ_out ,by="link_id")

  final_out<-final_out %>%
    dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),0,.))) %>%
    dplyr::select(-tidyselect::any_of("pour_point_id"))

  data.table::fwrite(x=final_out,
                     file=file.path(temp_dir,out_filename),
                     buffMB = 128L,
                     nThread = 1,
                     showProgress = F)

  utils::zip(zip_loc,
             file.path(temp_dir,out_filename),
             flags = '-r9Xjq'
  )


  suppressWarnings(file.remove(list.files(temp_dir,full.names = T,recursive = T)))


  return(final_out)

}

#' @export
#' @importFrom carrier crate
#' @importFrom data.table fwrite rbindlist data.table fread
#' @importFrom DBI dbConnect dbDisconnect dbAppendTable
#' @importFrom dplyr select bind_rows distinct filter tbl collect left_join bind_cols mutate summarize pull rename_with
#' @importFrom exactextractr exact_extract
#' @importFrom furrr future_pmap furrr_options
#' @importFrom future plan tweak multisession future futureOf resolved
#' @importFrom progressr with_progress progress progressor
#' @importFrom purrr map pmap map2
#' @importFrom RSQLite SQLite
#' @importFrom sf st_as_sf write_sf read_sf
#' @importFrom stats setNames
#' @importFrom terra rast subset extract cellFromRowColCombine readValues readStop
#' @importFrom tibble as_tibble
#' @importFrom tidyselect any_of
#' @importFrom carrier crate
#' @importFrom data.table fwrite rbindlist data.table fread
#' @importFrom DBI dbConnect dbDisconnect dbAppendTable
#' @importFrom dplyr select bind_rows distinct filter tbl collect left_join bind_cols mutate summarize pull rename_with
#' @importFrom exactextractr exact_extract
#' @importFrom furrr future_pmap furrr_options
#' @importFrom future plan tweak multisession future futureOf resolved
#' @importFrom progressr with_progress progress progressor
#' @importFrom purrr map pmap map2
#' @importFrom RSQLite SQLite
#' @importFrom sf st_as_sf write_sf read_sf
#' @importFrom stats setNames
#' @importFrom terra rast subset extract cellFromRowColCombine readValues readStop
#' @importFrom tibble as_tibble
#' @importFrom tidyselect any_of
#' @keywords internal
#' @noRd

parallel_layer_processing <- function(n_cores,
                                      attr_db_loc,
                                      polygons,
                                      n_per_cycle,
                                      rasts,
                                      cols,
                                      temp_dir,
                                      tbl_nm,
                                      sub_nm,
                                      link_id_nm=c("subb_link_id","catch_link_id"),
                                      progress=T
) {
  #browser() # Here setup readStart(), readValues(), and readStop() #https://github.com/rspatial/terra/issues/251
  options(scipen = 999)
  options(future.rng.onMisuse = "ignore")

  temp_dir_O<-temp_dir
  temp_dir<-file.path(temp_dir,basename(tempfile()))
  dir.create(temp_dir)


  link_id_nm<-match.arg(link_id_nm,c("subb_link_id","catch_link_id"),several.ok = F)

  loi_cols<-cols
  loi_rasts_exists<-rasts
  all_subb_v<-polygons

  n_cores_2<-n_cores

  if (n_cores>1) {
    n_cores_2<-n_cores_2-1
    oplan <- future::plan(list(future::tweak(future::multisession, workers = 2), future::tweak(future::multisession, workers = n_cores_2)))
    on.exit(future::plan(oplan), add = TRUE)
  }


  all_subb<-sf::st_as_sf(all_subb_v) %>%
    dplyr::select(link_id)
  fp<-file.path(temp_dir,paste0(tbl_nm,"_",basename(tempfile()),".shp"))
  ot<-sf::write_sf(all_subb,fp,overwrite=T)

  all_subb$core<-rep(1:(n_cores_2),length.out=nrow(all_subb))
  splt<-split(all_subb,all_subb$core)
  splt<-lapply(splt,function(x){
    x$split<-rep(1:ceiling(nrow(x)/n_per_cycle),length.out=nrow(x))
    x<-tibble::as_tibble(x) %>%
      dplyr::select(link_id,core,split)
    split(x,x$split)
  })

  total_outs<-sum(unlist(purrr::map(splt,length)))

  out<-future::future(
    #purrr::pmap(
    furrr::future_pmap(
      list(x=splt,
           core_numb=names(splt),
           loi_rasts_exists=list(loi_rasts_exists),
           loi_cols=list(loi_cols),
           temp_dir=list(temp_dir),
           link_id_nm=list(link_id_nm),
           sub_nm=list(sub_nm),
           fp=list(fp),
           attr_db_loc=list(attr_db_loc)
           # cell_fp=list(cell_fp),
           # cell_tbl=list(cell_tbl)
      ),
      .options = furrr::furrr_options(globals = FALSE),
      carrier::crate(
        function(x,
                 core_numb,
                 loi_rasts_exists,
                 loi_cols,
                 temp_dir,
                 link_id_nm,
                 sub_nm,
                 fp,
                 attr_db_loc
        ){
          #browser()

          options(scipen = 999)
          `%>%` <- magrittr::`%>%`

          loi_rasts<-purrr::map(loi_rasts_exists,terra::rast)
          loi_rasts_comb<-terra::rast(loi_rasts)
          names(loi_rasts_comb)<-unlist(sapply(loi_rasts,names))
          loi_rasts_comb<-terra::subset(loi_rasts_comb,loi_cols)

          # terra::readStart(loi_rasts_comb)
          #
          splt<-x

          xx<-dplyr::bind_rows(x) %>%
            dplyr::select(link_id) %>%
            dplyr::distinct()

          if (F){
            poly<-sf::read_sf(fp) %>%
              dplyr::filter(link_id %in% xx$link_id)
          } else {
            poly<-NULL
          }

          out<-purrr::pmap(list(xx=splt,
                                core_numb=names(splt),
                                loi_rasts_comb=list(loi_rasts_comb),
                                temp_dir=list(temp_dir),
                                sub_nm=list(sub_nm),
                                link_id_nm=list(link_id_nm),
                                fp=list(fp),
                                poly=list(poly),
                                attr_db_loc=list(attr_db_loc)
          ),
          carrier::crate(
            function(xx,
                     core_numb,
                     loi_rasts_comb,
                     temp_dir,
                     sub_nm,
                     link_id_nm,
                     fp,
                     poly,
                     attr_db_loc
                     #cell_tbl_sub
            ){
              #browser()
              options(scipen = 999)
              `%>%` <- magrittr::`%>%`

              if (T){ #this is using exact() with cell numbers
                con_attr2<-DBI::dbConnect(RSQLite::SQLite(),attr_db_loc)

                cell_tbl_sub<-try(stop(""),silent=T)

                while(inherits(cell_tbl_sub,'try-error')){
                  Sys.sleep(1)

                  cell_tbl_sub<-try({
                    if (link_id_nm=="subb_link_id"){
                      dplyr::tbl(con_attr2,"link_id_cellstats") %>%
                        dplyr::filter(subb_link_id %in% local(as.character(xx$link_id))) %>%
                        dplyr::select(subb_link_id,cell_number) %>%
                        dplyr::collect()
                    } else {
                      dplyr::tbl(con_attr2,"us_flowpaths") %>%
                        dplyr::filter(pour_point_id %in% local(as.character(xx$link_id))) %>%
                        dplyr::left_join(
                          dplyr::tbl(con_attr2,"link_id_cellstats"),
                          by=c("origin_link_id"="subb_link_id")
                        ) %>%
                        dplyr::select(catch_link_id=pour_point_id,cell_number) %>%
                        dplyr::distinct() %>%
                        dplyr::collect()
                    }
                  },silent=T)

                }

                DBI::dbDisconnect(con_attr2)

                out<-cell_tbl_sub %>%
                  dplyr::bind_cols(
                    terra::extract(
                      loi_rasts_comb,
                      cell_tbl_sub$cell_number
                    )
                  ) %>%
                  data.table::fwrite(file=file.path(temp_dir,paste0(sub_nm,"_s_",xx$core[[1]],"_",xx$split[[1]],".csv")),
                                     buffMB = 128L,
                                     nThread = 1,
                                     showProgress = F
                  )
              }

              if (F){ #this is using exact_extract() with a polygon

                out<-exactextractr::exact_extract(
                  loi_rasts_comb,
                  poly,
                  weights=NULL,
                  include_cell=T,
                  fun=NULL,
                  #include_cols="link_id",
                  progress=F
                ) %>%
                  purrr::map2(split(poly$link_id,poly$link_id),~dplyr::mutate(.x,link_id=.y)) %>%
                  data.table::rbindlist() %>%
                  dplyr::select(-coverage_fraction) %>%
                  stats::setNames(c(link_id_nm,names(loi_rasts_comb),"cell_number")) %>%
                  data.table::fwrite(file=file.path(temp_dir,paste0(sub_nm,"_s_",xx$core[[1]],"_",xx$split[[1]],".csv")),
                                     buffMB = 128L,
                                     nThread = 1,
                                     showProgress = F
                  )
              }

              if (F) { #this is using readValues() per subbasin
                con_attr2<-DBI::dbConnect(RSQLite::SQLite(),attr_db_loc)

                cell_tbl_sub<-try(stop(""),silent=T)

                while(inherits(cell_tbl_sub,'try-error')){
                  Sys.sleep(1)

                  cell_tbl_sub<-try({
                    if (link_id_nm=="subb_link_id"){
                      dplyr::tbl(con_attr2,"link_id_cellstats") %>%
                        dplyr::filter(subb_link_id %in% local(as.character(xx$link_id))) %>%
                        dplyr::select(link_id=subb_link_id,cell_number,row,col) %>%
                        dplyr::collect()
                    } else {
                      dplyr::tbl(con_attr2,"us_flowpaths") %>%
                        dplyr::filter(pour_point_id %in% local(as.character(xx$link_id))) %>%
                        dplyr::left_join(
                          dplyr::tbl(con_attr2,"link_id_cellstats"),
                          by=c("origin_link_id"="subb_link_id")
                        ) %>%
                        dplyr::select(link_id=pour_point_id,cell_number,row,col) %>%
                        dplyr::distinct() %>%
                        dplyr::collect()
                    }
                  },silent=T)

                }


                DBI::dbDisconnect(con_attr2)

                out<-purrr::map(xx$link_id,
                                function(xxx){
                                  #browser()

                                  target_cell_range<-cell_tbl_sub %>%
                                    dplyr::filter(link_id %in% local(as.character(xxx))) %>%
                                    dplyr::summarize(
                                      row_start=min(row),
                                      row_end=max(row),
                                      n_row=max(row)-min(row),
                                      col_start=min(col),
                                      col_end=max(col),
                                      n_col=max(col)-min(col)
                                    )

                                  target_cell_numbs<-cell_tbl_sub %>%
                                    dplyr::filter(link_id %in% local(as.character(xxx))) %>%
                                    dplyr::pull(cell_number)


                                  out_cell_nums<-terra::cellFromRowColCombine(loi_rasts_comb,
                                                                              row=target_cell_range$row_start:target_cell_range$row_end,
                                                                              col=target_cell_range$col_start:target_cell_range$col_end)

                                  out<-data.table::data.table(
                                    terra::readValues(loi_rasts_comb,
                                                      row=target_cell_range$row_start,
                                                      nrows=target_cell_range$n_row+1,
                                                      col=target_cell_range$col_start,
                                                      ncols=target_cell_range$n_col+1,
                                                      dataframe=T
                                    ))

                                  # if (nrow(out) != length(out_cell_nums)){
                                  #   browser()
                                  # }

                                  out<-out%>%
                                    dplyr::mutate(cell_number=out_cell_nums) %>%
                                    dplyr::filter(cell_number %in% target_cell_numbs) %>%
                                    dplyr::mutate(abcdefg = xxx) %>%
                                    dplyr::rename_with(.cols=tidyselect::any_of("abcdefg"),~paste0(link_id_nm))

                                  return(out)
                                }) %>%
                  dplyr::bind_rows()

                data.table::fwrite(out,
                                   file=file.path(temp_dir,paste0(sub_nm,"_s_",xx$core[[1]],"_",xx$split[[1]],".csv")),
                                   buffMB = 128L,
                                   nThread = 1,
                                   showProgress = F
                )
                #data.table::fwrite(out,file=file.path(temp_dir,paste0(sub_nm,"_s_",xx$core[[1]],"_",xx$split[[1]],".csv")))

              }

              if (file.exists(file.path(temp_dir,paste0(sub_nm,xx$core[[1]],"_",xx$split[[1]],".csv"))))
                file.remove(file.path(temp_dir,paste0(sub_nm,xx$core[[1]],"_",xx$split[[1]],".csv")))

              file.rename(
                file.path(temp_dir,paste0(sub_nm,"_s_",xx$core[[1]],"_",xx$split[[1]],".csv")),
                file.path(temp_dir,paste0(sub_nm,xx$core[[1]],"_",xx$split[[1]],".csv"))
              )

              try(rm(out),silent=T)
              gg<-gc()
              return(NA)

            }
          ))
          #})


          terra::readStop(loi_rasts_comb)


          return(NA)
        })

    ))

  #browser()
  future_out <- future::futureOf(out)

  total_procs<-0


  progressr::with_progress(enable=progress,{
    p <- progressr::progressor(steps = total_outs)


    while(!future::resolved(future_out)) {
      Sys.sleep(1)

      fl_attr<-list.files(temp_dir,sub_nm,full.names = T)
      fl_attr<-fl_attr[grepl(".csv",fl_attr)]
      fl_attr<-fl_attr[!grepl(paste0(sub_nm,"_s_"),fl_attr)]

      if (length(fl_attr)>0) {
        con_attr_sub<-DBI::dbConnect(RSQLite::SQLite(),attr_db_loc,cache_size=1000000)

        template<-dplyr::tbl(con_attr_sub,tbl_nm) %>%
          dplyr::collect(n=1) %>%
          sapply(class)

        df<-purrr::map(fl_attr,~try(data.table::fread(.,
                                                      select=template,
                                                      nThread=1,
                                                      showProgress=F,),silent = T))
        #df<-purrr::map(fl_attr,~try(utils::read.csv(.),silent = T))
        fl_attr<-fl_attr[!sapply(df,function(x) inherits(x,"try-error"))]
        df<-df[!sapply(df,function(x) inherits(x,"try-error"))]

        df<-df %>%
          dplyr::bind_rows()

        if (nrow(df)>0){

          ot<-try(stop(""),silent=T)

          while(inherits(ot,"try-error")){
            Sys.sleep(1)
            ot<-try({
              DBI::dbAppendTable(conn=con_attr_sub,
                                 name=tbl_nm,
                                 value=df)
            },silent=T)
          }

        }

        DBI::dbDisconnect(con_attr_sub)

        fr<-suppressMessages(file.remove(fl_attr))
        for (i in seq_along(fr)){
          p()
        }

        try(rm(df),silent=T)
        gg<-gc()
        total_procs<-total_procs+length(fl_attr)
      }
    }

  })

  Sys.sleep(3)
  #browser()

  if (length(out$result$conditions)>0){
    err<-out$result$conditions[[1]]$condition
    if (inherits(err,"error")){
      stop(err)
    }
  }

  fl_attr<-list.files(temp_dir,sub_nm,full.names = T)
  fl_attr<-fl_attr[grepl(".csv",fl_attr)]
  fl_attr<-fl_attr[!grepl(paste0(sub_nm,"_s_"),fl_attr)]

  if (length(fl_attr)>0) {

    if (length(fl_attr)>0) {
      con_attr_sub<-DBI::dbConnect(RSQLite::SQLite(),attr_db_loc,cache_size=1000000)

      template<-tbl(con_attr_sub,tbl_nm) %>%
        collect(n=1) %>%
        sapply(class)

      df<-purrr::map(fl_attr,~try(data.table::fread(.,
                                                    select=template,
                                                    nThread=1,
                                                    showProgress=F,),silent = T))
      #df<-purrr::map(fl_attr,~try(utils::read.csv(.),silent = T))
      fl_attr<-fl_attr[!sapply(df,function(x) inherits(x,"try-error"))]
      df<-df[!sapply(df,function(x) inherits(x,"try-error"))]

      df<-df %>%
        dplyr::bind_rows()

      if (nrow(df)>0){

        ot<-try(stop(""),silent=T)

        while(inherits(ot,"try-error")){
          Sys.sleep(1)
          ot<-try({
            DBI::dbAppendTable(conn=con_attr_sub,
                               name=tbl_nm,
                               value=df)
          },silent=T)
        }

      }

      DBI::dbDisconnect(con_attr_sub)

      fr<-suppressMessages(file.remove(fl_attr))
      try(rm(df),silent=T)
      gg<-gc()
      total_procs<-total_procs+length(fl_attr)

    }

  }

  fl_attr<-list.files(temp_dir,sub_nm,full.names = T)
  fl_attr<-fl_attr[grepl(".csv",fl_attr)]

  if (any(grepl(paste0(sub_nm,"_s_"),fl_attr))) stop("Some intermediate files could not be read into database")

  if (total_procs<total_outs) {
    stop("An error occured, not all attributes added to database")
  }

  fr<-file.remove(list.files(temp_dir,paste0(gsub(".shp","",basename(fp))),full.names = T))

  unlink(temp_dir,recursive = T,force=T)

  #t1<-DBI::dbExecute(con_attr_sub,"PRAGMA journal_mode = OFF")
  #DBI::dbDisconnect(con_attr_sub)


  return(NULL)

}
