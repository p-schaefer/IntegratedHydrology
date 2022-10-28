
#' Attributes stream segments/sampling points with layers of interest (loi)
#'
#' @param input output from `process_hydrology()` (if `process_loi()` was not run on `process_hydrology()`, `loi_file` must be specified)
#' @param loi_file filepath of `process_loi()` output (optional, will overwrite data in `process_hydrology()` output if present).
#' @param spec table containing which sampling points (and/or stream segments) to attribute. Must contain a column with the same `site_id_col` used in `process_hydrology()`, and `loi` column containing a named list of loi `variable_names`, and associated `loi_numeric_stats` for each. See example.
#' @param weighting_scheme character. One or more weighting schemes: c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS")
#' @param loi_numeric_stats character. One or more of c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max", "sum", "cell_count"). Those without distwtd_ are simple "lumped" statistics.
#' @param inv_function function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.
#' @param remove_region character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/lu.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster}. Regions to remove when summarizing the attributes (e.g., remove lake from catchment)
#' @param dw_dir character. File path for stored `hydroweight::hydroweight()` outputs, if separately calculated. Note file names must match format of \code{paste0(`site_id_col`,"_inv_distances.zip")}.
#' @param all_reaches logical. If \code{TRUE}, attributes are calculated for all reaches (sampling points are ignored). Warning, can be very slow.
#' @param OS_combine logical. Should target_O and target_S be merged as targets for iEucS, iFLS, and/or HAiFLS? Use \code{TRUE} or \code{FALSE}. This allows cells surrounding \code{target_O} to flow directly into \code{target_O} rather than be forced through \code{target_S}.
#' @param target_streamseg logical. If \code{TRUE}, `target_O` is considered the entire stream segment, else `target_O` is just the most downstream sampling point
#' @param catch_buffer numeric. Tolerance values used for \code{sf::st_buffer} in meters in `get_catchment()`.
#' @param return_products logical. If \code{TRUE}, a list containing all geospatial analysis products. If \code{FALSE}, folder path to resulting .zip file.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#'

#' @importFrom carrier crate
#' @importFrom data.table fread
#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom dplyr collect tbl mutate across na_if left_join select filter if_any pull distinct group_by ungroup rename arrange desc case_when sym summarise
#' @importFrom furrr future_map2 furrr_options
#' @importFrom future nbrOfWorkers availableCores
#' @importFrom hydroweight process_input hydroweight hydroweight_attributes
#' @importFrom progressr with_progress progressor
#' @importFrom purrr map pmap map_chr reduce
#' @importFrom RSQLite SQLite
#' @importFrom sf read_sf st_union st_length st_area st_as_sf st_geometry st_buffer
#' @importFrom stats setNames
#' @importFrom terra writeRaster rast crs crop vect wrap
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr nest unnest
#' @importFrom tidyselect any_of ends_with everything
#' @importFrom units set_units
#' @importFrom utils unzip head zip
#' @export


attrib_points<-function(
    input,
    loi_file=NULL,
    spec=NULL,
    all_reaches=F,
    weighting_scheme = c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS"),
    loi_numeric_stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max", "sum", "cell_count"),
    OS_combine=F,
    dw_dir=NULL,
    target_streamseg=F,
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    catch_buffer=0,
    return_products=F,
    remove_region=NULL,
    temp_dir=NULL,
    verbose=F
){
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")
  options(dplyr.summarise.inform = FALSE)

  if (!is.null(spec) && !inherits(spec,"data.frame")) stop("'spec' must be a data frame")
  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (return_products) message("Size of 'return_products' may be very large and result in slow calculation time")

  match.arg(weighting_scheme,several.ok = T)
  match.arg(loi_numeric_stats,several.ok = T)

  zip_loc<-input$outfile
  loi_loc<-loi_file
  if (is.null(loi_loc)) loi_loc<-zip_loc

  temp_dir<-file.path(gsub(basename(zip_loc),"",zip_loc),basename(tempfile()))
  if (!dir.exists(temp_dir)) dir.create(temp_dir)

  if (!is.null(dw_dir)){
    dw_fl<-list.files(dw_dir,full.names = T,recursive=T)
  } else {
    dw_fl<-NULL
  }

  fl<-utils::unzip(list=T,zip_loc)
  fl_loi<-utils::unzip(list=T,loi_loc)

  # Get site name column ----------------------------------------------------
  if (any(grepl("snapped_points",fl)) & !all_reaches){
    site_id_col<-c("link_id",paste0(data.table::fread(cmd=paste("unzip -p ",zip_loc,"site_id_col.csv"))))

    if (!is.null(spec)) {
      site_id_col<-site_id_col[site_id_col %in% colnames(spec)]
    } else {
      site_id_col<-"link_id"
    }

  } else {
    site_id_col<-"link_id"
  }

  db_fp<-input$db_loc
  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)
  stream_links<-dplyr::collect(dplyr::tbl(con,"stream_links")) %>%
    dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
    dplyr::mutate(dplyr::across(tidyselect::any_of(site_id_col),dplyr::na_if,""))
  DBI::dbDisconnect(con)

  all_points<-sf::read_sf(file.path("/vsizip",zip_loc,"stream_links.shp"))%>%
    dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>% #stream_links.shp
    dplyr::left_join(stream_links, by = c("link_id"))

  # all_points<-read_sf(file.path("/vsizip",zip_loc,"stream_links.shp")) %>%
  #   left_join(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv")) %>%
  #               mutate(across(any_of(site_id_col),na_if,"")),
  #             by="link_id")

  all_catch<-sf::read_sf(file.path("/vsizip",zip_loc,"Catchment_poly.shp"))

  # Setup remove_region -----------------------------------------------------
  remove_region<-hydroweight::process_input(remove_region,input_name="remove_region")
  if (!is.null(remove_region)){
    if (inherits(remove_region,"SpatRaster")) {
      terra::writeRaster(remove_region,file.path(temp_dir,"remove_region.tif"))
      remove_region<-file.path(temp_dir,"remove_region.tif")
    }
    if (inherits(remove_region,"SpatVector")) {
      terra::writeRaster(remove_region,file.path(temp_dir,"remove_region.shp"))
      remove_region<-file.path(temp_dir,"remove_region.shp")
    }
  }

  #browser()

  # Setup loi  --------------------------------------------------------------

  utils::unzip(loi_loc,file="loi_meta.rds",exdir = temp_dir)

  loi_meta<-readRDS(file.path(temp_dir,"loi_meta.rds"))

  loi_rasts_exists<-purrr::map(loi_meta,~unlist(purrr::map(.,~.$output_filename)))
  loi_rasts_names<-unlist(purrr::map(loi_meta,~unlist(purrr::map(.,~.$lyr_variables))))

  loi_rasts<-purrr::map(loi_rasts_exists,terra::rast)
  loi_rasts_comb<-terra::rast(loi_rasts)
  names(loi_rasts_comb)<-unlist(sapply(loi_rasts,names))

  loi_rasts_exists_names<-loi_rasts_names
  loi_rasts_exists_names<-purrr::map(loi_rasts_exists_names,~purrr::map(.,~stats::setNames(as.list(.),.)) %>% unlist(recursive=T))
  loi_rasts_exists_names<-purrr::map(loi_rasts_exists_names,~purrr::map(.,~loi_numeric_stats))

  if (F){
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
  # loi_rasts_exists_names<-map(loi_rasts,names)
  # loi_rasts_exists_names<-map(loi_rasts_exists_names,~map(.,~setNames(as.list(.),.)) %>% unlist(recursive=T))
  # loi_rasts_exists_names<-map(loi_rasts_exists_names,~map(.,~loi_numeric_stats))
  #
  # loi_rasts_names<-map(loi_rasts,names) %>% unlist()
  # names(loi_rasts_names)<-loi_rasts_names

  target_crs<-terra::crs(loi_rasts[[1]])

  # Setup spec table if missing ---------------------------------------------
  if (is.null(spec)){
    spec<-tibble::tibble(uid=all_points %>%
                   dplyr::select(tidyselect::any_of(site_id_col)) %>%
                   dplyr::filter(!dplyr::if_any(tidyselect::any_of(site_id_col),is.na)) %>%
                   dplyr::pull(1)
    ) %>%
      stats::setNames(site_id_col) %>%
      dplyr::mutate(loi=list(stats::setNames(unlist(loi_rasts_exists_names,recursive = F,use.names = F),loi_rasts_names)))

    if (all_reaches){
      spec<-spec %>%
        dplyr::mutate(link_id=as.character(floor(as.numeric(link_id)))) %>%
        dplyr::distinct()
    }

    spec <- spec %>%
      dplyr::mutate(dplyr::across(c(tidyselect::any_of(site_id_col),tidyselect::any_of("link_id")),as.character))

  }

  if (!any(colnames(spec) %in% "loi")) stop("'spec' must have a column named 'loi'")
  if (!any(colnames(spec) %in% site_id_col)) stop(paste0("'spec' must have a column named '",site_id_col,"'"))

  incorrect_loi<-unique(unlist(sapply(spec$loi,names)))
  incorrect_loi<-incorrect_loi[!incorrect_loi %in% loi_rasts_names]
  if (length(incorrect_loi)>0) stop(paste0("'spec$loi' contains names not in specified loi:", paste0(incorrect_loi,collapse = ", ")))

  if (!any(colnames(spec) == "link_id")){
    spec<-spec %>%
      dplyr::left_join(all_points %>%
                  tibble::as_tibble() %>%
                  dplyr::select(tidyselect::any_of(site_id_col),link_id) %>%
                  dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)),
                by=site_id_col
      )
  }

  # loi_cols<-unlist(spec$loi,recursive=F)
  # loi_cols_group<-map_chr(loi_cols,~paste0(sort(.),collapse=""))
  # loi_cols_group<-split(names(loi_cols_group),loi_cols_group)
  # #loi_cols_group<-map(loi_cols_group,~)
  #
  # map(loi_rasts_exists_names,names) %>%
  #   map()

  # Assemble Output Table ---------------------------------------------------

  if (target_streamseg) {
    target_O<-sf::read_sf(file.path("/vsizip",zip_loc,"stream_lines.shp")) %>%
      dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
      dplyr::left_join(spec %>% dplyr::select(link_id,tidyselect::any_of(site_id_col)),
                by="link_id")

    if (all_reaches){
      target_O<-target_O %>%
        dplyr::mutate(link_id=as.character(floor(as.numeric(link_id)))) %>%
        dplyr::select(link_id) %>%
        dplyr::group_by(link_id) %>%
        dplyr::mutate(geometry=sf::st_union(geometry)) %>%
        dplyr::ungroup()

    }

    # if (any(!spec[[site_id_col]] %in% target_O$link_id)){
    #
    #   ms<-spec[[site_id_col]][!spec[[site_id_col]] %in% target_O$link_id]
    #
    #   extra_o<-all_points %>%
    #     filter(link_id %in% ms) %>%
    #     select(any_of(colnames(target_O)))
    #
    #   target_O<-bind_rows(target_O,extra_o) %>% distinct()
    # }


  } else {
    target_O<-all_points %>%
      dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character))
  }

  target_O<-target_O %>%
    dplyr::filter(link_id %in% spec$link_id)

  #browser()

  # Generate Catchments
  out<-spec %>%
    dplyr::select(tidyselect::any_of(site_id_col),loi)%>%
    dplyr::left_join(all_catch %>%
                dplyr::mutate(link_id=as.character(link_id)) %>%
                dplyr::left_join(all_points %>%
                            tibble::as_tibble() %>%
                            dplyr::mutate(link_id=as.character(link_id)) %>%
                            dplyr::mutate(dplyr::across(tidyselect::any_of(site_id_col),as.character)) %>%
                            dplyr::select(link_id,tidyselect::any_of(site_id_col)),
                          by="link_id") %>%
                dplyr::select(tidyselect::any_of(site_id_col))
              ,by=site_id_col) %>%
    stats::setNames(c("UID","loi","clip_region")) %>%
    dplyr::left_join(target_O %>%
                dplyr::select(tidyselect::any_of(site_id_col)) %>%
                stats::setNames(c("UID","geometry")) %>%
                dplyr::rename(target_O=geometry),
              by="UID")

  if (target_streamseg){
    out<-out%>%
      dplyr::filter(as.numeric(sf::st_length(target_O))>0)
  }


  if (verbose) print("Generating Stream Level Attributes")
  # Get Stream Level Weights
  hw_streams<-hydroweight::hydroweight(hydroweight_dir=temp_dir,
                                       target_S = file.path("/vsizip",zip_loc,"dem_streams_d8.tif"),
                                       target_uid = 'ALL',
                                       OS_combine = FALSE,
                                       dem=file.path("/vsizip",zip_loc,"dem_final.tif"),
                                       flow_accum = file.path("/vsizip",zip_loc,"dem_accum_d8.tif"),
                                       weighting_scheme = weighting_scheme[grepl("lumped|FLS",weighting_scheme)],
                                       inv_function = inv_function,
                                       clean_tempfiles=T,
                                       return_products=F)


  n_cores<-future::nbrOfWorkers()
  if (is.infinite(n_cores)) n_cores<-future::availableCores(logical = F)

  #browser()
  if (verbose) print("Calculating Attributes")
  progressr::with_progress(enable=T,{
  p <- progressr::progressor(steps = nrow(out))

  out2<-out %>%
    dplyr::arrange(dplyr::desc(sf::st_area(clip_region))) %>%
    dplyr::mutate(t_dir=temp_dir,
           dw_dir=ifelse(is.null(dw_dir),NA,dw_dir),
           target_crs=target_crs,
           loi=loi,
           weighting_scheme=rep(list(weighting_scheme),nrow(out)),
           OS_combine=OS_combine,
           site_id_col=site_id_col,
           loi_rasts_exists=rep(list(unlist(loi_rasts_exists)),nrow(out)),
           inv_fun=rep(list(inv_function),nrow(out)),
           loi_rasts_exists_names=rep(list(loi_rasts_names),nrow(out)),
           buffer=catch_buffer,
           return_products=return_products,
           stream_weights=hw_streams,
           p=rep(list(p),nrow(out))
    ) %>%
    dplyr::mutate(core=rep(1:n_cores,length.out=nrow(.))) %>%
    dplyr::group_by(core) %>%
    tidyr::nest() %>%
    dplyr::ungroup() %>%
    dplyr::mutate(data2=purrr::map(data, function(x){
      #browser()
      new_temp_dir<-file.path(x$t_dir[[1]],basename(tempfile()))
      dir.create(new_temp_dir)

      utils::unzip(
        gsub(paste0("/",basename(x$loi_rasts_exists[[1]][[1]])),"",gsub("/vsizip/","",x$loi_rasts_exists[[1]][[1]])),
        files=sapply(x$loi_rasts_exists[[1]], basename),
        exdir = new_temp_dir
      )

      utils::unzip(
        zip_loc,
        files=c("dem_final.tif","dem_accum_d8.tif"),#"dem_streams_d8.tif"
        exdir = new_temp_dir
      )

      utils::unzip(x$stream_weights[[1]],
                   exdir = new_temp_dir)

      stream_w<-file.path(new_temp_dir,utils::unzip(list=T,x$stream_weights[[1]])$Name)

      return(list(stream_w=stream_w,
                  rast_load=as.list(file.path(new_temp_dir,sapply(x$loi_rasts_exists[[1]], basename))) %>%
                    stats::setNames(gsub("\\.tif$","",sapply(x$loi_rasts_exists[[1]], basename))),
                  flow_accum=file.path(new_temp_dir,"dem_accum_d8.tif"),
                  dem=file.path(new_temp_dir,"dem_final.tif"),
                  ts=file.path(new_temp_dir,"dem_streams_d8.tif"),
                  new_temp_dir=new_temp_dir
      )
      )

    }))

  out2<-out2 %>%
    dplyr::mutate(attrib=furrr::future_map2(data2,data,
                              .options = furrr::furrr_options(globals = FALSE),
                              carrier::crate(function(data2,x){
                                # mutate(attrib=map2(data2,data,carrier::crate(function(data2,x){
                                #browser()
                                options(scipen = 999)
                                `%>%` <- magrittr::`%>%`

                                loi_rasts<-purrr::map(data2$rast_load,terra::rast) %>%
                                  stats::setNames(names(x$loi_rasts_exists_names[[1]]))

                                dem<-terra::rast(data2$dem)
                                #ts<-terra::rast(data2$ts)
                                ts<-NA
                                flow_accum<-terra::rast(data2$flow_accum)
                                stream_weights<-purrr::map(data2$stream_w,terra::rast)
                                names(stream_weights)<-sapply(stream_weights,names)

                                purrr::pmap(list(UID=x$UID,
                                                 tar_O=split(x$target_O,1:length(x$target_O)),
                                                 ts=rep(list(ts),nrow(x)),
                                                 clip_region=split(x$clip_region,1:length(x$clip_region)),
                                                 t_dir=rep(data2$new_temp_dir[[1]],nrow(x)),
                                                 target_crs=x$target_crs,
                                                 dw_dir=x$dw_dir,
                                                 loi_cols=x$loi,
                                                 weighting_s=x$weighting_scheme,
                                                 OS_comb=x$OS_combine,
                                                 site_id_c=x$site_id_col,
                                                 loi_rasts=rep(list(loi_rasts),nrow(x)),
                                                 loi_rasts_exists_nm=x$loi_rasts_exists_names,
                                                 buff=x$buffer,
                                                 return_p=x$return_products,
                                                 stream_weights=rep(list(stream_weights),nrow(x)),
                                                 flow_accum=rep(list(flow_accum),nrow(x)),
                                                 dem=rep(list(dem),nrow(x)),
                                                 inv_fun=x$inv_fun,
                                                 p=x$p
                                ),
                                function(UID,
                                         tar_O,
                                         ts,
                                         clip_region,
                                         t_dir,
                                         target_crs,
                                         loi_cols,
                                         weighting_s,
                                         OS_comb,
                                         site_id_c,
                                         loi_rasts,
                                         loi_rasts_exists_nm,
                                         buff,
                                         inv_fun,
                                         return_p,
                                         stream_weights,
                                         flow_accum,
                                         dem,
                                         dw_dir,
                                         p){
                                  #browser()
                                  options(scipen = 999)

                                  attrib_fn<-carrier::crate(function(uid,
                                                                     tar_O,
                                                                     ts,
                                                                     clip_region,
                                                                     loi_cols,
                                                                     t_dir,
                                                                     tar_crs,
                                                                     weighting_s,
                                                                     OS_comb,
                                                                     inv_fun,
                                                                     site_id_c,
                                                                     loi_rasts,
                                                                     loi_rasts_exists_nm,
                                                                     buff,
                                                                     return_p,
                                                                     stream_weights,
                                                                     dem,
                                                                     flow_accum,
                                                                     dw_dir,
                                                                     p) {
                                    #browser()
                                    options(scipen = 999)

                                    `%>%` <- magrittr::`%>%`

                                    save_file<-file.path(t_dir,paste0(uid,"_inv_distances.zip"))

                                    if (!is.null(dw_dir) && !is.na(dw_dir)){
                                      dw_zip<-file.path(dw_dir,paste0(uid,"_inv_distances.zip"))

                                      if (file.exists(file.exists(dw_zip))){
                                        file.copy(
                                          dw_zip,
                                          save_file
                                        )
                                      }
                                    } else {
                                      #browser()
                                      cr<-sf::st_as_sf(tibble::tibble(UID=uid,geometry=sf::st_geometry(clip_region[[1]])),crs=tar_crs)
                                      to<-sf::st_as_sf(tibble::tibble(UID=uid,geometry=sf::st_geometry(tar_O[[1]])),crs=tar_crs)

                                      if (length(weighting_s[grepl("FLO",weighting_s)])==0){
                                        hw<-file.path(t_dir,paste0(uid,"_inv_distances.zip"))
                                      } else {
                                        #browser() #This is broken
                                        hw_zip<-hydroweight::hydroweight(hydroweight_dir=t_dir,
                                                                         target_O = to,
                                                                         #target_S = ts,
                                                                         target_uid = uid,
                                                                         OS_combine = OS_comb,
                                                                         clip_region = sf::st_buffer(cr,units::set_units(buff,"m"),nQuadSegs = 1),
                                                                         dem=dem,
                                                                         flow_accum=flow_accum,
                                                                         weighting_scheme = weighting_s[grepl("FLO",weighting_s)],
                                                                         inv_function = inv_fun,
                                                                         clean_tempfiles=T,
                                                                         return_products = F,
                                                                         wrap_return_products=F,
                                                                         save_output=T)

                                        hw<-purrr::map(file.path("/vsizip",hw_zip,utils::unzip(list=T,hw_zip)$Name),terra::rast)
                                        names(hw)<-sapply(hw,names)

                                        #hw<-purrr::map(hw,terra::unwrap)

                                        # hw_fl<-unzip(list=T,hw)
                                        # hw_fl<-purrr::map(hw_fl,file.path("/vsizip",hw,.))
                                      }

                                      t_dir2<-file.path(t_dir,basename(tempfile()))
                                      dc<-dir.create(t_dir2)

                                      hw_all_strm<-stream_weights

                                      hw_all_strm<-purrr::map(hw_all_strm,~terra::crop(.,
                                                                                       snap="in",
                                                                                       terra::vect(sf::st_buffer(cr,units::set_units(buff,"m"),nQuadSegs = 1)),
                                                                                       mask=T
                                      ))

                                      # hw_all_strm<-purrr::map(hw_all_strm,~terra::writeRaster(.,file.path(t_dir2,paste0(names(.),".tif")),overwrite=T,gdal="COMPRESS=NONE"))
                                      # hw_all_strm<-purrr::map(hw_all_strm,~file.path(t_dir2,paste0(names(.),".tif")))
                                      #
                                      # utils::zip(hw,
                                      #            unlist(hw_all_strm),
                                      #            flags = '-r9Xjq'
                                      # )

                                      save_file<-c(hw,hw_all_strm)

                                    }

                                    if (return_p) {
                                      distance_weights<-save_file
                                      # fls<-utils::unzip(hw,list=T)
                                      # fls<-file.path("/vsizip",hw,fls$Name)
                                      # distance_weights<-lapply(fls,terra::rast)
                                      distance_weights<-lapply(distance_weights,terra::wrap)
                                    } else {
                                      distance_weights<-NULL
                                    }

                                    attr_out<-purrr::pmap(list(loi_nms=as.list(names(loi_rasts)) %>% stats::setNames(names(loi_rasts)),
                                                               loi_path=loi_rasts,
                                                               loi_lyr_nms=loi_rasts_exists_nm,
                                                               loi_cols=rep(list(loi_cols),length(loi_rasts))),
                                                          function(loi_nms,loi_path,loi_lyr_nms,loi_cols){
                                                            #browser()

                                                            # Add another loop for each unique set of stats

                                                            loi_cols_group<-purrr::map_chr(loi_cols,~paste0(sort(.),collapse=""))
                                                            if (loi_nms=="cat_rast") loi_cols_group[]<-"mean"
                                                            loi_cols_group<-split(names(loi_cols_group),loi_cols_group)


                                                            out<-purrr::map(loi_cols_group,function(gp){

                                                              if (grepl("num_rast",loi_nms)) {
                                                                loi_numeric_stats<-loi_cols[names(loi_cols) %in% gp][[1]]
                                                              } else {
                                                                loi_numeric_stats<-NULL
                                                              }

                                                              #browser()

                                                              cls<-names(loi_cols)[names(loi_cols) %in% gp &
                                                                                     names(loi_cols) %in% loi_lyr_nms]

                                                              if (length(cls)==0) return(NULL)

                                                              out<-hydroweight::hydroweight_attributes(
                                                                loi=loi_path,
                                                                loi_columns = cls,
                                                                loi_numeric=grepl("num_rast",loi_nms),
                                                                loi_numeric_stats = loi_numeric_stats,
                                                                #roi=cr,
                                                                roi_uid=uid,
                                                                roi_uid_col=site_id_c,
                                                                distance_weights=save_file,
                                                                #remove_region=remove_region,
                                                                return_products = return_p
                                                              )

                                                              if (return_p) names(out$return_products)<-paste0(names(out$return_products),"_",loi_nms)
                                                              return(out)
                                                            })

                                                            return(out)
                                                          })

                                    file.remove(hw_zip)

                                    out<-list(
                                      attr=unlist(purrr::map(attr_out,~purrr::map(.,~.$attribute_table)),recursive=F) %>%
                                        .[!sapply(.,is.null)] %>%
                                        purrr::reduce(dplyr::left_join,by=site_id_c) %>%
                                        dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~ifelse(is.na(.),"0",.))),
                                      distance_weights=distance_weights,
                                      weighted_attr=unlist(unlist(purrr::map(attr_out,~unlist(purrr::map(.,~.$return_products),recursive = F)),recursive = F),recursive = F)
                                    )

                                    p()

                                    return(out)
                                  })
                                  attrib_fn(uid=UID,
                                            tar_O=tar_O,
                                            ts=ts,
                                            clip_region=clip_region,
                                            loi_cols=loi_cols,
                                            t_dir=t_dir,
                                            tar_crs=target_crs,
                                            weighting_s=weighting_s,
                                            OS_comb=OS_comb,
                                            inv_fun=inv_fun,
                                            site_id_c=site_id_c,
                                            loi_rasts=loi_rasts,
                                            loi_rasts_exists_nm=loi_rasts_exists_nm,
                                            buff=buff,
                                            return_p=return_p,
                                            stream_weights=stream_weights,
                                            flow_accum=flow_accum,
                                            dem=dem,
                                            dw_dir=dw_dir,
                                            p=p)
                                })
                              })))

  out3<-out2%>%
    tidyr::unnest(attrib) %>%
    tidyr::unnest(data) %>%
    tibble::tibble(-core) %>%
    dplyr::arrange(UID)
  })

  #browser()
  out4<-out3 %>%
    #select(-target_O,-loi,-clip_region,-UID) %>%
    dplyr::select(attrib) %>%
    dplyr::mutate(attrib_out=purrr::map(attrib,~.$attr)) %>%
    dplyr::mutate(attrib_out=purrr::map(attrib_out,~dplyr::mutate(.,dplyr::across(c(tidyselect::everything(),-tidyselect::any_of(site_id_col)),as.numeric)))) %>%
    dplyr::mutate(distance_weights=purrr::map(attrib,~.$distance_weights)) %>%
    dplyr::mutate(weighted_attr=purrr::map(attrib,~.$weighted_attr)) %>%
    tidyr::unnest(attrib_out) %>%
    dplyr::select(tidyselect::any_of(site_id_col),tidyselect::any_of("distance_weights"),tidyselect::any_of("weighted_attr"),tidyselect::everything(),-attrib) %>%
    dplyr::mutate(dplyr::across(c(tidyselect::everything(),-tidyselect::any_of(site_id_col),-tidyselect::any_of("distance_weights"),-tidyselect::any_of("weighted_attr")),as.numeric)) %>%
    dplyr::mutate(dplyr::across(tidyselect::ends_with("_prop"),~dplyr::case_when(is.na(.) | is.nan(.) ~ 0, T ~ .))) %>%
    dplyr::group_by(!!dplyr::sym(site_id_col)) %>%
    dplyr::summarise(dplyr::across(tidyselect::everything(),utils::head,1)) %>%
    dplyr::ungroup()
  #dplyr::distinct()


  zip_loc<-input$outfile

  out_file<-zip_loc

  if (verbose) print("Generating Output")

  saveRDS(out4,file.path(temp_dir,"point_attributes.rds"))

  utils::zip(out_file,
      file.path(temp_dir,"point_attributes.rds"),
      flags = '-r9Xjq'
  )

  output<-input

  output<-c(
    list(attrib=out4),
    output
  )

  suppressWarnings(file.remove(list.files(temp_dir,recursive = T,full.names = T)))

  return(out4)
}


