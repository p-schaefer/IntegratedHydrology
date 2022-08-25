
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
#' @param buffer numeric. Amount to buffer the catchment (in meters) when calculating `hydroweight::hydroweight()`. `hydroweight::hydroweight()` sometimes misses point `target_O` values unless buffered. Doesn't affect attribute values
#' @param tolerance numeric. Tolerance values used for \code{sf::st_snap} in meters in `get_catchment()`.
#' @param catch_buffer numeric. Tolerance values used for \code{sf::st_buffer} in meters in `get_catchment()`.
#' @param return_products logical. If \code{TRUE}, a list containing all geospatial analysis products. If \code{FALSE}, folder path to resulting .zip file.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'
#' @examples

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
    buffer=30,
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
    tolerance=0.000001,
    catch_buffer=0.001,
    return_products=F,
    remove_region=NULL,
    temp_dir=NULL,
    verbose=F
){

  if (!is.null(spec) && !inherits(spec,"data.frame")) stop("'spec' must be a data frame")
  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (return_products) warning("Size of 'return_products' may be very large and result in slow calculation time")

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

  fl<-unzip(list=T,zip_loc)
  fl_loi<-unzip(list=T,loi_loc)
  all_points<-read_sf(file.path("/vsizip",zip_loc,"stream_links.shp"))

  # Get site name column ----------------------------------------------------
  if (any(grepl("snapped_points",fl)) & !all_reaches){
    site_id_col<-colnames(read_sf(file.path("/vsizip",zip_loc,"snapped_points.shp")))[[1]]
  } else {
    site_id_col<-"link_id"
  }

  # Setup remove_region -----------------------------------------------------
  remove_region<-hydroweight::process_input(remove_region,input_name="remove_region")
  if (!is.null(remove_region)){
    if (inherits(remove_region,"SpatRaster")) {
      writeRaster(remove_region,file.path(temp_dir,"remove_region.tif"))
      remove_region<-file.path(temp_dir,"remove_region.tif")
    }
    if (inherits(remove_region,"SpatVector")) {
      writeRaster(remove_region,file.path(temp_dir,"remove_region.shp"))
      remove_region<-file.path(temp_dir,"remove_region.shp")
    }
  }

  #browser()

  # Setup loi  --------------------------------------------------------------
  loi_rasts_exists<-c("num_rast.tif","cat_rast.tif")
  if (!any(loi_rasts_exists %in% fl_loi$Name)) stop("No 'loi' present in input, please run 'process_loi()' first, or specify location of process_loi() ouput")
  loi_rasts_exists<-fl_loi$Name[grepl("num_rast|cat_rast",fl_loi$Name)]
  loi_rasts_exists<-map(loi_rasts_exists,~file.path("/vsizip",loi_loc,.))
  names(loi_rasts_exists)<-gsub("\\.tif","",sapply(loi_rasts_exists,basename))

  loi_rasts_exists_names<-map(loi_rasts_exists,~rast(.) %>% names())
  loi_rasts_exists_names<-map(loi_rasts_exists_names,~map(.,~setNames(as.list(.),.)) %>% unlist(recursive=T))
  loi_rasts_exists_names<-map(loi_rasts_exists_names,~map(.,~loi_numeric_stats))

  loi_rasts<-map(loi_rasts_exists,rast)
  loi_rasts_names<-map(loi_rasts,names) %>% unlist()
  names(loi_rasts_names)<-loi_rasts_names

  target_crs<-crs(loi_rasts[[1]])

  # Setup spec table if missing ---------------------------------------------
  if (is.null(spec)){
    spec<-tibble(uid=all_points %>%
                   select(any_of(site_id_col)) %>%
                   filter(!if_any(site_id_col,is.na)) %>%
                   pull(1)
    ) %>%
      setNames(site_id_col) %>%
      mutate(loi=list(setNames(unlist(loi_rasts_exists_names,recursive = F,use.names = F),loi_rasts_names)))

    if (all_reaches){
      spec<-spec %>%
        mutate(link_id=floor(link_id)) %>%
        distinct()
    }
  }

  if (!any(colnames(spec) %in% "loi")) stop("'spec' must have a column named 'loi'")
  if (!any(colnames(spec) %in% site_id_col)) stop(paste0("'spec' must have a column named '",site_id_col,"'"))

  incorrect_loi<-unique(unlist(sapply(spec$loi,names)))
  incorrect_loi<-incorrect_loi[!incorrect_loi %in% loi_rasts_names]
  if (length(incorrect_loi)>0) stop(paste0("'spec$loi' contains names not in specified loi:", paste0(incorrect_loi,collapse = ", ")))

  # loi_cols<-unlist(spec$loi,recursive=F)
  # loi_cols_group<-map_chr(loi_cols,~paste0(sort(.),collapse=""))
  # loi_cols_group<-split(names(loi_cols_group),loi_cols_group)
  # #loi_cols_group<-map(loi_cols_group,~)
  #
  # map(loi_rasts_exists_names,names) %>%
  #   map()

  # Assemble Output Table ---------------------------------------------------

  if (target_streamseg) {
    target_O<-read_sf(file.path("/vsizip",zip_loc,"stream_lines.shp"))

    if (all_reaches){
      target_O<-target_O %>%
        mutate(link_id=floor(link_id)) %>%
        select(link_id) %>%
        distinct()
    }
  } else {
    target_O<-all_points
  }

  with_progress({
    p <- progressor(steps = nrow(spec))

    out<-spec %>%
      setNames(c("UID","loi")) %>%
      rowwise() %>%
      left_join(target_O %>%
                  select(any_of(site_id_col)) %>%
                  setNames(c("UID","geometry")) %>%
                  rename(target_O=geometry)) %>%
      mutate(clip_region=future_map(UID,~get_catchment(input=input,
                                                       site_id_col=site_id_col,
                                                       target_points=.,
                                                       tolerance =tolerance ,
                                                       buffer =catch_buffer ) %>%
                                      select(geometry) %>%
                                      rename(clip_region=geometry))) %>%
      unnest(clip_region) %>%
      mutate(attrib=future_pmap(list(uid=UID,
                                     target_O=target_O,
                                     clip_region=clip_region,
                                     loi_cols=loi),
                                function(uid,target_O,clip_region,loi_cols) {

                                  save_file<-file.path(temp_dir,paste0(uid,"_inv_distances.zip"))

                                  dw_zip<-file.path(dw_dir,paste0(uid,"_inv_distances.zip"))

                                  if (!is.null(dw_dir) && file.exists(dw_zip)){
                                    file.copy(
                                      dw_zip,
                                      save_file
                                    )
                                  } else {
                                    cr<-st_as_sf(tibble(UID=uid,geometry=st_geometry(clip_region)),crs=crs(target_crs))
                                    to<-st_as_sf(tibble(UID=uid,geometry=st_geometry(target_O)),crs=crs(target_crs))

                                    hw<-hydroweight::hydroweight(hydroweight_dir=temp_dir,
                                                                 target_O = to,
                                                                 target_S = file.path("/vsizip",zip_loc,"dem_streams_d8.tif"),
                                                                 target_uid = uid,
                                                                 OS_combine = OS_combine,
                                                                 clip_region = st_buffer(cr,units::set_units(buffer,"m"),nQuadSegs = 1),
                                                                 dem=file.path("/vsizip",zip_loc,"dem_final.tif"),
                                                                 flow_accum = file.path("/vsizip",zip_loc,"dem_accum_d8.tif"),
                                                                 weighting_scheme = weighting_scheme,
                                                                 inv_function = inv_function,
                                                                 clean_tempfiles=T,
                                                                 return_products=F)
                                  }

                                  if (return_products) {
                                    fls<-unzip(hw,list=T)
                                    fls<-file.path("/vsizip",hw,fls$Name)
                                    distance_weights<-lapply(fls,rast)
                                    distance_weights<-lapply(distance_weights,wrap)
                                    names(distance_weights)<-gsub("//.tif","",basename(fls))
                                  } else {
                                    distance_weights<-NULL
                                  }

                                  attr_out<-pmap(list(loi_nms=as.list(names(loi_rasts_exists)) %>% setNames(names(loi_rasts_exists)),
                                                      loi_path=loi_rasts_exists,
                                                      loi_lyr_nms=loi_rasts_exists_names),
                                                 function(loi_nms,loi_path,loi_lyr_nms){
                                                   #browser()

                                                   # Add another loop for each unique set of stats

                                                   loi_cols_group<-map_chr(loi_cols,~paste0(sort(.),collapse=""))
                                                   if (loi_nms=="cat_rast") loi_cols_group[]<-"mean"
                                                   loi_cols_group<-split(names(loi_cols_group),loi_cols_group)


                                                   out<-map(loi_cols_group,function(gp){

                                                     if (grepl("num_rast",loi_nms)) {
                                                       loi_numeric_stats<-loi_cols[names(loi_cols) %in% gp][[1]]
                                                     } else {
                                                       loi_numeric_stats<-NULL
                                                     }

                                                     #browser()
                                                     out<-hydroweight::hydroweight_attributes(
                                                       loi=loi_path,
                                                       loi_columns = names(loi_cols)[names(loi_cols) %in% gp &
                                                                                       names(loi_cols) %in% names(loi_lyr_nms)],
                                                       loi_numeric=grepl("num_rast",loi_nms),
                                                       loi_numeric_stats = loi_numeric_stats,
                                                       roi=cr,
                                                       roi_uid=uid,
                                                       roi_uid_col=site_id_col,
                                                       distance_weights=save_file,
                                                       remove_region=remove_region,
                                                       return_products = return_products
                                                     )

                                                     if (return_products) names(out$return_products)<-paste0(names(out$return_products),"_",loi_nms)
                                                     return(out)
                                                   })



                                                   return(out)
                                                 })

                                  out<-list(
                                    attr=unlist(map(attr_out,~map(.,~.$attribute_table)),recursive=F) %>% reduce(left_join,by=site_id_col),
                                    distance_weights=distance_weights,
                                    weighted_attr=unlist(unlist(map(attr_out,~unlist(map(.,~.$return_products),recursive = F)),recursive = F),recursive = F)
                                  )

                                  p()

                                  return(out)

                                })
      )
  })

  out<-out %>%
    select(-target_O,-loi,-clip_region,-UID) %>%
    mutate(attrib_out=map(attrib,~.$attr)) %>%
    mutate(distance_weights=map(attrib,~.$distance_weights)) %>%
    mutate(weighted_attr=map(attrib,~.$weighted_attr)) %>%
    unnest(attrib_out) %>%
    select(site_id_col,distance_weights,weighted_attr,everything(),-attrib)

  file.remove(list.files(temp_dir,recursive = T,full.names = T))

  return(out)
}
