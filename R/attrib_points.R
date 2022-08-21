
#' Title
#'
#' @param input
#' @param spec
#' @param weighting_scheme
#' @param loi_numeric_stats
#' @param OS_combine_points
#' @param OS_combine_subbasin
#' @param points_target_streamseg
#' @param subbasin_target_streamseg
#' @param inv_function
#' @param remove_region
#' @param temp_dir
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples

attrib_points<-function(
    input,
    spec=NULL,
    all_reaches=F,
    weighting_scheme = c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS"),
    loi_numeric_stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max", "sum", "cell_count"),
    OS_combine=F,
    target_streamseg=F,
    buffer=30,
    inv_function = function(x) {
      (x * 0.001 + 1)^-1
    },
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

  temp_dir<-file.path(gsub(basename(zip_loc),"",zip_loc),basename(tempfile()))
  if (!dir.exists(temp_dir)) dir.create(temp_dir)

  fl<-unzip(list=T,zip_loc)
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

  browser()

  # Setup loi  --------------------------------------------------------------
  loi_rasts_exists<-c("numeric_rasters.tif","cat_rasters.tif")
  if (!any(loi_rasts_exists %in% fl$Name)) stop("No 'loi' present in input, please run 'process_loi()' first")
  loi_rasts_exists<-fl$Name[grepl("numeric_rasters|cat_rasters",fl$Name)]
  loi_rasts_exists<-map(loi_rasts_exists,~file.path("/vsizip",zip_loc,.))
  names(loi_rasts_exists)<-gsub("\\.tif","",sapply(loi_rasts_exists,basename))

  loi_rasts_exists_names<-map(loi_rasts_exists,~rast(.) %>% names())
  loi_rasts_exists_names<-map(loi_rasts_exists_names,~map(.,~setNames(as.list(.),.)) %>% unlist(recursive=T))
  loi_rasts_exists_names<-map(loi_rasts_exists_names,~map(.,~c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max", "sum", "cell_count")))

  loi_rasts<-map(loi_rasts_exists,rast)
  loi_rasts_names<-map(loi,names) %>% unlist()
  names(loi_rasts_names)<-loi_rasts_names

  target_crs<-crs(loi_rasts[[1]])

  # Setup spec table if missing ---------------------------------------------
  if (is.null(spec)){
    spec<-tibble(uid=all_points %>%
                   select(site_id_col) %>%
                   filter(!if_any(site_id_col,is.na)) %>%
                   pull(1)
    ) %>%
      setNames(site_id_col) %>%
      mutate(loi=list(loi_rasts_exists_names))
  }

  if (!any(colnames(spec) %in% "loi")) stop("'spec' must have a column named 'loi'")
  if (!any(colnames(spec) %in% site_id_col)) stop(paste0("'spec' must have a column named '",site_id_col,"'"))

  incorrect_loi<-unique(unlist(sapply(spec$loi,names)))
  incorrect_loi<-incorrect_loi[!incorrect_loi %in% loi_rasts_names]
  if (length(incorrect_loi)>0) stop(paste0("'spec$loi' contains names not in specified loi:", paste0(incorrect_loi,collapse = ", ")))





  # Assemble Output Table ---------------------------------------------------

  if (target_streamseg) {
    target_O<-read_sf(file.path("/vsizip",zip_loc,"stream_lines.shp"))
  } else {
    target_O<-all_points
  }

  # with_progress({
  #   p <- progressor(steps = nrow(spec))

    out<-spec %>%
      setNames(c("UID","loi")) %>%
      rowwise() %>%
      left_join(target_O %>%
                  select(any_of(site_id_col)) %>%
                  setNames(c("UID","geometry")) %>%
                  rename(target_O=geometry)) %>%
      mutate(clip_region=future_map(UID,~get_catchment(input=input,
                                                       site_id_col=site_id_col,
                                                       target_sites=.) %>%
                                      select(geometry) %>%
                                      rename(clip_region=geometry))) %>%
      unnest(clip_region) %>%
      mutate(attrib=pmap(list(uid=UID,
                                     target_O=target_O,
                                     clip_region=clip_region,
                                     loi_cols=loi),
                                function(uid,target_O,clip_region,loi_cols) {
                                  browser()

                                  save_file<-file.path(temp_dir,paste0(uid,"_inv_distances.zip"))

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
                                                   browser()

                                                   # Add another loop for each unique set of stats

                                                   loi_cols_group<-map_chr(loi_cols,paste0,collapse="")
                                                   loi_cols_group<-split(names(loi_cols_group),loi_cols_group)

                                                   out<-map(loi_cols_group,function(gp){
                                                     browser()
                                                     out<-hydroweight::hydroweight_attributes(
                                                       loi=loi_path,
                                                       loi_columns = names(loi_cols)[names(loi_cols) %in% gp],
                                                       loi_numeric=grepl("numeric_rasters",loi_nms),
                                                       loi_numeric_stats = loi_cols[names(loi_cols) %in% gp][[1]],
                                                       roi=cr,
                                                       roi_uid=uid,
                                                       roi_uid_col=site_id_col,
                                                       distance_weights=save_file,
                                                       remove_region=remove_region,
                                                       return_products = return_products
                                                     )

                                                     names(out$return_products)<-paste0(names(out$return_products),"_",loi_nms)
                                                   })



                                                   return(out)
                                                 })

                                  out<-list(
                                    attr=map(attr_out,~.$attribute_table) %>% reduce(left_join,by=site_id_col),
                                    distance_weights=distance_weights,
                                    weighted_attr=unlist(unlist(map(attr_out,~.$return_products),recursive = F),recursive = F)
                                  )

                                  #p()

                                  return(out)

                                })
      )
  # })

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
