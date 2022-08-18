#' Process DEM by generating various terrain products
#'
#' #' \code{IntegratedHydrology::process_hydrology()} processes a DEM into various geospatial analysis products.
#'
#' @param dem character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/dem.tif"), \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of GeoTiFF type. Digital elevation model raster.
#' @param threshold numeric. Threshold in flow accumulation values for channelization.
#' @param save_dir character. File path to write resulting .zip file.
#' @param return_products logical. If \code{TRUE}, a list containing all geospatial analysis products. If \code{FALSE}, folder path to resulting .zip file.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical. If \code{FALSE}, the function will not print output prints.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'
#' @examples
#'
#' library(whitebox)
#' library(terra)
#' library(sf)
#' library(hydroweight)
#'
#' ## Import toy_dem from whitebox package
#' toy_file<-sample_dem_data()
#' toy_file <- system.file("extdata", "DEM.tif", package = "whitebox")
#' toy_dem <- rast(raster::raster(x = toy_file)) # reading the file from terra directly sometimes crashes R for some reason
#' crs(toy_dem) <- "epsg:3161"
#'
#' ## Generate hydroweight_dir as a temporary directory
#' save_dir <- tempdir()
#'
#' hydro_out<-process_hydrology(
#'   dem=toy_dem,
#'   threshold=1000,
#'   return_products=T,
#'   save_dir=save_dir,
#'   temp_dir=NULL,
#'   verbose=F
#' )
#'
#
process_hydrology<-function(
    dem,
    threshold,
    save_dir,
    return_products=F,
    temp_dir=NULL,
    verbose=F
) {

  require(sf)
  require(terra)
  require(whitebox)
  require(tidyverse)

  if (!is.integer(threshold)) stop("'threshold' must be an integer value")
  if (!dir.exists(save_dir)) dir.create(save_dir)
  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)

  wbt_options(exe_path=wbt_exe_path(),
              verbose=verbose,
              wd=temp_dir)

  dem<-hydroweight::process_input(dem,input_name="dem")
  if (!inherits(dem,"SpatRaster")) stop("dem must be a class 'SpatRaster'")
  target_crs<-crs(dem)

  writeRaster(dem,file.path(temp_dir,"dem_final.tif"),overwrite=T)


  # Process flow dirrection -------------------------------------------------
  ## Generate d8 pointer
  if (verbose) print("Generating d8 pointer")
  wbt_d8_pointer(
    dem="dem_final.tif",
    output="dem_d8.tif"
  )

  ## Generate d8 flow accumulation in units of cells
  if (verbose) print("Generating d8 flow accumulation")
  wbt_d8_flow_accumulation(
    input = "dem_d8.tif",
    output = "dem_accum_d8.tif",
    out_type = "cells",
    pntr=T
  )
  wbt_d8_flow_accumulation(
    input = "dem_d8.tif",
<<<<<<< HEAD
    output = "dem_accum_d8_sca.tif",
=======
    output = "dem_accum_d8.tif",
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
    out_type = "sca",
    pntr=T
  )


  ## Generate streams with a stream initiation threshold
  if (verbose) print("Extracting Streams")
  wbt_extract_streams(
    flow_accum = "dem_accum_d8.tif",
    output = "dem_streams_d8.tif",
    threshold = threshold
  )

  # Generate subbasin polygons ----------------------------------------------
  if (verbose) print("Generating subbasins")
  wbt_subbasins(
    d8_pntr="dem_d8.tif",
    streams="dem_streams_d8.tif",
    output="Subbasins.tif"
  )

  if (verbose) print("Converting subbasins to polygons")
  subb<-rast(file.path(temp_dir,"Subbasins.tif"))
  subb<-as.polygons(subb,dissolve = TRUE)
  subb<-st_as_sf(subb) %>%
    mutate(sbbsn_area=st_area(.))
  names(subb)[1]<-"link_id"

  write_sf(subb,file.path(temp_dir,"Subbasins_poly.shp"))

  # Attribute Stream network ------------------------------------------------
  if (verbose) print("Calculating stream link attributes")

  wbt_stream_link_identifier(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="link_id.tif"
  )

  wbt_stream_link_class(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="link_class.tif"
  )

  wbt_stream_link_length(
    d8_pntr= "dem_d8.tif",
    linkid="link_id.tif",
    output="link_length.tif"
  )

  wbt_stream_link_slope( # slope gradient in degrees
    d8_pntr= "dem_d8.tif",
    linkid="link_id.tif",
    dem= "dem_final.tif",
    output="link_slope.tif"
  )
  # # convert above to % change, %slope= tan(Angle in degrees*pi/180)*100

  wbt_tributary_identifier(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="trib_id.tif"
  )

  wbt_hack_stream_order(
    streams =  "dem_streams_d8.tif",
    d8_pntr =  "dem_d8.tif",
    output =  "StOrd_Hack.tif"
  )
  wbt_strahler_stream_order(
    streams =  "dem_streams_d8.tif",
    d8_pntr =  "dem_d8.tif",
    output =  "StOrd_Str.tif"
  )
  wbt_horton_stream_order(
    streams =  "dem_streams_d8.tif",
    d8_pntr =  "dem_d8.tif",
    output =  "StOrd_Hort.tif"
  )
  wbt_shreve_stream_magnitude(
    streams =  "dem_streams_d8.tif",
    d8_pntr =  "dem_d8.tif",
    output =  "StOrd_Shr.tif"
  )

  wbt_stream_slope_continuous(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    dem= "dem_final.tif",
    output="Stream_Slope.tif"
  )
  # # convert above to % change, %slope= tan(Angle in degrees*pi/180)*100

  wbt_length_of_upstream_channels(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="USChanLen_Total.tif"
  )

  wbt_farthest_channel_head(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="USChanLen_Furthest.tif"
  )

  ## Generate stream lines from raster ---------------
  wbt_raster_streams_to_vector(
    streams = "link_id.tif",
    d8_pntr = "dem_d8.tif",
    output = "strm_link_id.shp"
  )
  strm<-read_sf(file.path(temp_dir, "strm_link_id.shp"))
  st_crs(strm)<-crs(dem)
  colnames(strm)[2]<-"link_id"

  write_rds(strm,file.path(temp_dir, "strm_link_id.rds"))

  # Generate point attributions ---------------------------------------------
  if (verbose) print("Extracting stream link attributes")

  browser()
  id1<-rast(file.path(temp_dir,"link_id.tif"))
  id2<-as.points(id1)
  id21<-st_as_sf(id2)

  attr<-c(
    rast(file.path(temp_dir,"link_class.tif")),
    rast(file.path(temp_dir,"link_length.tif")),
    rast(file.path(temp_dir,"link_slope.tif")),
    rast(file.path(temp_dir,"trib_id.tif")),
    rast(file.path(temp_dir,"USChanLen_Total.tif")),
    rast(file.path(temp_dir,"USChanLen_Furthest.tif")),
    rast(file.path(temp_dir,"dem_final.tif")),
    rast(file.path(temp_dir,"StOrd_Hack.tif")),
    rast(file.path(temp_dir,"StOrd_Str.tif")),
    rast(file.path(temp_dir,"StOrd_Hort.tif")),
    rast(file.path(temp_dir,"StOrd_Shr.tif")),
    rast(file.path(temp_dir,"Stream_Slope.tif"))
  )


  id3<-terra::extract(attr,id2)

  id4<-id3 %>%
    mutate(link_type=case_when(
      link_class==1~"Exterior Link",
      link_class==2~"Interior Link",
      link_class==3~"Source Node (head water)",
      link_class==4~"Link Node",
      link_class==5~"Sink Node",
    )) %>%
    tibble() %>%
    select(ID,link_class,link_type,everything())

  names(id4)<-abbreviate(names(id4),10)

  final_points<-bind_cols(id21,id4)

  links<-final_points %>%
    group_by(link_id) %>%
    filter(USChnLn_Fr==max(USChnLn_Fr)) %>%
    ungroup()

  # Add columns: for us and ds link_id and trib_id ---------------------------------

  r1<-rast(file.path(temp_dir,"link_class.tif"))
  r2<-rast(file.path(temp_dir,"link_id.tif"))
  r3<-rast(file.path(temp_dir,"trib_id.tif"))

  st_r<-rast(file.path(temp_dir, "dem_streams_d8.tif"))

  d8_pntr=rast(file.path(temp_dir, "dem_d8.tif"))

  # Upstream ----------------------------------------------------------
  if (verbose) print("Identifying Upstream Links")

  nodes<-final_points %>%
    filter(link_class==4) %>%
    select(ID) %>%
    arrange(ID) %>%
    vect()

  cv <- cells(d8_pntr,nodes)

  cv1 <- tibble(target=cv[,2]) %>%
    bind_cols(terra::adjacent(d8_pntr, cells=.$target, directions="queen")) %>%
    rename(TL=...2,
           TM=...3,
           TR=...4,
           ML=...5,
           MR=...6,
           BL=...7,
           BM=...8,
           BR=...9) %>%
    gather("direction","cell_num",-target) %>%
    arrange(target) %>%
    mutate(on_stream=terra::extract(st_r,.$cell_num)$dem_streams_d8) %>%
    filter(on_stream==1) %>%
    mutate(link_id=terra::extract(r2,.$cell_num)$link_id) %>%
    mutate(trib_id=terra::extract(r3,.$cell_num)$trib_id) %>%
    mutate(target_link_id=terra::extract(r2,.$target)$link_id) %>%
    mutate(target_trib_id=terra::extract(r3,.$target)$trib_id)

  cv2<-cv1  %>%
    mutate(flow_dir=terra::extract(d8_pntr,.$cell_num)$dem_d8) %>%
    mutate(flow_in=case_when(
      direction=="TL" & flow_dir == 4 ~ T,
      direction=="TM" & flow_dir == 8 ~ T,
      direction=="TR" & flow_dir == 16 ~ T,
      direction=="ML" & flow_dir == 2 ~ T,
      direction=="MR" & flow_dir == 32 ~ T,
      direction=="BL" & flow_dir == 1 ~ T,
      direction=="BM" & flow_dir == 128 ~ T,
      direction=="BR" & flow_dir == 64 ~ T,
      T ~ F
    )) %>%
    group_by(target) %>%
    mutate(n_inflows=sum(flow_in)) %>%
    ungroup()

  final_us<-cv2 %>%
    filter(flow_in) %>%
    select(n_inflows,link_id,trib_id,target_link_id,target_trib_id) %>%
    group_by(target_link_id,target_trib_id) %>%
    mutate(nm=row_number()) %>%
    ungroup() %>%
    rename(ustrib_id = trib_id,
           uslink_id=link_id) %>%
    pivot_wider(id_cols = c(target_link_id,target_trib_id),
                names_from = c(nm),
                names_sep = "",
                values_from = c(ustrib_id,uslink_id)) %>%
    rename(link_id =target_link_id,
           trib_id=target_trib_id)

  names(final_us)<-abbreviate(names(final_us),10)


  # Downstream -----------------------------------------------------------
  if (verbose) print("Identifying Downstream Links")

  nodes<-final_points %>%
    group_by(link_id) %>%
    filter(USChnLn_Fr==max(USChnLn_Fr)) %>%
    ungroup() %>%
    select(ID) %>%
    arrange(ID) %>%
    vect()

  cv <- cells(d8_pntr,nodes)

  cv1 <- tibble(target=cv[,2]) %>%
    bind_cols(terra::adjacent(d8_pntr, cells=.$target, directions="queen")) %>%
    rename(TL=...2,
           TM=...3,
           TR=...4,
           ML=...5,
           MR=...6,
           BL=...7,
           BM=...8,
           BR=...9) %>%
    gather("direction","cell_num",-target) %>%
    arrange(target) %>%
    mutate(on_stream=terra::extract(st_r,.$cell_num)$dem_streams_d8) %>%
    filter(on_stream==1) %>%
    mutate(link_id=terra::extract(r2,.$cell_num)$link_id) %>%
    mutate(trib_id=terra::extract(r3,.$cell_num)$trib_id) %>%
    mutate(target_link_id=terra::extract(r2,.$target)$link_id) %>%
    mutate(target_trib_id=terra::extract(r3,.$target)$trib_id)

  cv2<-cv1  %>%
    mutate(flow_dir=terra::extract(d8_pntr,.$target)$dem_d8) %>%
    mutate(flow_in=case_when(
      direction=="TL" & on_stream & flow_dir == 64 ~ T,
      direction=="TM" & on_stream & flow_dir == 128 ~ T,
      direction=="TR" & on_stream & flow_dir == 1 ~ T,
      direction=="ML" & on_stream & flow_dir == 32 ~ T,
      direction=="MR" & on_stream & flow_dir == 2 ~ T,
      direction=="BL" & on_stream & flow_dir == 16 ~ T,
      direction=="BM" & on_stream & flow_dir == 8 ~ T,
      direction=="BR" & on_stream & flow_dir == 4 ~ T,
      T ~ F
    )) %>%
    group_by(target) %>%
    mutate(n_inflows=sum(flow_in)) %>%
    ungroup()

  final_ds<-cv2 %>%
    filter(flow_in) %>%
    select(n_inflows,link_id,trib_id,target_link_id,target_trib_id) %>%
    group_by(target_link_id,target_trib_id) %>%
    mutate(nm=row_number()) %>%
    ungroup() %>%
    rename(dstrib_id = trib_id,
           dslink_id=link_id) %>%
    pivot_wider(id_cols = c(target_link_id,target_trib_id),
                names_from = c(nm),
                names_sep = "",
                values_from = c(dstrib_id,dslink_id)) %>%
    rename(link_id =target_link_id,
           trib_id=target_trib_id)

  names(final_ds)<-abbreviate(names(final_ds),10)


  # Putting it all together -------------------------------------------------
<<<<<<< HEAD
#browser()
=======
browser()
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
  final_links<-links %>%
    full_join(final_us) %>%
    full_join(final_ds)

  # check both ID and link_id are unique
  check_link_id<-any(duplicated(final_links$link_id))
  check_id<-any(duplicated(final_links$ID))
  if (check_link_id | check_id) warning("Some link_id's and/or ID's are duplicated, this may indicate an issue with the upstream/downstream IDs")

  final_lines<-strm %>%
    left_join(final_points %>%
                as_tibble() %>%
                filter(grepl("Link",link_type)) %>%
                select(link_id,trib_id,link_lngth,link_slope,USChnLn_Tt,USChnLn_Fr,StOrd_Hack,StOrd_Str,StOrd_Hort,StOrd_Shr) %>%
                group_by(link_id,trib_id) %>%
                summarise(across(everything(),mean)) %>%
                distinct()) %>%
    left_join(final_links %>%
                as_tibble() %>%
                select(link_id,trib_id,starts_with("dstrib_id"),starts_with("dslink_id"),starts_with("ustrib_id"),starts_with("uslink_id"))) %>%
    left_join(subb %>%
                as_tibble() %>%
                select(link_id,sbbsn_area)
    )

  write_sf(final_links,file.path(temp_dir,"stream_links.shp"))
  write_sf(final_lines,file.path(temp_dir,"stream_lines.shp"))
  write_sf(final_points,file.path(temp_dir,"stream_points.shp"))


  # Trace Flow Paths --------------------------------------------------------
  if (verbose) print("Tracing Flow Paths")

  ds_flowpaths<-trace_ds_flowpath(final_links)
  us_flowpaths<-trace_us_flowpath(final_links)

  if (verbose) print("Calculating Pairwise Distances")
  pwise_dist<-pairwise_dist(ds_flowpaths=ds_flowpaths,
                            us_flowpaths=us_flowpaths,
                            subbasins=subb)

  saveRDS(ds_flowpaths,file.path(temp_dir,"ds_flowpaths.rds"))
  saveRDS(us_flowpaths,file.path(temp_dir,"us_flowpaths.rds"))

  write_csv(pwise_dist$long_pwise,file.path(temp_dir,"long_pwise.csv"))
  write_csv(pwise_dist$wide_prop_shared_catchment,file.path(temp_dir,"wide_prop_shared_catchment.csv"))
  write_csv(pwise_dist$wide_directed_path_length,file.path(temp_dir,"wide_directed_path_length.csv"))


  # Generate Output ---------------------------------------------------------
  if (verbose) print("Generating Output")

  dist_list_out<-list(
    "dem_final.tif",
    "dem_d8.tif",
    "dem_accum_d8.tif",
    "dem_streams_d8.tif",
    "ds_flowpaths.rds",
    "us_flowpaths.rds",
    "long_pwise.csv",
    "wide_prop_shared_catchment.csv",
    "wide_directed_path_length.csv"
  )

  dist_list_out<-c(
    dist_list_out,
    list.files(temp_dir,"Subbasins_poly"),
    list.files(temp_dir,"stream_links"),
    list.files(temp_dir,"stream_lines"),
    list.files(temp_dir,"stream_points")
  )

  dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

  out_file<-file.path(save_dir,paste0("Processed_Hydrology.zip"))
  if (file.exists(out_file)) {
    out_file<-file.path(save_dir,paste0(basename(tempfile()), "Processed_Hydrology.zip"))
    warning(paste0("Target .zip file already exists. Saving as: ",out_file))
  }

  zip(out_file,
      unlist(dist_list_out),
      flags = '-r9Xjq'
  )

  output<-list(
    outfile=out_file
  )

  if (return_products){
    output<-c(
      list(
        subbasins=subb,
        stream_lines=final_lines,
        links=final_links,
        points=final_points,
        ds_flowpaths=ds_flowpaths,
        us_flowpaths=ds_flowpaths,
        pwise_dist=pwise_dist
      ),
      output
    )
  }

  return(output)

  file.remove(list.files(temp_dir,full.names = T))
}
