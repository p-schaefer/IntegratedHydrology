
<<<<<<< HEAD
#' Title
#'
#' @param input
#' @param extra_attr
#' @param return_products
#' @param temp_dir
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
#'

attrib_streamline<-function(
    input,
    extra_attr=c(
      "link_slope",
      "cont_slope",
      "USChnLn_To",
      "Elevation",
      "StOrd_Hack",
      "StOrd_Str",
      "StOrd_Hort",
      "StOrd_Shr"
    ),
    return_products=F,
    temp_dir=NULL,
    verbose=F
=======

attrib_streamline<-function(
    input,
  extra_attrib=c("")
>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
) {
  require(sf)
  require(terra)
  require(whitebox)
  require(tidyverse)

<<<<<<< HEAD
  extra_attr<-match.arg(extra_attr,several.ok = T)

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  wbt_options(exe_path=wbt_exe_path(),
              verbose=verbose,
              wd=temp_dir)

  terra::terraOptions(verbose = verbose,
                      tempdir = temp_dir
  )

  zip_loc<-input$outfile
  fl<-unzip(list=T,zip_loc)

  unzip(zip_loc,
        c("dem_d8.tif","dem_streams_d8.tif","dem_final.tif"),
        exdir=temp_dir,
        overwrite=T,
        junkpaths=T)

  # dem_d8_streams<-rast(file.path("/vsizip",zip_loc,"dem_streams_d8.tif"))
  # writeRaster(dem_d8_streams,file.path(temp_dir,"dem_streams_d8.tif"),overwrite=T,gdal="COMPRESS=NONE")
  dem_final<-rast(file.path(temp_dir,"dem_final.tif"))
  names(dem_final)<-"Elevation"
  writeRaster(dem_final,file.path(temp_dir,"Elevation.tif"),overwrite=T,gdal="COMPRESS=NONE")

  subb<-st_as_sf(vect(file.path("/vsizip",zip_loc,"Subbasins_poly.shp")))

  # Attribute Stream network ------------------------------------------------

  # # Essential attributes
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
    output="link_lngth.tif"
  )

  wbt_tributary_identifier(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="trib_id.tif"
  )

  wbt_farthest_channel_head(
    d8_pntr= "dem_d8.tif",
    streams= "dem_streams_d8.tif",
    output="USChnLn_Fr"
  )

  # Extra Attributes --------------------------------------------------------
  if ("link_slope" %in% extra_attr)
    wbt_stream_link_slope( # slope gradient in degrees
      d8_pntr= "dem_d8.tif",
      linkid="link_id.tif",
      dem= "Elevation.tif",
      output="link_slope.tif"
    )
  # # convert above to % change, %slope= tan(Angle in degrees*pi/180)*100
  if ("cont_slope" %in% extra_attr)
    wbt_stream_slope_continuous(
      d8_pntr= "dem_d8.tif",
      streams= "dem_streams_d8.tif",
      dem= "Elevation.tif",
      output="cont_slope.tif"
    )
  # # convert above to % change, %slope= tan(Angle in degrees*pi/180)*100

  if ("StOrd_Hack" %in% extra_attr)
    wbt_hack_stream_order(
      streams =  "dem_streams_d8.tif",
      d8_pntr =  "dem_d8.tif",
      output =  "StOrd_Hack.tif"
    )
  if ("StOrd_Str" %in% extra_attr)
    wbt_strahler_stream_order(
      streams =  "dem_streams_d8.tif",
      d8_pntr =  "dem_d8.tif",
      output =  "StOrd_Str.tif"
    )

  if ("StOrd_Hort" %in% extra_attr)
    wbt_horton_stream_order(
      streams =  "dem_streams_d8.tif",
      d8_pntr =  "dem_d8.tif",
      output =  "StOrd_Hort.tif"
    )

  if ("StOrd_Shr" %in% extra_attr)
    wbt_shreve_stream_magnitude(
      streams =  "dem_streams_d8.tif",
      d8_pntr =  "dem_d8.tif",
      output =  "StOrd_Shr.tif"
    )

  if ("USChnLn_To" %in% extra_attr)
    wbt_length_of_upstream_channels(
      d8_pntr= "dem_d8.tif",
      streams= "dem_streams_d8.tif",
      output="USChnLn_To.tif"
    )



  ## Generate stream lines from raster ---------------
  wbt_raster_streams_to_vector(
    streams = "link_id.tif",
    d8_pntr = "dem_d8.tif",
    output = "strm_link_id.shp"
  )
  strm<-read_sf(file.path(temp_dir, "strm_link_id.shp"))
  st_crs(strm)<-crs(dem_final)
  colnames(strm)[2]<-"link_id"

  write_rds(strm,file.path(temp_dir, "strm_link_id.rds"))

  # Generate point attributions ---------------------------------------------
  if (verbose) print("Extracting stream link attributes")

  id1<-rast(file.path(temp_dir,"link_id.tif"))
  id2<-as.points(id1)
  id21<-st_as_sf(id2)

  attr_main<-c(
    file.path(temp_dir,"link_class.tif"),
    file.path(temp_dir,"link_lngth.tif"),
    file.path(temp_dir,"trib_id.tif"),
    file.path(temp_dir,"USChnLn_Fr.tif")
  )

  extra_attr<-file.path(temp_dir,paste0(extra_attr,".tif"))
  extra_attr<-extra_attr[file.exists(extra_attr)]

  extra_attr<-c(attr_main,extra_attr)
  names(extra_attr)<-gsub("\\.tif","",basename(extra_attr))

  attr<-lapply(extra_attr,function(x) rast(x))

  id3<-terra::extract(do.call(c,setNames(attr,NULL)),id2)

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

  r1<- attr$link_class #rast(file.path(temp_dir,"link_class.tif"))
  r2<- id1  #rast(file.path(temp_dir,"link_id.tif"))
  r3<- attr$trib_id #rast(file.path(temp_dir,"trib_id.tif"))

  st_r<-rast(file.path(temp_dir, "dem_streams_d8.tif"))
  d8_pntr<-rast(file.path(temp_dir, "dem_d8.tif"))

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
  #browser()
  final_links<-links %>%
    full_join(final_us) %>%
    full_join(final_ds)

  # check both ID and link_id are unique
  check_link_id<-any(duplicated(final_links$link_id))
  check_id<-any(duplicated(final_links$ID))
  if (check_link_id | check_id) warning("Some link_id's and/or ID's are duplicated, this may indicate an issue with the upstream/downstream IDs")

  #browser()
  final_lines<-strm %>%
    left_join(final_points %>%
                as_tibble() %>%
                filter(grepl("Link",link_type)) %>%
                select(link_id,trib_id,link_lngth,link_slope,any_of(names(attr))) %>%
                group_by(link_id,trib_id) %>%
                summarise(across(everything(),mean)) %>%
                distinct()) %>%
    left_join(final_links %>%
                as_tibble() %>%
                select(link_id,trib_id,starts_with("dstrib_id"),starts_with("dslink_id"),starts_with("ustrib_id"),starts_with("uslink_id"))) %>%
    left_join(subb %>%
                as_tibble() %>%
                select(link_id,sbbsn_area)
    ) %>%
    select(link_id,trib_id,link_lngth,sbbsn_area,USChnLn_Fr,everything())

  write_sf(final_links,file.path(temp_dir,"stream_links.shp"))
  write_sf(final_lines,file.path(temp_dir,"stream_lines.shp"))
  write_sf(final_points,file.path(temp_dir,"stream_points.shp"))

  # Generate Output ---------------------------------------------------------
  if (verbose) print("Generating Output")


  dist_list_out<-c(
    list.files(temp_dir,"stream_links"),
    list.files(temp_dir,"stream_lines"),
    list.files(temp_dir,"stream_points")
  )

  dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

  out_file<-zip_loc

  zip(out_file,
      unlist(dist_list_out),
      flags = '-r9Xjq'
  )

  output<-input

  if (return_products){
    output<-c(
      list(
        stream_lines=final_lines,
        links=final_links,
        points=final_points
      ),
      output
    )
  }

  file.remove(list.files(temp_dir,full.names = T))

  return(output)
=======


>>>>>>> c50188ded16dc12e032435fcda4ceaffb482418d
}
