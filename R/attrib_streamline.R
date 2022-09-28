
#' Generate and attribute stream line and points along the lines
#'
#' @param input resulting object from `generate_subbasins()`
#' @param extra_attr character. One or more of c("link_slope", "cont_slope", "USChnLn_To", "Elevation", "StOrd_Hack", "StOrd_Str", "StOrd_Hort", "StOrd_Shr"). Optional attributes to add to stream vector outputs.
#' @param points character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/points.shp"), or any GIS data object that will be converted to spatial points. Points representing sampling locations.
#' @param site_id_col character. Variable name in `points` that corresponds to unique site identifiers. This column will be included in all vector geospatial analysis products. Note, if multiple points have the same `site_id_col`, their centroid will be used and returned; if multiple points overlap after snapping, only the first is used.
#' @param snap_distance integer. Maximum distance which points will be snapped to stream lines in map units
#' @param break_on_noSnap logical. Should the function stop if any points are not snapped to a stream segment (i.e., are beyon snap_distance)
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param compress logical. Should output rasters be compressed, slower but more space efficient.
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
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
    points=NULL,
    site_id_col=NULL,
    snap_distance=10L,
    break_on_noSnap=T,
    return_products=F,
    temp_dir=NULL,
    compress=F,
    verbose=F
) {
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")

  if (!is.null(points) & is.null(site_id_col)) stop("`site_id_col` can not be NULL if `points` is provided")

  extra_attr<-match.arg(extra_attr,several.ok = T)

  if (!is.null(points) & !is.numeric(snap_distance)) stop("'snap_distance' must be an integer value")
  if (!is.null(points) & !is.logical(break_on_noSnap)) stop("'break_on_noSnap' must be logical")

  if (!is.null(points)) {
    if (!is.character(site_id_col)) stop("'site_id_col' must be a single character")
    if (length(site_id_col)>1 | length(site_id_col)==0) stop("length 'site_id_col' must be 1")
    if (site_id_col=="link_id") stop("'site_id_col' cannot be 'link_id'")
  } else {
    site_id_col<-"link_id"
  }

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")
  if (!is.logical(compress)) stop("'compress' must be logical")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  wbt_options(exe_path=wbt_exe_path(),
              verbose=verbose,
              wd=temp_dir,
              compress_rasters =compress)

  terra::terraOptions(verbose = verbose,
                      tempdir = temp_dir
  )

  gdal_arg<-NULL
  if (compress){
    gdal_arg<-"COMPRESS=NONE"
  }

  options(scipen = 999)
  options(dplyr.summarise.inform = FALSE)

  zip_loc<-input$outfile
  fl<-unzip(list=T,zip_loc)

  unzip(zip_loc,
        c("dem_d8.tif","dem_streams_d8.tif","dem_final.tif"),
        exdir=temp_dir,
        overwrite=T,
        junkpaths=T)

  dem_final<-rast(file.path(temp_dir,"dem_final.tif"))
  names(dem_final)<-"Elevation"
  writeRaster(dem_final,file.path(temp_dir,"Elevation.tif"),overwrite=T,gdal=gdal_arg)

  #subb<-st_as_sf(vect(file.path("/vsizip",zip_loc,"Subbasins_poly.shp")))

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
  strm<-read_sf(file.path(temp_dir, "strm_link_id.shp")) %>%
    select(STRM_VAL)
  st_crs(strm)<-crs(dem_final)
  colnames(strm)[1]<-"link_id"

  #saveRDS(strm,file.path(temp_dir, "strm_link_id.rds"))

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

  id3<-suppressWarnings(suppressMessages(terra::extract(do.call(c,setNames(attr,NULL)),id2)))

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


  # Instert Points if present -----------------------------------------------

  if (!is.null(points)){
    points<-hydroweight::process_input(points,target = vect(final_points),input_name="points")

    write_sf(st_as_sf(points) %>%
               select(any_of(site_id_col),everything()),
             file.path(temp_dir,"original_points.shp"))

    points<-st_as_sf(points) %>%
      group_by(!!rlang::sym(site_id_col)) %>%
      summarize(across(geometry,st_centroid)) %>%
      ungroup()

    if (!site_id_col %in% names(points)) stop("'site_id_col' must be a variable name in 'points'")

    print("Snapping Points")
    snapped_points<-points %>%
      st_join(final_points %>%
                filter(!link_type %in% c("Link Node","Sink Node","Source Node (head water)")) %>%
                select(ID,link_id),
              join=nngeo::st_nn,
              maxdist = snap_distance,
              progress =T
      ) %>%
      as_tibble() %>%
      select(-geometry) %>%
      left_join(final_points %>%
                  filter(!link_type %in% c("Link Node","Sink Node","Source Node (head water)")) %>%
                  select(ID,link_id),
                by = c("ID", "link_id")) %>%
      st_as_sf()


    if (break_on_noSnap){
      if (any(is.na(snapped_points$link_id))) stop(paste0("The following points could not be snapped: ", paste0(snapped_points[[site_id_col]][is.na(snapped_points$link_id)],collapse = ", ") ))
    }

    if (any(is.na(snapped_points$link_id))) warning(paste0("The following points could not be snapped and were not included: ", paste0(snapped_points[[site_id_col]][is.na(snapped_points$link_id)],collapse = ", ") ))

    snapped_points<-snapped_points %>%
      filter(!is.na(link_id)) %>%
      select(any_of(site_id_col),everything())

    new_final_points<-final_points
    new_final_points$link_class[new_final_points$ID %in% snapped_points$ID]<-6
    new_final_points$link_type[new_final_points$ID %in% snapped_points$ID]<-"Sample Point"
    new_final_points<-new_final_points %>%
      left_join(
        snapped_points %>% as_tibble() %>% select(ID,any_of(site_id_col)),
        by="ID"
      ) %>%
      filter(!is.na(!!sym(site_id_col))) %>%
      group_by(link_id) %>%
      arrange(link_id,USChnLn_Fr) %>%
      mutate(link_id_new=row_number()) %>%
      mutate(link_id_new=formatC(link_id_new,width=nchar(max(link_id_new)),format="d",flag=0)) %>%
      mutate(link_id_new=as.numeric(paste0(paste0(link_id),".",link_id_new)))

    final_points<-final_points %>%
      filter(!ID %in% new_final_points$ID) %>%
      bind_rows(new_final_points) %>%
      arrange(link_id,desc(USChnLn_Fr)) %>%
      group_by(link_id) %>%
      mutate(link_id_new=case_when(
        row_number()==1 & is.na(link_id_new) ~ link_id,
        T ~ link_id_new
      )) %>%
      fill(link_id_new,.direction = "down") %>%
      select(-link_id) %>%
      rename(link_id=link_id_new) %>%
      select(link_id,everything()) %>%
      ungroup() %>%
      arrange(ID)

    write_sf(final_points %>%
               select(link_id,any_of(site_id_col)) %>%
               filter(!is.na(link_id))
               ,file.path(temp_dir,"snapped_points.shp"))


    # Make new stream raster --------------------------------------------------

    strm<-final_points %>%
      # select(link_id,USChnLn_Fr) %>%
      # group_by(link_id) %>%
      # arrange(USChnLn_Fr) %>%
      # select(-USChnLn_Fr) %>%
      # summarise(do_union=F) %>%
      # st_cast("LINESTRING")
      select(link_id) %>%
      arrange(link_id) %>%
      vect() %>%
      rasterize(y=dem_final,field="link_id") %>%
      #hydroweight::process_input(target=dem_final,input_name="New Stream Raster",resample_type="near") %>%
      writeRaster(file.path(temp_dir,"new_stream_layer.tif"),overwrite=T,gdal="COMPRESS=NONE")

    wbt_raster_streams_to_vector(
      streams = file.path(temp_dir,"new_stream_layer.tif"),
      d8_pntr= file.path(temp_dir,"dem_d8.tif"),
      output = file.path(temp_dir,"new_stream_layer.shp")
    )

    # for some reason wbt_raster_streams_to_vector() rounds numbers weirdly

    un_ID<-unique(final_points$link_id)
    strm<-read_sf(file.path(temp_dir,"new_stream_layer.shp")) %>%
      select(STRM_VAL) %>%
      rename(link_id=STRM_VAL) %>%
      rowwise() %>%
      mutate(link_id=un_ID[which.min(abs(link_id-un_ID))]) %>%
      ungroup()
    st_crs(strm)<-crs(dem_final)

    write_sf(strm,file.path(temp_dir,"stream_lines.shp"))
    saveRDS(strm,file.path(temp_dir, "strm_link_id.rds"))

  }

  # Add columns: for us and ds link_id and trib_id ---------------------------------

  st_r<-rast(file.path(temp_dir, "dem_streams_d8.tif"))
  d8_pntr<-rast(file.path(temp_dir, "dem_d8.tif"))

  # Upstream ----------------------------------------------------------
  if (verbose) print("Identifying Upstream Links")
  #browser()

  nodes<-final_points %>%
    group_by(link_id) %>%
    filter(USChnLn_Fr==min(USChnLn_Fr)) %>%
    ungroup() %>%
    # filter(link_class %in% c(4,5)) %>%
    select(ID) %>%
    arrange(ID) %>%
    vect()

  cv <- cells(d8_pntr,nodes)
  all_cell <- cells(d8_pntr,vect(final_points))
  all_cell <- as_tibble(mutate(final_points,cell=all_cell[,2]))

  cv1 <- tibble(target=cv[,2]) %>%
    bind_cols(terra::adjacent(d8_pntr, cells=.$target, directions="queen") %>%
                data.frame() %>%
                setNames(c("TL","TM","TR","ML","MR","BL","BM","BR"))) %>%
    gather("direction","cell_num",-target) %>%
    arrange(target) %>%
    filter(!is.na(cell_num),!is.nan(cell_num)) %>%
    mutate(on_stream=data.frame(terra::extract(st_r,xyFromCell(st_r,.$cell_num)))$dem_streams_d8) %>%
    filter(on_stream==1) %>%
    left_join(
      all_cell %>% select(link_id,cell),
      by=c("cell_num"="cell")
    )%>%
    left_join(
      all_cell %>% select(trib_id,cell),
      by=c("cell_num"="cell")
    )%>%
    left_join(
      all_cell %>% select(target_link_id=link_id,cell),
      by=c("target"="cell")
    )%>%
    left_join(
      all_cell %>% select(target_trib_id=trib_id,cell),
      by=c("target"="cell")
    )

  cv2<-cv1  %>%
    mutate(flow_dir=terra::extract(d8_pntr,xyFromCell(d8_pntr,.$cell_num))$dem_d8) %>%
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
    bind_cols(terra::adjacent(d8_pntr, cells=.$target, directions="queen")%>%
                data.frame() %>%
                setNames(c("TL","TM","TR","ML","MR","BL","BM","BR"))) %>%
    gather("direction","cell_num",-target) %>%
    arrange(target) %>%
    mutate(on_stream=terra::extract(st_r,xyFromCell(st_r,.$cell_num))$dem_streams_d8) %>%
    filter(on_stream==1) %>%
    left_join(
      all_cell %>% select(link_id,cell),
      by=c("cell_num"="cell")
    )%>%
    left_join(
      all_cell %>% select(trib_id,cell),
      by=c("cell_num"="cell")
    )%>%
    left_join(
      all_cell %>% select(target_link_id=link_id,cell),
      by=c("target"="cell")
    )%>%
    left_join(
      all_cell %>% select(target_trib_id=trib_id,cell),
      by=c("target"="cell")
    )

  cv2<-cv1  %>%
    mutate(flow_dir=terra::extract(d8_pntr,xyFromCell(d8_pntr,.$target))$dem_d8) %>%
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
  links<-final_points %>%
    group_by(link_id) %>%
    filter(USChnLn_Fr==max(USChnLn_Fr)) %>%
    ungroup()

  final_links<-links %>%
    full_join(final_us, by = c("link_id", "trib_id")) %>%
    full_join(final_ds, by = c("link_id", "trib_id")) #

  # check both ID and link_id are unique
  check_link_id<-any(duplicated(final_links$link_id))
  check_id<-any(duplicated(final_links$ID))
  if (check_link_id | check_id) warning("Some link_id's and/or ID's are duplicated, this may indicate an issue with the upstream/downstream IDs")


  write_sf(final_links %>% select(link_id),file.path(temp_dir,"stream_links.shp"))
  write_sf(final_points %>% select(link_id),file.path(temp_dir,"stream_points.shp"))

  data.table::fwrite(final_links %>% as_tibble() %>% select(-geometry),file.path(temp_dir,"stream_links.csv"))
  data.table::fwrite(final_points %>% as_tibble() %>% select(-geometry),file.path(temp_dir,"stream_points.csv"))
  data.table::fwrite(tibble(site_id_col=site_id_col),file.path(temp_dir,"site_id_col.csv"))

  # Generate Output ---------------------------------------------------------
  if (verbose) print("Generating Output")


  dist_list_out<-c(
    list.files(temp_dir,"site_id_col"),
    list.files(temp_dir,"snapped_points"),
    list.files(temp_dir,"original_points"),
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

  output<-input[!names(input) %in% c("stream_lines",
                                     "links",
                                     "points",
                                     "snapped_points",
                                     "original_points"
                                     )]

  if (return_products){
    output<-c(
      list(
        stream_lines=strm,
        links=final_links,
        points=final_points,
        snapped_points=snapped_points,
        original_points=points
      ),
      output
    )
  }

  file.remove(list.files(temp_dir,full.names = T,recursive = T))

  return(output)
}
