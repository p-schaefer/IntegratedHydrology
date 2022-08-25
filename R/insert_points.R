
#' Inserts sampling points into vector layers, splitting lines, points, and catchments at insertion points
#'
#' @param input resulting object from `attrib_streamline()`
#' @param points character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/points.shp"), or any GIS data object that will be converted to spatial points. Points representing sampling locations.
#' @param site_id_col character. Variable name in `points` that corresponds to unique site identifiers. This column will be included in all vector geospatial analysis products. Note, if multiple points have the same `site_id_col`, their centroid will be used and returned; if multiple points overlap after snapping, only the first is used.
#' @param snap_distance integer. Maximum distance which points will be snapped to stream lines in map units
#' @param return_products logical. If \code{TRUE}, a list containing all geospatial analysis products. If \code{FALSE}, folder path to resulting .zip file.
#' @param temp_dir character. File path for intermediate products; these are deleted once the function runs successfully.
#' @param verbose logical. If \code{FALSE}, the function will not print output prints.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'
#' @examples
#'

insert_points<-function(
    input,
    points,
    site_id_col,
    snap_distance,
    return_products=F,
    temp_dir=NULL,
    verbose=F
) {
  # require(sf)
  # require(terra)
  # require(whitebox)
  # require(tidyverse)

  if (!is.integer(snap_distance)) stop("'snap_distance' must be an integer value")

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")
  if (!is.character(site_id_col)) stop("'site_id_col' must be a single character")
  if (length(site_id_col)>1 | length(site_id_col)==0) stop("length 'site_id_col' must be 1")
  if (site_id_col=="link_id") stop("'site_id_col' cannot be 'link_id'")
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
  fl<-unzip(list=T,zip_loc,overwrite = T,junkpaths =F)

  if (verbose) print("Extracting Data")
  unzip(zip_loc,
        c("dem_streams_d8.tif","dem_d8.tif",fl$Name[grepl("stream_links|stream_points|Subbasins_poly|stream_lines",fl$Name)]),
        exdir=temp_dir,
        overwrite=T,
        junkpaths=T)

  stream_links<-read_sf(file.path(temp_dir,"stream_links.shp"))
  points<-hydroweight::process_input(points,target = stream_links,input_name="points")

  points<-st_as_sf(points) %>%
    group_by(!!rlang::sym(site_id_col)) %>%
    summarize(across(geometry,st_centroid)) %>%
    ungroup()

  if (!site_id_col %in% names(points)) stop("'site_id_col' must be a variable name in 'points'")


  stream_points<-read_sf(file.path(temp_dir,"stream_points.shp"))
  Subbasins_poly<-read_sf(file.path(temp_dir,"Subbasins_poly.shp"))
  stream_lines<-read_sf(file.path(temp_dir,"stream_lines.shp"))

  print("Snapping Points")

  snapped_points<-points %>%
    st_join(stream_points %>%
              filter(link_class %in% c(1,2)) %>%
              select(ID,link_id),
            join=nngeo::st_nn,
            maxdist = snap_distance,
            progress =T#,
            #parallel = future::nbrOfWorkers()
            ) %>%
    as_tibble() %>%
    select(-geometry) %>%
    left_join(stream_points %>%
                filter(link_class %in% c(1,2)) %>%
                select(ID,link_id),
              by = c("ID", "link_id")) %>%
    st_as_sf() #%>%
    #select(-ID)

  if (any(is.na(snapped_points$link_id))) warning(paste0("The following points could not be snapped and were not included: ", paste0(snapped_points[[site_id_col]][is.na(snapped_points$link_id)],collapse = ", ") ))

  snapped_points<-snapped_points %>%
    filter(!is.na(link_id))

  write_sf(snapped_points,file.path(temp_dir,"snapped_points.shp"))


  # Split subbasins by sampling points ----------------------------------
  new_subbasins_fn<-function(l_id,
                             pnt,
                             p,
                             stream_links=file.path(temp_dir,"stream_links.shp"),
                             Subbasins_poly=file.path(temp_dir,"Subbasins_poly.shp"),
                             temp_dir=temp_dir){
    #browser()

    stream_links<-hydroweight::process_input(stream_links,"stream_links")
    Subbasins_poly<-hydroweight::process_input(Subbasins_poly,"Subbasins_poly")
    stream_links<-st_as_sf(stream_links)
    Subbasins_poly<-st_as_sf(Subbasins_poly)

    pnt_file<-file.path(temp_dir,paste0("Tempsite_",l_id,".shp"))

    pnt<-pnt %>%
      mutate(across(any_of(site_id_col),as.character)) %>%
      bind_rows(
        stream_links %>%
          filter(link_id==l_id) %>%
          select(geometry) %>%
          mutate(UID=paste0("pour_point_A1B2C3_",l_id)) %>%
          setNames(c("geometry",site_id_col))
      ) %>%
      mutate(x=st_coordinates(geometry)[,1],
             y=st_coordinates(geometry)[,2]) %>%
      group_by(x,y) %>%
      summarize(across(everything(),head,1)) %>%
      select(-x,-y) %>%
      ungroup()

    write_sf(pnt,pnt_file)

    cr<-Subbasins_poly %>% filter(link_id == l_id) %>% vect()

    rast(file.path(temp_dir,"dem_d8.tif")) %>%
      crop(y=cr,
           mask=T,
           overwrite=T,
           filename=file.path(temp_dir,paste0("d8_temp_",l_id,".tif"))) %>%
      writeRaster(
        filename=file.path(temp_dir,paste0("d8_",l_id,".tif")),
        overwrite=T
      )

    wbt_unnest_basins(
      d8_pntr=file.path(temp_dir,paste0("d8_",l_id,".tif")),
      pour_pts=pnt_file,
      output=file.path(temp_dir,paste0("Catch_",l_id,".tif"))
    )

    catch_fls<-list.files(temp_dir,pattern=paste0("Catch_",l_id,"_"))
    catch_rast<-rast(file.path(temp_dir,catch_fls))
    catch_rast[is.na(catch_rast)]<-0
    catch_rast[catch_rast>0]<-1
    catch_rast<-app(catch_rast,sum)

    catch_rast[catch_rast==0]<-NA

    #browser()
    catch_poly<-catch_rast %>%
      as.polygons() %>%
      st_as_sf() %>%
      st_join(pnt) %>%
      arrange(sum) %>%
      mutate(link_id=row_number()-1) %>%
      mutate(link_id=formatC(link_id,width=nchar(max(link_id)),format="d",flag=0)) %>%
      mutate(link_id=as.numeric(paste0(l_id,".",link_id))) %>%
      mutate(sbbsn_area=st_area(.)) %>%
      select(link_id,sbbsn_area,any_of(site_id_col), geometry) %>%
      st_intersection(Subbasins_poly %>% filter(link_id == l_id) %>% select(geometry))

    catch_poly[grepl("pour_point_A1B2C3_",catch_poly%>% as_tibble() %>% select(any_of(site_id_col)) %>% unlist()),site_id_col]<-NA_character_

    p()

    flrm<-unique(c(list.files(temp_dir,pattern=paste0("_",l_id,"."),full.names = T),
                   list.files(temp_dir,pattern=paste0("_",l_id,"_"),full.names = T)
    ))

    file.remove(flrm)
    return(catch_poly)
  }

  # Split points by sampling points ----------------------------------
  new_points_fn<-function(x,
                          p,
                          stream_points=file.path(temp_dir,"stream_points.shp"),
                          temp_dir=temp_dir) {

    stream_points<-hydroweight::process_input(stream_points,"stream_points")
    stream_points<-st_as_sf(stream_points)

    out<-x %>%
      select(-sbbsn_area) %>%
      rename(link_id_new=link_id) %>%
      st_buffer(-units::as_units(0.5,"m"),nQuadSegs = 1)

    out2<-stream_points %>%
      filter(link_id == min(out$link_id_new)) %>%
      select(geometry) %>%
      st_join(out,join=st_nearest_feature)

    out2[grepl("pour_point_A1B2C3_",out2%>% as_tibble() %>% select(any_of(site_id_col)) %>% unlist()),site_id_col]<-NA_character_

    p()
    return(out2)
  }

  # Split lines by sampling points ----------------------------------
  new_lines_fn<-function(x,
                         p,
                         stream_lines=file.path(temp_dir,"stream_lines.shp"),
                         temp_dir=temp_dir
  ) {

    stream_lines<-hydroweight::process_input(stream_lines,"stream_lines")
    stream_lines<-st_as_sf(stream_lines)

    trg_strm<-stream_lines %>%
      filter(link_id == min(x$link_id))

    out <- x %>%
      select(-sbbsn_area) %>%
      rename(link_id_new=link_id) %>%
      st_buffer(-units::as_units(0.5,"m"),nQuadSegs = 1)

    out2 <- st_intersection(trg_strm,out)
    missing_pieces<-st_difference(trg_strm %>% select(geometry),x %>% select(geometry) %>% st_union())

    missing_pieces<-st_join(missing_pieces,out2,join=st_nearest_feature)

    out3<-bind_rows(
      out2 %>% select(link_id_new),
      missing_pieces %>% select(link_id_new)
    ) %>%
      group_by(link_id_new) %>%
      summarize() %>%
      select(geometry) %>%
      st_join(out2,join=st_nearest_feature) %>%
      mutate(dslink_id1=ifelse(link_id_new==min(link_id_new),dslink_id1,lag(link_id_new))) %>%
      mutate(uslink_id1=ifelse(link_id_new==max(link_id_new),uslink_id1,lead(link_id_new))) %>%
      mutate(across(c(starts_with("uslink_id"),-uslink_id1),~ifelse(link_id_new==max(link_id_new),.,NA))) %>%
      select(-link_id) %>%
      rename(link_id=link_id_new) %>%
      select(any_of(colnames(trg_strm)),everything()) %>%
      mutate(across(c(contains('uslink_id'),contains('dslink_id')),~ifelse(.==0,NA_real_,.))) %>%
      mutate(across(c(everything(),-geometry),~ifelse(is.nan(.),NA,.)))

    out3[grepl("pour_point_A1B2C3_",out3%>% as_tibble() %>% select(any_of(site_id_col)) %>% unlist()),site_id_col]<-NA_character_

    p()
    return(out3)
  }

  # Split links by sampling points ----------------------------------
  new_links_fn<-function(pnts,
                         lns,
                         p,
                         stream_links=file.path(temp_dir,"stream_links.shp"),
                         temp_dir=temp_dir){

    stream_links<-hydroweight::process_input(stream_links,"stream_links")
    stream_links<-st_as_sf(stream_links)

    if (nrow(lns)==1) {
      replace_target<-stream_links %>%
        filter(link_id %in% min(lns$link_id)) %>%
        st_join(pnts %>% mutate(across(any_of(site_id_col),as.character)))

      return(replace_target)
    }

    lns_target<-lns %>%
      as_tibble() %>%
      filter(!if_any(any_of(site_id_col),is.na))

    replace_target<-stream_links %>%
      filter(link_id %in% min(lns$link_id)) %>%
      mutate(uslink_id1=min(lns$link_id[-c(1)])) %>%
      mutate(across(c(starts_with("uslink_id"),-uslink_id1),~NA)) %>%
      mutate(uslink_id1=ifelse(is.infinite(abs(uslink_id1)),NA,uslink_id1)) %>%
      st_join(pnts %>% mutate(across(any_of(site_id_col),as.character)))

    add_links<-pnts %>%
      mutate(across(any_of(site_id_col),as.character)) %>%
      bind_rows(
        stream_links %>%
          filter(link_id==min(lns$link_id)) %>%
          select(geometry) %>%
          mutate(UID=paste0("pour_point_A1B2C3_",min(lns$link_id))) %>%
          setNames(c("geometry",site_id_col))
      ) %>%
      st_join(lns, join=st_nearest_feature) %>%
      arrange(link_id) %>%
      mutate(link_id_new=link_id) %>%
      mutate(link_id=min(link_id)) %>%
      select(link_id,link_id_new,any_of(site_id_col)) %>%
      left_join(stream_links %>% as_tibble() %>% select(-geometry)) %>%
      mutate(dslink_id1=ifelse(link_id_new==min(link_id_new),dslink_id1,lag(link_id_new))) %>%
      mutate(uslink_id1=ifelse(link_id_new==max(link_id_new),uslink_id1,lead(link_id_new))) %>%
      mutate(across(c(starts_with("uslink_id"),-uslink_id1),~ifelse(link_id_new==max(link_id_new),.,NA))) %>%
      #filter(link_id_new %in% max(link_id_new)) %>%
      select(-link_id) %>%
      rename(link_id=link_id_new) %>%
      st_join(pnts %>% mutate(across(any_of(site_id_col),as.character))) %>%
      filter(!if_any(any_of(site_id_col),is.na)) %>%
      distinct() %>%
      filter(!is.na(link_id)) %>%
      filter(st_geometry(geometry)!=st_geometry(replace_target))

    out3<-bind_rows(replace_target,add_links) %>%
      distinct() %>%
      mutate(across(c(contains('uslink_id'),contains('dslink_id')),~ifelse(.==0,NA_real_,.))) %>%
      mutate(across(c(everything(),-geometry),~ifelse(is.nan(.),NA,.))) %>%
      distinct()

    out3[grepl("pour_point_A1B2C3_",out3%>% as_tibble() %>% select(any_of(site_id_col)) %>% unlist()),site_id_col]<-NA_character_

    p()
    return(out3)
  }

  if (verbose) print("Inserting Points")
  # Adjusts flow tracers in segment containing points -----------------------


  new_data<-snapped_points %>%
    select(any_of(site_id_col),link_id) %>%
    group_by(link_id) %>%
    nest()

  with_progress({
    print("Splitting Subbasins")
    p <- progressor(steps = nrow(new_data))

    new_data<-new_data  %>%
      mutate(new_subbasins=future_map2(link_id,data,~new_subbasins_fn(l_id=.x,
                                                                      pnt=.y,
                                                                      p=p,
                                                                      stream_links=file.path(temp_dir,"stream_links.shp"),
                                                                      Subbasins_poly=file.path(temp_dir,"Subbasins_poly.shp"),
                                                                      temp_dir=temp_dir)))
  })

  with_progress({
    print("Splitting Points")
    p <- progressor(steps = nrow(new_data))

    new_data<-new_data %>%
      mutate(new_points=future_map(new_subbasins,~new_points_fn(x=.,
                                                                p=p,
                                                                stream_points=file.path(temp_dir,"stream_points.shp"),
                                                                temp_dir=temp_dir)))
  })

  with_progress({
    print("Splitting Lines")
    p <- progressor(steps = nrow(new_data))

    new_data<-new_data %>%
      mutate(new_lines=future_map(new_subbasins, ~new_lines_fn(x=.,
                                                               p=p,
                                                               stream_lines=file.path(temp_dir,"stream_lines.shp"),
                                                               temp_dir=temp_dir)))
  })

  #browser()
  with_progress({
    print("Splitting Links")
    p <- progressor(steps = nrow(new_data))

    new_data<-new_data %>%
      mutate(new_links=future_map2(data,new_lines,~new_links_fn(pnts=.x,
                                                                lns=.y,
                                                                p=p,
                                                                stream_links=file.path(temp_dir,"stream_links.shp"),
                                                                temp_dir=temp_dir)))
  })

  # browser()
  # Add to master data
  new_lines<-new_data %>%
    ungroup() %>%
    select(new_lines) %>%
    unnest(new_lines) %>%
    st_as_sf()

  stream_lines<-stream_lines %>%
    filter(!link_id %in% new_lines$link_id) %>%
    bind_rows(new_lines)

  new_links<-new_data %>%
    ungroup() %>%
    select(new_links) %>%
    filter(map_lgl(new_links,~nrow(.)>0)) %>%
    unnest(new_links) %>%
    st_as_sf()

  stream_links<-stream_links %>%
    filter(!link_id %in% new_links$link_id) %>%
    bind_rows(new_links)

  # Fix flow tracers upstream of sampling points -------------------------------

  new_data<-new_data %>%
    mutate(new_lines=future_map(new_lines, function(lns) {

      out<-lns %>%
        bind_rows(
          stream_lines %>%
            filter(link_id %in% (lns %>%
                                   as_tibble() %>%
                                   filter(link_id==max(link_id)) %>%
                                   select(starts_with("uslink_id")) %>%
                                   unlist() %>%
                                   .[!is.na(.)])
            ) %>%
            mutate(dslink_id1=max(lns$link_id))
        ) %>%
        mutate(across(c(contains('uslink_id'),contains('dslink_id')),~ifelse(.==0,NA_real_,.))) %>%
        mutate(across(c(everything(),-geometry),~ifelse(is.nan(.),NA,.)))

      out[grepl("pour_point_A1B2C3_",out%>% as_tibble() %>% select(any_of(site_id_col)) %>% unlist()),site_id_col]<-NA_character_

      return(out)

    })) %>%
    mutate(new_links=future_map2(data,new_lines, function(pnts,lns) {
      lns_target<-lns %>%
        as_tibble() %>%
        filter(!if_any(any_of(site_id_col),is.na))

      us_IDs<- lns_target%>% # This needs to be a separate loop
        filter(link_id==max(link_id)) %>%
        select(starts_with("uslink_id")) %>%
        unlist() %>%
        .[!is.na(.)]

      replace_us<-stream_links %>% # these will replace existing points
        filter(link_id %in% us_IDs) %>%
        mutate(dslink_id1=max(lns_target$link_id)) %>%
        mutate(across(c(contains('uslink_id'),contains('dslink_id')),~ifelse(.==0,NA_real_,.))) %>%
        mutate(across(c(everything(),-geometry),~ifelse(is.nan(.),NA,.)))

      replace_us[grepl("pour_point_A1B2C3_",replace_us%>% as_tibble() %>% select(any_of(site_id_col)) %>% unlist()),site_id_col]<-NA_character_

      return(replace_us)
    }))

  if (verbose) print("Generating Output")

  # lines replace target lines, and upstream lines

  new_lines<-new_data %>%
    ungroup() %>%
    select(new_lines) %>%
    unnest(new_lines) %>%
    st_as_sf()

  stream_lines<-stream_lines %>%
    filter(!link_id %in% new_lines$link_id) %>%
    bind_rows(new_lines)

  write_sf(stream_lines,file.path(temp_dir,"stream_lines.shp"))

  # points replace along target

  new_points<-new_data %>%
    ungroup() %>%
    select(new_points) %>%
    unnest(new_points) %>%
    st_as_sf()

  stream_points<-stream_points %>%
    filter(!ID %in% new_points$ID) %>%
    bind_rows(new_points)

  write_sf(stream_points,file.path(temp_dir,"stream_points.shp"))

  # links, add and replace

  new_links<-new_data %>%
    ungroup() %>%
    select(new_links) %>%
    filter(map_lgl(new_links,~nrow(.)>0)) %>%
    unnest(new_links) %>%
    st_as_sf()

  stream_links<-stream_links %>%
    filter(!link_id %in% new_links$link_id) %>%
    bind_rows(new_links)

  write_sf(stream_links,file.path(temp_dir,"stream_links.shp"))

  snapped_points<-snapped_points %>%
    rename(orig_link_id=link_id) %>%
    st_join(stream_lines %>% select(link_id))

  # polygons replace at target

  new_poly<-new_data %>%
    ungroup() %>%
    select(new_subbasins) %>%
    unnest(new_subbasins) %>%
    st_as_sf()

  are_collect<-which(st_geometry_type(new_poly)=="GEOMETRYCOLLECTION")
  if (length(are_collect)>0){
    new_poly2<-new_poly[are_collect,] %>%
      st_collection_extract("POLYGON") %>%
      group_by(link_id) %>%
      summarize(across(c(everything(),-geometry),head,1))

    new_poly<-new_poly %>%
      filter(!link_id %in% new_poly2$link_id) %>%
      bind_rows(new_poly2)

  }

  Subbasins_poly<-Subbasins_poly %>%
    filter(!link_id %in% new_poly$link_id) %>%
    bind_rows(new_poly %>% mutate(sbbsn_area=as.numeric(sbbsn_area))) %>%
    mutate(link_id=as.character(link_id)) %>%
    select(everything(),geometry)

  write_sf(Subbasins_poly,file.path(temp_dir,"Subbasins_poly.shp"))

  # Generate Output ---------------------------------------------------------

  dist_list_out<-c(
    list.files(temp_dir,"snapped_points"),
    list.files(temp_dir,"Subbasins_poly"),
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
    output$subbasins<-Subbasins_poly
    output$stream_lines<-stream_lines
    output$links<-stream_links
    output$points<-stream_points
    output$snapped_points<-snapped_points
  }

  file.remove(list.files(temp_dir,full.names = T,recursive = T))

  return(output)

}
