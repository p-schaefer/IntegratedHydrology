
#' Title
#'
#' @param input
#' @param points
#' @param snap_distance
#' @param return_products
#' @param temp_dir
#' @param verbose
#'
#' @return
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
  require(sf)
  require(terra)
  require(whitebox)
  require(tidyverse)

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

  points<-hydroweight::process_input(points,input_name="points")
  points<-st_as_sf(points)

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

  if (verbose) print("Extracting Data")
  unzip(zip_loc,
        c("dem_streams_d8.tif","dem_d8.tif",fl$Name[grepl("stream_links|stream_points|Subbasins_poly|stream_lines",fl$Name)]),
        exdir=temp_dir,
        overwrite=T,
        junkpaths=T)

  stream_links<-read_sf(file.path(temp_dir,"stream_links.shp"))
  stream_points<-read_sf(file.path(temp_dir,"stream_points.shp"))
  Subbasins_poly<-read_sf(file.path(temp_dir,"Subbasins_poly.shp"))
  stream_lines<-read_sf(file.path(temp_dir,"stream_lines.shp"))

  write_sf(st_as_sf(points),file.path(temp_dir,"sample_points.shp"))

  if (verbose) print("Snapping Points")
  wbt_jenson_snap_pour_points(
    pour_pts = "sample_points.shp",
    streams = "dem_streams_d8.tif",
    snap_dist=snap_distance,
    output="snapped_sample_points.shp"
  )

  snapped_points<-read_sf(file.path(temp_dir,"snapped_sample_points.shp")) %>%
    st_join(stream_lines %>% select(link_id))

  write_sf(snapped_points,file.path(temp_dir,"snapped_points.shp"))

  # Add a catch for points that overlap links exactly

  if (verbose) print("Inserting Points")
  new_data<-snapped_points %>%
    select(any_of(site_id_col),link_id) %>%
    group_by(link_id) %>%
    nest() %>%
    mutate(new_subbasins=map2(link_id,data,function(l_id,pnt){

      # Split subbasins by sampling points ----------------------------------

      pnt_file<-file.path(temp_dir,paste0("Tempsite_",l_id,".shp"))

      pnt<-pnt %>%
        mutate(across(any_of(site_id_col),as.character)) %>%
        bind_rows(
          stream_links %>%
            filter(link_id==l_id) %>%
            select(geometry) %>%
            mutate(UID=paste0("pour_point_",l_id)) %>%
            setNames(c("geometry",site_id_col))
        )

      write_sf(pnt,pnt_file)

      rast(file.path(temp_dir,"dem_d8.tif")) %>%
        crop(y=Subbasins_poly %>% filter(link_id == l_id) %>% vect(),
             mask=T,overwrite=T) %>%
        writeRaster(
          file.path(temp_dir,paste0("d8_",l_id,".tif")),
          overwrite=T
        )

      wbt_unnest_basins(
        d8_pntr=paste0("d8_",l_id,".tif"),
        pour_pts=pnt_file,
        output=paste0("Catch_",l_id,".tif")
      )

      catch_fls<-list.files(temp_dir,pattern=paste0("Catch_",l_id))
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

      return(catch_poly)
    })) %>%
    mutate(new_points=map(new_subbasins, function(x) {
      # Split points by sampling points ----------------------------------

      x %>%
        select(-sbbsn_area) %>%
        rename(link_id_new=link_id) %>%
        st_buffer(-units::as_units(0.5,"m"),nQuadSegs = 1) %>%
        split(.$link_id_new) %>%
        map(st_intersection,stream_points) %>%
        map(select,-link_id ) %>%
        map(rename, link_id=link_id_new) %>%
        bind_rows()
    }
    )) %>%
    mutate(new_lines=map(new_subbasins, function(x) {
      # Split lines by sampling points ----------------------------------

      out<-x %>%
        select(-sbbsn_area) %>%
        rename(link_id_new=link_id) %>%
        st_buffer(-units::as_units(0.5,"m"),nQuadSegs = 1) %>%
        split(.$link_id_new) %>%
        map(st_intersection,stream_lines %>% filter(link_id == min(x$link_id))) %>%
        map(select,-link_id ) %>%
        bind_rows() %>%
        mutate(dslink_id1=ifelse(link_id_new==min(link_id_new),dslink_id1,lag(link_id_new))) %>%
        mutate(uslink_id1=ifelse(link_id_new==max(link_id_new),uslink_id1,lead(link_id_new))) %>%
        mutate(across(c(starts_with("uslink_id"),-uslink_id1),~ifelse(link_id_new==max(link_id_new),.,NA))) %>%
        rename(link_id=link_id_new)

      out<-out %>%
        bind_rows(
          stream_lines %>%
            filter(link_id %in% (out %>%
                                   as_tibble() %>%
                                   filter(link_id==max(link_id)) %>%
                                   select(starts_with("uslink_id")) %>%
                                   unlist() %>%
                                   .[!is.na(.)])
            ) %>%
            mutate(dslink_id1=max(x$link_id))
        )

      return(out)
    }
    )) %>%
    mutate(new_links=map2(data,new_lines, function(pnts,lns){
      # Split links by sampling points ----------------------------------

      lns_target<-lns %>%
        as_tibble() %>%
        filter(!if_any(any_of(site_id_col),is.na))

      us_IDs<- lns_target%>%
        filter(link_id==max(link_id)) %>%
        select(starts_with("uslink_id")) %>%
        unlist() %>%
        .[!is.na(.)]

      replace_us<-stream_links %>% # these will replace existing points
        filter(link_id %in% us_IDs) %>%
        mutate(dslink_id1=max(lns_target$link_id))

      replace_target<-stream_links %>%
        filter(link_id %in% min(lns_target$link_id)) %>%
        mutate(uslink_id1=max(lns_target$link_id)) %>%
        mutate(across(c(starts_with("uslink_id"),-uslink_id1),~NA))

      add_links<-pnts %>%
        mutate(across(any_of(site_id_col),as.character)) %>%
        bind_rows(
          stream_links %>%
            filter(link_id==min(lns_target$link_id)) %>%
            select(geometry) %>%
            mutate(UID=paste0("pour_point_",min(lns_target$link_id))) %>%
            setNames(c("geometry",site_id_col))
        ) %>%
        left_join(lns %>% as_tibble() %>% select(-geometry)) %>%
        arrange(link_id) %>%
        mutate(link_id_new=link_id) %>%
        mutate(link_id=min(link_id)) %>%
        select(link_id,link_id_new,any_of(site_id_col)) %>%
        left_join(stream_links %>% as_tibble() %>% select(-geometry)) %>%
        mutate(dslink_id1=ifelse(link_id_new==min(link_id_new),dslink_id1,lag(link_id_new))) %>%
        mutate(uslink_id1=ifelse(link_id_new==max(link_id_new),uslink_id1,lead(link_id_new))) %>%
        mutate(across(c(starts_with("uslink_id"),-uslink_id1),~ifelse(link_id_new==max(link_id_new),.,NA))) %>%
        filter(link_id_new %in% max(link_id_new)) %>%
        select(-link_id) %>%
        rename(link_id=link_id_new)

      return(bind_rows(replace_us,replace_target,add_links))
    }

    )
    )

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
    unnest(new_links) %>%
    st_as_sf()

  stream_links<-stream_links %>%
    filter(!link_id %in% new_links$link_id) %>%
    bind_rows(new_links)

  write_sf(stream_links,file.path(temp_dir,"stream_links.shp"))

  # polygons replace at target

  new_poly<-new_data %>%
    ungroup() %>%
    select(new_subbasins) %>%
    unnest(new_subbasins) %>%
    st_as_sf()

  Subbasins_poly<-Subbasins_poly %>%
    filter(!link_id %in% new_poly$link_id) %>%
    bind_rows(new_poly %>% mutate(sbbsn_area=as.numeric(sbbsn_area)))

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

  file.remove(list.files(temp_dir,full.names = T))

  return(output)

}
