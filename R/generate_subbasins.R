
#' Extracts polygon subbasins from 'process_flowdir()'
#'
#' @param input resulting object from `process_flowdir()`
#' @param points character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/points.shp"), or any GIS data object that will be converted to spatial points. Points representing sampling locations.
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param verbose logical.
#'
#' @seealso [whitebox::wbt_subbasins]
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'


generate_subbasins<-function(
    input,
    points,
    return_products=F,
    temp_dir=NULL,
    verbose=F
) {
  # require(sf)
  # require(terra)
  # require(whitebox)
  # require(tidyverse)
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")

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

  site_id_col<-paste0(data.table::fread(cmd=paste("unzip -p ",zip_loc,"site_id_col.csv")))

  unzip(zip_loc,
        c("dem_d8.tif","dem_streams_d8.tif"),
        exdir=temp_dir,
        overwrite=T,
        junkpaths=T)

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

  stream_links<-read_sf(file.path("/vsizip",zip_loc,"stream_links.shp")) %>%
    left_join(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv")),
              by="link_id")

  # Split subbasins at sampling points --------------------------------------

  if (!is.null(points)){

    with_progress(enable=T,{
      print("Splitting Subbasins")

      new_data<-stream_links %>%
        filter(floor(link_id) %in% floor(link_id[!is.na(!!sym(site_id_col))])) %>%
        mutate(link_id_base=floor(link_id)) %>%
        select(link_id_base,link_id,any_of(site_id_col)) %>%
        as_tibble() %>%
        rename(point=geometry) %>%
        group_by(link_id_base) %>%
        mutate(temp_dir=temp_dir,
               target_crs=list(st_crs(subb))) %>%
        nest() %>%
        ungroup() %>%
        left_join(
          subb %>% select(-sbbsn_area) %>% as_tibble() %>% rename(subb_poly=geometry),
          by=c("link_id_base"="link_id")
        )

      p <- progressor(steps = nrow(new_data))

      new_data<-new_data %>%
        mutate(p=rep(list(p),nrow(.))) %>%
        mutate(new_subb=future_pmap(list(data=data,
                                         link_id=link_id_base,
                                         subb_poly=subb_poly,
                                         temp_dir=temp_dir,
                                         p=p),
                                    carrier::crate(function(data,link_id,subb_poly,temp_dir,target_crs,p){
                                      #print(link_id)
                                      `%>%` <- magrittr::`%>%`

                                      if (nrow(data)==1){
                                        catch_poly<-data %>%
                                          dplyr::select(link_id) %>%
                                          dplyr::mutate(geometry=sf::st_geometry(subb_poly)) %>%
                                          sf::st_as_sf(crs = data$target_crs[[1]]) %>%
                                          dplyr::mutate(sbbsn_area=sf::st_area(.)) %>%
                                          dplyr::select(link_id,sbbsn_area, geometry)

                                        p()

                                        return(catch_poly)
                                      }


                                      pnt_file<-file.path(temp_dir,paste0("Tempsite_",link_id,".shp"))


                                      sf::write_sf(data %>% dplyr::select(link_id,point) %>% sf::st_as_sf(),pnt_file)

                                      cr<-terra::vect(subb_poly)

                                      t1<-terra::rast(file.path(temp_dir,"dem_d8.tif")) %>%
                                        terra::crop(y=cr,
                                                    mask=T) %>%
                                        terra::writeRaster(
                                          filename=file.path(temp_dir,paste0("d8_",link_id,".tif")),
                                          overwrite=T
                                        )

                                      whitebox::wbt_unnest_basins(
                                        d8_pntr=paste0("d8_",link_id,".tif"),
                                        pour_pts=paste0("Tempsite_",link_id,".shp"),
                                        output=paste0("Catch_",link_id,".tif")
                                      )


                                      catch_fls<-list.files(temp_dir,pattern=paste0("Catch_",link_id,"_"))
                                      catch_rast<-terra::rast(file.path(temp_dir,catch_fls))
                                      catch_rast[is.na(catch_rast)]<-0
                                      catch_rast[catch_rast>0]<-1
                                      catch_rast<-terra::app(catch_rast,sum)

                                      catch_rast[catch_rast==0]<-NA

                                      catch_poly<-catch_rast %>%
                                        terra::as.polygons() %>%
                                        sf::st_as_sf() %>%
                                        sf::st_join(data %>% sf::st_as_sf()) %>%
                                        dplyr::mutate(sbbsn_area=sf::st_area(.)) %>%
                                        dplyr::select(link_id,sbbsn_area, geometry)

                                      flrm<-unique(c(list.files(temp_dir,pattern=paste0("_",link_id,"."),full.names = T),
                                                     list.files(temp_dir,pattern=paste0("_",link_id,"_"),full.names = T)
                                      ))

                                      file.remove(flrm)

                                      p()

                                      return(catch_poly)

                                    })))

      subb2<-new_data %>%
        select(new_subb) %>%
        unnest(cols = c(new_subb)) %>%
        st_as_sf()

      subb<-subb %>%
        filter(!link_id %in% new_data$link_id_base) %>%
        bind_rows(subb2) %>%
        arrange(link_id)

      write_sf(subb,file.path(temp_dir,"Subbasins_poly.shp"))

    })
  }

  # Prepare Output ----------------------------------------------------------

  final_links<-stream_links %>%
    left_join(subb %>%
                as_tibble() %>%
                select(link_id,sbbsn_area),
              by = c("link_id"))

  write_sf(final_links %>% select(link_id) ,file.path(temp_dir,"stream_links.shp"))
  data.table::fwrite(final_links %>% as_tibble() %>% select(-geometry),file.path(temp_dir,"stream_links.csv"))

  dist_list_out<-list(
    list.files(temp_dir,"Subbasins_poly"),
    list.files(temp_dir,"stream_links")
  )

  dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

  out_file<-zip_loc

  if (verbose) print("Generating Output")

  zip(out_file,
      unlist(dist_list_out),
      flags = '-r9Xjq'
  )

  output<-input[!names(input) %in% c("subbasins","links")]

  if (return_products){
    output<-c(
      list(subbasins=subb,
           links=final_links),
      output
    )
  }
  file.remove(list.files(temp_dir,full.names = T,recursive=T))

  return(output)
}
