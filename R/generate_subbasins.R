
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


#' @importFrom carrier crate
#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom dplyr collect tbl mutate across na_if left_join filter select rename group_by ungroup bind_rows arrange copy_to
#' @importFrom furrr future_pmap furrr_options
#' @importFrom rlang sym
#' @importFrom RSQLite SQLite
#' @importFrom sf st_as_sf st_area write_sf read_sf st_crs st_geometry st_join
#' @importFrom terra terraOptions writeRaster rast as.polygons vect crop app
#' @importFrom tibble as_tibble
#' @importFrom tidyr nest unnest
#' @importFrom tidyselect any_of
#' @importFrom whitebox wbt_options wbt_exe_path wbt_subbasins wbt_unnest_basins

generate_subbasins<-function(
    input,
    points,
    return_products=F,
    temp_dir=NULL,
    verbose=F
) {
  if (!inherits(input,"ihydro")) stop("'input' must be of class('ihydro')")
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  whitebox::wbt_options(exe_path=whitebox::wbt_exe_path(),
              verbose=verbose>2,
              wd=temp_dir)

  terra::terraOptions(verbose = verbose>3,
                      tempdir = temp_dir
  )

  db_fp<-input$outfile

  #fl<-unzip(list=T,zip_loc)

  con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)

  site_id_col<-dplyr::collect(dplyr::tbl(con,"site_id_col"))$site_id_col

  for (i in c("dem_d8","dem_streams_d8")){
    terra::writeRaster(
      terra::rast(db_fp,lyrs=i),
      file.path(temp_dir,paste0(i,".tif")),overwrite=T
    )
  }

  # site_id_col<-paste0(data.table::fread(cmd=paste("unzip -p ",zip_loc,"site_id_col.csv")))
  #
  # unzip(zip_loc,
  #       c("dem_d8.tif","dem_streams_d8.tif"),
  #       exdir=temp_dir,
  #       overwrite=T,
  #       junkpaths=T)

  # Generate subbasin polygons ----------------------------------------------
  if (verbose) message("Generating subbasins")
  whitebox::wbt_subbasins(
    d8_pntr="dem_d8.tif",
    streams="dem_streams_d8.tif",
    output="Subbasins.tif"
  )

  # Only keep data with an associated stream line
  # SOmetimes you can get a sink at the edge of a DEM with no assicaited stream

  #lns<-sf::read_sf(db_fp,"stream_lines")


  if (verbose) message("Converting subbasins to polygons")
  subb<-terra::rast(file.path(temp_dir,"Subbasins.tif"))
  subb<-terra::as.polygons(subb,dissolve = TRUE)
  subb<-sf::st_as_sf(subb) %>%
    dplyr::mutate(sbbsn_area=sf::st_area(.))
  names(subb)[1]<-"link_id"

  sf::write_sf(subb,# %>% dplyr::filter(link_id %in% lns$link_id)
               file.path(temp_dir,"Subbasins_poly.shp")
               )

  # Split subbasins at sampling points --------------------------------------
  #db_fp<-input$db_loc
  #con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)

  if (!is.null(points)){
    #browser()

    stream_links<-dplyr::collect(dplyr::tbl(con,"stream_links_attr")) %>%
      dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(site_id_col),~dplyr::na_if(.,"")))

    stream_links<-sf::read_sf(db_fp,layer="stream_links")%>%
      dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>% #stream_links.shp
      dplyr::left_join(stream_links,
                       by = c("link_id")) %>%
      dplyr::mutate(link_id=as.numeric(link_id))

    message("Splitting Subbasins")

    new_data<-stream_links %>%
      dplyr::filter(floor(link_id) %in% floor(link_id[!is.na(!!rlang::sym(site_id_col))])) %>%
      dplyr::mutate(link_id_base=floor(link_id)) %>%
      dplyr::select(link_id_base,link_id,any_of(site_id_col)) %>%
      tibble::as_tibble() %>%
      dplyr::rename(point=geom) %>%
      dplyr::group_by(link_id_base) %>%
      dplyr::mutate(temp_dir=temp_dir,
             target_crs=list(sf::st_crs(subb))) %>%
      tidyr::nest() %>%
      dplyr::ungroup() %>%
      dplyr::left_join(
        subb %>% dplyr::select(-sbbsn_area) %>% tibble::as_tibble() %>% dplyr::rename(subb_poly=geometry),
        by=c("link_id_base"="link_id")
      )

    #browser()
    p <- progressor(steps = nrow(new_data))
    with_progress(enable=T,{

      new_data<-new_data %>%
        dplyr::mutate(p=rep(list(p),nrow(.))) %>%
        dplyr::mutate(new_subb=furrr::future_pmap(list(data=data, #
                                         link_id=link_id_base,
                                         subb_poly=subb_poly,
                                         temp_dir=temp_dir,
                                         p=p),
                                    .options = furrr::furrr_options(globals = FALSE),
                                    carrier::crate(function(data,link_id,subb_poly,temp_dir,target_crs,p){
                                      #browser()
                                      options(scipen = 999)
                                      `%>%` <- magrittr::`%>%`

                                      if (nrow(data)==1){
                                        catch_poly<-data %>%
                                          dplyr::select(link_id) %>%
                                          dplyr::mutate(geom=sf::st_geometry(subb_poly)) %>%
                                          sf::st_as_sf(crs = data$target_crs[[1]]) %>%
                                          dplyr::mutate(sbbsn_area=sf::st_area(.)) %>%
                                          dplyr::select(link_id,sbbsn_area, geom)

                                        p()

                                        return(catch_poly)
                                      }


                                      pnt_file<-file.path(temp_dir,paste0("Tempsite_",link_id,".shp"))


                                      sf::write_sf(data %>% dplyr::select(link_id,point) %>% sf::st_as_sf(),pnt_file)

                                      cr<-terra::vect(subb_poly)

                                      t1<-terra::rast(file.path(temp_dir,"dem_d8.tif")) %>%
                                        terra::crop(y=cr,
                                                    mask=T,
                                                    snap="in",
                                                    touches=F,
                                                    filename=file.path(temp_dir,paste0("d8_int_",link_id,".tif")),
                                                    overwrite=T) %>%
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

                                      p()

                                      return(catch_poly)

                                    })))
    })

    subb2<-new_data %>%
      select(new_subb) %>%
      tidyr::unnest(cols = c(new_subb)) %>%
      sf::st_as_sf() %>%
      dplyr::mutate(sbbsn_area=sf::st_area(.))

    subb<-subb %>%
      dplyr::filter(!link_id %in% new_data$link_id_base) %>%
      dplyr::bind_rows(subb2) %>%
      dplyr::arrange(link_id)

    sf::write_sf(subb,file.path(temp_dir,"Subbasins_poly.shp"))

    all_stream_links<-dplyr::collect(dplyr::tbl(con,"stream_links_attr")) %>%
      dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(site_id_col),~dplyr::na_if(.,"")))

    all_stream_links<-sf::read_sf(db_fp,layer="stream_links") %>%
      dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
      dplyr::left_join(all_stream_links, by="link_id")

    final_links<-all_stream_links %>%
      dplyr::left_join(subb %>%
                  tibble::as_tibble() %>%
                  dplyr::select(link_id,sbbsn_area) %>%
                  dplyr::mutate(dplyr::across(link_id,as.character)),
                by = c("link_id"))

  } else {
    stream_links<-dplyr::collect(dplyr::tbl(con,"stream_links_attr")) %>%
      dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
      dplyr::mutate(dplyr::across(tidyselect::any_of(site_id_col),~dplyr::na_if(.,"")))

    stream_links<-sf::read_sf(db_fp,layer="stream_links") %>%
      dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
      dplyr::left_join(stream_links,
                by="link_id")

    final_links<-stream_links %>%
      dplyr::left_join(subb %>%
                  tibble::as_tibble() %>%
                  dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
                  dplyr::select(link_id,sbbsn_area),
                by = c("link_id"))
  }

  # Prepare Output ----------------------------------------------------------

  sf::write_sf(final_links %>% dplyr::select(link_id) ,file.path(temp_dir,"stream_links.shp"))

  #db_fp<-input$db_loc

  ot<-final_links %>%
    tibble::as_tibble() %>%
    dplyr::select(-tidyselect::any_of("geom"),-tidyselect::any_of("geometry")) %>%
    dplyr::copy_to(df=.,
            con,
            "stream_links_attr",
            overwrite =T,
            temporary =F,
            indexes=c("link_id","trib_id"),
            analyze=T,
            in_transaction=T)

  # Write cell numbers to database ------------------------------------------

  # if (verbose) message("Generating cell number tables")
  #
  # decimalplaces <- function(x) {
  #   if (abs(x - round(x)) > .Machine$double.eps^0.5) {
  #     nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  #   } else {
  #     return(0)
  #   }
  # }
  #
  # n_dec<-max(sapply(subb$link_id,decimalplaces))
  # all_subb_rast<-terra::rasterize(subb,
  #                                 terra::rast(file.path("/vsizip",zip_loc,"dem_final.tif")),
  #                                 field="link_id")
  #
  # all_subb_out<-data.table::data.table(
  #   subb_link_id=terra::values(all_subb_rast),
  #   cell_number=1:terra::ncell(all_subb_rast)
  # ) %>%
  #   setNames(c("subb_link_id","cell_number")) %>%
  #   mutate(row=rowFromCell(all_subb_rast,cell_number)) %>%
  #   mutate(col=colFromCell(all_subb_rast,cell_number)) %>%
  #   mutate(subb_link_id=round(subb_link_id,n_dec)) %>%
  #   mutate(subb_link_id=as.character(subb_link_id)) %>%
  #   copy_to(df=.,
  #           con,
  #           "link_id_cellstats",
  #           overwrite =T,
  #           temporary =F,
  #           indexes=c("subb_link_id","cell_number"),
  #           analyze=T,
  #           in_transaction=T)


  # DBI::dbExecute(con,"CREATE INDEX inx_stream_links ON stream_links (link_id,trib_id)")


  #data.table::fwrite(final_links %>% as_tibble() %>% select(-geometry),file.path(temp_dir,"stream_links.csv"))

  # dist_list_out<-list(
  #   list.files(temp_dir,"Subbasins_poly"),
  #   list.files(temp_dir,"stream_links")
  # )
  #
  # dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))
  #
  # out_file<-zip_loc
  #
  # if (verbose) message("Generating Output")
  #
  # zip(out_file,
  #     unlist(dist_list_out),
  #     flags = '-r9Xjq'
  # )

  for (i in c("Subbasins_poly",
              "stream_links")){
    if (file.exists(file.path(temp_dir,paste0(i,".shp")))){
      t1<-sf::write_sf(
        sf::read_sf(file.path(temp_dir,paste0(i,".shp"))), #%>%
          #dplyr::filter(dplyr::if_any(.cols=tidyselect::any_of("link_id"),~.x %in% lns$link_id)),
        db_fp,
        layer=i,
        append = T,
        delete_layer=T#,
        #layer_options = "OVERWRITE=true"
      )
    }
  }

  output<-input[!names(input) %in% c("subbasins","links")]

  if (return_products){
    output<-c(
      list(subbasins=subb,
           links=final_links
      ),
      output
    )
  }

  DBI::dbDisconnect(con)

  suppressWarnings(file.remove(list.files(temp_dir,full.names = T,recursive=T)))

  class(output)<-"ihydro"
  return(output)
}
