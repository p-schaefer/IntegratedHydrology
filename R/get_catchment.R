
#' Generate upstream catchment areas for target points
#'
#' @param input output of `process_hydrology()` or one containing `generate_subbasins()` and `trace_flowpaths()`
#' @param site_id_col \code{NULL} or character. If Character, must match `site_id_col` from `insert_points()`
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param target_points character. Unique site identifier(s) from `site_id_col`, or `link_id` values corresponding to specific reaches.
#'
#' @return polygon of upstream catchments at `target_points`
#' @export
#'

#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom dplyr collect tbl mutate across na_if filter if_any rename group_by ungroup select
#' @importFrom furrr future_map
#' @importFrom progressr with_progress progressor
#' @importFrom RSQLite SQLite
#' @importFrom sf read_sf st_union st_as_sf
#' @importFrom tidyr nest unnest
#' @importFrom tidyselect any_of everything
#' @importFrom whitebox wbt_options wbt_exe_path

get_catchment<-function(
    input,
    site_id_col=NULL,
    target_points,
    temp_dir=NULL
) {
  options(future.rng.onMisuse="ignore")
  options(scipen = 999)

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  tdir<-file.path(temp_dir,basename(tempfile()))
  dir.create(tdir)

  options(dplyr.summarise.inform = FALSE,future.rng.onMisuse="ignore")

  whitebox::wbt_options(exe_path=whitebox::wbt_exe_path(),verbose=F)

  target_points<-as.character(target_points)

  if (is.null(site_id_col) || site_id_col=="link_id") site_id_col<-"link_id"

  if (!is.character(site_id_col)) stop("'site_id_col' must be a single character")
  if (length(site_id_col)>1 | length(site_id_col)==0) stop("length 'site_id_col' must be 1")

  zip_loc<-input$outfile
  db_loc<-input$db_loc

  subb<-sf::read_sf(file.path("/vsizip",zip_loc,"Subbasins_poly.shp"))

  con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
  points<-dplyr::collect(dplyr::tbl(con,"stream_links")) %>%
    dplyr::mutate(dplyr::across(c(link_id,tidyselect::any_of(site_id_col)),as.character)) %>%
    dplyr::mutate(dplyr::across(tidyselect::any_of(site_id_col),~dplyr::na_if(.,"")))
  DBI::dbDisconnect(con)

  # points<-as_tibble(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv"))) %>%
  #   mutate(across(any_of(site_id_col),na_if,""))

  if (!site_id_col %in% names(points)) stop("'site_id_col' must be a variable name in 'points'")

  sites<-points %>%
    select(tidyselect::any_of(site_id_col),link_id) %>%
    dplyr::filter(!dplyr::if_any(tidyselect::any_of(site_id_col),is.na)) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(),paste0)) %>%
    dplyr::filter(dplyr::if_any(tidyselect::any_of(site_id_col), ~.x %in% target_points))

  missing_sites<-target_points[!target_points %in% sites[[1]]]
  if (length(missing_sites)>0) stop(paste0("'target_points' not present in 'points' layer: ",paste0(missing_sites,collapse = ", ")))


  us_fp_fun<-function(link_id_in,db_loc=db_loc){
    con <- DBI::dbConnect(RSQLite::SQLite(), db_loc)
    out<-dplyr::tbl(con,"us_flowpaths") %>%
      dplyr::filter(pour_point_id %in% link_id_in) %>%
      dplyr::rename(link_id=origin_link_id) %>%
      dplyr::collect() %>%
      # DBI::dbGetQuery(con, paste0("SELECT source_id,link_id FROM us_flowpaths WHERE source_id IN (",paste0(link_id,collapse = ","),")")) %>%
      dplyr::group_by(pour_point_id) %>%
      tidyr::nest() %>%
      dplyr::ungroup()

    out2<-out$data
    names(out2)<-out$pour_point_id

    out2<-out2[link_id_in]

    DBI::dbDisconnect(con)
    return(out2)
  }

#browser()
  progressr::with_progress(enable=T,{
    p <- progressr::progressor(steps = nrow(sites))

    out<-sites %>%
      dplyr::mutate(us_flowpaths=us_fp_fun(link_id,db_loc=db_loc)) %>%
      dplyr::mutate(subb=furrr::future_map(us_flowpaths,function(x) {
        suppressPackageStartupMessages(library(sf))
        subb %>%
          dplyr::filter(link_id %in% x$link_id)
      })) %>%
      dplyr::mutate(geometry=furrr::future_map(subb,function(x) {#
        #browser()
        out<-dplyr::select(x,geometry) %>%
          sf::st_union() #%>%
          # st_as_sf() %>%
          # rename(geometry=x) %>%
          # mutate(link_id=y) %>%
          # select(link_id,geometry)
        p()

        #out<-write_sf(out,file.path(tdir,paste0(y,"_Catchment_poly.shp")))
        return(out)
      })) %>%
      tidyr::unnest(geometry) %>%
      dplyr::select(-us_flowpaths) %>%
      sf::st_as_sf()
  })

  suppressWarnings(unlink(tdir,force = T,recursive = T))

  # out<-list.files(tdir,"_Catchment_poly.shp",full.names=T) %>%
  #   future_map(read_sf) %>%
  #   bind_rows()

  return(out)

}
