
#' Generate upstream catchment areas for target points
#'
#' @param input output of `process_hydrology()` or one containing `generate_subbasins()` and `trace_flowpaths()`
#' @param site_id_col \code{NULL} or character. If Character, must match `site_id_col` from `insert_points()`
#' @param target_points character. Unique site identifier(s) from `site_id_col`, or `link_id` values corresponding to specific reaches.
#' @param tolerance numeric. Tolerance values used for \code{sf::st_snap} in meters
#' @param buffer numeric. Tolerance values used for \code{sf::st_buffer} in meters
#'
#' @return polygon of upstream catchments at `target_points`
#' @export
#'

get_catchment<-function(
    input,
    site_id_col=NULL,
    target_points,
    tolerance=0.000001,
    buffer=0.001
) {
  options(future.rng.onMisuse="ignore")
  options(scipen = 999)

  tdir<-tempfile()
  dir.create(tdir)

  options(dplyr.summarise.inform = FALSE,future.rng.onMisuse="ignore")

  wbt_options(exe_path=wbt_exe_path(),verbose=F)

  target_points<-as.character(target_points)

  if (is.null(site_id_col) || site_id_col=="link_id") site_id_col<-"link_id"

  if (!is.character(site_id_col)) stop("'site_id_col' must be a single character")
  if (length(site_id_col)>1 | length(site_id_col)==0) stop("length 'site_id_col' must be 1")

  zip_loc<-input$outfile

  subb<-read_sf(file.path("/vsizip",zip_loc,"Subbasins_poly.shp"))
  unzip(zip_loc,files =c("Subbasins_poly.shp","Subbasins_poly.shx","Subbasins_poly.prj","Subbasins_poly.dbf"),exdir=tdir)
  points<-read_sf(file.path("/vsizip",zip_loc,"stream_links.shp"))

  if (!site_id_col %in% names(points)) stop("'site_id_col' must be a variable name in 'points'")

  sites<-points %>%
    as_tibble() %>%
    select(any_of(site_id_col),link_id) %>%
    filter(!if_any(any_of(site_id_col),is.na)) %>%
    mutate(across(everything(),paste0)) %>%
    filter(if_any(any_of(site_id_col), ~.x %in% target_points))

  missing_sites<-target_points[!target_points %in% sites[[1]]]
  if (length(missing_sites)>0) stop(paste0("'target_points' not present in 'points' layer: ",paste0(missing_sites,collapse = ", ")))

  conn<-unz(zip_loc,"us_flowpaths.rds")
  us_flowpaths<-readRDS(gzcon(conn))
  close(conn)

  # geo_fn<-function(x,subb=subb,buffer=buffer,tolerance=tolerance) {
  #   filter(subb,if_any(contains('link_id'), ~.x %in% x$link_id)) %>%
  #     select(geometry) %>%
  #     # st_buffer(units::set_units(buffer,"m"),nQuadSegs = 1) %>%
  #     # st_snap(x=.,y=., tolerance = units::set_units(tolerance,"m")) %>%
  #     st_union()
  # }

  with_progress(enable=T,{
    p <- progressor(steps = nrow(sites))

    out<-sites %>%
      mutate(us_flowpaths=us_flowpaths[link_id]) %>%
      filter(!map_lgl(us_flowpaths,is.null)) %>%
      filter(map_dbl(us_flowpaths,nrow)>0) %>%
      #mutate(sub=rep(list(subb),nrow(.))) %>%
      mutate(subb=future_map(us_flowpaths,function(x) {
        suppressPackageStartupMessages(library(sf))
        subb %>%
        #y %>%
          #read_sf(file.path("/vsizip",zip_loc,"Subbasins_poly.shp")) %>%
          filter(link_id %in% x$link_id)
      })) %>%
      mutate(geometry=future_map(subb,function(x) {
        out<-select(x,geometry) %>%
          # st_buffer(units::set_units(buffer,"m"),nQuadSegs = 1) %>%
          # st_snap(x=.,y=., tolerance = units::set_units(tolerance,"m")) %>%
          st_union()
        p()
        return(out)
      })) %>%
      unnest(geometry) %>%
      select(-us_flowpaths) %>%
      st_as_sf()
  })

  unlink(tdir,force = T,recursive = T)

  return(out)

}
