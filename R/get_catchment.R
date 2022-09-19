
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
    mutate(across(everything(),as.character)) %>%
    filter(if_any(any_of(site_id_col), ~.x %in% target_points))

  missing_sites<-target_points[!target_points %in% sites[[1]]]
  if (length(missing_sites)>0) stop(paste0("'target_points' not present in 'points' layer: ",paste0(missing_sites,collapse = ", ")))

  conn<-unz(zip_loc,"us_flowpaths.rds")
  us_flowpaths<-readRDS(gzcon(conn))
  close(conn)

  geo_fn<-function(x,subb=subb,buffer=buffer,tolerance=tolerance) {
    # tdir<-tempfile()
    # tc<-dir.create(tdir)
    # tfile<-file.path(tdir,paste0(basename(tempfile()),".shp"))
    #
    # ploys<-filter(subb,if_any(contains('link_id'), ~.x %in% x$link_id)) %>%
    #   select(geometry) %>%
    #   write_sf(tfile)
    #
    # wbt_clean_vector(
    #   input=tfile,
    #   output=tfile
    # )
    #
    # wbt_dissolve(
    #   input=tfile,
    #   output=tfile
    # )
    #
    # out<-read_sf(tfile)
    #
    # #unlink(tdir,force = T,recursive = T)
    #
    # return(out)

    filter(subb,if_any(contains('link_id'), ~.x %in% x$link_id)) %>%
      select(geometry) %>%
      # st_buffer(units::set_units(buffer,"m"),nQuadSegs = 1) %>%
      # st_snap(x=.,y=., tolerance = units::set_units(tolerance,"m")) %>%
      st_union()
  }

  #browser()

  # s1<-as.list(rep(list(subb),length(sites$link_id)))
  # names(s1)<-sites$link_id
  #
  # future_map(us_flowpaths[sites$link_id],function(x,y) {
  #   sb1<-read_sf(file.path(tdir,"Subbasins_poly.shp"))
  #
  #   filter(sb1,link_id %in% y$link_id)
  # })


  out<-sites %>%
    mutate(us_flowpaths=us_flowpaths[link_id]) %>%
    filter(!map_lgl(us_flowpaths,is.null)) %>%
    filter(map_dbl(us_flowpaths,nrow)>0) %>%
    mutate(subb=future_map(us_flowpaths,function(x) {
      read_sf(file.path(tdir,"Subbasins_poly.shp")) %>%
        filter(link_id %in% x$link_id)
    })) %>%
    mutate(geometry=future_map(subb,function(x) {
      select(x,geometry) %>%
        # st_buffer(units::set_units(buffer,"m"),nQuadSegs = 1) %>%
        # st_snap(x=.,y=., tolerance = units::set_units(tolerance,"m")) %>%
        st_union()
    })) %>%
    # mutate(subb=as.list(rep(list(subb),length(link_id))),
    #        buffer=buffer,
    #        tolerance=tolerance) %>%
    #   mutate(subb=future_map2(subb,
    #                           us_flowpaths,
    #                           function(x,y) {
    #                             filter(x,link_id %in% y$link_id)
    #                             #browser()
    #                           }
    #                           # ~ filter(.x,link_id %in% .y$link_id) %>%
    #                           #   #filter(if_any(contains('link_id'), function(x) ~x %in% .y$link_id)) %>%
  #                           #   select(geometry)
  #   )) %>%
  #   mutate(geometry=future_pmap(
  #     list(
  #       us_flowpaths=us_flowpaths,
  #       subb=subb,
  #       buffer=buffer,
  #       tolerance=tolerance
  #     ),
  #     function(us_flowpaths,
  #              subb,
  #              buffer,
  #              tolerance)
  #       geo_fn(x=us_flowpaths,subb=subb) #,buffer=buffer,tolerance=tolerance
  #   )) %>%
  unnest(geometry) %>%
    select(-us_flowpaths) %>%
    st_as_sf()

  unlink(tdir,force = T,recursive = T)

  return(out)

}
