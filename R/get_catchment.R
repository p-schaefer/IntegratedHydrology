
#' Title
#'
#' @param input
#' @param site_id_col
#' @param target_sites
#' @param tolerance
#'
#' @return
#' @export
#'
#' @examples
#'
#'
get_catchment<-function(
    input,
    site_id_col=NULL,
    target_sites,
    tolerance=0.000001
) {

  target_sites<-as.character(target_sites)

  if (is.null(site_id_col) || site_id_col=="link_id") site_id_col<-"link_id"

  if (!is.character(site_id_col)) stop("'site_id_col' must be a single character")
  if (length(site_id_col)>1 | length(site_id_col)==0) stop("length 'site_id_col' must be 1")

  zip_loc<-input$outfile

  subb<-read_sf(file.path("/vsizip",zip_loc,"Subbasins_poly.shp"))
  points<-read_sf(file.path("/vsizip",zip_loc,"stream_links.shp"))

  if (!site_id_col %in% names(points)) stop("'site_id_col' must be a variable name in 'points'")

  sites<-points %>%
    as_tibble() %>%
    select(any_of(site_id_col),link_id) %>%
    filter(!if_any(site_id_col,is.na)) %>%
    mutate(across(everything(),as.character))

  sites<-sites[sites[[site_id_col]] %in% target_sites,]

  missing_sites<-target_sites[!target_sites %in% sites[[1]]]
  if (length(missing_sites)>0) stop(paste0("'target_sites' not present in 'points' layer: ",paste0(missing_sites,collapse = ", ")))

  us_flowpaths<-readRDS(gzcon(unz(zip_loc,"us_flowpaths.rds")))

  # if (site_id_col=="link_id") {
  #   out_file<-tibble(
  #     link_id=sites$link_id[match(sites[[1]],target_sites,nomatch = 0)]
  #   )
  # } else {
  #   out_file<-tibble(
  #     uid=target_sites,
  #     link_id=sites$link_id[match(sites[[1]],target_sites,nomatch = 0)]
  #   ) %>%
  #     setNames(c(site_id_col,"link_id"))
  # }


  sites %>%
    mutate(us_flowpaths=us_flowpaths[link_id]) %>%
    mutate(geometry=future_map(us_flowpaths,~filter(subb,link_id %in% .$link_id) %>%
                                 st_buffer(0.001,nQuadSegs = 1)%>%
                                 st_snap(x=.,y=., tolerance = tolerance) %>%
                                 st_union())
    ) %>%
    unnest(geometry) %>%
    select(-us_flowpaths) %>%
    st_as_sf()

}
