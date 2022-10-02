
#' Generate upstream catchment areas for target points
#'
#' @param input output of `process_hydrology()` or one containing `generate_subbasins()` and `trace_flowpaths()`
#' @param site_id_col \code{NULL} or character. If Character, must match `site_id_col` from `insert_points()`
#' @param target_points character. Unique site identifier(s) from `site_id_col`, or `link_id` values corresponding to specific reaches.
#'
#' @return polygon of upstream catchments at `target_points`
#' @export
#'

get_catchment<-function(
    input,
    site_id_col=NULL,
    target_points
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

  points<-as_tibble(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv"))) %>%
    mutate(across(any_of(site_id_col),na_if,""))

  if (!site_id_col %in% names(points)) stop("'site_id_col' must be a variable name in 'points'")

  sites<-points %>%
    select(any_of(site_id_col),link_id) %>%
    filter(!if_any(any_of(site_id_col),is.na)) %>%
    mutate(across(everything(),paste0)) %>%
    filter(if_any(any_of(site_id_col), ~.x %in% target_points))

  missing_sites<-target_points[!target_points %in% sites[[1]]]
  if (length(missing_sites)>0) stop(paste0("'target_points' not present in 'points' layer: ",paste0(missing_sites,collapse = ", ")))


  unzip(zip_loc,files =c("flowpaths_out.db"),exdir=tdir)
  db_fp<-file.path(tdir,"flowpaths_out.db")

  us_fp_fun<-function(link_id,db_fp=db_fp){
    con <- DBI::dbConnect(RSQLite::SQLite(), db_fp)
    out<-DBI::dbGetQuery(con, paste0("SELECT * FROM us_flowpaths WHERE source_id IN (",paste0(link_id,collapse = ","),")")) %>%
      group_by(source_id) %>%
      nest() %>%
      ungroup()

    out2<-out$data
    names(out2)<-out$source_id

    out2<-out2[link_id]

    DBI::dbDisconnect(con)
    return(out2)
  }


  with_progress(enable=T,{
    p <- progressor(steps = nrow(sites))

    out<-sites %>%
      mutate(us_flowpaths=us_fp_fun(link_id,db_fp=db_fp)) %>%
      mutate(subb=future_map(us_flowpaths,function(x) {
        suppressPackageStartupMessages(library(sf))
        subb %>%
          filter(link_id %in% x$link_id)
      })) %>%
      mutate(geometry=future_map(subb,function(x) {
        out<-select(x,geometry) %>%
          st_union()
        p()
        return(out)
      })) %>%
      unnest(geometry) %>%
      select(-us_flowpaths) %>%
      st_as_sf()
  })

  suppressWarnings(unlink(tdir,force = T,recursive = T))

  return(out)

}
