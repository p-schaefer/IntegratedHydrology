
#' @export
trace_ds_flowpath<-function(
    stream_lines
) {
  require(tidyverse)

  stream_lines<-as_tibble(stream_lines)
  stream_lines$link_id<-as.character(stream_lines$link_id)

  unique_link_id<-unique(stream_lines$link_id)
  unique_link_id<-split(unique_link_id,unique_link_id)

  pb <- txtProgressBar(min = 0,
                       max = length(unique_link_id),
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,
                       char = "=")

  unique_link_id<-lapply(unique_link_id,function(x){
    out<-x
    repeat {
      ds_id <- stream_lines %>%
        filter(link_id %in% out) %>%
        select(starts_with("dslink_id")) %>%
        unlist() %>%
        as.character() %>%
        unique()

      if (all(ds_id %in% out)) break

      out<-unique(c(out,ds_id))
    }

    out<-out[!is.na(out)]
    out<-tibble(link_id=out,
                link_lngth=stream_lines$link_lngth[match(out,stream_lines$link_id)])

    setTxtProgressBar(pb, which(names(unique_link_id)==out$link_id[[1]]))

    return(out)
  })

  return(unique_link_id)

}

#' @export
trace_us_flowpath<-function(
    stream_lines
) {
  require(tidyverse)

  stream_lines<-as_tibble(stream_lines)
  stream_lines$link_id<-as.character(stream_lines$link_id)

  unique_link_id<-unique(stream_lines$link_id)
  unique_link_id<-split(unique_link_id,unique_link_id)

  pb <- txtProgressBar(min = 0,
                       max = length(unique_link_id),
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,
                       char = "=")

  unique_link_id<-lapply(unique_link_id,function(x){
    out<-x
    repeat {
      ds_id <- stream_lines %>%
        filter(link_id %in% out) %>%
        select(starts_with("uslink_id")) %>%
        unlist() %>%
        as.character() %>%
        unique()

      if (all(ds_id %in% out)) break

      out<-unique(c(out,ds_id))
    }

    out<-out[!is.na(out)]
    out<-tibble(link_id=out,
                link_lngth=stream_lines$link_lngth[match(out,stream_lines$link_id)],
                sbbsn_area=stream_lines$sbbsn_area[match(out,stream_lines$link_id)])

    setTxtProgressBar(pb, which(names(unique_link_id)==out$link_id[[1]]))

    return(out)
  })

  return(unique_link_id)
}

#' @export
pairwise_dist<-function(
    ds_flowpaths=NULL,
    us_flowpaths=NULL,
    stream_lines=NULL,
    subbasins
) {
  require(tidyverse)

  if (is.null(pairwise_dist) & is.null(ds_flowpaths)) stop("Either 'ds_flowpaths' or 'stream_lines' must be provided")
  if (is.null(pairwise_dist) & is.null(us_flowpaths)) stop("Either 'us_flowpaths' or 'stream_lines' must be provided")

  if (is.null(ds_flowpaths)) ds_flowpaths<-trace_ds_flowpath(stream_lines)
  if (is.null(us_flowpaths)) us_flowpaths<-trace_us_flowpath(stream_lines)
  us_catchment_areas<-lapply(us_flowpaths,function(x) sum(x$sbbsn_area,na.rm=T))

  if (any(duplicated(names(ds_flowpaths)))) stop("'ds_flowpaths' cannot contain duplicate names")
  if (any(duplicated(names(us_flowpaths)))) stop("'us_flowpaths' cannot contain duplicate names")
  if (any(!names(ds_flowpaths) %in% subbasins$link_id)) stop("'ds_flowpaths' contains link_ds not present in 'subbasins'")
  if (any(!names(us_flowpaths) %in% subbasins$link_id)) stop("'us_flowpaths' contains link_ds not present in 'subbasins'")
  if (any(!names(ds_flowpaths) %in% names(us_flowpaths))) stop("'ds_flowpaths' contains link_ds not present in 'us_flowpaths'")
  if (any(!names(us_flowpaths) %in% names(ds_flowpaths))) stop("'us_flowpaths' contains link_ds not present in 'ds_flowpaths'")

  out_tbl<-tibble(
    origin=rep(names(ds_flowpaths),each=length(ds_flowpaths)),
    destination=rep(names(ds_flowpaths),length.out=length(ds_flowpaths)*length(ds_flowpaths))) %>%
    mutate(
      directed_path_length=ds_flowpaths[origin]
    ) %>%
    mutate(directed_path_length=map2_dbl(directed_path_length,destination,
                                         ~mutate(.x,is_target=link_id==.y) %>%
                                           #.[-c(1),] %>% #The first row is always the origin
                                           filter(!is.na(link_lngth)) %>%
                                           mutate(link_lngth=cumsum(link_lngth)) %>%
                                           filter(is_target) %>%
                                           pull(link_lngth) %>%
                                           ifelse(length(.)==0,0,.)
    )) %>%
    mutate(
      origin_catchment=unlist(us_catchment_areas[origin]),
      destination_catchment=unlist(us_catchment_areas[destination])
    ) %>%
    mutate(prop_shared_catchment=case_when(
      directed_path_length>0 ~ origin_catchment/destination_catchment,
      T ~ 0
    )) %>%
    select(-origin_catchment,-destination_catchment)

  out_prop<-out_tbl %>%
    select(-directed_path_length) %>%
    mutate(origin=paste0("prop_link_id_",origin)) %>%
    rename(link_id=destination) %>%
    pivot_wider(names_from=origin,values_from=prop_shared_catchment)

  out_dist<-out_tbl %>%
    select(-prop_shared_catchment) %>%
    mutate(origin=paste0("dist_link_id_",origin)) %>%
    rename(link_id=destination) %>%
    pivot_wider(names_from=origin,values_from=directed_path_length)

  return(list(
    long_pwise=out_tbl,
    wide_prop_shared_catchment=out_prop,
    wide_directed_path_length=out_dist
  ))

}

