
#' Generate upstream and downstream flow direction lists originating from each subbasin (and sampling point)
#'
#' @param input resulting object from `attrib_streamline()` or `insert_points()`
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export
#'
#' @examples

trace_flowpaths<-function(
    input,
    return_products=F,
    temp_dir=NULL,
    verbose=F
){

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

  options(dplyr.summarise.inform = FALSE)

  zip_loc<-input$outfile
  fl<-unzip(list=T,zip_loc)

  final_links<-read_sf(file.path("/vsizip",zip_loc,"stream_links.shp"))

  browser()

  fp<-trace_flowpath_fn(input=final_links,verbose=verbose)
  ds_flowpaths<-trace_ds_flowpath(input=final_links,verbose=verbose)
  us_flowpaths<-trace_us_flowpath(input=final_links,verbose=verbose)

  saveRDS(ds_flowpaths,file.path(temp_dir,"ds_flowpaths.rds"))
  saveRDS(us_flowpaths,file.path(temp_dir,"us_flowpaths.rds"))

  dist_list_out<-list(
    "ds_flowpaths.rds",
    "us_flowpaths.rds"
  )

  dist_list_out<-lapply(dist_list_out,function(x) file.path(temp_dir,x))

  out_file<-zip_loc

  if (verbose) print("Generating Output")

  zip(out_file,
      unlist(dist_list_out),
      flags = '-r9Xjq'
  )

  output<-input

  if (return_products){
    output<-c(
      list(ds_flowpaths=ds_flowpaths,
           us_flowpaths=us_flowpaths),
      output
    )
  }
  file.remove(list.files(temp_dir,full.names = T))

  return(output)

}

#' @export
trace_flowpath_fn<-function(
    input,
    verbose=F
) {

  input_tib<-input %>%
    as_tibble() %>%
    select(link_id,trib_id,link_lngth,sbbsn_area,USChnLn_Fr,
           starts_with("uslink_id"),starts_with("dslink_id"),starts_with("dstrib_id"),starts_with("ustrib_id")) %>%
    filter(!is.na(link_id)) %>%
    mutate(link_id=as.character(link_id)) %>%
    pivot_longer(cols=c(starts_with("uslink_id"),starts_with("ustrib_id"))) %>%
    filter(case_when(
      name %in% c("uslink_id1","ustrib_id1") ~ T,
      is.na(value) ~ F,
      T ~ T
    )) %>%
    mutate(name=ifelse(grepl("trib",name),"ustrib_id1","uslink_id1")) %>%
    pivot_wider(names_from=name,values_from=value,values_fn=list) %>%
    unnest(cols=c(uslink_id1,ustrib_id1)) %>%
    # pivot_longer(cols=starts_with("uslink_id"),values_to="uslink_id1") %>%
    # #filter(name!="uslink_id1" & is.na(uslink_id1)) %>%
    # select(-name) %>%
    # pivot_longer(cols=starts_with("ustrib_id"),values_to="ustrib_id1") %>%
    # #filter(name!="ustrib_id1" & is.na(ustrib_id1)) %>%
    # select(-name) %>%
    # mutate(ustrib_id1=ifelse(is.na(uslink_id1),NA_real_,ustrib_id1)) %>%
    distinct()

  unique_link_id<-input_tib %>%
    filter(!is.na(link_id)) %>%
    #filter(if_all(starts_with("uslink_id"),is.na)) %>%
    pull(link_id) %>%
    unique()

  unique_link_id<-split(unique_link_id,unique_link_id)
  print("Generating Downstream flowpaths")

  safe_in<-function(x,y){
    as.logical(match(y[!is.na(y)],x[!is.na(x)],nomatch=0))
  }

  id_fn<-function(x,input_tib,p){
    #browser()
    out<-x

    x_trib<-input_tib %>%
      filter(link_id==x) %>%
      select(link_id,trib_id,USChnLn_Fr,dslink_id1,dstrib_id1,uslink_id1,ustrib_id1)

    if (all(!is.na(x_trib$uslink_id1))) {
      ds_out<-NULL
    } else {
      out_samp_trib<-input_tib %>% #Everything below x on same trib
        filter(trib_id==x_trib$trib_id[1] & USChnLn_Fr>=x_trib$USChnLn_Fr[1]) %>%
        select(link_id,trib_id,USChnLn_Fr,dslink_id1,dstrib_id1)

      while(any(!safe_in(out_samp_trib$trib_id,out_samp_trib$dstrib_id1))){

        missing_tribs<-out_samp_trib %>%
          filter(!dstrib_id1 %in% trib_id )%>%
          select(dstrib_id1,USChnLn_Fr) %>%
          distinct()%>%
          filter(!is.na(dstrib_id1))

        for (i in 1:nrow(missing_tribs)) {
          out_samp_trib<-out_samp_trib %>%
            bind_rows(
              input_tib %>%
                filter(trib_id==missing_tribs$dstrib_id1[i] & USChnLn_Fr>=missing_tribs$USChnLn_Fr[i]) %>%
                select(link_id,trib_id,USChnLn_Fr,dslink_id1,dstrib_id1)
            )
        }
      }

      ds_out<-out_samp_trib %>%
        select(-starts_with("dslink_id"),-starts_with("dstrib_id")) %>%
        left_join(input_tib %>% dplyr::select(link_id,trib_id,link_lngth,sbbsn_area),by=c("link_id","trib_id")) %>%
        distinct() %>%
        filter(!is.na(link_id))
    }


    out_samp_trib<-input_tib %>%
      filter(trib_id==x_trib$trib_id[1] & USChnLn_Fr<=x_trib$USChnLn_Fr[1]) %>%
      select(link_id,trib_id,USChnLn_Fr,uslink_id1,ustrib_id1)

    while(any(!safe_in(out_samp_trib$trib_id,out_samp_trib$ustrib_id1))){

      missing_tribs<-out_samp_trib %>%
        filter(!ustrib_id1 %in% trib_id) %>%
        select(ustrib_id1) %>%
        distinct() %>%
        filter(!is.na(ustrib_id1))

      for (i in 1:nrow(missing_tribs)) {
        out_samp_trib<-out_samp_trib %>%
          bind_rows(
            input_tib %>%
              filter(trib_id==missing_tribs$ustrib_id1[i]) %>%
              select(link_id,trib_id,USChnLn_Fr,uslink_id1,ustrib_id1)
          )
      }

    }

    us_out<-out_samp_trib %>%
      select(-starts_with("uslink_id"),-starts_with("ustrib_id")) %>%
      left_join(input_tib %>% dplyr::select(link_id,trib_id,link_lngth,sbbsn_area),by=c("link_id","trib_id")) %>%
      distinct() %>%
      filter(!is.na(link_id))

    p()
    return(list(us_out=us_out,ds_out=ds_out))
  }

  input_tib_list<-as.list(rep(list(input_tib),length(unique_link_id)))

  browser()
  with_progress({
    p <- progressor(steps = length(unique_link_id))

    unique_link_id<-future_map2(unique_link_id,input_tib_list,~id_fn(x=.x,input_tib=.y,p=p))
  })

  # Get remaining distances by unnesting the list
  final_out<-unique_link_id

  for (i in final_out){
    i<-i %>%
      group_by(trib_id) %>%
      arrange(USChnLn_Fr)
    nrw<-nrow(i)
    new_entries<-lapply(2:nrw,function(x) i[seq(x,nrw),])
    names(new_entries)<-sapply(new_entries,function(x) head(x$link_id,1))

    keep_entries<-new_entries[!names(new_entries) %in% names(final_out)]

    final_out<-c(final_out,keep_entries)
  }

  final_out<-final_out[order(names(final_out))]
  final_out<-final_out[!is.na(names(final_out))]

  return(final_out)

}

#' @export
trace_ds_flowpath<-function(
    input,
    verbose=F
) {

  input_tib<-input %>%
    as_tibble() %>%
    select(link_id,trib_id,link_lngth,sbbsn_area,USChnLn_Fr,
           starts_with("uslink_id"),starts_with("dslink_id"),starts_with("dstrib_id"),starts_with("ustrib_id")) %>%
    filter(!is.na(link_id)) %>%
    mutate(link_id=as.character(link_id))

  unique_link_id<-input_tib %>%
    filter(!is.na(link_id)) %>%
    filter(if_all(starts_with("uslink_id"),is.na)) %>%
    pull(link_id) %>%
    unique()

  unique_link_id<-split(unique_link_id,unique_link_id)
  print("Generating Downstream flowpaths")

  safe_in<-function(x,y){
    as.logical(match(y[!is.na(y)],x[!is.na(x)],nomatch=0))
  }

  id_fn<-function(x,input_tib,p){
    #browser()
    out<-x

    x_trib<-input_tib %>%
      filter(link_id==x) %>%
      select(link_id,trib_id,USChnLn_Fr)

    out_samp_trib<-input_tib %>% #Everything below x on same trib
      filter(trib_id==x_trib$trib_id & USChnLn_Fr>=x_trib$USChnLn_Fr) %>%
      select(link_id,trib_id,USChnLn_Fr,dslink_id1,dstrib_id1)

    while(any(!safe_in(out_samp_trib$trib_id,out_samp_trib$dstrib_id1))){

      missing_tribs<-out_samp_trib %>%
        filter(!dstrib_id1 %in% trib_id )%>%
        select(dstrib_id1,USChnLn_Fr) %>%
        distinct()

      for (i in 1:nrow(missing_tribs)) {
        out_samp_trib<-out_samp_trib %>%
          bind_rows(
            input_tib %>%
              filter(trib_id==missing_tribs$dstrib_id1[i] & USChnLn_Fr>=missing_tribs$USChnLn_Fr[i]) %>%
              select(link_id,trib_id,USChnLn_Fr,dslink_id1,dstrib_id1)
          )

      }


    }

    out<-out_samp_trib %>%
      select(-starts_with("dslink_id"),-starts_with("dstrib_id")) %>%
      left_join(input_tib %>% dplyr::select(link_id,trib_id,link_lngth,sbbsn_area),by=c("link_id","trib_id")) %>%
      distinct() %>%
      filter(!is.na(link_id))

    p()
    return(out)
  }

  input_tib_list<-as.list(rep(list(input_tib),length(unique_link_id)))

  browser()
  with_progress({
    p <- progressor(steps = length(unique_link_id))

    unique_link_id<-future_map2(unique_link_id,input_tib_list,~id_fn(x=.x,input_tib=.y,p=p))
  })

  # Get remaining distances by unnesting the list
  final_out<-unique_link_id

  for (i in final_out){
    i<-i %>%
      group_by(trib_id) %>%
      arrange(USChnLn_Fr)
    nrw<-nrow(i)
    new_entries<-lapply(2:nrw,function(x) i[seq(x,nrw),])
    names(new_entries)<-sapply(new_entries,function(x) head(x$link_id,1))

    keep_entries<-new_entries[!names(new_entries) %in% names(final_out)]

    final_out<-c(final_out,keep_entries)
  }

  final_out<-final_out[order(names(final_out))]
  final_out<-final_out[!is.na(names(final_out))]

  return(final_out)

}

#' @export
trace_us_flowpath<-function(
    input,
    verbose=F
) {

  input_tib<-input %>%
    as_tibble() %>%
    select(link_id,trib_id,link_lngth,sbbsn_area,USChnLn_Fr,
           starts_with("dslink_id"),starts_with("dstrib_id"),starts_with("uslink_id"),starts_with("ustrib_id")) %>%
    filter(!is.na(link_id)) %>%
    mutate(link_id=as.character(link_id)) %>%
    pivot_longer(cols=starts_with("uslink_id"),values_to="uslink_id1") %>%
    select(-name) %>%
    pivot_longer(cols=starts_with("ustrib_id"),values_to="ustrib_id1") %>%
    select(-name) %>%
    filter(!is.na(uslink_id1),!is.na(ustrib_id1))

  unique_link_id<-input_tib %>%
    filter(!is.na(link_id)) %>%
    # filter(if_all(starts_with("dslink_id"),is.na)) %>% # this won't work, cant get by unnesting
    pull(link_id) %>%
    unique()

  unique_link_id<-split(unique_link_id,unique_link_id)
  print("Generating Upstream flowpaths")

  safe_in<-function(x,y){
    as.logical(match(y[!is.na(y)],x[!is.na(x)],nomatch=0))
  }

  id_fn<-function(x,input_tib,p){

    out<-x

    x_trib<-input_tib %>%
      filter(link_id==x) %>%
      select(link_id,trib_id,USChnLn_Fr)

    out_samp_trib<-input_tib %>% #Everything below x on same trib
      filter(trib_id==x_trib$trib_id & USChnLn_Fr<=x_trib$USChnLn_Fr) %>%
      select(link_id,trib_id,USChnLn_Fr,uslink_id1,ustrib_id1)

    while(any(!safe_in(out_samp_trib$trib_id,out_samp_trib$ustrib_id1))){

      missing_tribs<-out_samp_trib %>%
        filter(!ustrib_id1 %in% trib_id) %>%
        select(ustrib_id1) %>%
        distinct()

      for (i in 1:nrow(missing_tribs)) {
        out_samp_trib<-out_samp_trib %>%
          bind_rows(
            input_tib %>%
              filter(trib_id==missing_tribs$ustrib_id1[i]) %>%
              select(link_id,trib_id,USChnLn_Fr,uslink_id1,ustrib_id1)
          )
      }

    }

    out<-out_samp_trib %>%
      select(-starts_with("uslink_id"),-starts_with("ustrib_id")) %>%
      left_join(input_tib %>% dplyr::select(link_id,trib_id,link_lngth,sbbsn_area),by=c("link_id","trib_id")) %>%
      distinct() %>%
      filter(!is.na(link_id))


    p()
    return(out)
  }

  input_tib_list<-as.list(rep(list(input_tib),length(unique_link_id)))

  browser()

  with_progress({
    p <- progressor(steps = length(unique_link_id))

    unique_link_id<-future_map2(unique_link_id,input_tib_list,~id_fn(x=.x,input_tib=.y,p=p))
  })

  final_out<-unique_link_id

  final_out<-final_out[order(names(final_out))]
  final_out<-final_out[!is.na(names(final_out))]

  return(final_out)
}

