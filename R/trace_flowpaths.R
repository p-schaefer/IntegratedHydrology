
#' Generate upstream and downstream flow direction lists originating from each subbasin (and sampling point)
#'
#' @param input resulting object from `attrib_streamline()` or `insert_points()`
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and resulting GIS products. If \code{FALSE}, file path only.
#' @param temp_dir character. File path for temporary file storage, If \code{NULL}, `tempfile()` will be used
#' @param verbose logical.
#'
#' @return If \code{return_products = TRUE}, all geospatial analysis products are returned. If \code{return_products = FALSE}, folder path to resulting .zip file.
#' @export

trace_flowpaths<-function(
    input,
    return_products=F,
    temp_dir=NULL,
    verbose=F
){
  options(future.rng.onMisuse="ignore")
  options(scipen = 999)

  if (!is.logical(return_products)) stop("'return_products' must be logical")
  if (!is.logical(verbose)) stop("'verbose' must be logical")

  if (is.null(temp_dir)) temp_dir<-tempfile()
  if (!dir.exists(temp_dir)) dir.create(temp_dir)
  temp_dir<-normalizePath(temp_dir)

  options(dplyr.summarise.inform = FALSE)

  zip_loc<-input$outfile

  final_links<-as_tibble(data.table::fread(cmd=paste("unzip -p ",zip_loc,"stream_links.csv"))) %>%
    mutate(across(any_of(site_id_col),na_if,""))

  #browser()

  fp<-trace_flowpath_fn(input=final_links,verbose=verbose)
  ds_flowpaths<-fp$final_out_ds
  us_flowpaths<-fp$final_out_us

  saveRDS(ds_flowpaths,file.path(temp_dir,"ds_flowpaths.rds"),compress = F)
  saveRDS(us_flowpaths,file.path(temp_dir,"us_flowpaths.rds"),compress = F)

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

  output<-input[!names(input) %in% c("catchment_poly")]

  all_catch<-get_catchment(input = output,
                           target_points = final_links[["link_id"]]) %>%
    select(link_id)

  write_sf(all_catch %>% select(link_id),file.path(temp_dir,"Catchment_poly.shp"))

  zip(out_file,
      list.files(temp_dir,"Catchment_poly",full.names = T),
      flags = '-r9Xjq'
  )

  if (return_products){
    output<-c(
      list(ds_flowpaths=ds_flowpaths,
           us_flowpaths=us_flowpaths,
           catchments=all_catch %>% select(link_id)
           ),
      output
    )
  }
  suppressWarnings(file.remove(list.files(temp_dir,full.names = T)))

  return(output)

}

#' @export
trace_flowpath_fn<-function(
    input,
    verbose=F
) {
  options(scipen = 999)
  options(future.rng.onMisuse="ignore")

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
    distinct()

  unique_link_id<-input_tib %>%
    filter(!is.na(link_id)) %>%
    pull(link_id) %>%
    unique()

  unique_link_id<-suppressWarnings(split(unique_link_id,rep(1:future::nbrOfWorkers(),length.out=length(unique_link_id))))
  unique_link_id<-map(unique_link_id,~setNames(.,.))
  names(unique_link_id)<-NULL
  #unique_link_id<-split(unique_link_id,unique_link_id)
  print("Generating Flowpaths")

  safe_in<-function(x,y){
    as.logical(match(y[!is.na(y)],x[!is.na(x)],nomatch=0))
  }

  id_fn<-function(y,input_tib,p){

    map(y,function(x){
      x_trib<-input_tib %>%
        filter(link_id==x) %>%
        select(link_id,trib_id,USChnLn_Fr,dslink_id1,dstrib_id1,uslink_id1,ustrib_id1)

      out<-x

      # Downstream portion
      if (all(!is.na(x_trib$uslink_id1))) {
        ds_out<-NULL
      } else {

        #print(x)
        repeat {
          ds_id <- input_tib %>%
            filter(link_id %in% out) %>%
            dplyr::select(starts_with("dslink_id")) %>%
            unlist() %>%
            as.character() %>%
            unique() %>%
            .[!is.na(.)]

          if (all(ds_id %in% out)) break

          out<-unique(c(out,ds_id))
        }

        out<-out[!is.na(out)]
        ds_out<-tibble(link_id=out) %>%
          left_join(input_tib %>% dplyr::select(link_id,trib_id,link_lngth,sbbsn_area),by="link_id") %>%
          distinct() %>%
          filter(!is.na(link_id))

      }

      # Upstream portion
      out<-x
      repeat {
        us_id <- input_tib %>%
          filter(link_id %in% out) %>%
          dplyr::select(starts_with("uslink_id")) %>%
          unlist() %>%
          as.character() %>%
          unique() %>%
          .[!is.na(.)]

        if (all(us_id %in% out)) break

        out<-unique(c(out,us_id))
      }

      out<-out[!is.na(out)]

      us_out<-tibble(link_id=out)%>%
        left_join(input_tib %>% dplyr::select(link_id,trib_id,link_lngth,sbbsn_area),by="link_id") %>%
        distinct() %>%
        filter(!is.na(link_id))

      p()
      return(list(us_out=us_out,ds_out=ds_out))
    })
  }

  input_tib_list<-as.list(rep(list(input_tib),length(unique_link_id)))

  with_progress(enable=T,{
    p <- progressor(steps = length(unlist(unique_link_id)))

    unique_link_id<-future_map2(unique_link_id,input_tib_list,~id_fn(y=.x,input_tib=.y,p=p)) %>%
      unlist(recursive=F)
  })

  # Get remaining DS distances by unnesting the list
  final_out_ds<-map(unique_link_id,~.$ds_out)
  final_out_ds<-final_out_ds[!sapply(final_out_ds,is.null)]
  final_out_us<-map(unique_link_id,~.$us_out)

  for (i in final_out_ds){
    nrw<-nrow(i)
    new_entries<-lapply(2:nrw,function(x) i[seq(x,nrw),])
    names(new_entries)<-sapply(new_entries,function(x) head(x$link_id,1))

    keep_entries<-new_entries[!names(new_entries) %in% names(final_out_ds)]

    final_out_ds<-c(final_out_ds,keep_entries)
  }

  final_out_ds<-final_out_ds[order(names(final_out_ds))]
  final_out_ds<-final_out_ds[!is.na(names(final_out_ds))]

  final_out_us<-final_out_us[order(names(final_out_us))]
  final_out_us<-final_out_us[!is.na(names(final_out_us))]

  return(list(final_out_us=final_out_us,final_out_ds=final_out_ds))

}

