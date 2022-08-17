IntegratedHydrology: Integrated hydrology tools for environmental
science
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Contents

-   [1.0 Introduction](#10-introduction)
-   [2.0 System setup and
    installation](#20-system-setup-and-installation)
-   [3.0 Prepare DEM and Sampling Points for for
    analysis](#30-inverse-distance-weighted-rasters-using-hydroweight)
-   [3.1 Generate toy terrain dataset and sampling
    points](#31-generate-toy-terrain-dataset)
-   [3.2 Process DEM with `process_hydrology()`](#32-generate-targets)

## 1.0 Introduction

Add Introduction!

[Back to top](#contents)

## 2.0 System setup and installation

*WhiteboxTools* and *whitebox* are required for
***IntegratedHydrology***. See
[whiteboxR](https://github.com/giswqs/whiteboxR) or below for
installation.

``` r
## Follow instructions for whitebox installation accordingly
## devtools::install_github("giswqs/whiteboxR") # For development version
## whitebox is now available on CRAN
#install.packages("whitebox")

library(whitebox)

if (F){
  install_whitebox()
  # Possible warning message:
  # ------------------------------------------------------------------------
  # Could not find WhiteboxTools!
  # ------------------------------------------------------------------------
  #
  # Your next step is to download and install the WhiteboxTools binary:
  #     > whitebox::install_whitebox()
  #
  # If you have WhiteboxTools installed already run `wbt_init(exe_path=...)`':
  #    > wbt_init(exe_path='/home/user/path/to/whitebox_tools')
  #
  # For whitebox package documentation, ask for help:
  #    > ??whitebox
  #
  # For more information visit https://giswqs.github.io/whiteboxR/
  #
  # ------------------------------------------------------------------------
}
```

[Back to top](#contents)

## 3.0 Prepare DEM and Sampling Points for for analysis

### 3.1 Generate toy terrain dataset and sampling points

We begin by bringing in our toy digital elevation model and using it to
generate terrain products. Burning stream networks into the DEM is not
yet available within IntegratedHydrology, as there are several caveats
associated with this process. See
[here](https://proceedings.esri.com/library/userconf/proc99/proceed/papers/pap802/p802.htm#Trois)

``` r
## Load libraries
library(dplyr)
library(hydroweight)
library(terra)
library(sf)
library(viridis)
library(whitebox)
library(mapview)
library(IntegratedHydrology)

## Generate save_dir as a temporary directory
save_dir <- tempdir()

## Import toy_dem from whitebox package
toy_file<-sample_dem_data()
toy_file <- system.file("extdata", "DEM.tif", package = "whitebox")
toy_dem <- rast(raster::raster(x = toy_file)) # reading the file from terra directly sometimes crashes R for some reason
crs(toy_dem) <- "epsg:3161"
```

[Back to top](#contents)

### 3.2 Process DEM with `process_hydrology()`

``` r
hydro_out<-process_hydrology(
  dem=toy_dem,
  threshold=1000L,
  return_products=T,
  save_dir=save_dir,
  temp_dir=NULL, # location to store temporary files, may require significant space if processing a large DEM 
  verbose=F
)

mapview(hydro_out$subbasins,zcol="link_id",legend=F)+
  mapview(hydro_out$stream_lines,zcol="link_id",legend=F)+
  mapview(hydro_out$links,zcol="link_id",legend=F)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

[Back to top](#contents)
