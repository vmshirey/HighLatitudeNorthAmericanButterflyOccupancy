logit <- function(x) log((x)/(1-x))

make.data <- function(tree=NULL, scale, imputeThres=1){
  
  # Read in the species trait data
  sp_traits <- read.csv("../data/taxa/species_range_clim.csv") %>%
    dplyr::mutate(SPID=row_number())
  nsp <- nrow(sp_traits)
  
  # Process spatial scale
  if(scale==50){
    # The spatial scale is 50x50km
    grid <- readRDS("../output/data/grid_50.rds")
    grid_range <- readRDS("../output/data/grid_50_range.rds")
    grid_occur <- readRDS("../output/data/grid_50_occur.rds") %>%
      sf::st_drop_geometry()
    grid_temp <- readRDS("../output/data/temp50.rds")
    grid_precip <- readRDS("../output/data/precip50.rds")
    gridarea <- readRDS("../output/data/grid_50_area.rds") %>% as.vector()
    
    message("Finished loading data...")
  }else if(scale==100){
    # The spatial scale is 100x100km
    grid <- readRDS("../output/data/grid_100.rds")
    grid_range <- readRDS("../output/data/grid_100_range.rds")
    grid_occur <- readRDS("../output/data/grid_100_occur.rds") %>%
      sf::st_drop_geometry()
    grid_temp <- readRDS("../output/data/temp100.rds")
    grid_precip <- readRDS("../output/data/precip100.rds")
    gridarea <- readRDS("../output/data/grid_100_area.rds") %>% as.vector()
    
    message("Finished loading data...")
  }else if(scale==200){
    # The spatial scale in 200x200km
    grid <- readRDS("../output/data/grid_200.rds")
    grid_range <- readRDS("../output/data/grid_200_range.rds")
    grid_occur <- readRDS("../output/data/grid_200_occur.rds") %>%
      sf::st_drop_geometry()
    grid_temp <- readRDS("../output/data/temp200.rds")
    grid_precip <- readRDS("../output/data/precip200.rds")
    gridarea <- readRDS("../output/data/grid_200_area.rds") %>% as.vector()
    
    message("Finished loading data...")
  }else{
    message("Invalid spatial resolution supplied. Must be either 50, 100, or 200.")
    return(NA)
  }
  
  # Get valid indices with climate data
  temp.ind <- grid_temp %>%
    na.omit() %>%
    pull(GID)
  
  # Time periods occur every 5 years, each year is a visit.
  nyr=10
  nvisit=5
  
  # Process the gridded occurrence information
  grid_occur <- grid_occur %>%
    dplyr::mutate(era=(year-year%%5)) %>%
    dplyr::mutate(year=(year-era)+1) %>%
    dplyr::mutate(era=(era-1965)/5) %>%
    left_join(sp_traits, by=c("SPID"="SPID")) %>%
    dplyr::select(SPID, GID, era, year) %>%
    dplyr::filter(SPID!=(nsp+1)) %>%
    dplyr::filter(era<=nyr, era>0) %>%
    unique()
  
  # Process the climate rasters
  temp_layers <- data.frame(layers=names(grid_temp)) %>%
    dplyr::mutate(year=as.numeric(stringr::str_extract(layers, "\\d{4}"))) %>%
    dplyr::mutate(era=(year-year%%5)) %>%
    dplyr::mutate(year=(year-era)+1) %>%
    dplyr::mutate(era=(era-1965)/5) %>%
    dplyr::mutate(id=row_number())
  
  temp_01 <- grid_temp[dplyr::filter(temp_layers,era==1)$id] %>%
    rowMeans()
  temp_02 <- grid_temp[dplyr::filter(temp_layers,era==2)$id] %>%
    rowMeans()
  temp_03 <- grid_temp[dplyr::filter(temp_layers,era==3)$id] %>%
    rowMeans()
  temp_04 <- grid_temp[dplyr::filter(temp_layers,era==4)$id] %>%
    rowMeans()
  temp_05 <- grid_temp[dplyr::filter(temp_layers,era==5)$id] %>%
    rowMeans()
  temp_06 <- grid_temp[dplyr::filter(temp_layers,era==6)$id] %>%
    rowMeans()
  temp_07 <- grid_temp[dplyr::filter(temp_layers,era==7)$id] %>%
    rowMeans()
  temp_08 <- grid_temp[dplyr::filter(temp_layers,era==8)$id] %>%
    rowMeans()
  temp_09 <- grid_temp[dplyr::filter(temp_layers,era==9)$id] %>%
    rowMeans()
  temp_10 <- grid_temp[dplyr::filter(temp_layers,era==10)$id] %>%
    rowMeans()
  
  temp <- cbind(temp_01, temp_02, temp_03, temp_04, temp_05, temp_06, temp_07,
                temp_08, temp_09, temp_10)
  
  precip_layers <- data.frame(layers=names(grid_precip)) %>%
    dplyr::mutate(year=as.numeric(stringr::str_extract(layers, "\\d{4}"))) %>%
    dplyr::mutate(era=(year-year%%5)) %>%
    dplyr::mutate(year=(year-era)+1) %>%
    dplyr::mutate(era=(era-1965)/5) %>%
    dplyr::mutate(id=row_number())
  
  precip_01 <- grid_precip[dplyr::filter(precip_layers,era==1)$id] %>%
    rowMeans()
  precip_02 <- grid_precip[dplyr::filter(precip_layers,era==2)$id] %>%
    rowMeans()
  precip_03 <- grid_precip[dplyr::filter(precip_layers,era==3)$id] %>%
    rowMeans()
  precip_04 <- grid_precip[dplyr::filter(precip_layers,era==4)$id] %>%
    rowMeans()
  precip_05 <- grid_precip[dplyr::filter(precip_layers,era==5)$id] %>%
    rowMeans()
  precip_06 <- grid_precip[dplyr::filter(precip_layers,era==6)$id] %>%
    rowMeans()
  precip_07 <- grid_precip[dplyr::filter(precip_layers,era==7)$id] %>%
    rowMeans()
  precip_08 <- grid_precip[dplyr::filter(precip_layers,era==8)$id] %>%
    rowMeans()
  precip_09 <- grid_precip[dplyr::filter(precip_layers,era==9)$id] %>%
    rowMeans()
  precip_10 <- grid_precip[dplyr::filter(precip_layers,era==10)$id] %>%
    rowMeans()
  
  precip <- cbind(precip_01, precip_02, precip_03, precip_04, precip_05, precip_06, precip_07,
                  precip_08, precip_09, precip_10)
  
  # Create an array for detection data
  X <- array(0, dim=c(nsp=nsp,
                      nsite=nrow(grid),
                      nyr=nyr,
                      nvisit=nvisit))
  
  for(i in 1:nrow(grid_occur)){
    X[grid_occur$SPID[i],
      grid_occur$GID[i],
      grid_occur$era[i],
      grid_occur$year[i]] <- 1
  }
  Z <- apply(X, c(1,2,3), sum)
  
  # Process the imputation of non-detection data
  if(imputeThres>=5){
    message("Warning: number of detected species used to impute non-detections >= 5, may result in sparse data...")
  }
  # keep sites where at least 1 species were seen
  site.keep <- which(apply(X, c(2), sum, na.rm=TRUE) > 0)
  site.keep <- site.keep[site.keep %in% temp.ind]
  
  vis.arr=array(1, dim=c(nsp=nrow(sp_traits),
                         nsite=nrow(grid),
                         nyr=nyr,
                         nvisit=nvisit))
  
  X <- X[,site.keep,,,drop=FALSE]
  Z <- Z[,site.keep,]
  grid_range <- grid_range[site.keep,,drop=FALSE]
  vis.arr <- vis.arr[,site.keep,,,drop=FALSE]
  nsp <- nsp
  nsite <- length(site.keep)
  temp <- temp[site.keep,]
  precip <- precip[site.keep,]
  gridarea <- gridarea[site.keep]
  
  message("Data loaded and prepped for master indexing...")
  
  # Generate a number of unique species detected array
  nsp.detected <- apply(X, 2:4, sum, na.rm = TRUE)
  
  # Sub-function used to generate a master, vectorized index for computational efficiency.
  get.indices <- function(sp){
    vis.arr2 <- vis.arr[sp,,,]
    
    vis.arr2[TRUE] <- 1
    vis.arr2[nsp.detected < imputeThres] <- 0
    
    # Restrict inference to species ranges
    vis.arr2[!grid_range[,sp],,] <- 0
    tmp <- which(vis.arr2==1, arr.ind=TRUE)
    indices <- cbind(rep(sp, nrow(tmp)), tmp)
    
    return(indices)
  }
  master.index <- do.call(rbind, lapply(1:nsp, get.indices))
  colnames(master.index) <- c("sp", "site", "yr", "visit")
  master.index <- unique(master.index)
  
  # Reformat data for exporting
  my.data <- list(X=X[master.index])
  
  my.constants <- list(
    nsp=dim(X)['nsp'],
    nsite=dim(X)['nsite'],
    nyr=dim(X)['nyr'],
    nind=nrow(master.index),
    yrv=master.index[,'yr'],
    sitev=master.index[,'site'],
    spv=master.index[,'sp'],
    temp=temp,
    temp.trait=scale(sp_traits$ave_temp2)[1:nsp],
    precip=precip,
    precip.trait=scale(sp_traits$ave_precip2)[1:nsp],
    gridarea=gridarea,
    VCOV=tree,
    ID=diag(1, nrow=nsp, ncol=nsp),
    X=X[master.index])
  
  range.matrix <- matrix(FALSE, nrow=nsp, ncol=nsite)
  for(sp in 1:nsp){
    range.matrix[sp,which(grid_range[,sp]==1)] <- TRUE
    for(site in 1:nsite){
      if(range.matrix[sp,site]==TRUE){
        range.matrix[sp,site]=site
      }  
    }
  }
  range.matrix[range.matrix==0] <- NA
  
  my.info <- list(kept.sites=site.keep,
                  range.list=range.matrix)
  
  message("Master indices generated, prepping inits...")
  Zst <- array(0, dim=c(my.constants$nsite, 
                        my.constants$nyr, 
                        my.constants$nsp))
  
  for(site in 1:nsite){
    for(yr in 1:nyr){
      for(sp in 1:nsp){
        Zst[site, yr, sp] <- Z[sp, site, yr]
      }
    }
  }
  Zst <- replace(Zst, Zst>1, 1)
  
  my.inits <- list(Z=Zst)
  
  my.params <- c("mu.p.0",
                 "sigma.p.0.sp",
                 "sigma.p.0.site",
                 "p.yr",
                 "p.0.sp",
                 "p.0.site",
                 "psi.area",
                 "sigma.psi.site",
                 "psi.site",
                 "psi.0.sp",
                 "beta.temp",
                 "beta.trait.temp",
                 "beta.precip",
                 "beta.trait.precip",
                 "lambda.temp",
                 "lambda.precip",
                 "sigma.psi.temp",
                 "sigma.psi.precip",
                 "psi.beta.temp.sp",
                 "psi.beta.precip.sp")
  
  # Return model-ready data
  return(list(my.constants=my.constants, my.data=my.data, my.inits=my.inits, my.info=my.info,
              my.params=my.params))
  
}