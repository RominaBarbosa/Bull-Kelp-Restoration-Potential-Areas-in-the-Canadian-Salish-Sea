###
### Function to stack rasters from the same path

stack_rasters_path_funtion<- function(path= layers_path, terra_class= "Y"){
  files_rasters <- list.files(path = path, pattern = "\\.tif$", full.names = TRUE)
  
  rasters<- stack(files_rasters[1])
  filename <- strsplit(tools::file_path_sans_ext(files_rasters[1]), "/")[[1]][length(strsplit(tools::file_path_sans_ext(files_rasters[1]), "/")[[1]])]
  names(rasters)<- filename
  
  
  for (i in 2:length(files_rasters)) {
    raster_i<- raster(files_rasters[i])
    filename <- strsplit(tools::file_path_sans_ext(files_rasters[i]), "/")[[1]][length(strsplit(tools::file_path_sans_ext(files_rasters[1]), "/")[[1]])]
    names(raster_i)<- filename
    rasters<- stack(rasters, raster_i)
  }
  
  if(terra_class== "Y"){
    rasters<- terra::rast(rasters)
  }
  return(rasters)
}


