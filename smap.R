

setwd("/Users/Anika/Downloads/")

library(smapr)
library(sp)
library(raster)

# from: https://cran.r-project.org/web/packages/smapr/vignettes/smapr-intro.html

set_smap_credentials("anika.petach", "IBfhs2009")

available_data <- find_smap(id = 'SPL4SMAU', dates = '2018-06-01', version = 4) #root zone soil moisture

# download data
local_files <- download_smap(available_data, overwrite = FALSE, verbose = FALSE)

local_files$name[1:2]

# exploring data
list_smap(local_files[1, ]) #list all of the data in a file (hdf5 bundles data)

# extract data
sm_raster <- extract_smap(local_files, '/Analysis_Data/sm_rootzone_analysis') #extract all of the data in the data frame local_files, generating a RasterBrick with one layer per file

sm_raster

# visualize data
plot(sm_raster)

# crop extent
us_extent <- extent(c(-155, -67, 24, 62))
us_extent <- as(us_extent, "SpatialPolygons")
sp::proj4string(us_extent) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
us_extent

# check projection is same
proj_us_extent <- spTransform(us_extent, crs(sm_raster))

#proj_us_extent <- projectRaster(sm_raster, crs(us_extent))

# match extents
us_soil_moisture <- crop(sm_raster, proj_us_extent)
plot(us_soil_moisture[[1]])

myCRS <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")

# project to lat long
us_soil_projed <- projectRaster(us_soil_moisture, crs=myCRS)
plot(us_soil_projed[[1]])

writeRaster(us_soil_projed[[1]], "/Users/Anika/Documents/GradSchool/FIA_CompetitionModel/output/smap_proj.grd")




