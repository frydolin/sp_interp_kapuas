###### SPATIAL INTERPOLATION OF PRECIPITATION DATA ON KAPUAS AND SUBCATCHMENTS #####

#### SET UP ####
library("sp")
library("maptools")
library("raster")
library("gstat")
library("hydroTSM")

   projection=CRS("+proj=utm +zone=49 +datum=WGS84 +units=m +no_defs")

#### LOAD FILES ####
   subcatch_shp <-readShapePoly(fn="input/subcatchments/subcatchments.shp", IDvar="catchment", proj4string=projection)
#    kapuas_shp  <-readShapePoly(fn="input/kapuas-basin//kapuas-basin.shp", IDvar="DN", proj4string=projection)
   stations<-readShapePoints("input/stationmap/stationmap.shp")
   stations=as.data.frame(stations)

   m_df       <-read.csv("input/monthly_means.csv", row.names=1, header = TRUE)
   
#### ANALYSIS ####
# subset the data
   dates=as.Date(row.names(m_df))
dates[[2]]
   m_dfx=m_df[which(as.character(dates)=="2001-01-01"):
                   which(as.character(dates)=="2012-12-01"),]
   datesx=dates[which(as.character(dates)=="2001-01-01"):
                    which(as.character(dates)=="2012-12-01")]

####   IDW ####
# Monthly interpolated map of rainfall

# Block
sek_tay.idw.b <- hydrokrige(
   x.ts=m_df, dates=dates,
   from="2001-01-01", to="2012-12-01",
   x.gis=stations, sname="ID", X="long", Y="lat", 
   subcatchments=subcatch_shp, 
   type="block",               
   cell.size= 0.05, grid.type="regular", p4s=projection,
   write2disk=TRUE, out.fmt="csv", fname="output/subcatch_m_idw.csv")

# CELL
#produces a large file! 2X MB
sek_tay.idw.c=list()
for (i in 1:nrow(ms_dfx)){
   sek_tay.idw.c[i] <- hydrokrige(
      x.ts=unlist(ms_dfx[i,]), 
      x.gis=stations, X="long", Y="lat", sname="ID", 
      subcatchments=sek_tay_shp,
      type="cells",               
      cell.size= 0.05, p4s=projection, grid.type="regular",
      stations.plot=TRUE, stations.offset=c(100,100),
      ColorRamp= "Precipitation",   
      main= paste("IDW Precipitation in",i)) }

## Convert to raster brick for easier file handling
sek_tay_brick=lapply(sek_tay.idw.c, raster)
sek_tay_brick=lapply(sek_tay_brick, function(x) x$var1.pred)
sek_tay_brick=brick(sek_tay_brick)

###

   kapuas.idw <- hydrokrige(
      x.ts=ys_mean[-6], x.gis=stations, catchment.name="all",
      X="x", Y="y", sname="ID",
      type= "cells",
      subcatchments=kapuas_shp,
      cell.size= 2000, 
      p4s=projection, grid.type="regular",stations.plot=TRUE, stations.offset=c(100,100),
      ColorRamp= "Precipitation",   
      main= "IDW Precipitation of Kapuas")

# #### ORDINARY KRIGING ####
#    sekayam.ok <- hydrokrige(x.ts=ys_mean, x.gis=stations, catchment.name="all",
#                              X="x", Y="y", sname="ID",
#                              type= "cells", formula=value~1,
#                              subcatchments=sekayam_shp,
#                              cell.size= 1000, 
#                              p4s=projection, grid.type="regular",stations.plot=TRUE, stations.offset=c(100000,100000),
#                              ColorRamp= "Precipitation",   
#                              main= "OK Precipitation of Sekayam")
#    kapuas.ok <- hydrokrige(x.ts=ys_mean[-6], x.gis=stations, catchment.name="all",
#                          X="x", Y="y", sname="ID",
#                          type= "cells", formula=value~1,
#                          subcatchments=kapuas_shp,
#                          cell.size= 2000, 
#                          p4s=projection, grid.type="regular",stations.plot=TRUE, stations.offset=c(100,100),
#                          ColorRamp= "Precipitation",   
#                          main= "OK Precipitation of Kapuas")
#### END ####