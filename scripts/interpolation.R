###### SPATIAL INTERPOLATION OF PRECIPITATION DATA ON KAPUAS AND SUBCATCHMENTS #####

#### SET UP ####
library(foreign)
library(sp)
library(maptools)
library(hydroTSM)

projection=CRS("+proj=utm +zone=49 +datum=WGS84 +units=m +no_defs")

#### LOAD FILES ####
   sekayam_shp <-readShapePoly(fn="input/sekayam_wgs84utm49N/sekayam-subcatchment.shp", IDvar="DN", proj4string=projection)
   kapuas_shp  <-readShapePoly(fn="input/kapuas_wgs84utm49N/kauas-catchment.shp", IDvar="DN", proj4string=projection)
   stations    <-read.dbf("input/stationmap_wgs84utm49N/stationmap.dbf")
   ys_df       <-read.csv2("input/yearly_sums.csv", row.names=1, header = TRUE, sep=";")
   
#### ANALYSIS ####
   # for now we use long term annual mean sums
   ys_mean=colMeans(ys_df, na.rm=TRUE)

####   IDW ####
   sekayam.idw <- hydrokrige(
               x.ts=ys_mean, x.gis=stations, catchment.name="all",
               X="x", Y="y", sname="ID",
               type= "cells",
               subcatchments=sekayam_shp,
               cell.size= 1000, 
               p4s=projection, grid.type="regular",stations.plot=TRUE, stations.offset=c(100,100),
               ColorRamp= "Precipitation",   
               main= "IDW Precipitation of Sekayam")

   kapuas.idw <- hydrokrige(
      x.ts=ys_mean[-6], x.gis=stations, catchment.name="all",
      X="x", Y="y", sname="ID",
      type= "cells",
      subcatchments=kapuas_shp,
      cell.size= 2000, 
      p4s=projection, grid.type="regular",stations.plot=TRUE, stations.offset=c(100,100),
      ColorRamp= "Precipitation",   
      main= "IDW Precipitation of Kapuas")

#### ORDINARY KRIGING ####
   sekayam.ok <- hydrokrige(x.ts=ys_mean, x.gis=stations, catchment.name="all",
                             X="x", Y="y", sname="ID",
                             type= "cells", formula=value~1,
                             subcatchments=sekayam_shp,
                             cell.size= 1000, 
                             p4s=projection, grid.type="regular",stations.plot=TRUE, stations.offset=c(100000,100000),
                             ColorRamp= "Precipitation",   
                             main= "OK Precipitation of Sekayam")
   kapuas.ok <- hydrokrige(x.ts=ys_mean[-6], x.gis=stations, catchment.name="all",
                         X="x", Y="y", sname="ID",
                         type= "cells", formula=value~1,
                         subcatchments=kapuas_shp,
                         cell.size= 2000, 
                         p4s=projection, grid.type="regular",stations.plot=TRUE, stations.offset=c(100,100),
                         ColorRamp= "Precipitation",   
                         main= "OK Precipitation of Kapuas")
#### END ####