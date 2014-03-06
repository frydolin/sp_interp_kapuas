###### SPATIAL INTERPOLATION OF PRECIPITATION DATA ON KAPUAS AND SUBCATCHMENTS #####

#### SET UP ####
## load libraries
   library("sp")
   library("maptools")
   library("raster")
   library("gstat")
   library("hydroTSM")

## set up time locale to get English time format
   Sys.setlocale("LC_TIME", "en_US.UTF-8") 

## define projection
   projection=CRS("+proj=utm +zone=49 +datum=WGS84 +units=m +no_defs")
### END SET UP ###

#### LOAD FILES ####
   subcatch_shp <-readShapePoly(fn="input/subcatchments/subcatchments.shp", IDvar="catchment", proj4string=projection)
#    kapuas_shp  <-readShapePoly(fn="input/kapuas-basin//kapuas-basin.shp", IDvar="DN", proj4string=projection)
   stations_shp<-readShapePoints("input/stationmap/stationmap.shp")
   stations=as.data.frame(stations_shp)
   
   d_df  <-read.csv("input/daily_data.csv", row.names=1, header = TRUE)
   m_df  <-read.csv("input/monthly_means.csv", row.names=1, header = TRUE)
### END LOAD FILES ###

####   IDW ####
### Block
# Daily
subcatch.d.idw.b <- hydrokrige(
   x.ts=d_df, 
   dates=as.Date(row.names(d_df)),from="2001-01-01", to="2012-12-31",
   x.gis=stations, sname="ID", X="long", Y="lat", 
   subcatchments=subcatch_shp, 
   type="block",               
   cell.size= 0.05, grid.type="regular", p4s=projection,
   plot=FALSE,
   write2disk=TRUE, out.fmt="csv", fname="output/subcatch_d_idw.csv")

# Monthly
subcatch.m.idw.b <- hydrokrige(
   x.ts=m_df, 
   dates=as.Date(row.names(m_df)),from="2001-01-01", to="2012-12-01",
   x.gis=stations, sname="ID", X="long", Y="lat", 
   subcatchments=subcatch_shp, 
   type="block",               
   cell.size= 0.05, grid.type="regular", p4s=projection,
   plot=FALSE,
   write2disk=TRUE, out.fmt="csv", fname="output/subcatch_m_idw.csv")

### CELL
# For monthly data
# subset the data
   dates=as.Date(row.names(m_df))
   m_dfx=m_df[which(as.character(dates)=="2001-01-01"):which(as.character(dates)=="2012-12-01"),]
   datesx=dates[which(as.character(dates)=="2001-01-01"):which(as.character(dates)=="2012-12-01")]

subcatch.m.idw.c=list()
for (j in 1:nrow(m_dfx)){
   subcatch.m.idw.c[j] <- hydrokrige(
      x.ts=unlist(m_dfx[j,]), 
      x.gis=stations, X="long", Y="lat", sname="ID", 
      subcatchments=subcatch_shp,
      type="cells",               
      cell.size= 0.05, p4s=projection, grid.type="regular",
      plot=FALSE)
   }
   rm(j)
### END IDW ###
#### CONVERSION ####
## Convert to raster brick for easier file handling
   m.idw_brick=lapply(subcatch.m.idw.c, raster)
   m.idw_brick=lapply(m.idw_brick, function(x) x$var1.pred)
   m.idw_brick=brick(m.idw_brick)
   m.idw_brick=setZ(m.idw_brick, z=datesx, name="time")
   library("rts")
   m.idw_rts=rts(m.idw_brick,time=datesx)
###
#### AGGREGATION ####
   # By month
   mon.fac <- format.Date(datesx,format="%m")
   mon.fac <- factor(mon.fac)
   idw.bymonth=zApply(m.idw_brick, by=mon.fac, fun=mean, name='months')
   names(idw.bymonth)=format.Date(datesx,format="%b")[1:12]
   # Yearly means
   idw.yearly=apply.yearly(m.idw_rts, mean)
   names(idw.yearly@raster)=as.character(c(2001:2012))
   idw.ov.av=period.apply(idw.yearly, 12, mean) # because it's 12 years 
### END AGGREGATION ###
#### VISUALIZATION ####
   library("rasterVis")
   #IDW BY MONTH
   levelplot(idw.bymonth, par.settings=rast.theme, xlab="longitude", ylab="latitude")+ layer(sp.polygons(subcatch_shp, lwd=0.8, col="#333333"))+ layer(sp.points(stations_shp, col="red"))
   # IDW FOR THE LAST 12 YEARS
   levelplot(idw.yearly@raster, par.settings=rast.theme, xlab="longitude", ylab="latitude")+ layer(sp.polygons(subcatch_shp, lwd=0.8, col="#333333"))+ layer(sp.points(stations_shp, col="red"))

   # LONG TERM MEAN
   levelplot(idw.ov.av@raster, par.settings=rast.theme, xlab="longitude", ylab="latitude") + layer(sp.polygons(subcatch_shp, lwd=0.8, col="#333333"))+ layer(sp.points(stations_shp, col="red")) + layer(sp.pointLabel(stations_shp, label=stations_shp$ID),theme=label.theme)
### END VISUALIZATION ###
#### END ####