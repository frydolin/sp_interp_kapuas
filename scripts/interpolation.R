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
   projection=CRS("+proj=longlat +datum=WGS84 +no_defs")
### END SET UP ###

#### LOAD FILES ####
   subcatch_shp <-readShapePoly(fn="input/subcatchments/subcatchments.shp", IDvar="catchment", proj4string=projection)
    kapuas_shp  <-readShapePoly(fn="input/kapuas-basin//kapuas-basin.shp", IDvar="DN", proj4string=projection)
   stations_shp<-readShapePoints("input/stationmap/stationmap.shp")
   stations_shp=stations_shp[c(-8,-9,-13),] # exclude SGU04, SGU19, STG03
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
      rm(dates)

idw.cell=function(x, catchment){
   idw.c=list()
for (j in 1:nrow(x)){
      idw.c[j] <- hydrokrige(
      x.ts=unlist(x[j,]), 
      x.gis=stations, X="long", Y="lat", sname="ID", 
      subcatchments=catchment,
      type="cells",               
      cell.size= 0.05, p4s=projection, grid.type="regular",
      plot=FALSE)
   }
   rm(j)
   return(idw.c)   
}
subcatch.m.idw.c=idw.cell(x=m_dfx, catchment=subcatch_shp)
kapuas.m.idw.c=idw.cell(x=m_dfx, catchment=kapuas_shp)
### END IDW ###
#### CONVERSION ####
## Convert to raster brick for easier file handling
   make.idw.brick=function(x, time){
   idw_brick=lapply(x, raster)
   idw_brick=lapply(idw_brick, function(x) x$var1.pred)
   idw_brick=brick(idw_brick)
   idw_brick=setZ(idw_brick, z=time, name="time")
   }
   subcatch.m.idw_brick=make.idw.brick(subcatch.m.idw.c, time=datesx)
   kapuas.m.idw_brick=make.idw.brick(kapuas.m.idw.c, time=datesx)
library("rts")
   subcatch.m.idw_rts=rts( subcatch.m.idw_brick,time=datesx)
   kapuas.m.idw_rts=rts(kapuas.m.idw_brick,time=datesx)
rm(subcatch.m.idw.c,kapuas.m.idw.c)
###
#### AGGREGATION ####
   # By month
   mon.fac <- format.Date(datesx,format="%m")
   mon.fac <- factor(mon.fac)
   subcatch.idw.bymonth=zApply(subcatch.m.idw_brick, by=mon.fac, fun=mean, name='months')
      names(subcatch.idw.bymonth)=format.Date(datesx,format="%b")[1:12]
   save(subcatch.idw.bymonth, file="output/subcatch.idw.bymonth")
   kapuas.idw.bymonth=zApply(kapuas.m.idw_brick, by=mon.fac, fun=mean, name='months')
      names(kapuas.idw.bymonth)=format.Date(datesx,format="%b")[1:12]
   # Yearly means
   subcatch.idw.yearly=apply.yearly(subcatch.m.idw_rts, mean)
   names(subcatch.idw.yearly@raster)=as.character(c(2001:2012))
   subcatch.idw.ov.av=period.apply(subcatch.m.idw_rts, 144, mean) # because it's 144month
   save(subcatch.idw.ov.av, file="output/subcatch.idw.ov.av")
   kapuas.idw.yearly=apply.yearly(kapuas.m.idw_rts, mean)
   names(kapuas.idw.yearly@raster)=as.character(c(2001:2012))
   kapuas.idw.ov.av=period.apply(kapuas.idw.yearly, 12, mean) # because it's 12 years 
### END AGGREGATION ###
#### VISUALIZATION ####
   library("rasterVis")
   source("scripts/graphic_pars.R")
   #IDW BY MONTH
      name="output/subcatch_idw_bymonth.png"
      png(filename=name, pointsize = 11, width=17, height=19, units="cm", res=300)
   levelplot(subcatch.idw.bymonth, par.settings=rast.theme, xlab="longitude", ylab="latitude")+ layer(sp.polygons(obj=subcatch_shp, lwd=0.5, col="#333333"))+ layer(sp.points(stations_shp, col="red", cex=0.5))
   dev.off()

      name="output/kapuas_idw_bymonth.png"
      png(filename=name, pointsize = 11, width=17, height=19, units="cm", res=300)
   levelplot(kapuas.idw.bymonth, par.settings=rast.theme, xlab="longitude", ylab="latitude")+ layer(sp.polygons(obj=subcatch_shp, lwd=0.5, col="#777777"))+layer(sp.polygons(obj=kapuas_shp, lwd=0.5, col="#333333"))+ layer(sp.points(stations_shp, col="red", cex=0.5))
   dev.off()

   # IDW FOR THE LAST 12 YEARS
      name="output/subcatch_idw_years.png"
      png(filename=name, pointsize = 11, width=17, height=19, units="cm", res=300)
   levelplot(subcatch.idw.yearly@raster, par.settings=rast.theme, xlab="longitude", ylab="latitude")+ layer(sp.polygons(subcatch_shp, lwd=0.5, col="#333333"))+ layer(sp.points(stations_shp, col="red", cex=0.5))
   dev.off()

      name="output/kapuas_idw_years.png"
      png(filename=name, pointsize = 11, width=17, height=19, units="cm", res=300)
   levelplot(kapuas.idw.yearly@raster, par.settings=rast.theme, xlab="longitude", ylab="latitude")+ layer(sp.polygons(obj=subcatch_shp, lwd=0.5, col="#777777"))+layer(sp.polygons(kapuas_shp, lwd=0.5, col="#333333"))+ layer(sp.points(stations_shp, col="red", cex=0.5))
   dev.off()

   # LONG TERM MEAN
      name="output/subcatch_idw_lt_av.png"
      png(filename=name, pointsize = 11, width=16, height=16, units="cm", res=300)
   levelplot(subcatch.idw.ov.av@raster, par.settings=rast.theme, xlab="longitude", ylab="latitude") + layer(sp.polygons(subcatch_shp, lwd=0.8, col="#333333"))+ layer(sp.points(stations_shp, col="red")) + layer(sp.pointLabel(stations_shp, label=stations_shp$ID),theme=label.theme)
   dev.off()

      name="output/kapuas_idw_lt_av.png"
      png(filename=name, pointsize = 11, width=16, height=16, units="cm", res=300)
   levelplot(kapuas.idw.ov.av@raster, par.settings=rast.theme, xlab="longitude", ylab="latitude") + layer(sp.polygons(obj=subcatch_shp, lwd=0.5, col="#777777"))+layer(sp.polygons(kapuas_shp, lwd=0.8, col="#333333"))+ layer(sp.points(stations_shp, col="red")) + layer(sp.pointLabel(stations_shp, label=stations_shp$ID),theme=label.theme)
   dev.off()
### END VISUALIZATION ###
#### END ####
