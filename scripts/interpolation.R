###### SPATIO-TEMPORAL RAINFALL PATTERNS IN KAPUAS BASIN ######
	### SPATIAL INTERPOLATION OF GAUGE DATA ###
	
## interpolation.R:
## loads, interpolates and plots the results

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
   stations_shp=stations_shp[c(-8,-9,-13),] # exclude inhomogenous stattions SGU04, SGU19, STG03
   stations=as.data.frame(stations_shp)
   
   d_df  <-read.csv("input/daily_data.csv", row.names=1, header = TRUE)
   m_df  <-read.csv("input/monthly_means.csv", row.names=1, header = TRUE)
### END LOAD FILES ###

#### IDW INTERPOLATION ####
### Block
	# only for the subcatchments
	# Daily
	subcatch.d.idw.b <- hydrokrige(
		 x.ts=d_df, 
		 dates=as.Date(row.names(d_df)),from="2001-01-01", to="2012-12-31",
		 x.gis=stations, sname="ID", X="long", Y="lat", 
		 subcatchments=subcatch_shp, 
		 type="block",               
		 cell.size= 0.05, grid.type="regular", p4s=projection,
		 plot=FALSE,
		 write2disk=TRUE, out.fmt="csv", fname="output/subcatch_d_idw.csv") # output written directly do disk

	# Monthly
	subcatch.m.idw.b <- hydrokrige(
		 x.ts=m_df, 
		 dates=as.Date(row.names(m_df)),from="2001-01-01", to="2012-12-01",
		 x.gis=stations, sname="ID", X="long", Y="lat", 
		 subcatchments=subcatch_shp, 
		 type="block",               
		 cell.size= 0.05, grid.type="regular", p4s=projection,
		 plot=FALSE,
		 write2disk=TRUE, out.fmt="csv", fname="output/subcatch_m_idw.csv") # output written directly do disk

### CELL
	# For monthly data
   # subset the data
      dates=as.Date(row.names(m_df))
      m_dfx=m_df[which(as.character(dates)=="2001-01-01"):which(as.character(dates)=="2012-12-01"),]
      datesx=dates[which(as.character(dates)=="2001-01-01"):which(as.character(dates)=="2012-12-01")]
      rm(dates)
	# define custom wrapper function:
	idw.cell=function(x, catchment){
		 idw.c=list()
	for (j in 1:nrow(x)){ # hydrokrige is run for each date seperately
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
	# execute the interpolation:
	subcatch.m.idw.c=idw.cell(x=m_dfx, catchment=subcatch_shp)
	kapuas.m.idw.c=idw.cell(x=m_dfx, catchment=kapuas_shp)
### END IDW ###

#### FORMAT CONVERSION ####
# Convert to raster brick for easier file handling
# and set time slot to make it a time series
   make.idw.brick=function(x, time){
   idw_brick=lapply(x, raster)
   idw_brick=lapply(idw_brick, function(x) x$var1.pred)
   idw_brick=brick(idw_brick)
   idw_brick=setZ(idw_brick, z=time, name="time")
   }
   subcatch.m.idw_brick=make.idw.brick(subcatch.m.idw.c, time=datesx)
   kapuas.m.idw_brick=make.idw.brick(kapuas.m.idw.c, time=datesx)
   
# Convert to raster time series format
	library("rts")
   subcatch.m.idw_rts=rts( subcatch.m.idw_brick,time=datesx)
   kapuas.m.idw_rts=rts(kapuas.m.idw_brick,time=datesx)
rm(subcatch.m.idw.c,kapuas.m.idw.c)
###

#### AGGREGATION ####
   # By month (i.e. mean in all Januaries, etc.)
   mon.fac <- format.Date(datesx,format="%m")
   mon.fac <- factor(mon.fac)
   
   subcatch.idw.bymonth=zApply(subcatch.m.idw_brick, by=mon.fac, fun=mean, name='months')
		 names(subcatch.idw.bymonth)=format.Date(datesx,format="%b")[1:12]
		 save(subcatch.idw.bymonth, file="output/subcatch.idw.bymonth")	# save to disk as it is needed later
		 
   kapuas.idw.bymonth=zApply(kapuas.m.idw_brick, by=mon.fac, fun=mean, name='months')
   	names(kapuas.idw.bymonth)=format.Date(datesx,format="%b")[1:12]
   	
   # Yearly means
   subcatch.idw.yearly=apply.yearly(subcatch.m.idw_rts, mean)
   	names(subcatch.idw.yearly@raster)=as.character(c(2001:2012))   	
   kapuas.idw.yearly=apply.yearly(kapuas.m.idw_rts, mean)
   	names(kapuas.idw.yearly@raster)=as.character(c(2001:2012))

	# Long term average 2001-2012
   subcatch.idw.ov.av=period.apply(subcatch.m.idw_rts, 144, mean) # because it's 144month
   	save(subcatch.idw.ov.av, file="output/subcatch.idw.ov.av")   	# save to disk as it is needed later
   kapuas.idw.ov.av=period.apply(kapuas.idw.yearly, 12, mean) # because it's 12 years 
### END AGGREGATION ###

#### VISUALIZATION ####
   library("rasterVis")
   source("scripts/graphic_pars.R")
   
   #IDW BY MONTH
   subcatch.idw.bymonth=extend(x=subcatch.idw.bymonth, y=(extent(subcatch.idw.bymonth)+0.1)) #enlarge area shown
      name="output/subcatch_idw_bymonth.svg"
   svg(filename=name, pointsize = 11, width=16/2.54, height=16/2.54, family="Lato")
   levelplot(subcatch.idw.bymonth, par.settings=rast.theme, xlab="longitude", ylab="latitude")+ layer(sp.polygons(obj=subcatch_shp, lwd=0.5, col="#333333"))+ layer(sp.points(stations_shp, col="red", cex=0.5))
   dev.off()

   kapuas.idw.bymonth=extend(x=kapuas.idw.bymonth, y=(extent(kapuas.idw.bymonth)+0.2)) #enlarge area shown
      name="output/kapuas_idw_bymonth.png"
      png(filename=name, pointsize = 11, width=17, height=19, units="cm", res=300)
   levelplot(kapuas.idw.bymonth, par.settings=rast.theme, xlab="longitude", ylab="latitude")+ layer(sp.polygons(obj=subcatch_shp, lwd=0.5, col="#777777"))+layer(sp.polygons(obj=kapuas_shp, lwd=0.5, col="#333333"))+ layer(sp.points(stations_shp, col="red", cex=0.5))
   dev.off()

   # FOR THE LAST 12 YEARS
   subcatch.idw.yearly@raster=extend(x=subcatch.idw.yearly@raster, y=(extent(subcatch.idw.yearly@raster)+0.1)) #enlarge area shown
      name="output/subcatch_idw_years.svg"
   svg(filename=name, pointsize = 11, width=16/2.54, height=16/2.54, family="Lato")
   levelplot(subcatch.idw.yearly@raster, par.settings=rast.theme, xlab="longitude", ylab="latitude")+ layer(sp.polygons(subcatch_shp, lwd=0.5, col="#333333"))+ layer(sp.points(stations_shp, col="red", cex=0.5))
   dev.off()

   kapuas.idw.yearly@raster=extend(x=kapuas.idw.yearly@raster, y=(extent(kapuas.idw.yearly@raster)+0.2)) #enlarge area shown
      name="output/kapuas_idw_years.png"
      png(filename=name, pointsize = 11, width=17, height=19, units="cm", res=300)
   levelplot(kapuas.idw.yearly@raster, par.settings=rast.theme, xlab="longitude", ylab="latitude")+ layer(sp.polygons(obj=subcatch_shp, lwd=0.5, col="#777777"))+layer(sp.polygons(kapuas_shp, lwd=0.5, col="#333333"))+ layer(sp.points(stations_shp, col="red", cex=0.5))
   dev.off()

   # LONG TERM MEAN
   subcatch.idw.ov.av=extend(x=subcatch.idw.ov.av@raster, y=(extent(subcatch.idw.ov.av@raster)+0.1)) #enlarge area shown
      name="output/subcatch_idw_lt_av.svg"
      svg(filename=name, pointsize = 11, width=16/2.54, height=16/2.54, family="Lato")
   levelplot(subcatch.idw.ov.av, par.settings=rast.theme, xlab="longitude", ylab="latitude") + layer(sp.polygons(subcatch_shp, lwd=0.8, col="#222222"))+layer(sp.polygons(kapuas_shp, lwd=0.8, col="#666666"))+ layer(sp.points(stations_shp, col="green")) + layer(sp.pointLabel(stations_shp, label=stations_shp$ID),theme=label.theme)
   dev.off()

   kapuas.idw.ov.av=extend(x=kapuas.idw.ov.av@raster, y=(extent(kapuas.idw.ov.av@raster)+0.2)) #enlarge area shown
      name="output/kapuas_idw_lt_av.svg"
			svg(filename=name, pointsize = 11, width=20/2.54, height=16/2.54, family="Lato")
   levelplot(kapuas.idw.ov.av, par.settings=rast.theme, xlab="longitude", ylab="latitude") + layer(sp.polygons(obj=subcatch_shp, lwd=0.5, col="#666666"))+layer(sp.polygons(kapuas_shp, lwd=0.8, col="#222222"))+ layer(sp.points(stations_shp, col="green")) + layer(sp.pointLabel(stations_shp, label=stations_shp$ID),theme=label.theme)
   dev.off()
### END VISUALIZATION ###

###### END interpolation.R ######
