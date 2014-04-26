###### SPATIO-TEMPORAL RAINFALL PATTERNS IN KAPUAS BASIN ######
	### SPATIAL INTERPOLATION OF GAUGE DATA ###

## graphic_pars.R
## sets graphic paramters and creates color scheme
## called from within interpolation.R

	col.fun <- colorRampPalette(c("#af5752","peachpuff","royalblue","#0f1c66"))
	rast.cols=col.fun(20)
	rast.theme<-rasterTheme(region=rast.cols, cex=0.8)
	label.theme <- list(add.text=list(cex=0.8,fontfamily='Lato', col="violetred"))

###### END graphic_pars.R ######

