###### SPATIO-TEMPORAL RAINFALL PATTERNS IN KAPUAS BASIN ######
	### SPATIAL INTERPOLATION OF GAUGE DATA ###

## graphic_pars.R
## sets graphic paramters and creates color scheme
## called from within interpolation.R

# pal <- function(col, border = "transparent", ...)
# {
#    n <- length(col)
#    plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
#         axes = FALSE, xlab = "", ylab = "", ...)
#    rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
# }
col.fun <- colorRampPalette(c("#af5752","peachpuff","royalblue","#0f1c66"))
rast.cols=col.fun(20)
rast.theme<-rasterTheme(region=rast.cols, cex=0.8)
label.theme <- list(add.text=list(cex=0.8,fontfamily='Lato', col="violetred"))
colors()

###### END graphic_pars.R ######
