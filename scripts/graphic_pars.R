

## graphic_pars.R
# pal <- function(col, border = "transparent", ...)
# {
#    n <- length(col)
#    plot(0, 0, type="n", xlim = c(0, 1), ylim = c(0, 1),
#         axes = FALSE, xlab = "", ylab = "", ...)
#    rect(0:(n-1)/n, 0, 1:n/n, 1, col = col, border = border)
# }

col.fun <- colorRampPalette(c("#D7635D","peachpuff","royalblue","#023b95"))
rast.cols=col.fun(10)
rast.theme<-rasterTheme(region=rast.cols)
label.theme <- list(add.text=list(cex=0.9,fontfamily='Lato'))

