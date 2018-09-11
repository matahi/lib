# mean +/- deviation
data_summary <- function(x) {
   m <- mean(x,na.rm=T)
   ymin <- m - sd(x, na.rm=T)
   ymax <- m + sd(x, na.rm=T)
   return(c(y=m,ymin=ymin,ymax=ymax))
}

# quantile summary
quantile_summary <- function(x) {
   m <- median(x,na.rm=T)
   range.info <- quantile(x, probs = c(0.25, 0.75), na.rm=T)
   return(c(y=m,ymin=range.info[[1]],ymax=range.info[[2]]))
}
