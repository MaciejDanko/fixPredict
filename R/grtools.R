
#' Plotting confidence intervals
#' 
#' @param x,fit x and y argument of the fit
#' @param ci a matrix (\code{length(x)} x \code{2}) with upper and lower confidence intervals for the fit.
#' @param first a logical indicating if this is a first plot.
#' @param xlab,ylab,xlim,ylim,col,lwd,lty,xaxt,yaxt graphical parameters passed to 
#' the \code{\link{plot}{graphics}} and the \code{\link{line}{graphics}} functions; see also \code{\link{par}{graphics}}.
#' @param border,density,angle graphical parameters passed the \code{\link{polygon}{graphics}}.
#' @param lwd.ci,lty.ci lty and lwd eqivalents passed to the \code{\link{polygon}{graphics}}.
#' @param af \code{alpha.f} of the \code{col} passed to the \code{\link{polygon}{graphics}}.
#' @param ... other graphics parmeters.
#' @export
#' @seealso fixPredict
CIplot.ci <- function(x, fit, ci, first=FALSE,xlab='',ylab='',xlim=NULL,ylim=NULL,
                      col=1,af=0.3,border=NA,lwd=1.5,lty=1,density=40,angle=45,xaxt='n',yaxt='n',
                      lwd.ci=1,lty.ci=1,...){
  y <- fit
  Y <- c(ci[,1],rev(ci[,2]))
  X <- c(x, rev(x))
  if (first) graphics::plot(x, y, type='l', xlab=xlab,xlim=xlim,ylab=ylab,ylim=ylim,col=col,xaxt=xaxt,yaxt=yaxt,lty=lty,...) 
  graphics::polygon(X,Y,col=grDevices::adjustcolor(col,alpha.f = af), border=border, density=density,angle=angle,lty=lty.ci,lwd=lwd.ci,...)
  graphics::lines(x, y, col=col, lwd=lwd,lty=lty,...)
}


