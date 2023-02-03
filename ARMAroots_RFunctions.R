##  Plot the AR and MA roots in one complex plain
##	Required package: forecast
##________________________________________________
##  plot_roots() produce a regular R plot
##  input is a fit object from arima(), Arima() or auto.arima()
  

plot_roots = function(obj,...){
  q <- p <- 0
  # AR component
  if (length(obj$model$phi) > 0) {
    test <- abs(obj$model$phi) > 1e-08
    if (any(test)) {
      p <- max(which(test))
    }
  }
  # MA component
  if (length(obj$model$theta) > 0) {
    test <- abs(obj$model$theta) > 1e-08
    if (any(test)) {
      q <- max(which(test))
    }
  }
  if((p == 0) && (q == 0)) warning("No roots to plot")
  if((p > 0) && (q == 0)){
        plot(obj, type = "ar", main = "AR roots")
  }
  if((p == 0) && (q > 0)){
        plot(obj, type = "ma", main = "MA roots")
  }
  if((p > 0) && (q > 0)){
    arvec <- obj$model$phi
    arroots = polyroot(c(1, -arvec))
    mavec <- obj$model$theta
    maroots = polyroot(c(1, mavec))
    oldpar <- par(pty = "s")
    on.exit(par(oldpar))
    plot(c(-1, 1), c(-1, 1), xlab = "Real", ylab = "Imaginary",
        type = "n", bty = "n", xaxt = "n", yaxt = "n", main = "AR and MA roots", ...)
    axis(1, at = c(-1, 0, 1), line = 0.5, tck = -0.025)
    axis(2, at = c(-1, 0, 1), labels = c("-i", "0", "i"), line = 0.5, tck = -0.025)
    circx <- seq(-1, 1, l = 501)
    circy <- sqrt(1 - circx ^ 2)
    lines(c(circx, circx), c(circy, -circy), col = "gray")
    lines(c(-2, 2), c(0, 0), col = "gray")
    lines(c(0, 0), c(-2, 2), col = "gray")
    col1 = "#CB181D"; col2 =  "#2171B5"
    inside = abs(arroots) > 1
    points(1/arroots[inside], pch = 19, col = col1, cex = 1.25)
    if (sum(!inside) > 0) {
        points(1 /arroots[!inside], pch = 3, col = col1, lwd = 1.5, cex = 1.25)
    }
    inside = abs(maroots) > 1
    points(1/maroots[inside], pch = 19, col = col2, cex = 1.25)
    if (sum(!inside) > 0) {
        points(1/maroots[!inside], pch = 3, col = col2, lwd = 1.5, cex = 1.25)
    }
    legend("topright", legend = c("AR", "MA"), pch = c(19,19), col = c(col1,col2), pt.cex = 1.25)
  }
}

##_______________________________________________________
##
##  Plot the AR and MA roots in one complex plain as  plot.roots()
##  autoplot.roots() produce a ggplot2 plot
##  input: a fit object from arima(), Arima() or auto.arima()
  

autoplot_roots <-function(obj,...){
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is needed for this function to work. Install it via install.packages(\"ggplot2\")", call. = FALSE)
    }
      q <- p <- 0
  # AR component
  if (length(obj$model$phi) > 0) {
    test <- abs(obj$model$phi) > 1e-08
    if (any(test)) {
      p <- max(which(test))
    }
  }
  # MA component
  if (length(obj$model$theta) > 0) {
    test <- abs(obj$model$theta) > 1e-08
    if (any(test)) {
      q <- max(which(test))
    }
  }
  if((p == 0) && (q == 0)) warning("No roots to plot")
  if((p > 0) && (q == 0)){
        g = autoplot(obj, type = "ar", main = "AR roots")
  }
  if((p == 0) && (q > 0)){
        g = autoplot(obj, type = "ma", main = "MA roots")
  }
  if((p > 0) && (q > 0)){
        arvec <- obj$model$phi
        arroots = polyroot(c(1, -arvec))
        mavec <- obj$model$theta
        maroots = polyroot(c(1, mavec))
        arData <- maData <- NULL
        allRoots <- data.frame(roots = numeric(0), type = character(0))
        allRoots <-rbind(allRoots, data.frame(roots = arroots, type = "AR"))
        allRoots<-rbind(allRoots, data.frame(roots = maroots, type = "MA"))
        allRoots$Real<-Re(1/allRoots$roots)
        allRoots$Imaginary<-Im(1/allRoots$roots)
        allRoots$UnitCircle <- factor(ifelse((abs(allRoots$roots) > 1), "Within", "Outside"))
        g <- ggplot2::ggplot(ggplot2::aes_(x = ~Real, y = ~Imaginary, colour = ~type), data = allRoots)
        g<-g + ggplot2::scale_colour_manual(values = c("#CB181D", "#2171B5"))
        g <- g + ggplot2::coord_fixed(ratio = 1)
        g <- g + ggplot2::annotate("path", x = cos(seq(0, 2 * pi, length.out = 100)), y = sin(seq(0, 2 * pi, length.out = 100)))
        g <- g + ggplot2::geom_vline(xintercept = 0)
        g <- g + ggplot2::geom_hline(yintercept = 0)
        g <- g + ggplot2::ggtitle(label ="AR and MA roots")
        g <- g + ggplot2::xlab("Real") 
        g <- g + ggplot2::ylab("Imaginary")
        g <- g + ggplot2::geom_point(size = 4)
        
  }
  g
}












