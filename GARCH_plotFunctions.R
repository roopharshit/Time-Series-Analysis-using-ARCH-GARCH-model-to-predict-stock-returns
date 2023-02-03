##  R functions for plotting GARCH
#########################################################################

##	plotPI(): plot PI and outliers
##________________________________________________________________________
##		input 	
##			iidWN: returned object of arfimaforecast()
##			garch: returned object of ugarchforecast()
##		Optional
##			col.0: color of iid WN PI	col.out.0: color of iid WN outliers
##			col.g: color of GARCH PI	col.out.g: color of GARCH outliers
##			border: a length 2 vector, border color of PIs
##			col.fore: color of forecasts
##			col.ts: color of actual returns
##
###########################################################################

plot_PI = function(iidWN, garch, col.0 = "#F4CAE4", col.g = "#C6DBEF", border = c("#807DBA","#807DBA"),
 col.ts = "#737373", col.fore = "#7CFC00", col.out.0 = "#D94801", col.out.g = "#2171B5"){
	model = garch@model
	
	## sample size, validation set size and time stamps
	n = model$modeldata$T
	n.fore = model$n.start
	Ynh = model$modeldata$data[(n+1):(n+n.fore)]
	Time  = model$modeldata$index[(n+1):(n+n.fore)]
	
	## model, distribution and coefficients 
	vmodel = model$modeldesc$vmodel
	dist = model$modeldesc$distribution
	coef = model$pars[,"Level"]
	
	## forecasts with GARCH Noises
	foret = garch@forecast$seriesFor["T+1",]	#point forecast
	sigmat = garch@forecast$sigmaFor["T+1",]	#sigma forecast
	cri =  qdist(dist, p = c(.025, .975), skew = coef["skew"], shape = coef["shape"])
	lo.s = foret + cri[1]*sigmat 	# lower bd of PI
	hi.s = foret + cri[2]*sigmat 	# upper bd of PI
	
	## iid Noise
	dist0 = iidWN@model$modeldesc$distribution
	model0 = iidWN@model
	coef0 = model0$pars[,"Level"]	
	cri0 = qdist(dist0, p = c(.025, .975), sigma = coef0["sigma"], skew = coef0["skew"], shape = coef0["shape"])
	forec = iidWN@forecast$seriesFor["T+1",]
	lo.c = forec + cri0[1]
	hi.c = forec + cri0[2]
	
	par(las = 1, mgp = c(1.75,.75,0))
	ymin = min(Ynh,lo.s,lo.c)
	ymax = max(Ynh,hi.s,hi.c)		
	plot(Time, Ynh, type = "n", ylab = c("daily return"), xlab = "time", ylim = c(ymin, ymax), axes = FALSE)
	box(col = "#BDBDBD")
	x.at = Axis(Time, side = 1)
	invisible(sapply(x.at, function(u) abline(v = u, lty = "dotted", col = "#BDBDBD")))
	y.at = axTicks(side = 2)
	axis(2,col = "#BDBDBD", col.ticks = "#000000")
	invisible(sapply(y.at, function(u) abline(h = u, lty = "dotted", col = "#BDBDBD")))
	
	## PI for iid
	colc = col2rgb(col.0)/255
	colc = rgb(colc[1],colc[2], colc[3], alpha = 0.60)
	polygon(c(Time, rev(Time)), c(lo.c, rev(hi.c)),col = colc ,border =  border[1], lwd = 0.5)
	
	## PI for GARCH
	cols = as.vector(col2rgb(col.g))/255   
	cols = rgb(cols[1],cols[2], cols[3], alpha = 0.60)
	polygon(c(Time, rev(Time)), c(lo.s, rev(hi.s)), col = cols, border = border[2], lwd = 0.5)
	
	## Time Series and forecasts
	lines(Time, Ynh, lwd = 0.5, col = col.ts)
	lines(Time, foret, lwd = 0.5, col = col.fore)
	
	 ## outliers of iid PIs
	 points(Time[Ynh > hi.c | Ynh < lo.c], Ynh[Ynh > hi.c | Ynh < lo.c], cex = 0.85,  col = col.out.0,  pch = 16)
	 points(Time[Ynh > hi.s | Ynh < lo.s], Ynh[Ynh > hi.s | Ynh < lo.s], cex = 0.85, lwd = 1.5, col = col.out.g)
	
	legend("topleft", inset = c(0,-0.05), xpd = T, title = c("i.i.d WN"),  legend = c("95% PI ","Outlier"), 
	col = c(border[1],col.out.0),  pt.cex = c(2.5,1), pch = c(22,16), 
	pt.lwd = c(0.5,2), pt.bg = c(colc,"#F16913"), bty = "n")
	
	legend("topleft", inset = c(0.17,-0.05), xpd = T, title = c("GARCH"),  legend = c("95% PI ", "Outlier"), 
	col = c(border[2], col.out.g), pt.cex = c(2.5,1), pch = c(22,1), 
	pt.lwd = c(0.5,2), pt.bg = c(cols, "#2171B5" ), bty = "n")
	
	legend("topright",inset = c(0,0), legend = c("forecast"), lty = 1, col = col.fore, seg.len = 4, bty = "n", lwd = 2)

}

##  Plot the AR and MA roots in one complex plain
##________________________________________________
##  plot_roots() produce a regular R plot
##  input: coef(fit) ## coefficients from a fit


plot_Roots = function(coef,...){
  col1 = "#CB181D"; col2 =  "#2171B5"
  q <- p <- 0
  # AR component
  phi = coef[substr(names(coef),1,2) == "ar"]
  if (length(phi) > 0) {
    test <- abs(phi) > 1e-08
    if (any(test)) {
      p <- max(which(test))
    }
  }
  # MA component
  theta = coef[substr(names(coef),1,2) == "ma"]
  if (length(theta) > 0) {
    test <- abs(theta) > 1e-08
    if (any(test)) {
      q <- max(which(test))
    }
  }
  if((p == 0) && (q == 0)) warning("No roots to plot")
  if((p > 0) && (q == 0)){       
    arvec <- phi
    arroots = polyroot(c(1, -arvec))
    oldpar <- par(pty = "s")
    on.exit(par(oldpar))
    plot(c(-1, 1), c(-1, 1), xlab = "Real", ylab = "Imaginary",
        type = "n", bty = "n", xaxt = "n", yaxt = "n", main = "AR roots", ...)
    axis(1, at = c(-1, 0, 1), line = 0.5, tck = -0.025)
    axis(2, at = c(-1, 0, 1), labels = c("-i", "0", "i"), line = 0.5, tck = -0.025)
    circx <- seq(-1, 1, l = 501)
    circy <- sqrt(1 - circx ^ 2)
    lines(c(circx, circx), c(circy, -circy), col = "gray")
    lines(c(-2, 2), c(0, 0), col = "gray")
    lines(c(0, 0), c(-2, 2), col = "gray")
    inside = abs(arroots) > 1
    points(1/arroots[inside], pch = 19, col = col1, cex = 1.25)
    if (sum(!inside) > 0) {
        points(1 /arroots[!inside], pch = 3, col = col1, lwd = 1.5, cex = 1.25)
    }
    legend("topright", legend = c("AR"), pch = c(19), col = c(col1), pt.cex = 1.25)
  }
  if((p == 0) && (q > 0)){
    mavec <- theta
    maroots = polyroot(c(1, mavec))
    oldpar <- par(pty = "s")
    on.exit(par(oldpar))
    plot(c(-1, 1), c(-1, 1), xlab = "Real", ylab = "Imaginary",
        type = "n", bty = "n", xaxt = "n", yaxt = "n", main = "MA roots", ...)
    axis(1, at = c(-1, 0, 1), line = 0.5, tck = -0.025)
    axis(2, at = c(-1, 0, 1), labels = c("-i", "0", "i"), line = 0.5, tck = -0.025)
    circx <- seq(-1, 1, l = 501)
    circy <- sqrt(1 - circx ^ 2)
    lines(c(circx, circx), c(circy, -circy), col = "gray")
    lines(c(-2, 2), c(0, 0), col = "gray")
    lines(c(0, 0), c(-2, 2), col = "gray")
    inside = abs(maroots) > 1
    points(1/maroots[inside], pch = 19, col = col2, cex = 1.25)
    if (sum(!inside) > 0) {
        points(1/maroots[!inside], pch = 3, col = col2, lwd = 1.5, cex = 1.25)
    }
    legend("topright", legend = c("MA"), pch = c(19), col = c(col2), pt.cex = 1.25)
  }
  if((p > 0) && (q > 0)){
    arvec <- phi
    arroots = polyroot(c(1, -arvec))
    mavec <- theta
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

##  Plot the AR and MA roots in one complex plain as  plot.roots()
##_______________________________________________________
##
##  autoplot_roots() produce a ggplot2 plot
##  input: coef(fit) # coefficients from a fit 
  

autoplot_Roots <-function(coef,...){
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 is needed for this function to work. Install it via install.packages(\"ggplot2\")", call. = FALSE)
    }
      q <- p <- 0
  # AR component
  phi = coef[substr(names(coef),1,2) == "ar"]
  if (length(phi) > 0) {
    test <- abs(phi) > 1e-08
    if (any(test)) {
      p <- max(which(test))
    }
  }
  # MA component
  theta = coef[substr(names(coef),1,2) == "ma"]
  if (length(theta) > 0) {
    test <- abs(theta) > 1e-08
    if (any(test)) {
      q <- max(which(test))
    }
  }
  if((p == 0) && (q == 0)) warning("No roots to plot")
  if((p > 0) && (q == 0)){
          arvec <- phi
          arroots = polyroot(c(1, -arvec))
          arData <- maData <- NULL
          allRoots <- data.frame(roots = numeric(0), type = character(0))
          allRoots <-rbind(allRoots, data.frame(roots = arroots, type = "AR"))
          allRoots$Real<-Re(1/allRoots$roots)
          allRoots$Imaginary<-Im(1/allRoots$roots)
          allRoots$UnitCircle <- factor(ifelse((abs(allRoots$roots) > 1), "Within", "Outside"))
          g <- ggplot2::ggplot(ggplot2::aes_(x = ~Real, y = ~Imaginary, colour = ~type), data = allRoots)
          g<-g + ggplot2::scale_colour_manual(values = c("#CB181D"))
          g <- g + ggplot2::coord_fixed(ratio = 1)
          g <- g + ggplot2::annotate("path", x = cos(seq(0, 2 * pi, length.out = 100)), y = sin(seq(0, 2 * pi, length.out = 100)))
          g <- g + ggplot2::geom_vline(xintercept = 0)
          g <- g + ggplot2::geom_hline(yintercept = 0)
          g <- g + ggplot2::ggtitle(label ="AR roots")
          g <- g + ggplot2::xlab("Real") 
          g <- g + ggplot2::ylab("Imaginary")
          g <- g + ggplot2::geom_point(size = 4)
  }
  if((p == 0) && (q > 0)){
        mavec <- theta
        maroots = polyroot(c(1, mavec))
        arData <- maData <- NULL
        allRoots <- data.frame(roots = numeric(0), type = character(0))
        allRoots<-rbind(allRoots, data.frame(roots = maroots, type = "MA"))
        allRoots$Real<-Re(1/allRoots$roots)
        allRoots$Imaginary<-Im(1/allRoots$roots)
        allRoots$UnitCircle <- factor(ifelse((abs(allRoots$roots) > 1), "Within", "Outside"))
        g <- ggplot2::ggplot(ggplot2::aes_(x = ~Real, y = ~Imaginary, colour = ~type), data = allRoots)
        g<-g + ggplot2::scale_colour_manual(values = c("#2171B5"))
        g <- g + ggplot2::coord_fixed(ratio = 1)
        g <- g + ggplot2::annotate("path", x = cos(seq(0, 2 * pi, length.out = 100)), y = sin(seq(0, 2 * pi, length.out = 100)))
        g <- g + ggplot2::geom_vline(xintercept = 0)
        g <- g + ggplot2::geom_hline(yintercept = 0)
        g <- g + ggplot2::ggtitle(label ="MA roots")
        g <- g + ggplot2::xlab("Real") 
        g <- g + ggplot2::ylab("Imaginary")
        g <- g + ggplot2::geom_point(size = 4)
  }
  if((p > 0) && (q > 0)){
        arvec <- phi
        arroots = polyroot(c(1, -arvec))
        mavec <- theta
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












