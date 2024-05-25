# Copied and slightly altered functions from ggridges package
## Necessary to hijack plotting of quantile lines on density plots and instead
## plot 95% HPD Credible Intervals and optionally posterior modes
### (from `coda::HPDinterval()` and `MCMCglmm::posterior.mode()`)
postPlot <- function(posterior, bw = "nrd", #TODO make separate prior/posterior bws?
                     plotHist = TRUE, histbreaks = 100, histcol = NULL,
                     prior = NULL, prange = c("prior", "posterior"),
                     main, sub, ylim,
                     denslwd = 6, denslty = "solid", denscol = "black",
                     hpdcol = "grey70",
                     meanpch = 23, meanbg = "white", meancol = "red",
                     modepch = 3, modebg = NA, modecol = "blue",
                     priorlwd = 3, priorlty = "solid", priorcol = "green",
                     at1 = NULL, at2 = NULL,
                     labels1 = NULL, labels2 = NULL,
                     plot = TRUE, ...){
  
  cl <- match.call()
  # bw.nrd() is like `bwf()` in `coda::densplot()` 
  poDens <- stats::density(posterior, bw = bw, kernel = "gaussian", n = 2^13)
  bw <- poDens$bw
  main.title <- if(missing(main)) "" else main
  sub.title <- if(missing(sub)) "" else sub
  histout <- if(plotHist){
    graphics::hist(posterior, breaks = histbreaks, plot = FALSE)
  } else list(breaks = NULL, counts = NULL, density = NULL, mids = NULL)
  maxYhistpoDens <- max(c(poDens$y, histout$density))
  if(missing(ylim)){
    ylimit <- c(0, maxYhistpoDens)
  } else ylimit <- ylim
  
  # Adjust density at boundaries
  ## reflect density estimates back across a boundary
  constraint <- "unbounded"
  po <- posterior
  if(max(po) <=1 && 1 - max(po) < 2 * bw){
    if(min(po) >= 0 && min(po) < 2 * bw){
      constraint <- "zero-one"
      po <- c(po, -po, 2 - po)
    } else{
      if(min(po) < 0 && 1 + min(po) < 2 * bw){
        if(min(po) >= -1){  #<-- otherwise 'unbounded' happens to have max near 1
          constraint <- "correlation"
          po <- c(po, -2 - po, 2 - po)
        }
      } else{
        constraint <- "upbound1"  #<-- most likely correlation converge on 1
        po <- c(po, 2-po)
      }
    }
  } else{
    if(min(po) >= 0 && min(po) < 2 * bw){
      constraint <- "positive"
      po <- c(po, -po)
    }
    if(min(po) >= -1 && 1 + min(po) < 2 * bw){
      constraint <- "loboundNeg1"  #<-- most likely correlation converge on -1
      po <- c(po, -2 - po)
    }
  }
  poDens <- stats::density(po, kernel = "gaussian", width = 4 * bw, n = 2^13)
  if(constraint == "zero-one"){
    poDens$y <- 3 * poDens$y[poDens$x >= 0 & poDens$x <=1]
    poDens$x <- poDens$x[poDens$x >= 0 & poDens$x <= 1]
  }
  if(constraint == "correlation"){
    poDens$y <- 3 * poDens$y[poDens$x >= -1 & poDens$x <=1]
    poDens$x <- poDens$x[poDens$x >= -1 & poDens$x <= 1]
  }
  if(constraint == "upbound1"){  #<-- most likely correlation converge on 1
    poDens$y <- 2 * poDens$y[poDens$x <=1]
    poDens$x <- poDens$x[poDens$x <= 1]
  }
  if(constraint == "loboundNeg1"){  #<-- most likely correlation converge on -1
    poDens$y <- 2 * poDens$y[poDens$x >= -1]
    poDens$x <- poDens$x[poDens$x >= -1]
  }
  if(constraint == "positive"){
    poDens$y <- 2 * poDens$y[poDens$x >= 0]
    poDens$x <- poDens$x[poDens$x >= 0]
  }
  
  if(plot){
    # set axis label off a little
    ## do before place actual axis
    par(mgp = par("mgp") + c(1, 0, 0))
    plot(poDens, type = "n", axes = FALSE,
         main = main.title, sub = sub.title, ylim = ylimit, ...)
    if(plotHist) graphics::hist(posterior, breaks = histbreaks, col = histcol,
                                freq = FALSE,
                                add = TRUE, ...)
    graphics::lines(poDens, lwd = denslwd, col = denscol, ...)
    yaxmax <- ifelse(is.null(at2), par("yaxp")[2], max(at2))
    lineylim <- ifelse(yaxmax >= maxYhistpoDens, yaxmax, maxYhistpoDens)
    porange <- range(posterior)
    hpd <- coda::HPDinterval(posterior)
    pomean <- mean(posterior)
    pomode <- MCMCglmm::posterior.mode(posterior)
    
    # plot range of all samples 
    axis(1, at = porange, labels = FALSE, tcl = 0, lwd = 2.5, lend = "square",
         line = 0.25)
    # Plot HPD, mean, and median
    ## get the y-coordinates for axis lines
    poRngAsAx <- graphics::grconvertY(graphics::grconvertY(0, from = "npc",
                                                           to = "lines") - 0.25,
                                      from = "lines", to = "user")
    xpd_old <- par("xpd")    #<-- capture old xpd value
    par(xpd = TRUE)
    # Plot HPD credible interval
    graphics::rect(xleft = hpd[1], ybottom = 1.5 * poRngAsAx,
                   xright = hpd[2], ytop = 0.5 * poRngAsAx,
                   col = hpdcol, lwd = 0.001)
    ## POINTS
    graphics::points(x = c(pomean, pomode), y = rep(poRngAsAx, 2),
                     pch = c(meanpch, modepch),
                     bg = c(meanbg, modebg), col = c(meancol, modecol),
                     cex = 1.75, lwd = 3) 
    par(xpd = xpd_old)   #<-- reset xpd
    
    #    axis(1, at = at1, labels = labels1, line = par("mgp")[3] + 1)
    axis(1, at = at1, labels = labels1, line = par("mgp")[3] + 1)
    axis(2, at = at2, labels = labels2)
  }
  
  #######
  # Prior
  #######
  if(!is.null(prior)){
    if(is.list(prior)){
      if(names(prior)[1L] != "x" || names(prior)[2L] != "y"){
        stop("prior must be either an object of class `mcmc`, a `list` containing a density, or a function that returns either of these")
      }
      prDens <- prior
    }
    
    if(coda::is.mcmc(prior)){
      ## range of prior/posterior
      if(missing(prange)) prange <- "prior"
      if(!startsWith(prange, "prior") && !startsWith(prange, "posterior")){
        stop("argument to 'prange' must start with the full word of either 'prior' or 'posterior'")
      }
      if(startsWith(prange, "prior")){
        prraN <- as.numeric(strsplit(prange, "prior")[[1L]][2L])
        prra <- range(prior)
      }
      if(startsWith(prange, "posterior")){
        prraN <- as.numeric(strsplit(prange, "posterior")[[1L]][2L])
        prra <- range(posterior)
      }
      if(is.na(prraN)) prraN <- 1
      # General formula for finding bounds
      ## Lower: 2*bound - prior (then select all >= bound)
      ## Upper: 2*bound - prior (then select all <= bound)
      if(constraint == "zero-one"){
        pr <- prior[prior >= 0 & prior <= 1]
        pr <- c(pr, -pr, 2 - pr)
        prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
        prDens$y <- 3 * prDens$y[prDens$x >= 0 & prDens$x <=1]
        prDens$x <- prDens$x[prDens$x >= 0 & prDens$x <= 1]
      }
      if(constraint == "upbound1" | constraint == "loboundNeg1"){
        if(min(prior) >= -1 && 1 + min(prior) < 2 * bw &&
           max(prior) <= 1 && 1 - max(prior) < 2 * bw) constraint <- "correlation"
      }
      if(constraint == "correlation"){
        pr <- prior[prior >= -1 & prior <= 1]
        pr <- c(pr, -2 - pr, 2 - pr)
        prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
        prDens$y <- 3 * prDens$y[prDens$x >= -1 & prDens$x <=1]
        prDens$x <- prDens$x[prDens$x >= -1 & prDens$x <= 1]
      }
      if(constraint == "upbound1"){  #<-- most likely correlation converge on 1
        pr <- prior[prior <= 1]
        pr <- c(pr, 2 - pr)
        prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
        prDens$y <- 2 * prDens$y[prDens$x <=1]
        prDens$x <- prDens$x[prDens$x <= 1]
      }
      if(constraint == "loboundNeg1"){  #<-- most likely correlation converge on -1
        pr <- prior[prior >= -1]
        pr <- c(pr, -2 - pr)
        prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
        prDens$y <- 2 * prDens$y[prDens$x >= -1]
        prDens$x <- prDens$x[prDens$x >= -1]
      }
      if(constraint == "positive"){
        #FIXME issues when prior range is used on very flat parameter expanded prior
        ## Density is high near zero, so usually above posterior plot and off figure
        #XXX NOTE if problem occurs, consider calculate own `bw` and not use
        ## `width` with posterior bw (see commented code below)
        if(max(prior) > prra[2L]*prraN){
          pr <- prior[prior >= 0 & prior <= prra[2L]*prraN]
          pr <- c(pr, -pr, 2*prra[2L]*prraN - pr)
          #          prDens <- stats::density(pr, bw = "nrd", kernel = "gaussian", n = 2^13)
          prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
          prDens$y <- 3 * prDens$y[prDens$x >= 0 & prDens$x <= prra[2L]*prraN]
          prDens$x <- prDens$x[prDens$x >= 0 & prDens$x <= prra[2L]*prraN]
        } else{
          pr <- prior[prior >= 0]
          pr <- c(pr, -pr)
          #            prDens <- stats::density(pr, bw = "nrd", kernel = "gaussian", n = 2^13)
          prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
          prDens$y <- 2 * prDens$y[prDens$x >= 0]
          prDens$x <- prDens$x[prDens$x >= 0]
        }
      }
      if(constraint == "unbounded"){
        #XXX NOTE if problem occurs, consider calculate own `bw` and not use
        ## `width` with posterior bw (see commented code below)
        if(prra[1L] < 0) prra[1L] <- prra[1L] * prraN
        if(prra[1L] >= 0) prra[1L] <- prra[1L] * (1/prraN)
        if(prra[2L] < 0) prra[2L] <- prra[2L] * (1/prraN)
        if(prra[2L] >= 0) prra[2L] <- prra[2L] * prraN
        pr <- prior[prior >= prra[1L] & prior <= prra[2L]]
        pr <- c(pr, 2*prra[1L] - pr, 2*prra[2L] - pr)
        #        prDens <- stats::density(pr, bw = "nrd", kernel = "gaussian", n = 2^13)
        prDens <- stats::density(pr, kernel = "gaussian", width = 4 * bw, n = 2^13)
        prDens$y <- 3 * prDens$y[prDens$x >= prra[1L] & prDens$x <= prra[2L]]
        prDens$x <- prDens$x[prDens$x >= prra[1L] & prDens$x <= prra[2L]]
      }
    }
    
    if(!coda::is.mcmc(prior) && !is.list(prior)){
      warning("prior is not a `list` or `mcmc` object - no prior added to plot")
    } else{
      # Add prior density to plot
      if(plot){
        # to keep asymptotically increasing priors from 'running away'
        ##TODO check the following restriction for lots of different priors
        prDensSub <- which(prDens$y <= lineylim)
        graphics::lines(prDens$y[prDensSub] ~ prDens$x[prDensSub],
                        col = priorcol, lty = priorlty, lwd = priorlwd, ...)
      }
    }
  } else prDens <- NULL
  
  return(invisible(list(call = cl,
                        postDensity = list(constraint = constraint,
                                           density = poDens,
                                           histogram = histout),
                        priorDensity = list(constraint = constraint,
                                            density = prDens)
  )))
  
}





stat_density_ridges_HPDCrI <- function(mapping = NULL, data = NULL, geom = "density_ridges",
                     position = "identity", na.rm = FALSE, show.legend = NA,
                     inherit.aes = TRUE, bandwidth = NULL, from = NULL, to = NULL,
                     jittered_points = FALSE, quantile_lines = FALSE, calc_ecdf = FALSE, quantiles = 4,
                     quantile_fun = quantile, n = 512, ...)
{
  layer(
    stat = StatDensityRidgesHPDCrI,
    data = data,
    mapping = mapping,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(bandwidth = bandwidth,
                  from = from,
                  to = to,
                  calc_ecdf = calc_ecdf,
                  quantiles = quantiles,
                  jittered_points = jittered_points,
                  quantile_lines = quantile_lines,
                  quantile_fun = quantile_fun,
                  n = n,
                  na.rm = na.rm, ...)
  )
}


################################################################################
StatDensityRidgesHPDCrI <- ggproto("StatDensityRidgesHPDCrI", Stat,
  required_aes = "x",

  default_aes = aes(height = ..density..),

  calc_panel_params = function(data, params) {
    if (is.null(params$bandwidth)) {
      xdata <- na.omit(data.frame(x=data$x, group=data$group))
      xs <- split(xdata$x, xdata$group)
      xs_mask <- vapply(xs, length, numeric(1)) > 1
      bws <- vapply(xs[xs_mask], bw.nrd0, numeric(1))
      bw <- mean(bws, na.rm = TRUE)
      message("Picking joint bandwidth of ", signif(bw, 3))

      params$bandwidth <- bw
    }

    if (is.null(params$from)) {
      params$from <- min(data$x, na.rm=TRUE) - 3 * params$bandwidth
    }

    if (is.null(params$to)) {
      params$to <- max(data$x, na.rm=TRUE) + 3 * params$bandwidth
    }

    data.frame(
      bandwidth = params$bandwidth,
      from = params$from,
      to = params$to
    )
  },

  setup_params = function(self, data, params) {
    # calculate bandwidth, min, and max for each panel separately
    panels <- split(data, data$PANEL)
    pardata <- lapply(panels, self$calc_panel_params, params)
    pardata <- ggridges:::reduce(pardata, rbind)

    if (length(params$quantiles) > 1 &&
        (max(params$quantiles, na.rm = TRUE) > 1 || min(params$quantiles, na.rm = TRUE) < 0)) {
      stop('invalid quantiles used: c(', paste0(params$quantiles, collapse = ','), ') must be within [0, 1] range')
    }

    params$bandwidth <- pardata$bandwidth
    params$from <- pardata$from
    params$to <- pardata$to
    params
  },

  compute_group = function(data, scales, from, to, bandwidth = 1,
                           calc_ecdf = FALSE, jittered_points = FALSE, quantile_lines = FALSE,
                           quantiles = 4, quantile_fun = quantile, n = 512) {
    # ignore too small groups
    if(nrow(data) < 3) return(data.frame())

    if (is.null(calc_ecdf)) calc_ecdf <- FALSE
    if (is.null(jittered_points)) jittered_points <- FALSE
    if (is.null(quantile_lines)) quantile_lines <- FALSE

    # when quantile lines are requested, we also calculate ecdf
    # this simplifies things for now; in principle, could disentangle
    # the two
    if (quantile_lines) calc_ecdf <- TRUE

    panel <- unique(data$PANEL)
    if (length(panel) > 1) {
      stop("Error: more than one panel in compute group; something's wrong.")
    }
    panel_id <- as.numeric(panel)

    d <- stats::density(
      data$x,
      bw = bandwidth[panel_id], from = from[panel_id], to = to[panel_id], na.rm = TRUE,
      n = n
    )

    # calculate maximum density for scaling
    maxdens <- max(d$y, na.rm = TRUE)

    # make interpolating function for density line
    densf <- approxfun(d$x, d$y, rule = 2)

    # calculate jittered original points if requested
    if (jittered_points) {
      df_jittered <- data.frame(
        x = data$x,
        # actual jittering is handled in the position argument
        density = densf(data$x),
        ndensity = densf(data$x) / maxdens,
        datatype = "point", stringsAsFactors = FALSE)

      # see if we need to carry over other point data
      # capture all data columns starting with "point", as those are relevant for point aesthetics
      df_points <- data[grepl("point_", names(data))]

      # uncomment following line to switch off carrying over data
      #df_points <- data.frame()

      if (ncol(df_points) == 0) {
        df_points <- NULL
        df_points_dummy <- NULL
      }
      else {
        # combine additional points data into results dataframe
        df_jittered <- cbind(df_jittered, df_points)
        # make a row of dummy data to merge with the other dataframes
        df_points_dummy <- na.omit(df_points)[1, , drop = FALSE]
      }
    } else {
      df_jittered <- NULL
      df_points_dummy <- NULL
    }

    # calculate quantiles, needed for both quantile lines and ecdf
    if ((length(quantiles)==1) && (all(quantiles >= 1))) {
      if (quantiles > 1) {
        probs <- seq(0, 1, length.out = quantiles + 1)[2:quantiles]
      }
      else {
        probs <- NA
      }
    } else {
      probs <- quantiles
      probs[probs < 0 | probs > 1] <- NA
    }

    qx <- na.omit(quantile_fun(data$x, probs = probs))

    # if requested, add data frame for quantile lines
    df_quantiles <- NULL

    if (quantile_lines && length(qx) > 0) {
      qy <- densf(qx)
      df_quantiles <- data.frame(
        x = qx,
        density = qy,
        ndensity = qy / maxdens,
        datatype = "vline",
        stringsAsFactors = FALSE
      )
      if (!is.null(df_points_dummy)){
        # add in dummy points data if necessary
        df_quantiles <- data.frame(df_quantiles, as.list(df_points_dummy))
      }
    }

    # combine the quantiles and jittered points data frames into one, the non-density frame
    df_nondens <- rbind(df_quantiles, df_jittered)

    if (calc_ecdf) {
      n <- length(d$x)
      ecdf <- c(0, cumsum(d$y[1:(n-1)]*(d$x[2:n]-d$x[1:(n-1)])))
      ecdf_fun <- approxfun(d$x, ecdf, rule = 2)
      ntile <- findInterval(d$x, qx, left.open = TRUE) + 1 # if make changes here, make them also below

      if (!is.null(df_nondens)) {
        # we add data for ecdf and quantiles back to all other data points
        df_nondens <- data.frame(
          df_nondens,
          ecdf = ecdf_fun(df_nondens$x),
          quantile = findInterval(df_nondens$x, qx, left.open = TRUE) + 1
        )
      }

      df_density <- data.frame(
        x = d$x,
        density = d$y,
        ndensity = d$y / maxdens,
        ecdf = ecdf,
        quantile = ntile,
        datatype = "ridgeline",
        stringsAsFactors = FALSE
      )
    }
    else {
      df_density <- data.frame(
        x = d$x,
        density = d$y,
        ndensity = d$y / maxdens,
        datatype = "ridgeline",
        stringsAsFactors = FALSE
      )
    }

    if (!is.null(df_points_dummy)){
      # add in dummy points data if necessary
      df_density <- data.frame(df_density, as.list(df_points_dummy))
    }

    # now combine everything and turn quantiles into factor
    df_final <- rbind(df_density, df_nondens)
    if ("quantile" %in% names(df_final)) {
      df_final$quantile <- factor(df_final$quantile)
    }

    df_final
  }
)


################################################################################
# user defined functions for argument "quantile_fun" in geom_density_ridges(2)
## 95% HPD Credible Intervals and posterior mode
HPDpstMode_fun <- function(x, probs){
  qx <- c(HPDinterval(as.mcmc(x), prob = probs)[1, , drop = TRUE],
    posterior.mode(as.mcmc(x)))
  qx <- qx[c(1,3,2)]
  names(qx)[2] <- "mode"
 qx
}    
## Just 95% HPD Credible Intervals
HPD_fun <- function(x, probs){
  HPDinterval(as.mcmc(x), prob = probs)[1, , drop = TRUE]
}

