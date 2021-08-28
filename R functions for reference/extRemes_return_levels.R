> extRemes::rl.fevd
function (x, period = 100) 
{
  pmods <- x$par.models
  if (is.fixedfevd(x)) {
    p <- x$results$par
    if (any(names(p) == "log.scale")) {
      id <- names(p) == "log.scale"
      p[id] <- exp(p[id])
      names(p)[id] <- "scale"
    }
    if (!is.null(x$threshold)) 
      lam <- mean(datagrabber(x)[, 1] > x$threshold)
    else lam <- NULL
    if ("shape" %in% names(p)) {
      return(rlevd(period = period, loc = p["location"], 
                   scale = p["scale"], shape = p["shape"], threshold = x$threshold, 
                   type = x$type, npy = x$npy, rate = lam))
    }
    else return(rlevd(period = period, loc = p["location"], 
                      scale = p["scale"], shape = 0, threshold = x$threshold, 
                      type = x$type, npy = x$npy, rate = lam))
  }
  else return(erlevd(x, period = period))
}
<bytecode: 0x7f91aacc47c8>
  <environment: namespace:extRemes>

  
  
##################################################### 

> extRemes::rlevd
function (period, loc = 0, scale = 1, shape = 0, threshold = 0, 
          type = c("GEV", "GP", "PP", "Gumbel", "Frechet", "Weibull", 
                   "Exponential", "Beta", "Pareto"), npy = 365.25, rate = 0.01) 
{
  if (any(period <= 1)) 
    stop("rlevd: invalid period argument.  Must be greater than 1.")
  type <- match.arg(type)
  type <- tolower(type)
  if (missing(loc)) 
    loc <- 0
  else if (is.null(loc)) 
    loc <- 0
  if (is.element(type, c("gumbel", "weibull", "frechet"))) {
    if (type == "gumbel" && shape != 0) {
      warning("rlevd: shape is not zero, but type is Gumbel.  Re-setting shape parameter to zero.")
      shape <- 0
      type <- "gev"
    }
    else if (type == "gumbel") 
      type <- "gev"
    else if (type == "frechet" && shape <= 0) {
      if (shape == 0) {
        warning("rlevd: shape is zero, but type is Frechet!  Re-setting type to Gumbel.")
        shape <- 0
      }
      else {
        warning("rlevd: type is Frechet, but shape < 0.  Negating shape to force it to be Frechet.")
        shape <- -shape
      }
      type <- "gev"
    }
    else if (type == "frechet") 
      type <- "gev"
    else if (type == "weibull" && shape >= 0) {
      if (shape == 0) {
        warning("rlevd: shape is zero, but type is Weibull!  Re-setting type to Gumbel.")
        shape <- 0
      }
      else {
        warning("rlevd: type is Weibull, but shape > 0.  Negating shape to force it to be Weibull.")
        shape <- -shape
      }
      type <- "gev"
    }
    else if (type == "weibull") 
      type <- "gev"
  }
  if (is.element(type, c("beta", "pareto", "exponential"))) {
    if (type == "exponential" && shape != 0) {
      warning("rlevd: shape is not zero, but type is Exponential.  Re-setting shape parameter to zero.")
      shape <- 0
      type <- "gp"
    }
    else if (type == "exponential") 
      type <- "gp"
    else if (type == "beta" && shape >= 0) {
      if (shape == 0) {
        warning("rlevd: shape is zero, but type is Beta!  Re-setting type to Exponential.")
        shape <- 0
      }
      else {
        warning("rlevd: type is Beta, but shape > 0.  Negating shape to force it to be Beta.")
        shape <- -shape
      }
      type <- "gp"
    }
    else if (type == "beta") 
      type <- "gp"
    else if (type == "pareto" && shape <= 0) {
      if (shape == 0) {
        warning("rlevd: shape is zero, but type is Pareto!  Re-setting type to Exponential.")
        shape <- 0
      }
      else {
        warning("rlevd: type is Pareto, but shape < 0.  Negating shape to force it to be Pareto.")
        shape <- -shape
      }
      type <- "gp"
    }
    else if (type == "pareto") 
      type <- "gp"
  }
  if (is.element(type, c("gev", "pp"))) {
    p <- 1 - 1/period
    res <- qevd(p = p, loc = loc, scale = scale, shape = shape, 
                type = "GEV")
  }
  else if (type == "gp") {
    m <- period * npy * rate
    if (shape == 0) 
      res <- threshold + scale * log(m)
    else res <- threshold + (scale/shape) * (m^shape - 1)
  }
  names(res) <- as.character(period)
  return(res)
}
<bytecode: 0x7f91aad51408>
  <environment: namespace:extRemes>