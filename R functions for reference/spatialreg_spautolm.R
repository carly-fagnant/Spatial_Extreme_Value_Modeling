
> spatialreg::spautolm
function (formula, data = list(), listw, weights, na.action, 
          family = "SAR", method = "eigen", verbose = NULL, trs = NULL, 
          interval = NULL, zero.policy = NULL, tol.solve = .Machine$double.eps, 
          llprof = NULL, control = list()) 
{
  timings <- list()
  .ptime_start <- proc.time()
  con <- list(tol.opt = .Machine$double.eps^(2/3), fdHess = NULL, 
              optimHess = FALSE, optimHessMethod = "optimHess", Imult = 2, 
              cheb_q = 5, MC_p = 16, MC_m = 30, super = NULL, spamPivot = "MMD", 
              in_coef = 0.1, type = "MC", correct = TRUE, trunc = TRUE, 
              SE_method = "LU", nrho = 200, interpn = 2000, small_asy = TRUE, 
              small = 1500, SElndet = NULL, LU_order = FALSE, pre_eig = NULL)
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  if (!inherits(listw, "listw")) 
    stop("No neighbourhood list")
  if (is.null(verbose)) 
    verbose <- get("verbose", envir = .spatialregOptions)
  stopifnot(is.logical(verbose))
  if (is.null(zero.policy)) 
    zero.policy <- get("zeroPolicy", envir = .spatialregOptions)
  stopifnot(is.logical(zero.policy))
  if (family == "SMA" && method != "eigen") 
    stop("SMA only for eigen method")
  if (method == "spam" || method == "spam_update") 
    stop("spam not supported as method")
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "na.action"), 
             names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  na.act <- attr(mf, "na.action")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy = zero.policy)
    if (!is.null(con$pre_eig)) {
      warning("NAs found, precomputed eigenvalues ignored")
      con$pre_eig <- NULL
    }
  }
  Y <- model.extract(mf, "response")
  if (any(is.na(Y))) 
    stop("NAs in dependent variable")
  X <- model.matrix(mt, mf)
  if (any(is.na(X))) 
    stop("NAs in independent variable")
  n <- nrow(X)
  if (n != length(listw$neighbours)) 
    stop("Input data and neighbourhood list have different dimensions")
  weights <- as.vector(model.extract(mf, "weights"))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (is.null(weights)) 
    weights <- rep(as.numeric(1), n)
  if (any(is.na(weights))) 
    stop("NAs in weights")
  if (any(weights < 0)) 
    stop("negative weights")
  lm.base <- lm(Y ~ X - 1, weights = weights)
  aliased <- is.na(coefficients(lm.base))
  cn <- names(aliased)
  names(aliased) <- substr(cn, 2, nchar(cn))
  if (any(aliased)) {
    nacoef <- which(aliased)
    X <- X[, -nacoef]
  }
  can.sim <- FALSE
  if (listw$style %in% c("W", "S")) 
    can.sim <- can.be.simmed(listw)
  sum_lw <- sum(log(weights))
  env <- new.env()
  assign("Y", Y, envir = env)
  assign("X", X, envir = env)
  assign("n", n, envir = env)
  assign("weights", weights, envir = env)
  assign("can.sim", can.sim, envir = env)
  assign("family", family, envir = env)
  assign("method", method, envir = env)
  assign("verbose", verbose, envir = env)
  assign("listw", listw, envir = env)
  assign("sum_lw", sum_lw, envir = env)
  W <- as(listw, "CsparseMatrix")
  if (family == "CAR") 
    if (!isTRUE(all.equal(W, t(W)))) 
      warning("Non-symmetric spatial weights in CAR model")
  assign("W", W, envir = env)
  I <- as_dsCMatrix_I(n)
  assign("I", I, envir = env)
  Sweights <- as(as(Diagonal(x = weights), "symmetricMatrix"), 
                 "CsparseMatrix")
  assign("Sweights", Sweights, envir = env)
  timings[["set_up"]] <- proc.time() - .ptime_start
  .ptime_start <- proc.time()
  if (verbose) 
    cat(paste("\nJacobian calculated using "))
  interval <- jacobianSetup(method, env, con, pre_eig = con$pre_eig, 
                            trs = trs, interval = interval)
  assign("interval", interval, envir = env)
  if (family == "SMA") 
    interval <- -rev(interval)
  nm <- paste(method, "set_up", sep = "_")
  timings[[nm]] <- proc.time() - .ptime_start
  .ptime_start <- proc.time()
  if (!is.null(llprof)) {
    if (length(llprof) == 1L) 
      llprof <- seq(interval[1], interval[2], length.out = llprof)
    ll_prof <- numeric(length(llprof))
    for (i in seq(along = llprof)) ll_prof[i] <- .opt.fit(llprof[i], 
                                                          env = env, tol.solve = tol.solve)
    nm <- paste(method, "profile", sep = "_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()
  }
  opt <- optimize(.opt.fit, interval = interval, maximum = TRUE, 
                  tol = con$tol.opt, env = env, tol.solve = tol.solve)
  lambda <- opt$maximum
  if (isTRUE(all.equal(lambda, interval[1])) || isTRUE(all.equal(lambda, 
                                                                 interval[2]))) 
    warning("lambda on interval bound - results should not be used")
  names(lambda) <- "lambda"
  LL <- opt$objective
  nm <- paste(method, "opt", sep = "_")
  timings[[nm]] <- proc.time() - .ptime_start
  .ptime_start <- proc.time()
  fit <- .SPAR.fit(lambda = lambda, env, out = TRUE, tol.solve = tol.solve)
  fit$signal_trend <- drop(X %*% fit$coefficients)
  fit$signal_stochastic <- drop(lambda * W %*% (Y - fit$signal_trend))
  fit$fitted.values <- fit$signal_trend + fit$signal_stochastic
  fit$residuals <- drop(Y - fit$fitted.values)
  LL0 <- .opt.fit(lambda = 0, env, tol.solve = tol.solve)
  LLNullLlm <- logLik(lm(Y ~ 1, weights = weights))
  nm <- paste(method, "output", sep = "_")
  timings[[nm]] <- proc.time() - .ptime_start
  .ptime_start <- proc.time()
  do_asy <- FALSE
  if (is.null(con$fdHess)) {
    con$fdHess <- !do_asy
    fdHess <- NULL
  }
  stopifnot(is.logical(con$fdHess))
  lambda.se <- NULL
  if (con$fdHess) {
    coefs <- c(lambda, fit$coefficients)
    fdHess <- getVcovmat(coefs, env, tol.solve = tol.solve, 
                         optim = con$optimHess, optimM = con$optimHessMethod)
    lambda.se <- sqrt(fdHess[1, 1])
  }
  timings[["fdHess"]] <- proc.time() - .ptime_start
  rm(env)
  GC <- gc()
  res <- list(fit = fit, lambda = lambda, LL = LL, LL0 = LL0, 
              call = match.call(), parameters = (ncol(X) + 2), aliased = aliased, 
              method = method, family = family, zero.policy = zero.policy, 
              weights = weights, interval = interval, trs = trs, timings = do.call("rbind", 
                                                                                   timings)[, c(1, 3)], LLNullLlm = LLNullLlm, fdHess = fdHess, 
              lambda.se = lambda.se, X = X, Y = Y)
  if (!is.null(na.act)) 
    res$na.action <- na.act
  if (is.null(llprof)) 
    res$llprof <- llprof
  else {
    res$llprof <- list(lambda = llprof, ll = ll_prof)
  }
  if (zero.policy) {
    zero.regs <- attr(listw$neighbours, "region.id")[which(card(listw$neighbours) == 
                                                             0)]
    if (length(zero.regs) > 0L) 
      attr(res, "zero.regs") <- zero.regs
  }
  class(res) <- "Spautolm"
  res
}
<bytecode: 0x7f918ef7f2a8>
  <environment: namespace:spatialreg>