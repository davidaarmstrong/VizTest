## TODO: fix simulation method to take quantile (credible interval) or HPD to calculate the overlaps. 

globalVariables(c("bound_end", "bound_start", "est", "label", "lwr", "stim_end", "stim_start", "upr", "vbl"))

#' Calculate Correspondence Between Pairwise Test and CI Overlaps
#'
#' @param obj A model object (or any object) where `coef()` and `vcov()` return estimates of coefficients and sampling variability.
#' @param test_level The type I error rate of the pairwise tests.
#' @param range_levels The range of confidence levels to try.
#' @param level_increment Step size of increase between the values of `range_levels`.
#' @param adjust Multiplicity adjustment to use when calculating the p-values for normal theory pairwise tests.
#' @param cifun For simulation results, the method used to calculate the confidence/credible interval either "quantile" (default) or "hdi" for highest density region. 
#' @param include_intercept Logical indicating whether the intercept should be included in the tests, defaults to `FALSE`.
#' @param include_zero Should univariate tests at zero be included, defaults to `TRUE`.
#' @param ... Other arguments, currently not implemented.
#' @export
#' @importFrom stats coef vcov qt pt p.adjust
#' @importFrom utils combn
#'
viztest <- function(obj,
                    test_level = 0.05,
                    range_levels = c(.25, .99),
                    level_increment = 0.01,
                    adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                    cifun = c("quantile", "hdi"), 
                    include_intercept = FALSE,
                    include_zero = TRUE,
                    ...){
  UseMethod("viztest")
}

#' Default method for viztest function
#'
#' @param obj A model object (or any object) where `coef()` and `vcov()` return estimates of coefficients and sampling variability.
#' @param test_level The type I error rate of the pairwise tests.
#' @param range_levels The range of confidence levels to try.
#' @param level_increment Step size of increase between the values of `range_levels`.
#' @param adjust Multiplicity adjustment to use when calculating the p-values for the pairwise tests.
#' @param include_intercept Logical indicating whether the intercept should be included in the tests, defaults to `FALSE`.
#' @param include_zero Should univariate tests at zero be included, defaults to `TRUE`.
#' @param ... Other arguments, currently not implemented.
#' @method viztest default
#' @importFrom dplyr tibble bind_rows
#' @export
viztest.default <- function(obj,
                     test_level = 0.05,
                     range_levels = c(.25, .99),
                     level_increment = 0.01,
                     adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                     include_intercept = FALSE,
                     include_zero = TRUE,
                     ...){
  adj <- match.arg(adjust)
  lev_seq <- seq(range_levels[1], range_levels[2], by=level_increment)
  resdf <- Inf
  if(!is.null(obj$df.residual) & inherits(obj, "lm")){
    resdf <- obj$df.residual
  }
  bhat <- coef(obj)
  if(is.null(names(bhat)))names(bhat) <- 1:length(bhat)
  V <- vcov(obj)
  if(!include_intercept){
    w_int <- grep("ntercept", names(bhat))
    if(length(w_int) > 0){
      bhat <- bhat[-w_int]
      V <- V[-w_int, -w_int]
    }
  }
  if(include_zero)bhat <- c(bhat, zero=0)
  o <- order(bhat, decreasing=TRUE)
  bhat <- bhat[o]
  if(include_zero)V <- cbind(rbind(V, 0), 0)
  V <- V[o, o]

  combs <- combn(length(bhat), 2)
  D <- matrix(0, nrow=length(bhat), ncol=ncol(combs))
  D[cbind(combs[1,], 1:ncol(combs))] <- 1
  D[cbind(combs[2,], 1:ncol(combs))] <- -1
  diffs <- bhat %*% D
  se_diffs <- sqrt(diag(t(D) %*% V %*% D))
  p_diff <- pt(diffs/se_diffs, resdf, lower.tail=FALSE)
  p_diff <- p.adjust(p_diff, method=adj)
  s <- p_diff < test_level
  L <- sapply(lev_seq, \(l)bhat - qt(1-(1-l)/2, resdf)*sqrt(diag(V)))
  U <- sapply(lev_seq, \(l)bhat + qt(1-(1-l)/2, resdf)*sqrt(diag(V)))
  s_star <- L[combs[1,], ] >= U[combs[2,], ]
  smat <- array(s, dim=dim(s_star))
  if("zero" %in% rownames(L)){
    w <- which(rownames(L) == "zero")
    new_L <- L[-w, ]
    new_U <- U[-w, ]
    out <- which(combs[1,] == w | combs[2,] == w)
    new_combs <- combs[, -out]
    new_combs[which(new_combs > w, arr.ind = TRUE)] <- new_combs[which(new_combs > w, arr.ind = TRUE)] - 1
    new_s <- s[-out]
  }else{
    new_L <- L
    new_U <- U
    new_combs <- combs
    new_s <- s
  }
  diff_sig <- new_L[new_combs[1,new_s],, drop=FALSE] - new_U[new_combs[2,new_s],, drop=FALSE]
  diff_insig <- new_U[new_combs[2,!new_s],, drop=FALSE] - new_L[new_combs[1,!new_s],, drop=FALSE]
  diff_sig[which(diff_sig <= 0, arr.ind=TRUE)] <- NA
  diff_insig[which(diff_insig <= 0, arr.ind=TRUE)] <- NA
  d_sig <- suppressWarnings(apply(diff_sig, 2, min, na.rm=TRUE))
  d_insig <- suppressWarnings(apply(diff_insig, 2, min, na.rm=TRUE))
  d_sig <- ifelse(is.finite(d_sig), d_sig, 0)
  d_insig <- ifelse(is.finite(d_insig), d_insig, 0)
  easiness <- d_sig*d_insig
  res <- data.frame(level = lev_seq,
                    psame = apply(s_star, 2, \(x)mean(x == s)),
                    pdiff = mean(s),
                    easy = easiness)
  est_data <- tibble(vbl = names(bhat), 
                     est = bhat, 
                     se = sqrt(diag(V)))
  res <- list(tab = res,
              pw_test = s,
              ci_tests = s_star,
              combs = combs,
              param_names = names(bhat),
              L = L,
              U = U, 
              est = est_data)
  class(res) <- "viztest"
  return(res)
}
#' Simulation method for viztest function
#' 
#' @param obj A model object (or any object) where `coef()` and `vcov()` return estimates of coefficients and sampling variability.
#' @param test_level The type I error rate of the pairwise tests.
#' @param range_levels The range of confidence levels to try.
#' @param level_increment Step size of increase between the values of `range_levels`.
#' @param adjust Multiplicity adjustment to use when calculating the p-values for the pairwise tests.
#' @param cifun For simulation results, the method used to calculate the confidence/credible interval either "quantile" (default) or "hdi" for highest density region. 
#' @param include_intercept Logical indicating whether the intercept should be included in the tests, defaults to `FALSE`.
#' @param include_zero Should univariate tests at zero be included, defaults to `TRUE`.
#' @param ... Other arguments, currently not implemented.
#' @importFrom stats quantile
#' @importFrom HDInterval hdi
#' @method viztest vtsim
#' @export
viztest.vtsim <- function(obj,
                          test_level = 0.05,
                          range_levels = c(.25, .99),
                          level_increment = 0.01,
                          adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                          cifun = c("quantile", "hdi"), 
                          include_intercept = FALSE,
                          include_zero = TRUE,
                          ...){
  cif <- match.arg(cifun)
  est <- obj$est
  if(include_zero){
    est <- cbind(est, zero=rep(0, nrow(est)))
  }
  cm <- colMeans(est, na.rm=TRUE)
  o <- order(cm, decreasing=TRUE)
  est <- est[,o]
  combs <- combn(ncol(est), 2)
  D <- matrix(0, nrow=ncol(est), ncol=ncol(combs))
  D[cbind(combs[1,], 1:ncol(combs))] <- 1
  D[cbind(combs[2,], 1:ncol(combs))] <- -1
  diffs <- est %*% D
  pvals <- apply(diffs, 2, \(x)mean(x > 0))
  pvals <- ifelse(pvals > .5, 1-pvals, pvals)
  s <- pvals < test_level
  lev_seq <- seq(range_levels[1], range_levels[2], by=level_increment)
  if(cif == "quantile"){
    L <- sapply(lev_seq, \(l)apply(est, 2, quantile, probs=((1-l)/2)))
    U <- sapply(lev_seq, \(l)apply(est, 2, quantile, probs=(1-(1-l)/2)))
  }else{
    LU <- lapply(lev_seq, \(l)apply(est, 2, hdi, credMass = l))
    L <- sapply(LU, \(x)x[1,])
    U <- sapply(LU, \(x)x[2,])
  }
  s_star <- L[combs[1,], , drop=FALSE] >= U[combs[2,], , drop=FALSE]
  smat <- array(s, dim=dim(s_star))
  if(include_zero){
    w <- which(apply(est, 2, \(x)all(x == 0)))
    new_L <- L[-w, ]
    new_U <- U[-w, ]
    out <- which(combs[1,] == w | combs[2,] == w)
    new_combs <- combs[, -out]
    new_combs[which(new_combs > w, arr.ind = TRUE)] <- new_combs[which(new_combs > w, arr.ind = TRUE)] - 1
    new_s <- s[-out]
  }else{
    new_L <- L
    new_U <- U
    new_combs <- combs
    new_s <- s
  }
  diff_sig <- new_L[new_combs[1,new_s],, drop=FALSE] - new_U[new_combs[2,new_s],, drop=FALSE]
  diff_insig <- new_U[new_combs[2,!new_s],, drop=FALSE] - new_L[new_combs[1,!new_s],, drop=FALSE]
  diff_sig[which(diff_sig <= 0, arr.ind=TRUE)] <- NA
  diff_insig[which(diff_insig <= 0, arr.ind=TRUE)] <- NA
  d_sig <- suppressWarnings(apply(diff_sig, 2, min, na.rm=TRUE))
  d_insig <- suppressWarnings(apply(diff_insig, 2, min, na.rm=TRUE))
  d_sig <- ifelse(is.finite(d_sig), d_sig, 0)
  d_insig <- ifelse(is.finite(d_insig), d_insig, 0)
  easiness <- d_sig*d_insig
  res <- data.frame(level = lev_seq,
                    psame = apply(s_star, 2, \(x)mean(x == s)),
                    pdiff = mean(s),
                    easy = easiness)
  res <- list(tab = res,
              pw_test = s,
              ci_tests = s_star,
              combs = combs,
              param_names = colnames(est),
              L = L,
              U = U)
  class(res) <- "viztest"
  return(res)
}

#' Print Method for viztest Objects
#'
#' Prints the results from the viztest function
#' @param x An object of class `viztest`.
#' @param best Logical indicating whether the results should be filtered to include only the best level(s) or include all levels
#' @param missed_tests Logical indicating whether the tests not represented by the optimal visual testing intervals should be displayed
#' @param level Which level should be used as the optimal one.  If `NULL`, the easiest optimal level will be used.  Easiness is measured by the sum of the overlap
#' in confidence intervals for insignificant tests plus the distance between the lower and upper bound for tests that are significant.
#' @param ... Other arguments, currently not implemented.
#' @export
#' @method print viztest
print.viztest <- function(x, ..., best=TRUE, missed_tests=TRUE, level=NULL){
  cat("\nCorrespondents of PW Tests with CI Tests\n")
  tmp <- x$tab[which(x$tab$psame == max(x$tab$psame)), ]
  lowest <- tmp[1,]
  highest <- tmp[nrow(tmp), ]
  middle <- tmp[floor(median(1:nrow(tmp))), ]
  easiest <- tmp[which.max(tmp$easy), ]
  b_levs <- bind_rows(lowest, middle, highest, easiest) %>% 
    mutate(method = c("Lowest", "Middle", "Highest", "Easiest"))
  if(best){
    print(b_levs)
  }else{
    print(x$tab)
  }
  if(missed_tests){
    if(is.null(level)){
      ## missed tests for lowest level
      l1 <- b_levs$level[1]
      w11 <- which(round(x$tab$level, 10) == round(l1, 10))
      mt1 <- data.frame(bigger = x$param_names[x$combs[1,]],
                       smaller = x$param_names[x$combs[2,]],
                       pw_test = ifelse(x$pw_test, "Sig", "Insig"),
                       ci_olap = ifelse(x$ci_tests[,w11], "No", "Yes"))
      w21 <- which((mt1$pw_test == "Sig" & mt1$ci_olap == "Yes") |
                    (mt1$pw_test == "Insig" & mt1$ci_olap == "No"))
      
      if(length(w21) > 0){
        cat("\nMissed Tests for Lowest Level (n=", length(w21), " of ", length(x$pw_test), ")\n", sep="")
        print(mt1[w21, ])
      }
      
      ## missed tests for middle level
      l2 <- b_levs$level[2]
      w12 <- which(round(x$tab$level, 10) == round(l2, 10))
      mt2 <- data.frame(bigger = x$param_names[x$combs[1,]],
                        smaller = x$param_names[x$combs[2,]],
                        pw_test = ifelse(x$pw_test, "Sig", "Insig"),
                        ci_olap = ifelse(x$ci_tests[,w12], "No", "Yes"))
      w22 <- which((mt2$pw_test == "Sig" & mt2$ci_olap == "Yes") |
                     (mt2$pw_test == "Insig" & mt2$ci_olap == "No"))
      if(length(w22) > 0){
        cat("\nMissed Tests for Lowest Level (n=", length(w22), " of ", length(x$pw_test), ")\n", sep="")
        print(mt2[w22, ])
      }
      
      ## missed tests for highest level
      l3 <- b_levs$level[3]
      w13 <- which(round(x$tab$level, 10) == round(l3, 10))
      mt3 <- data.frame(bigger = x$param_names[x$combs[1,]],
                        smaller = x$param_names[x$combs[2,]],
                        pw_test = ifelse(x$pw_test, "Sig", "Insig"),
                        ci_olap = ifelse(x$ci_tests[,w13], "No", "Yes"))
      w23 <- which((mt3$pw_test == "Sig" & mt3$ci_olap == "Yes") |
                     (mt3$pw_test == "Insig" & mt3$ci_olap == "No"))
      if(length(w23) > 0){
        cat("\nMissed Tests for Lowest Level (n=", length(w23), " of ", length(x$pw_test), ")\n", sep="")
        print(mt3[w23, ])
      }
      
      ## missed tests for easiest level
      l4 <- b_levs$level[4]
      w14 <- which(round(x$tab$level, 10) == round(l4, 10))
      mt4 <- data.frame(bigger = x$param_names[x$combs[1,]],
                        smaller = x$param_names[x$combs[2,]],
                        pw_test = ifelse(x$pw_test, "Sig", "Insig"),
                        ci_olap = ifelse(x$ci_tests[,w14], "No", "Yes"))
      w24 <- which((mt4$pw_test == "Sig" & mt4$ci_olap == "Yes") |
                     (mt4$pw_test == "Insig" & mt4$ci_olap == "No"))
            tmp <- x$tab[which(x$tab$psame == max(x$tab$psame)), ]
      
      if(length(w24) > 0){
        cat("\nMissed Tests for Lowest Level (n=", length(w24), " of ", length(x$pw_test), ")\n", sep="")
        print(mt4[w24, ])
      }
      if(length(w21) == 0){
        cat("\nAll ", length(x$pw_test), " tests properly represented for by CI overlaps.\n")
      }
    }else{
      w <- which(round(x$tab$level, 10) == round(level, 10))
      mt <- data.frame(bigger = x$param_names[x$combs[1,]],
                       smaller = x$param_names[x$combs[2,]],
                       pw_test = ifelse(x$pw_test, "Sig", "Insig"),
                       ci_olap = ifelse(x$ci_tests[,w], "No", "Yes"))
      w2 <- which((mt$pw_test == "Sig" & mt$ci_olap == "Yes") |
                    (mt$pw_test == "Insig" & mt$ci_olap == "No"))
      if(length(w2) > 0){
        cat("\nMissed Tests (n=", length(w2), " of ", length(x$pw_test), ")\n", sep="")
        print(mt[w2, ])
      }else{
        cat("\nAll ", length(x$pw_test), " tests properly represented for by CI overlaps.\n")
      }
    }
  }
}

#'Coefficient Method for vtcustom Objects.
#'
#' Returns coefficients from an object of class `vtcustom`.  This allows users to build a
#' list with estimates and a variance covariance matrix.  If the object's class is set to `vtcustom`
#' and the object has an element `coef`, this coefficient method will return the vector of estimates.  
#' For maximum usability, the entries of the coefficient vector should be named. 
#' @param object An object of class `vtcustom`. 
#' @param ... Currently not implemented
#' @method coef vtcustom
#' @export
coef.vtcustom <- function(object, ...){
  object$coef
}

#' Variance-Covariance Method for vtcustom Objects.
#'
#' Returns the variance-covariance matrix from an object of class `vtcustom`.  This allows users to build a
#' list with estimates and a variance covariance matrix.  If the object's class is set to `vtcustom`
#' and the object has an element `vcov`, this coefficient method will return the variance-covariance matrix.  
#' The matrix may have row and column names, but they will have no effect on the result. 
#' @param object An object of class `vtcustom`. 
#' @param ... Currently not implemented
#' @method vcov vtcustom
#' @export
vcov.vtcustom <- function(object, ...){
  object$vcov
}

#' Make Reference Segments
#' 
#' Makes reference segments for the `plot.viztest`.  These segments run along the upper bound of a stimulus from the stimulus in question to the last comparison stimulus whose lower bound is smaller than the upper bound of the stimulus in question. 
#' @param .data Data to be used in the calculating segments.  The data should contain the variables `lwr` and `upr` for the lower and upper confidence bounds, respectively. 
#' @param vdt Visual difficulty threshold, see the details of `help(plot.viztest)` for an explanation of the parameter. 
#' @param ... Other arguments to be passed down.  Currently not implemented. 
#' @returns A data frame with variables `stim_start` and `stim_end` identifying the starting and ending stimuli for the segments as well as `bound_start` and `bound_end` which give the relevant value of the reference bound.  There is also a variable called `ambiguous` which indicates if any comparisons with that bound are ambiguous.
#' @importFrom dplyr filter
#' @export
make_segs <- function(.data, vdt = .02, ...){
  segs <- NULL
  for(i in 1:nrow(.data)){
    if(any(.data$lwr[i:nrow(.data)] < .data$upr[i])){
      segs <- rbind(segs, data.frame(stim_start=i, stim_end=(i-1) + max(which(.data$lwr[i:nrow(.data)] <= .data$upr[i])), bound_start=.data$upr[i], bound_end=.data$upr[i]))
    }
  }
  rg <- max(.data$upr, na.rm=TRUE) - min(.data$lwr, na.rm=TRUE)
  amb_thresh <- vdt*rg
  segs$ambiguous <- FALSE
  for(i in 1:(nrow(segs)-1)){
    if(any(abs(.data$lwr[(segs$stim_start[i]+1):nrow(.data)] - segs$bound_start[i]) < amb_thresh)){
      segs$ambiguous[i] <- TRUE
      segs$stim_end[i] <- max(which((abs(.data$lwr[(segs$stim_start[i]+1):nrow(.data)] - segs$bound_start[i]) < amb_thresh)==T))+i
    }
  } 
  segs
}


#' Plot Method for viztest Objects
#' 
#' Plots the output of viztest objects with optional reference lines 
#' @param x Object to be plotted, should be of class `viztest`
#' @param ref_lines Reference lines to be plotted - one of "all", "ambiguous", "none".  This could also be a vector of stimulus names to plot - they should be the same as the names of the estimates in `x$est`. See details for explanation. 
#' @param viz_diff_thresh Threshold for identifying visual difficulty, see details. 
#' @param make_plot Logical indicating whether the plot should be constructed or the data returned. 
#' @param level Level at which to plot the estimates.  Accepts both numeric entries or one of "ce", "max", "min", "median" - defaults to "ce", the cognitively easiest level.  
#' @param trans A function to transform the estimates and their confidence intervals like `plogis`.
#' @param ... Other arguments passed down.  Currently not implemented.
#' @details The `ref_lines` argument identifies what reference lines will be plotted in the figure.  For any particular stimulus, the reference lines run along the upper bound of the stimulus from the stimulus location to the most distant stimulus with overlapping confidence intervals.  
#' When `ref_lines = "all"`, all lines are plotted, though in displays with many stimuli, this can make for a messy graph.  When `"ref_lines = ambiguous"` is specified, then only the ones that help discriminate in cases where the result might be visually difficult to discern are plotted. 
#' A comparison is determined to be visually difficult if the upper bound of the stimulus in question is within `viz_diff_thresh` times the difference between the smallest lower bound and the largest upper bound.  If `ref_lines = "non"`, then none of the reference lines are plotted. 
#' Alternatively, you can specify the names of stimuli whose reference lines will be plotted.  These should be the same as the names in the data.  The `viztest()` function returns an object `est`, which contains the data that are used as input to this function.  The variable `vbl` in 
#' The `est` data frame contains the stimulus names. 
#' @method plot viztest
#' @importFrom dplyr left_join arrange `%>%` join_by
#' @importFrom ggplot2 ggplot geom_pointrange geom_segment aes
#' @export
plot.viztest <- function(x, ..., ref_lines="none", viz_diff_thresh = .02, make_plot=TRUE, level=c("ce","max","min","median"),trans=I){
  inp <- x$est
  tmp <- x$tab[which(x$tab$psame == max(x$tab$psame)), ]
  if(!is.numeric(level)){
    lvl <- match.arg(level)
    level <- switch(lvl,
                    "ce" = tmp[which(tmp$easy == max(tmp$easy)), ]$level,
                    "max" = tmp[which(tmp$level == max(tmp$level)), ]$level,
                    "min" = tmp[which(tmp$level == min(tmp$level)), ]$level,
                    "median" = tmp[which(round(tmp$level,2) == round(median(tmp$level),2)), ]$level)
  }
  w <- which(round(level, 10) == round(x$tab$level, 10))
  if(length(w) == 0)stop("level must be one in x$tab$level or one of ce, max, min, or median.\n")
  if(!(level %in% tmp$level))warning("chosen level outside of range of maximally representing CI overlaps. Visual tests may not be faithfull to pairwise test results!!!")
  inp$lwr <- x$L[,w]
  inp$upr <- x$U[,w]
  inp <- inp %>% arrange(est)
  inp <- inp %>% filter(vbl != "zero")
  segs <- make_segs(inp, vdt=viz_diff_thresh)
  segs$vbl <- rownames(segs)
  inp$label <- factor(1:nrow(inp), labels=inp$vbl)
  inp <- left_join(inp, segs, by=join_by(vbl))
  cols <- c("est","se","lwr","upr","bound_start","bound_end")
  inp[,cols] <- apply(inp[,cols],2,trans)
  if(any(inp$vbl == "zero"))inp <- inp[-which(inp$vbl == "zero"), ]
  if(!make_plot){
    res <- inp
  }else{
    g <- ggplot(inp) + 
      geom_pointrange(aes(x=est, xmin = lwr, xmax=upr, y=label), size=.01) 
    if("all" %in% ref_lines){
      g <- g + geom_segment(aes(x=bound_start, xend=bound_end, y=stim_start, yend=stim_end), color="gray75", linetype=3)
    }  
    if( "ambiguous" %in% ref_lines){
      g <- g + geom_segment(data = inp[which(inp$ambiguous), ], aes(x=bound_start, xend=bound_end, y=stim_start, yend=stim_end), color="gray75", linetype=3)
    }
    if(!any(c("ambiguous", "all", "none") %in% ref_lines)){
      incl <- which(ref_lines %in% inp$vbl)
      if(length(incl) == 0)stop("ref_lines should either be one of (all, ambiguous, or none) or a vector of names consistent with x$est$vbl.\n")
      g <- g + geom_segment(data = inp[incl, ], aes(x=bound_start, xend=bound_end, y=stim_start, yend=stim_end), color="gray75", linetype=3)
    }
    res <- g
  }
  return(res)
}
