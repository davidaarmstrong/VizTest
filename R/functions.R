#' Calculate Correspondence Between Pairwise Test and CI Overlaps
#'
#' @param obj A model object (or any object) where `coef()` and `vcov()` return estimates of coefficients and sampling variability.
#' @param test_level The type I error rate of the pairwise tests.
#' @param range_levels The range of confidence levels to try.
#' @param level_increment Step size of increase between the values of `range_levels`.
#' @param adjust Multiplicity adjustment to use when calculating the p-values for the pairwise tests.
#' @param include_intercept Logical indicating whether the intercept should be included in the tests, defaults to `FALSE`.
#' @param include_zero Should univariate tests at zero be included, defaults to `TRUE`.
#' @param easy_thresh Proportion of distance between the biggest upper bound and the smallest lower bound to use
#' as the threshold in the calculation of "easiness".
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
                            include_intercept = FALSE,
                            include_zero = TRUE,
                            easy_thresh=.05,
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
#' @param easy_thresh Proportion of distance between the biggest upper bound and the smallest lower bound to use
#' as the threshold in the calculation of "easiness".
#' @param ... Other arguments, currently not implemented.
#' @method viztest default
#' @export
viztest.default <- function(obj,
                     test_level = 0.05,
                     range_levels = c(.25, .99),
                     level_increment = 0.01,
                     adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                     include_intercept = FALSE,
                     include_zero = TRUE,
                     easy_thresh=.05,
                     ...){
  adj <- match.arg(adjust)
  lev_seq <- seq(range_levels[1], range_levels[2], by=level_increment)
  resdf <- Inf
  if(!is.null(obj$df.residual) & inherits(obj, "lm")){
    resdf <- obj$df.residual
  }
  bhat <- coef(obj)
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
  delta <- apply(U, 2, max) - apply(L, 2, min)
  thresh <- easy_thresh*delta
  s_star <- L[combs[1,], ] >= U[combs[2,], ]
  smat <- array(s, dim=dim(s_star))
  easiness <- sapply(1:ncol(L), \(m){
    e <- ifelse(s, L[combs[1,],m] - U[combs[2,],m], U[combs[2,],m] - L[combs[1,],m])
    e <- ifelse(abs(e) > thresh[m], sign(e)*thresh[m], e)
    sum(e)
  } )
  res <- data.frame(level = lev_seq,
                    psame = apply(s_star, 2, \(x)mean(x == s)),
                    pdiff = mean(s),
                    easy = easiness)
  res <- list(tab = res,
              pw_test = s,
              ci_tests = s_star,
              combs = combs,
              param_names = names(bhat),
              thresh = thresh,
              L = L,
              U = U)
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
#' @param include_intercept Logical indicating whether the intercept should be included in the tests, defaults to `FALSE`.
#' @param include_zero Should univariate tests at zero be included, defaults to `TRUE`.
#' @param easy_thresh Proportion of distance between the biggest upper bound and the smallest lower bound to use
#' as the threshold in the calculation of "easiness".
#' @param ... Other arguments, currently not implemented.
#' @importFrom stats quantile
#' @method viztest vtsim
#' @export
viztest.vtsim <- function(obj,
                          test_level = 0.05,
                          range_levels = c(.25, .99),
                          level_increment = 0.01,
                          adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                          include_intercept = FALSE,
                          include_zero = TRUE,
                          easy_thresh=.05,
                          ...){
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
  L <- sapply(lev_seq, \(l)apply(est, 2, quantile, probs=((1-l)/2)))
  U <- sapply(lev_seq, \(l)apply(est, 2, quantile, probs=(1-(1-l)/2)))
  delta <- apply(U, 2, max) - apply(L, 2, min)
  thresh <- easy_thresh*delta
  s_star <- L[combs[1,], ] >= U[combs[2,], ]
  smat <- array(s, dim=dim(s_star))
  easiness <- sapply(1:ncol(L), \(m){
    e <- ifelse(s, L[combs[1,],m] - U[combs[2,],m], U[combs[2,],m] - L[combs[1,],m])
    e <- ifelse(abs(e) > thresh[m], sign(e)*thresh[m], e)
    sum(e)
  } )
  res <- data.frame(level = lev_seq,
                    psame = apply(s_star, 2, \(x)mean(x == s)),
                    pdiff = mean(s),
                    easy = easiness)
  res <- list(tab = res,
              pw_test = s,
              ci_tests = s_star,
              combs = combs,
              param_names = colnames(est),
              thresh = thresh,
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
  if(best){
    print(x$tab[which(x$tab$psame == max(x$tab$psame)), ])
  }else{
    print(x$tab)
  }
  if(missed_tests){
    if(is.null(level)){
      tmp <- x$tab[which(x$tab$psame == max(x$tab$psame)), ]
      tmp <- tmp[which(tmp$easy == max(tmp$easy)), ]
      level <- tmp[1, "level"]
    }
    w <- which(round(x$tab$level, 10) == level)
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

