globalVariables(c("bound_end", "bound_start", "est", "label", "lwr", "stim_end", "stim_start", "upr", "vbl", "i", "j", "vij", "zb", "smaller", "larger", "ll", "us", "olap", "ambig", "lwr_add", "upr_add", "x", "hjust"))

#' Calculate Correspondence Between Pairwise Test and CI Overlaps
#' 
#' @description `viztest()` does a grid search over `range_levels` to find the confidence level(s) such that the (non-)overlaps in 
#' confidence intervals corresponds as closely as possible with the results of pairwise tests.  To the extent that 
#' a level is found that accounts for all pairwise tests, confidence bounds at this level can be added to coefficient or marginal 
#' effects plots to enable readers to reliably identify estimates that are statistically different from each other. 
#' 
#' @details The algorithm first calculates results of a set of pairwise tests. For objects with estimates and a variance-covariance matrix, 
#' normal theory tests are calculated.  Optionally, these tests can be subjected to a multiplicity adjustment.  In the case of simulation results, 
#' something akin to p-values are calculated by identifying the probability that one estimate is larger than another.  To mimic the way we use p-values 
#' in the frequentist case, we subtract the probability of difference from 1, such that smaller values indicate more confidence in the difference.  
#' The algorithm then performs a grid search over `range_levels` at increments of `level_increment`.  For each candidate level, the 
#' confidence intervals for all parameters are calculated.  For each pair of estimates, it identifies whether the confidence intervals 
#' (or credible intervals if the input is a matrix of Bayesian simulation draws) overlaps.  For each candidate level, it calculates the proportion of times where
#' differences are significant/credible and confidence/credible intervals do not overlap or differences are not significant/credible and the intervals do overlap.  
#' The main idea is to find the level(s) such that the (non-)overlaps perfectly correspond with whether the differences are significant.  
#' 
#' If such a level can be found, a visual inspection of confidence or credible intervals at that level will identify whether a pair of estimates is 
#' statistically different or not.  
#' 
#' While most of the parameters are straightforward, the `sig_diffs` argument must be specified such that the stimuli are in order from highest to lowest.  This is most 
#' easily done by using `make_diff_template()` to identify the appropriate order of the comparisons.  
#'
#' @param obj A model object (or any object) where `coef()` and `vcov()` return estimates of coefficients and sampling variability.
#' @param test_level The type I error rate of the pairwise tests.
#' @param range_levels The range of confidence levels to try.
#' @param level_increment Step size of increase between the values of `range_levels`.
#' @param adjust Multiplicity adjustment to use when calculating the p-values for normal theory pairwise tests.
#' @param cifun For simulation results, the method used to calculate the confidence/credible interval either "quantile" (default) or "hdi" for highest density region. 
#' @param include_intercept Logical indicating whether the intercept should be included in the tests, defaults to `FALSE`.
#' @param include_zero Should univariate tests at zero be included, defaults to `TRUE`.
#' @param sig_diffs An optional vector of values identify whether each pair of values is statistically different (1) or not (0).  See Details for more information on specifying this value; there is some added complexity here. 
#' @param tol Tolerance for evaluation of symmetry and positive definiteness. 
#' @param ... Other arguments, currently not implemented.
#' @export
#' 
#' @references David A. Armstrong II and William Poirier. "Decoupling Visualization and Testing when Presenting Confidence Intervals" Political Analysis <doi:10.1017/pan.2024.24>.
#' @importFrom stats coef vcov qt pt p.adjust
#' @importFrom utils combn
#' @importFrom multcomp glht adjusted
#' @importFrom dplyr left_join 
#' @returns A list (of class "viztest") with the following elements: 
#' 1. tab: a data frame with results from the grid search.  The data frame has four variables: `level` - is the confidence level used in the grid search; `psame` - the proportion of (non-)overlaps that match the 
#' normal theory tests; `pdiff` - the proportion of pairwise tests that are statistically significant; `easy` - the ease with which the comparisons are made. 
#' 2. pw_tests: A logical vector indicating which tests are significantly significant. 
#' 3. ci_tests: A logical vector indicating whether the confidence intervals are disjoint (`TRUE`) or overlap (`FALSE`). 
#' 4. combs: The pairwise combinations of stimuli used in the test.  Note, the stimuli are reordered from largest to smallest, so the numbers do not represent the position in the original ordering. 
#' 5. param_names: A vector of the names of the parameters reordered by size - largest to smallest. 
#' 6. L: The lower confidence bounds from the grid search. 
#' 7. U: The upper confidence bounds from the grid search. 
#' 8. est: A data frame with the variables `vbl` - the parameter name; `est` - the parameter estimate; `se` - the parameter standard error. 
#' 9. call: model call
#' @examples
#' data(mtcars)
#' mtcars$cyl <- as.factor(mtcars$cyl)
#' mtcars$hp <- scale(mtcars$hp)
#' mtcars$wt <- scale(mtcars$wt)
#' mod <- lm(qsec ~ hp + wt + cyl, data=mtcars)
#' viztest(mod)
#' 
viztest <- function(obj,
                    test_level = 0.05,
                    range_levels = c(.25, .99),
                    level_increment = 0.01,
                    adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                    cifun = c("quantile", "hdi"), 
                    include_intercept = FALSE,
                    include_zero = TRUE,
                    sig_diffs = NULL,
                    tol = 1e-08,
                    ...){
  UseMethod("viztest")
}

#' @method viztest default
#' @importFrom dplyr tibble bind_rows
#' @importFrom stats sd
#' @export
viztest.default <- function(obj,
                     test_level = 0.05,
                     range_levels = c(.25, .99),
                     level_increment = 0.01,
                     adjust = c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                     cifun = c("quantile", "hdi"), 
                     include_intercept = FALSE,
                     include_zero = TRUE,
                     sig_diffs = NULL,
                     tol = 1e-08,
                     ...){
  
  adj <- match.arg(adjust)
  lev_seq <- seq(range_levels[1], range_levels[2], by=level_increment)
  resdf <- Inf
  if(inherits(obj, "lm")){
    if(!is.null(obj$df.residual)){
      resdf <- obj$df.residual
    }  
  }
  bhat <- coef(obj)
  if(is.null(names(bhat)))names(bhat) <- 1:length(bhat)
  V <- vcov(obj)
  if(any(is.na(bhat))){
    na_coef <- which(is.na(bhat))
    if(inherits(V, "matrix")){
      na_var <- which(apply(V, 1, \(x)all(is.na(x))))
    }else{
      na_var <- which(is.na(V))
    }
    coef_var_agree <- try(all(na_var == na_coef), silent=TRUE)
    if(!inherits(coef_var_agree, "try-error")){
      bhat <- bhat[-na_coef]
      V <- V[-na_var, -na_var]
      message(paste0("Rank-deficient model detected - parameter name(s):", paste(names(na_coef), collapse=", "), " removing NAs from estimates/variances and proceeding.\n"))
    }else{
      stop("NAs found in estimates and/or variances, but they are not in compatible places, please solve the problem and try again.\n")
    }
  }
  if(!is_pd(V, tol=tol)){
    stop("Variance-covariance matrix is not positive definite.  Cannot proceed with viztest().\n")
  }
  if(!(length(bhat) == nrow(V) & is_sym(V, tol=tol))){
    stop("Length of coefficient vector does not match dimensions of variance-covariance matrix.  Cannot proceed with viztest().\n")
  }
  if(!include_intercept){
    w_int <- grep("ntercept", names(bhat))
    if(length(w_int) > 0){
      bhat <- bhat[-w_int]
      V <- V[-w_int, -w_int]
    }
  }
  if(include_zero){
    bhat <- c(bhat, zero=0)
    V <- cbind(rbind(V, 0), 0)
  }
  o <- order(bhat, decreasing=TRUE)
  bhat <- bhat[o]
  V <- V[o, o]

  combs <- combn(length(bhat), 2)
  if(is.null(sig_diffs)){
    D <- matrix(0, nrow=length(bhat), ncol=ncol(combs))
    D[cbind(combs[1,], 1:ncol(combs))] <- 1
    D[cbind(combs[2,], 1:ncol(combs))] <- -1
    diffs <- bhat %*% D
    se_diffs <- sqrt(diag(t(D) %*% V %*% D))
    p_diff <- 2*pt(diffs/se_diffs, resdf, lower.tail=FALSE)
    p_diff <- p.adjust(p_diff, method=adj)
    s <- p_diff < test_level
  }else{
    s <- as.logical(sig_diffs)
  }
  L <- sapply(lev_seq, \(l)bhat - qt(1-(1-l)/2, resdf)*sqrt(diag(V)))
  U <- sapply(lev_seq, \(l)bhat + qt(1-(1-l)/2, resdf)*sqrt(diag(V)))
  s_star <- L[combs[1,], ] >= U[combs[2,], ]
  smat <- array(s, dim=dim(s_star))
  # remove tests against zero to calculate "easiness"
  if("zero" %in% rownames(L)){
    w <- which(rownames(L) == "zero")
    new_L <- L[-w, ]
    new_U <- U[-w, ]
    out <- which(combs[1,] == w | combs[2,] == w)
    new_combs <- combs[, -out, drop = FALSE]
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
  easiness <- -abs(d_sig-d_insig)
  res <- data.frame(level = lev_seq,
                    psame = apply(s_star, 2, \(x)mean(x == s)),
                    pdiff = mean(s),
                    easy = easiness)
  est_data <- tibble(vbl = names(bhat), 
                     est = bhat, 
                     se = sqrt(diag(V)))
  if(adj == "none"){
    est_data$lwr_add <- bhat - qt(1-test_level/2, resdf)*sqrt(diag(V))
    est_data$upr_add <- bhat + qt(1-test_level/2, resdf)*sqrt(diag(V))
  }else{
    if("zero" %in% names(bhat)){
      tmp_bhat <- bhat
      tmp_V <- V
      zero_ind <- which(names(bhat) == "zero")
      tmp_bhat <- bhat[-zero_ind]
      tmp_V <- V[-zero_ind, -zero_ind]
    }else {
      tmp_bhat <- bhat
      tmp_V <- V
    }
    K <- diag(length(tmp_bhat))
    g <- glht(model = NULL, linfct = K, coef.=tmp_bhat, vcov.=tmp_V)
    smry <- summary(g, test=adjusted(adj))
    cis <- confint(smry)
    dat_add <- data.frame(vbl = names(tmp_bhat), 
                          lwr_add = cis$confint[,"lwr"], 
                          upr_add = cis$confint[,"upr"])
    est_data <- left_join(est_data, dat_add, by="vbl")
  }    
  res <- list(tab = res,
              pw_test = s,
              ci_tests = s_star,
              combs = combs,
              param_names = names(bhat),
              L = L,
              U = U, 
              est = est_data, 
              call = match.call(viztest, call = sys.call()))
  class(res) <- "viztest"
  attr(res, "adjusted") <- adj
  attr(res, "test_level") <- test_level 
  return(res)
}
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
                          sig_diffs = NULL,
                          tol = 1e-08, 
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
  if(is.null(sig_diffs)){
    D <- matrix(0, nrow=ncol(est), ncol=ncol(combs))
    D[cbind(combs[1,], 1:ncol(combs))] <- 1
    D[cbind(combs[2,], 1:ncol(combs))] <- -1
    diffs <- est %*% D
    pvals <- apply(diffs, 2, \(x)mean(x > 0))
    pvals <- 2*ifelse(pvals > .5, 1-pvals, pvals)
    s <- pvals < test_level
  }else{
    s <- sig_diffs
  }
  lev_seq <- seq(range_levels[1], range_levels[2], by=level_increment)
  if(cif == "quantile"){
    L <- sapply(lev_seq, \(l)apply(est, 2, quantile, probs=((1-l)/2)))
    U <- sapply(lev_seq, \(l)apply(est, 2, quantile, probs=(1-(1-l)/2)))
    tl_L <- apply(est, 2, quantile, probs=(1 - (1 - test_level/2)))
    tl_U <- apply(est, 2, quantile, probs=(1 - test_level / 2))
  }else{
    LU <- lapply(lev_seq, \(l)apply(est, 2, hdi, credMass = l))
    L <- sapply(LU, \(x)x[1,])
    U <- sapply(LU, \(x)x[2,])
    tl_LU <- apply(est, 2,  hdi, credMass = 1-test_level)
    tl_L <- tl_LU[1,]
    tl_U <- tl_LU[2,]
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
  easiness <- -abs(d_sig-d_insig)
  res <- data.frame(level = lev_seq,
                    psame = apply(s_star, 2, \(x)mean(x == s)),
                    pdiff = mean(s),
                    easy = easiness)
  cme <- colMeans(est)
  esd <- apply(est, 2, sd)
  est_data <- tibble(vbl = names(cme), 
                     est = cme, 
                     sd = esd)
  est_data$lwr_add <- tl_L
  est_data$upr_add <- tl_U
  res <- list(tab = res,
              pw_test = s,
              ci_tests = s_star,
              combs = combs,
              param_names = colnames(est),
              L = L,
              U = U, 
              est = est_data, 
              call = match.call(viztest, call = sys.call()))
  class(res) <- "viztest"
  attr(res, "adjusted") <- "none"
  attr(res, "test_level") <- test_level 
  return(res)
}

#' Print Method for viztest Objects
#'
#' @description Prints a summary of the results from the `viztest()` function.  
#' 
#' @details The results are printed in such a way that the range of optional levels is produced including the range along with two candidates for the 
#' best levels to use - middle and easiest.  
#'
#' Prints the results from the viztest function
#' @param x An object of class `viztest`.
#' @param best Logical indicating whether the results should be filtered to include only the best level(s) or include all levels
#' @param missed_tests Logical indicating whether the tests not represented by the optimal visual testing intervals should be displayed
#' @param level Which level should be used as the optimal one.  If `NULL`, the easiest optimal level will be used.  Easiness is measured by the sum of the overlap
#' in confidence intervals for insignificant tests plus the distance between the lower and upper bound for tests that are significant.
#' @param ... Other arguments, currently not implemented.
#' @returns Printed results that give the level(s) that correspond most closely with the pairwise test results.  The values returned are the smallest, 
#' largest, middle and easiest.  By default this function also reports the tests that are not captured by the (non-)overlaps in confidence intervals
#' when each different level is used. 
#' 
#' @export
#' @importFrom dplyr mutate
#' @importFrom stats median
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
        cat("\nMissed Tests for Middle Level (n=", length(w22), " of ", length(x$pw_test), ")\n", sep="")
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
        cat("\nMissed Tests for Highest Level (n=", length(w23), " of ", length(x$pw_test), ")\n", sep="")
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
        cat("\nMissed Tests for Easiest Level (n=", length(w24), " of ", length(x$pw_test), ")\n", sep="")
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

#' Make custom visual testing data
#' 
#' @description Makes custom visual testing objects that can be used as input to the `viztest()` function.  This is useful in the case
#' where `coef()` and `vcov()` do not function as expected on objects of interest, where the user wants to intervene with some 
#' modification to the usual estimates or (more likely) variance-covariance matrix or where normal theory tests may not be 
#' as useful (e.g., in the case of simulations of non-normal values).  The examples section below shows how this could be leveraged
#' to use a heteroskedasticity-consistent covariance matrix in the test rather than the one returned by `lm()`. 
#' 
#' @param estimates A vector of estimates if type is `"est_var"` and or a number of simulations by 
#' number of parameters matrix of simulated values if type is `"sim"`.  
#' @param variances In the case of independent estimates, a vector of variances of the same length 
#' as `estimates` if type is `"est_var"`.  These will be used as the diagonal elements in a variance-covariance matrix with 
#' zero covariances.  Alternatively, if type is `"est_var"`, this could be a variance-covariance matrix, with the same number 
#' of rows and columns as there are elements in the `estimates` vector.  If type is `"sim"`, variances should be `NULL`, but 
#' will be disregarded in any event. Also, note, these should be variances of the estimates (e.g., squared standard errors) and not 
#' raw variances from the data. 
#' @param type Indicates the type of input data either estimates with variances or a variance-covariance matrix or data from
#' a simulation. 
#' @param tol Tolerance for evaluation of symmetry and positive definiteness. 
#' @param ... Other arguments passed down, currently not implemented. 
#' @returns 
#' 1. If the input is a vector of parameter estimates and a variance-covariance matrix, then a list with estimates and a variance-covariance matrix of class `"vtcustom"` is returned.  In this case, the functionms `coef.vtcustom()` and `vcov.vtcustom()` are 
#' used to extract the coefficients and variance-covariance matrix in a way that will work with `viztest.default()`. 
#' 2. If the input is a matrix of simulation draws, an object of class `"vtsim"` that has a single element - the data giving the draws from the simulation is returned.  In this case, `viztest.vtsim()` does the relevant testing.  
#' @export
#' 
#' @examples
#' data(mtcars)
#' mtcars$cyl <- as.factor(mtcars$cyl)
#' mtcars$hp <- scale(mtcars$hp)
#' mtcars$wt <- scale(mtcars$wt)
#' mod <- lm(qsec ~ hp + wt + cyl, data=mtcars)
#' V <- sandwich::vcovHC(mod, "HC3")
#' vtdat <- make_vt_data(coef(mod), V)
#' viztest(vtdat, 
#'         test_level = .025, 
#'         include_intercept = FALSE, 
#'         include_zero = FALSE)
make_vt_data <- function(estimates, variances=NULL, type=c("est_var", "sim"), tol = 1e-08, ...){
  typ <- match.arg(type)
  if(typ == "est_var"){
    if(is.null(variances)){
      stop("When specifying type 'est_var', variances must also be specified.\n")
    }
    if(!inherits(variances, "matrix") & length(estimates) != length(variances)){
      stop("Variances must either be a matrix with rows and columns equal to those in estimates or a vector of the same length as estimates.\n")
    }
    if(inherits(variances, "matrix")){
      if (nrow(variances) != length(estimates) | ncol(variances) != length(estimates)){
        stop("Variance-covariance matrix must have the same number of rows and columns as there are elements in the estimates vector.\n")
      }
    }
    if(any(is.na(estimates))){
      na_coef <- which(is.na(estimates))
      if(inherits(variances, "matrix")){
        na_var <- which(apply(variances, 1, \(x)all(is.na(x))))
      }else{
        na_var <- which(is.na(variances))
      }
      coef_var_agree <- try(all(na_var == na_coef), silent=TRUE)
      if(!inherits(coef_var_agree, "try-error")){
        estimates <- estimates[-na_coef]
        variances <- variances[-na_var, -na_var]
        message(paste0("Rank-deficient model detected - parameter name(s):", paste(names(na_coef), collapse=", "), " removing NAs from estimates/variances and proceeding.\n"))
      }else{
        stop("NAs found in estimates and/or variances, but they are not in compatible places, please solve the problem and try again.\n")
      }
    }
    if(!inherits(variances, "matrix")){
      V <- diag(variances)
    }else{
      V <- variances
    }
    if(!is_pd(V, tol=tol)){
      stop("Variance-covariance matrix is not positive definite.  Cannot proceed with viztest().\n")
    }
    out <- structure(.Data = list(coef = estimates, vcov = V), class="vtcustom")
  }else{
    if(!is.null(variances)){
      message("When type is 'sim', variances are disregarded.\n")
    }
    out <- structure(.Data = list(est = estimates), class="vtsim")
  }
  return(out)
}


#' @method coef vtcustom
coef.vtcustom <- function(object, ...){
  object$coef
}

#' @method coef eff
coef.eff <- function(object, ...){
  object$fit
}

#' @method coef emmGrid
coef.emmGrid <- function(object, ...){
  nms <- apply(object@grid[,-ncol(object@grid), drop=FALSE], 1, paste, collapse=":")
  res <- object@bhat
  names(res) <- nms
  res
}


#' @method vcov vtcustom
vcov.vtcustom <- function(object, ...){
  object$vcov
}

#' @importFrom dplyr filter
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
#' @param add_test_level Add the (1-test level) confidence interval to the plot.  For this to work, you must have specified `add_test_level` in the call to `viztest()` so that the appropriate confidence intervals can be calculated.  
#' @param ref_lines Reference lines to be plotted - one of "all", "ambiguous", "none".  This could also be a vector of stimulus names to plot - they should be the same as the names of the estimates in `x$est`. See details for explanation. 
#' @param viz_diff_thresh Threshold for identifying visual difficulty, see details. 
#' @param make_plot Logical indicating whether the plot should be constructed or the data returned. 
#' @param level Level at which to plot the estimates.  Accepts both numeric entries or one of "ce", "max", "min", "median" - defaults to "ce", the cognitively easiest level.  
#' @param trans A function to transform the estimates and their confidence intervals like `plogis`.
#' @param est_point_args A list of arguments to be passed to `geom_point()` that plots the point estimates. 
#' @param opt_ci_args A list of arguments to be passed to `geom_linerange()` to plot the optimal visual testing intervals. 
#' @param test_ci_args A list of arguments to be passed to `geom_linerange()` to plot the (1-test level) confidence intervals.
#' @param ref_line_args A list of arguments to be passed to `geom_segment()` to plot the reference lines.
#' @param scale_linewidth_args A list of arguments to be passed to `scale_linewidth_manual()` to change the thickness of the confidence intervals.
#' @param scale_color_args A list of arguments to be passed to `scale_color_manual()` to change the default color of the different confidence intervals when `add_test_level=TRUE`.
#' @param overall_theme A theme function that will be passed to the `ggplot` call before `theme()`. Default is `theme_bw`.
#' @param theme_arg A list of arguments to be passed to `theme()` to modify the theme of the plot.
#' @param remove_caption Logical indicating whether caption should be removed.  By default, it is printed to alert the user. 
#' @param ... Other arguments passed down.  Currently not implemented.
#' @details The `ref_lines` argument identifies what reference lines will be plotted in the figure.  For any particular stimulus, the reference lines run along the upper bound of the stimulus from the stimulus location to the most distant stimulus with overlapping confidence intervals.  
#' When `ref_lines = "all"`, all lines are plotted, though in displays with many stimuli, this can make for a messy graph.  When `"ref_lines = ambiguous"` is specified, then only the ones that help discriminate in cases where the result might be visually difficult to discern are plotted. 
#' A comparison is determined to be visually difficult if the upper bound of the stimulus in question is within `viz_diff_thresh` times the difference between the smallest lower bound and the largest upper bound.  If `ref_lines = "non"`, then none of the reference lines are plotted. 
#' Alternatively, you can specify the names of stimuli whose reference lines will be plotted.  These should be the same as the names in the data.  The `viztest()` function returns an object `est`, which contains the data that are used as input to this function.  The variable `vbl` in 
#' The `est` data frame contains the stimulus names. 
#' @returns By default, a ggplot is returned.  If `make_plot = FALSE`, the data for the plot are returned, but the plot is not constructed.  If the data are returned, the following variables are in the dataset: 
#' * `vbl` - The name of the parameter. 
#' * `est` - The parameter estimate
#' * `se` - The standard error of the estimate
#' * `lwr`, `upr` - The inferential confidence bounds being used
#' * `lwr_add`, `upr_add` - The confidence intervals that come from `add_level`. 
#' * `label` - Factor giving the parameter names
#' * `stim_start`, `stim_end` - y-axis bounds of the reference line
#' * `bound_start`, `bound_end` - x-axis values for reference lines
#' * `ambiguous` - Logical vector indicating whether the comparison is considered "ambiguous". 
#' @method plot viztest
#' @importFrom dplyr left_join arrange `%>%` join_by
#' @importFrom ggplot2 ggplot geom_pointrange geom_segment aes labs geom_point geom_linerange unit theme_bw scale_color_manual scale_linewidth_manual theme margin
#' @importFrom ggtext element_textbox_simple
#' @examples
#' data(mtcars)
#' mod2 <- lm(mpg ~ as.factor(cyl) + vs + am + as.factor(gear), data = mtcars)
#' v <- viztest(mod2)
#' plot(v, ref_lines="ambiguous") + ggplot2::theme_classic()
#' 
#' @export
plot.viztest <- function(x, 
                         ..., 
                         add_test_level = TRUE, 
                         ref_lines="none", 
                         viz_diff_thresh = .02, 
                         make_plot=TRUE, 
                         level=c("ce","max","min","median"),
                         trans=I, 
                         est_point_args = list(color="black", size=2), 
                         opt_ci_args = list(),
                         test_ci_args = list(),
                         ref_line_args = list(color="gray75", linetype=3),
                         scale_linewidth_args = list(values=c(3.5, .5)),
                         scale_color_args = list(values = c("gray75", "black")),
                         overall_theme = theme_bw,
                         theme_arg = list(legend.position="top", 
                                          plot.caption = element_textbox_simple(width = unit(1, "npc"),  # Wraps to plot width
                                                                                halign = 0,
                                                                                margin = margin(1, 0, 0, 0,"lines"))),# Prevents overlap of caption and x-axis title
                         remove_caption=FALSE){
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
  inp[,c("est","lwr","upr","lwr_add","upr_add","bound_start","bound_end")] <- apply(inp[,c("est","lwr","upr","lwr_add","upr_add","bound_start","bound_end")],2,trans)
  if(any(inp$vbl == "zero"))inp <- inp[-which(inp$vbl == "zero"), ]
  if(!make_plot){
    res <- inp
  }else{
    if(add_test_level){
      opt_ci_args$mapping <- aes(x=est, xmin=lwr, xmax=upr, y=label, 
                                 colour=sprintf('Optimal Visual Intervals (%.1f%%)', level*100), 
                                 linewidth=sprintf('Optimal Visual Intervals (%.1f%%)', level*100))
      test_ci_args$mapping <- aes(xmin = lwr_add, xmax=upr_add, y=label, 
                                  colour=sprintf('Standard Confidence Intervals (%.1f%%)', (1-attr(x, "test_level"))*100), 
                                  linewidth=sprintf('Standard Confidence Intervals (%.1f%%)', (1-attr(x, "test_level"))*100))
      lab_args <- list(x="Estimate", y= "Parameter", colour="", linewidth="")
    }else{
      opt_ci_args$mapping <- aes(x=est, xmin=lwr, xmax=upr, y=label)
      lab_args <- list(x=sprintf("Estimate and Optimal Visual CI (%.1f%%)", 
                               level*100), 
                       y = "Parameter")
    }
    est_point_args$mapping <- aes(x=est, y=label)
    ref_line_args$mapping <- aes(x=bound_start, xend=bound_end, y=stim_start, yend=stim_end)
    g <- ggplot(inp) + 
      do.call(geom_linerange, opt_ci_args) +
      do.call(scale_linewidth_manual,scale_linewidth_args) +
      do.call(scale_color_manual,scale_color_args)+
      overall_theme() +
      do.call(theme,theme_arg)
    if(add_test_level & all(c("lwr_add", "upr_add") %in% colnames(inp))){
      g <- g + do.call(geom_linerange, test_ci_args)
    }
    g <- g + do.call(geom_point, est_point_args)
    if("all" %in% ref_lines){
      g <- g + do.call(geom_segment, ref_line_args)
    }  
    if( "ambiguous" %in% ref_lines){
      ref_line_args$data <- inp[which(inp$ambiguous), ]
      g <- g + do.call(geom_segment, ref_line_args)
    }
    if(!any(c("ambiguous", "all", "none") %in% ref_lines)){
      incl <- which(ref_lines %in% inp$vbl)
      if(length(incl) == 0)stop("ref_lines should either be one of (all, ambiguous, or none) or a vector of names consistent with x$est$vbl.\n")
      ref_line_args$data <- inp[incl, ]
      g <- g + do.call(geom_segment, ref_line_args)
    }
    if(!remove_caption)lab_args$caption <- "Confidence intervals have been adjusted to permit visual testing as per Armstrong and Poirier (2025) [doi:10.1017/pan.2024.24]"
    g <- g + do.call(labs, lab_args)
    res <- g
  }
  return(res)
}


#' Make Template for Pairwise Significance Input
#' 
#' Provides a template for producing a binary vector indicating whether each pair of 
#' estimates has a significant difference. 
#' 
#' @param estimates A vector of point estimates (ideally, a named vector). 
#' @param include_zero Logical indicating whether tests against zero should be included. 
#' @param include_intercept Logical indicating whether the intercept should be included. 
#' @param ... Other arguments passed down, currently not implemented. 
#' 
#' @details The `viztest()` function uses a normal difference of means test to identify
#' whether there is a significant difference or not.  While this test could be done 
#' with adjustments for multiplicity or robust standard errors of all different kinds, 
#' there may be times when the user would prefer to identify the significant differences 
#' manually.  The `viztest()` function internally reorders the estimates from largest to smallest
#' so this function does that and then prints the pairs that will correspond with the 
#' visual testing grid search being done by `viztest()`.  
#' 
#' Please note that the `include_zero` and `include_intercept` arguments should be set the same
#' here as they are in your call to `viztest()`.  If they are not, `viztest()` will stop because
#' the results from the comparison of confidence intervals will have different dimensions than the 
#' differences that are manually provides. 
#' 
#' @returns A two-column data frame containing the names of the larger and smaller parameters in the appropriate order. This can be 
#' used to identify the appropriate order in which to specify the `sig_diffs` argument to `viztest()`. 
#' @examples
#' make_diff_template(estimates = c(e1 = 2, e2 = 1, e3 = 3))
#' @export
make_diff_template <- function(estimates, include_zero=TRUE, include_intercept=FALSE, ...){
  if(is.null(names(estimates))){
    names(estimates) <- paste0("est_", 1:length(estimates))
  }
  if(include_zero)estimates <- c(estimates, zero=0)
  if(!include_intercept)estimates <- estimates[!grepl("ntercept", names(estimates))]
  est_num <- seq_along(estimates)
  o <- order(estimates, decreasing = TRUE)
  estimates <- estimates[o]
  est_num <- est_num[o]
  combs <- t(combn(names(estimates), 2))
  colnames(combs) <- c("Larger", "Smaller")
  as.data.frame(combs)
}

#' Calculate z-score for Confidence Interval Overlap
#' 
#' Calculates the z-score required such that confidence intervals do not overlap under the null hypothesis withe a specified probability. 
#' @param b A vector of estiamtes
#' @param v The variance-covariance matrix for `b`. 
#' @param alpha The desired probability at which the confidence intervals do not overlap under the null hypothesis.  
#' @param df Degrees of freedom for the t-distribution, defaults to `Inf` indicating a normal distribution. 
#' @param ... Other arguments passed down, currently not implemented.
#' 
#' @importFrom stats cov2cor qnorm pnorm pt qt sd
#' @importFrom dplyr summarise group_by row_number mutate select
#' @importFrom tidyr pivot_longer
#' 
#' @returns A list with two elements:
#' `ave_z`: A data frame with one row for each estimate in `b` and the following variables: 
#' * `vij`: observation number 
#' * `s_zb`: standard deviation of the z-scores across all pairs of intervals containing that estimate. 
#' * `min_zb`, `max_zb`: The minimum and maximum z-scores for the pairs of intervals containing that estimate.
#' * `zb`: The mean z-score for the pairs of intervals containing that estimate.
#' * `ci`: The confidence level corresponding to `zb`. 
#' `all_z`: A data frame with one row for each pair of estimates in `b` and the following variables:
#' * `i`, `j`: The indices of the two estimates in the pair.
#' * `s_i`, `s_j`: The standard errors of the two estimates in the pair.
#' * `theta`: The ratio of the standard errors of the two estimates.
#' * `rho`: The correlation between the two estimates.
#' * `zb`: The z-score for the pair of estimates.
#' * `ci` : The confidence level corresponding to `zb`.
#' * `olap_ave` The probability that the two intervals do not overlap under the null hypothesis. 
#' * `olap_84` The probability that two 84% confidence intervals for the estimates in the pair would not overlap under the null hypothesis. 
#' @references 
#' Harvey Goldstein and Michael J.R. Healy.  (1995) "The Graphical Presentation of A Collection of Means." Journal of the Royal Statistical Society, Series A 158(1): 175-177 <doi:10.2307/2983411>. 
#' David Afshartous and Richard A. Preston.  (2010) "Confidence Intervals for Dependent Data: Equating Non-overlap with Statistical Significance." Computational Statistics and Data Analysis 54: 2296-2305 <doi:10.1016/j.csda.2010.04.011>
#' @export
#' @examples
#' data(mtcars)
#' mod <- lm(mpg ~ wt + hp + disp + vs, data=mtcars)
#' gen_z(coef(mod), vcov(mod))
#' 
gen_z <- function(b, v, alpha=.05, df = Inf, ...){
  f2.theta.rho <- function(theta, rho) { theta/(theta^{2} +1 - 2*rho*theta)^(1/2) + (1/theta)/(1 + theta^(-2) - 2*rho*theta^{-1})^(1/2)}
  f2.alpha.theta.rho.t.beta <- function(alpha, theta, df, rho) {qt(alpha/2, df, lower.tail=FALSE)/f2.theta.rho(theta, rho) }
  se <- sqrt(diag(v))
  rho <- cov2cor(v)
  eg <- expand.grid(i=1:length(b), j=1:length(b)) 
  eg <- eg[which(eg[,1] < eg[,2]), ]
  eg$s_i <- se[eg[,1]]
  eg$s_j <- se[eg[,2]]
  eg$theta = eg$s_i/eg$s_j
  eg$rho <- rho[cbind(eg[,1], eg[,2])]
  eg$zb <- f2.alpha.theta.rho.t.beta(alpha, eg$theta, df, eg$rho)
  eg$ci <- 1-2*pt(eg$zb, df, lower.tail=FALSE)
  eg$olap_ave = 2*pnorm(mean(eg$zb)*f2.theta.rho(eg$theta, eg$rho), lower.tail=FALSE)
  eg$olap_84 = 2*pnorm(qnorm(.92)*f2.theta.rho(eg$theta, eg$rho), lower.tail=FALSE)
  
  tmp <- eg %>% select(i, j, zb) %>% 
    mutate(obs = row_number()) %>% 
    pivot_longer(c("i", "j"), names_to = "ij", values_to = "vij") %>% 
    group_by(vij) %>% 
    summarise(s_zb = sd(zb), 
              min_zb = min(zb), 
              max_zb = max(zb), 
              zb = mean(zb)
    ) %>% 
    mutate(ci = 1-2*pt(zb, df, lower.tail=FALSE))
  
  return(list(ave_z = tmp, all_z = eg))
}


#' Internal function to check for symmetry
#' @noRd
is_sym <- function(x, tol=1e-08){
  rnd <- abs(floor(log10(abs(tol))))
  x <- round(x,rnd)
  sqr <- nrow(x) == ncol(x) 
  if(sqr){
    lt_ut <- all(x[lower.tri(x)] == t(x)[lower.tri(x)])
  }else{
    lt_ut <- FALSE
  }
  lt_ut
}

#' Internal function to check for positive definiteness
#' @noRd
is_pd <- function(x, tol=1e-08){
  if(is_sym(x, tol = tol) & is.numeric(x)){
    ev <- eigen(x, only.values = TRUE)$values
    ev <- ifelse(abs(ev) < tol, 0, ev)
    !any(ev < 0)
  }else{
    FALSE
  }
}


#' Get Letters for Multiple Comparisons
#' 
#' Gets the letter matrix for a compact letter display.  This can be passed to the `letter_plot()` function from the `psre` package to produce plots of 
#' confidence intervals with a letter display.  
#' 
#' @param x An object that can be one of the following classes: an object of class `glht` produced by `glht()` from the `multcomp` package, 
#' an object of class `emmGrid` produced by the `emmeans` package, or a list with elements `est` - a vector of estimates and `var` - a variance-covariance matrix for the estimates.
#' @param ... Additional arguments passed down either to `cld` if the object is an `emmGrid` class or `summary.glht` if the object is a `glht` class.  
#' If `x` is a list of estimates and variances, it will be converted to an `emmGrid` object internally and the `emmGrid` method dispatched, 
#' `...` will be passed to `cld` in that case. Additional arguments for the `summary.glht()` include `test` and a host of others.  See the help file for `summary.glht` for examples. 
#' Additional arguments for `cld` include `adjust` for multiplicity corrections and others.  See `?emmeans:::cld.emmGrid` for options and details. 
#' 
#' @returns A logical index indicating which estimates are in which letter group. 
#' @export
#' 
#' @importFrom multcomp glht cld
get_letters <- function(x=NULL, ...){
  UseMethod("get_letters")
}

#' @importFrom emmeans emmobj
#' @method get_letters default
#' @export
get_letters.default <- function(x, ...){
  if(!all(c("est", "var") %in% names(x))){
    stop("x should be a list with elements 'est' - a vector of estimates and 'var' - a variance-covariance matrix for the estimates.\n")
  }
  est <- x$est
  nms <- names(est)
  if(is.null(nms)){
    nms <- paste0("b_", seq_along(est))
    names(est) <- nms
  }
  obj <- emmobj(est, x$var, levels=nms)
  get_letters(obj, ...)
}

#' @method get_letters glht 
#' @export
get_letters.glht <- function(x, ...){
  args <- list(...)
  args$object <- x
  s <- do.call(summary, args)
  cl <- cld(s)
  cl$mcletters$LetterMatrix
}

#' @method get_letters emmGrid
#' @export
get_letters.emmGrid <- function(x, ...){
  args <- list(...)
  args$Letters <- letters
  args$object <- x
  cl <- do.call(cld, args)
  grps <- cl$.group
  unletters <- unique(c(unlist(strsplit(grps, ""))))
  if(any(unletters == " ")){
    unletters <- unletters[-which(unletters == " ")]
  }
  unletters <- sort(unletters)
  out <- matrix(NA, nrow = nrow(cl), ncol=length(unletters))
  colnames(out) <- unletters
  rownames(out) <- cl[[1]]
  for(l in seq_along(unletters)){
    out[,l] <- grepl(unletters[l], grps)
  }
  out
}



#' Make Annotations for Significance Brackets
#' 
#' Makes a list of annotations for significance brakcets produced by the `geom_signif()` function from the `ggsignif` package.  The annotations are added for 
#' pairs of estimates whose confidence intervals overlap, but the estimates are nonetheless significantly different from each other.  
#' 
#' @param obj An object of class `viztest` produced by the `viztest()` function.
#' @param type Indicates whether annotations are produced for overlapping intervals that are significantly different from each other or not. The `"auto"` option will
#' find the type that produces the fewest annotations. If `type="discrepancies"`, annotations will be made for pairs of estimates whose test results do not correspond with the (non-)overlaps in the confidence intervals. 
#' @param tol Tolerance for determining whether intervals are close enough to be considered ambiguous.  This also plots significance flags for intervals that 
#' do not overlap, but the distance between them is smaller than the tolerance.  The default is zero, but increasing the value will potentially produce more significance flags. 
#' @param nudge A vector of the same length as the number of brackets.  This will nudge the y-position of the brakcet by the indicated amount.  This will be difficult to 
#' specify ahead of time, but can be specified to clean up a plot after an initial run. 
#' @param ... Other arguments, currently ignored. 
#' 
#' @importFrom stats confint
#' @export
#' @examples
#' data(chickwts)
#' chick_mod <- lm(weight ~ feed, data=chickwts)
#' library(marginaleffects)
#' chick_preds <- avg_predictions(chick_mod, variables="feed")
#' b <- coef(chick_preds)
#' names(b) <- chick_preds$feed
#' v <- vcov(chick_preds)
#' chick_vt_data <- make_vt_data(b, v)
#' chick_vt <- viztest(chick_vt_data, test_level = 0.0001, include_zero=FALSE)
#' chick_vt
#' 
#' make_annotations(chick_vt, type="discrepancies")
make_annotations <- function(obj, type = c("auto", "significant", "insignificant", "discrepancies"),  tol = 0, nudge=NULL, ...){
  typ <- match.arg(type)
  if(typ == "auto"){
    ns <- make_annotations(obj, type="insignificant", tol=tol, nudge=nudge, ...)
    s <- make_annotations(obj, type="significant", tol=tol, nudge=nudge, ...)
    if(length(s$annotations) > length(ns$annotations)){
      return(ns)
    }else{
      return(s)
    }    
  }
  lwr <- obj$est$lwr_add
  names(lwr) <- obj$est$vbl
  upr <- obj$est$upr_add
  names(upr) <- obj$est$vbl
  rng <- max(upr) - min(lwr)
  if(typ != "discrepancies"){
    tmp_combs <- data.frame(
      smaller = obj$param_names[obj$combs[2,]], 
      larger = obj$param_names[obj$combs[1,]] 
    ) %>% 
      mutate(us = upr[smaller], 
             ll = lwr[larger],
             ul = upr[larger], 
             olap = ll < us) 
    if(typ == "significant"){
    tmp_combs <- tmp_combs %>%
      mutate(s = obj$pw_test, 
             ambig = ll > us & (ll - us) <= tol) %>% 
      filter(s & (olap | ambig))
    }else{
      tmp_combs <- tmp_combs %>%
        mutate(s = !obj$pw_test) %>% 
        filter(s & olap)
    }
    if(is.null(nudge))nudge <- rep(0, nrow(tmp_combs))
    if(length(nudge) != nrow(tmp_combs)){
      warning("nudge must be NULL or same length as number of annotations, currently ignored")
      nudge <- rep(0, nrow(tmp_combs))
    } 
    annot_flag <- ifelse(typ == "significant", "*", "NS")
    list(
      annotations = rep(annot_flag, nrow(tmp_combs)), 
      y_position = tmp_combs$ul + .05*rng + nudge, 
      xmin = tmp_combs$smaller, 
      xmax = tmp_combs$larger
    )
  }else{
    tab <- obj$tab
    tab$row <- 1:nrow(tab)
    use_levs <- tab[which(tab$psame == max(tab$psame)), ]
    use_row <- use_levs[which.max(use_levs$easy), "row"]
    tests <- obj$pw_test
    olaps <- obj$ci_tests[,use_row]
    discrep <- which(olaps != tests)
    if(length(discrep) > 0){
      dcombs <- obj$combs[, discrep, drop=FALSE]
      tmp_combs <- data.frame(
        smaller = obj$param_names[obj$combs[2,discrep]], 
        larger = obj$param_names[obj$combs[1,discrep]] 
      ) %>% 
        mutate(ul = upr[larger]) 
      if(is.null(nudge))nudge <- rep(0, nrow(tmp_combs))
      if(length(nudge) != nrow(tmp_combs)){
        warning("nudge must be NULL or same length as number of annotations, nudge being ignored")
        nudge <- rep(0, nrow(tmp_combs))
      } 
      
      list(
        annotations = ifelse(tests[discrep], "*", "NS"), 
        y_position = tmp_combs$ul + .05*rng + nudge, 
        xmin = tmp_combs$smaller, 
        xmax = tmp_combs$larger
      )
    }else{
      stop("No discrepancies between pairwise tests and (non-)overlaps in inferential confidence intervals.\n")
    }
  }
  
}

#' Reorder a factor for a forest plot
#'
#' Orders levels of \code{x} by \code{bvar} (via \code{fn}), with "Summary"
#' optionally pinned to the bottom (or top).
#'
#' @param x        Character/factor column to turn into an ordered factor.
#' @param bvar     Numeric variable used for ordering (e.g. the point estimate).
#' @param fn       Aggregation function (default \code{mean}).
#' @param descending Logical; if \code{TRUE} largest value gets the top row.
#' @param summary_bottom Logical; if \code{TRUE} "Summary" is pinned to row 1
#'   (bottom of a ggplot y-axis).
#' @param ... Other arguments, currently unimplemented.
#' @importFrom stats aggregate
#' @export
reorder_forest <- function(x, bvar, fn = mean,
                           descending = TRUE,
                           summary_bottom = TRUE, ...) {
  ag <- aggregate(
    ifelse(descending, -1, 1) * bvar,   # fixed: was `byvar`
    list(x),
    fn,
    ...
  )
  names(ag) <- c("Group.1", "x")
  ag <- ag[order(ag$x), ]
  
  if (summary_bottom) {
    nms <- c("Summary", setdiff(ag$Group.1, "Summary"))
  } else {
    nms <- c(setdiff(ag$Group.1, "Summary"), "Summary")
  }
  ag <- ag[match(nms, ag$Group.1), ]
  factor(x, levels = ag$Group.1)
}


#' Forest plot points with precision-weighted squares and summary diamonds
#'
#' `geom_forestpoint()` draws the central markers in a forest plot:
#' non-summary rows are rendered as **squares** whose area is controlled by the
#' ggplot2 `size` aesthetic (typically proportional to precision), while summary
#' rows are rendered as **diamonds** whose width reflects the confidence interval
#' and whose height is derived from that width and bounded by row spacing.
#'
#' This geom is designed to work with standard ggplot2 scales (e.g.,
#' `scale_size_area()`) and pairs naturally with `geom_linerange()` for confidence
#' intervals and `geom_foreststripe()` for background striping.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param data The data to be displayed in this layer. If `NULL`, the data are
#'   inherited from the plot data as specified in the call to [ggplot2::ggplot()].
#' @param stat Statistical transformation to use. Defaults to `"identity"`.
#' @param position Position adjustment. Defaults to `"identity"`.
#' @param ... Additional arguments passed to the underlying `GeomForestPoint`
#'   ggproto object (e.g., color, fill, linewidth).
#' @param diamond_aspect Numeric. Controls how strongly the diamond height
#'   increases with confidence interval width. Larger values produce taller
#'   diamonds.
#' @param diamond_row_frac Numeric in (0, 1). Maximum fraction of the vertical
#'   row spacing that the diamond half-height may occupy, preventing overlap
#'   with adjacent rows.
#' @param diamond_min_frac Numeric in (0, 1). Minimum fraction of the row spacing
#'   used as the diamond half-height, ensuring visibility for very narrow
#'   confidence intervals.
#' @param na.rm Logical. If `TRUE`, silently removes missing values.
#' @param show.legend Logical or `NA`. Whether this layer should be included in
#'   the legend.
#' @param inherit.aes Logical. If `FALSE`, the layer does not inherit aesthetic
#'   mappings from the parent plot.
#'
#' @details
#' For non-summary rows, the square side length is derived from the ggplot2
#' `size` aesthetic (in mm units), so users can control point sizing using
#' standard size scales such as [ggplot2::scale_size_area()].
#'
#' For summary rows, the diamond width is determined by the supplied confidence
#' interval (`xmin`/`xmax`), and the height is computed as a bounded function of
#' that width. The diamond height **does not** use the `size` aesthetic.
#'
#' @return A ggplot2 layer object that can be added to a plot.
#'
#' @examples
#' # geom_forestpoint() is designed to be used as part of gg_forest().
#' # See ?gg_forest for a complete runnable example.
#'
#' @importFrom ggplot2 layer
#' @importFrom ggplot2 Geom
#' @importFrom grid grobTree polygonGrob rectGrob
#'
#' @export
geom_forestpoint <- function(mapping = NULL, data = NULL,
                             stat = "identity", position = "identity",
                             ...,
                             diamond_aspect   = 1.0,
                             diamond_row_frac = 0.40,
                             diamond_min_frac = 0.06,
                             na.rm       = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE) {
  ggplot2::layer(
    geom        = GeomForestPoint,
    mapping     = mapping,
    data        = data,
    stat        = stat,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      diamond_aspect   = diamond_aspect,
      diamond_row_frac = diamond_row_frac,
      diamond_min_frac = diamond_min_frac,
      na.rm            = na.rm,
      ...
    )
  )
}

#' @rdname geom_forestpoint
#' @format NULL
#' @usage NULL
#' @importFrom ggplot2 ggproto Geom aes draw_key_polygon alpha
#' @importFrom grid unit rectGrob polygonGrob grobTree gpar nullGrob
#' @export
GeomForestPoint <- ggplot2::ggproto(
  "GeomForestPoint", ggplot2::Geom,
  
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(
    xmin      = NA_real_,
    xmax      = NA_real_,
    is_summary = FALSE,
    size      = 1.5,
    colour    = "black",
    fill      = "white",
    alpha     = 1,
    linewidth = 0.5,
    linetype  = 1
  ),
  
  draw_key = ggplot2::draw_key_polygon,
  
  setup_data = function(data, params) {
    if (!("is_summary" %in% names(data))) data$is_summary <- FALSE
    data$is_summary <- as.logical(data$is_summary)
    data
  },
  
  draw_panel = function(data, panel_params, coord,
                        diamond_aspect   = 1.0,
                        diamond_row_frac = 0.40,
                        diamond_min_frac = 0.06,
                        na.rm = FALSE) {
    
    coords <- coord$transform(data, panel_params)
    keep   <- is.finite(coords$x) & is.finite(coords$y)
    coords <- coords[keep, , drop = FALSE]
    if (nrow(coords) == 0) return(grid::nullGrob())
    
    is_sum  <- coords$is_summary %in% TRUE
    studies <- coords[!is_sum, , drop = FALSE]
    sums    <- coords[ is_sum, , drop = FALSE]
    
    grobs <- list()
    
    # ── Squares (non-summary rows) ──────────────────────────────────────────
    if (nrow(studies) > 0) {
      side <- grid::unit(studies$size, "mm")
      grobs[[length(grobs) + 1L]] <- grid::rectGrob(
        x = studies$x,
        y = studies$y,
        width  = side,
        height = side,
        default.units = "native",
        gp = grid::gpar(
          col  = studies$colour,
          fill = ggplot2::alpha(studies$fill, studies$alpha),
          lwd  = studies$linewidth,
          lty  = studies$linetype
        )
      )
    }
    
    # ── Diamonds (summary rows) ─────────────────────────────────────────────
    if (nrow(sums) > 0) {
      xmin <- sums$xmin; xmax <- sums$xmax
      bad_min <- !is.finite(xmin); xmin[bad_min] <- sums$x[bad_min]
      bad_max <- !is.finite(xmax); xmax[bad_max] <- sums$x[bad_max]
      
      x_rng  <- panel_params$x.range
      x_span <- diff(x_rng)
      if (!is.finite(x_span) || x_span <= 0) x_span <- 1
      
      diamond_row_frac <- max(0, diamond_row_frac)
      diamond_min_frac <- max(0, min(diamond_min_frac, diamond_row_frac))
      
      for (i in seq_len(nrow(sums))) {
        y0    <- sums$y[i]
        dists <- abs(coords$y - y0); dists <- dists[dists > 0]
        dy_min <- if (length(dists) == 0) 1 else min(dists)
        
        dx      <- abs(xmax[i] - xmin[i])
        dx_norm <- dx / x_span
        h_frac  <- min(diamond_min_frac + diamond_aspect * dx_norm,
                       diamond_row_frac)
        h <- h_frac * dy_min
        
        grobs[[length(grobs) + 1L]] <- grid::polygonGrob(
          x = c(xmin[i], sums$x[i], xmax[i], sums$x[i]),
          y = c(y0, y0 + h, y0, y0 - h),
          default.units = "native",
          gp = grid::gpar(
            col  = sums$colour[i],
            fill = ggplot2::alpha(sums$fill[i], sums$alpha[i]),
            lwd  = sums$linewidth[i],
            lty  = sums$linetype[i]
          )
        )
      }
    }
    
    do.call(grid::grobTree, grobs)
  }
)


#' Render a forest-plot table panel as a ggplot2 layer
#'
#' `geom_foresttable()` builds a “table” panel for forest plots by expanding the
#' input data to long form (one row per *forest row × table column*) and drawing
#' formatted cell text at fixed x positions. Column headers are typically added
#' with `scale_x_foresttable()` (top axis tick labels), which makes the table
#' align cleanly with a forest plot when composing with \pkg{patchwork}.
#'
#' The `cols` argument controls which columns are printed. Formatting can be
#' customized via `fmt`, a named list of functions (one per column) that convert
#' cell values to character strings.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()]. If not
#'   supplied, the layer will map `x`, `label`, and `hjust` internally so users
#'   generally only need to map `y`.
#' @param data The data to be displayed in this layer. If `NULL`, the layer
#'   inherits the plot data. If a data frame is provided, it is expanded to
#'   long form immediately. If a function is provided, it is used as-is (advanced).
#' @param position Position adjustment. Defaults to `"identity"`.
#' @param ... Additional arguments passed to the underlying text geom
#'   (`GeomForestTableText`, which inherits from [ggplot2::GeomText]), such as
#'   `colour`, `family`, or `fontface`.
#' @param cols Character vector of column names to print in the table. Must be
#'   non-empty.
#' @param fmt Optional named list of formatting functions. Each function should
#'   take a single value and return a length-1 character (or coercible) result.
#'   Names should correspond to entries of `cols`.
#' @param fmt_default Default formatting function used for columns not present
#'   in `fmt`. Defaults to [base::as.character()].
#' @param col_gap Numeric. Spacing between adjacent table columns on the x axis.
#' @param col_align Character vector of alignments for each column in `cols`.
#'   Each entry must be one of `"left"`, `"center"`, or `"right"`. Length 1 is
#'   recycled to `length(cols)`.
#' @param col_nudge Numeric vector of per-column horizontal nudges (in x-axis
#'   units) applied to the column positions. Length 1 is recycled to
#'   `length(cols)`.
#' @param na.rm Logical. If `TRUE`, silently removes missing values.
#' @param show.legend Logical. Whether this layer should be included in legends.
#'   Defaults to `FALSE`.
#' @param inherit.aes Logical. If `FALSE`, the layer does not inherit aesthetic
#'   mappings from the parent plot.
#'
#' @details
#' `geom_foresttable()` uses ggplot2's feature where the `data` argument may be a
#' function: when `data = NULL`, this layer supplies a data function that receives
#' the plot data and expands it to long form while preserving all original
#' columns (including faceting variables and the y mapping). This makes the geom
#' compatible with faceting and grouping while keeping the user-facing API
#' simple.
#'
#' Column positions are computed as `seq_along(cols) * col_gap + col_nudge`.
#' Horizontal justification is determined by `col_align` (left/center/right).
#'
#' @return A ggplot2 layer object that can be added to a plot.
#'
#' @examples
#' # geom_foresttable() is designed to be used as part of gg_forest().
#' # See ?gg_forest for a complete runnable example.
#'
#' @seealso
#' scale_x_foresttable(), gg_forest()
#'
#' @importFrom ggplot2 layer aes ggproto GeomText
#' @importFrom stats setNames
#' @importFrom utils modifyList
#' @export
geom_foresttable <- function(mapping     = NULL,
                             data        = NULL,   # may be NULL → inherits
                             position    = "identity",
                             ...,
                             cols,
                             fmt         = NULL,
                             fmt_default = as.character,
                             col_gap     = 1,
                             col_align   = "left",
                             col_nudge   = 0, 
                             na.rm       = FALSE,
                             show.legend = FALSE,
                             inherit.aes = TRUE) {
  
  if (missing(cols) || length(cols) == 0)
    stop("geom_foresttable(): `cols` must be a non-empty character vector.")
  
  # ── alignment lookup ────────────────────────────────────────────────────────
  n_cols <- length(cols)
  if (length(col_align) == 1L) col_align <- rep(col_align, n_cols)
  if (length(col_align) != n_cols)
    stop("`col_align` must be length 1 or the same length as `cols`.")
  
  if (length(col_nudge) == 1L) col_nudge <- rep(col_nudge, n_cols)
  if (length(col_nudge) != n_cols)
    stop("`col_nudge` must be length 1 or the same length as `cols`.")
  names(col_nudge) <- cols
  
  hjust_by_col        <- ifelse(col_align == "right",  1,
                                ifelse(col_align == "center", 0.5, 0))
  names(hjust_by_col) <- cols
  
  # ── x positions ─────────────────────────────────────────────────────────────
  x_lines        <- (seq_along(cols)) * col_gap + col_nudge
  names(x_lines) <- cols
  
  # ── formatting ──────────────────────────────────────────────────────────────
  if (is.null(fmt)) fmt <- list()
  
  fmt_cell <- function(cn, v) {
    fn  <- if (!is.null(fmt[[cn]])) fmt[[cn]] else fmt_default
    out <- tryCatch(fn(v), error = function(e) "")
    if (length(out) != 1L) out <- out[1L]
    if (is.na(out))         out <- ""
    as.character(out)
  }
  
  # ── data function ────────────────────────────────────────────────────────────
  # ggplot2 calls the `data` argument as a function(plot_data) when it is a
  # function.  We use this to (a) keep *all* original columns, and (b) expand
  # to long form so GeomText gets one row per (row × column) cell.
  #
  # The closure captures cols, x_lines, hjust_by_col, fmt_cell.
  
  data_fn <- function(plot_data) {
    
    # Resolve which requested cols are actually present
    cols_here <- intersect(cols, names(plot_data))
    if (length(cols_here) == 0L) {
      warning("geom_foresttable: none of `cols` found in data. ",
              "Columns present: ", paste(names(plot_data), collapse = ", "))
      return(plot_data[0L, , drop = FALSE])
    }
    
    # Build long-form data: one row per (original row × displayed column)
    rows <- lapply(cols_here, function(cn) {
      d           <- plot_data                        # keep ALL columns (incl. y)
      d$x         <- unname(x_lines[cn])             # col_gap position + col_nudge
      d$hjust     <- unname(hjust_by_col[cn])
      d$.col_name <- cn
      d$label     <- vapply(plot_data[[cn]], function(v) fmt_cell(cn, v),
                            character(1L))
      d
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
  }
  
  # If the user supplied explicit data, wrap it; otherwise use NULL so ggplot2
  # inherits plot data — but we need to intercept it either way.
  # Strategy: if data is a plain data.frame, wrap it in data_fn logic directly.
  if (!is.null(data) && is.data.frame(data)) {
    layer_data <- data_fn(data)
  } else if (is.null(data)) {
    layer_data <- data_fn   # ggplot2 will call this with plot$data
  } else {
    layer_data <- data      # already a function supplied by user
  }
  
  # x, label, and hjust come from data_fn columns; inject into mapping so
  # GeomForestTableText can find them without the user having to spell them out.
  if (is.null(mapping)) mapping <- ggplot2::aes()
  if (!"x"     %in% names(mapping)) mapping <- modifyList(mapping, ggplot2::aes(x     = x))
  if (!"label" %in% names(mapping)) mapping <- modifyList(mapping, ggplot2::aes(label = label))
  if (!"hjust" %in% names(mapping)) mapping <- modifyList(mapping, ggplot2::aes(hjust = hjust))
  
  ggplot2::layer(
    stat        = "identity",          # data is already in final form
    geom        = GeomForestTableText,
    mapping     = mapping,
    data        = layer_data,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(na.rm = na.rm, ...)
  )
}

# ── GeomForestTableText ───────────────────────────────────────────────────────
# A thin wrapper around GeomText that declares x, y, label, hjust as required
# so ggplot2 does not complain about missing aesthetics.

#' @rdname geom_foresttable
#' @format NULL
#' @usage NULL
#' @export
GeomForestTableText <- ggplot2::ggproto(
  "GeomForestTableText", ggplot2::GeomText,
  
  required_aes = c("x", "y", "label"),
  
  default_aes = ggplot2::aes(
    colour    = "black",
    size      = 3.5,
    angle     = 0,
    hjust     = 0,
    vjust     = 0.5,
    alpha     = NA,
    family    = "",
    fontface  = 1,
    lineheight = 1.2
  )
)


#' X scale for forest tables with column headers on the top axis
#'
#' `scale_x_foresttable()` constructs a continuous x scale designed for use with
#' `geom_foresttable()`. It places table column headers on the **top axis tick
#' labels**, sets appropriate breaks and limits based on the requested columns,
#' and ensures consistent spacing so that table panels align cleanly with forest
#' plots when composed using \pkg{patchwork}.
#'
#' The scale uses numeric x positions internally (one per column), which allows
#' precise control over column spacing and alignment.
#'
#' @param cols Character vector of column names corresponding to the table
#'   columns. Determines the number, order, and spacing of x-axis breaks.
#' @param col_labels Optional character vector giving display labels for each
#'   column. If unnamed, must be the same length as `cols` and is assumed to be
#'   in the same order. If named, names are matched to `cols`.
#' @param col_gap Numeric. Spacing between adjacent columns on the x axis.
#' @param position Character. Position of the x axis; defaults to `"top"`.
#' @param expand Expansion applied to the x scale. Defaults to no expansion
#'   (`expansion(mult = c(0, 0))`).
#' @param limits Optional numeric vector of length 2 giving explicit x-axis
#'   limits. If `NULL`, limits are computed automatically to bracket all columns.
#' @param ... Additional arguments passed to [ggplot2::scale_x_continuous()].
#'
#' @details
#' This scale is typically paired with `geom_foresttable()` and a void theme
#' (`theme_void()`) to create a table-like panel whose headers are rendered as
#' axis labels. Because headers are true axis labels, ggplot2 can align the table
#' panel perfectly with a forest plot panel when combining plots.
#'
#' @return A ggplot2 scale object suitable for addition to a plot.
#'
#' @seealso
#' geom_foresttable(), gg_forest()
#'
#' @importFrom ggplot2 scale_x_continuous expansion
#' @importFrom stats setNames
#'
#' @export
scale_x_foresttable <- function(cols,
                                col_labels = NULL,
                                col_gap    = 1,
                                position   = "top",
                                expand     = ggplot2::expansion(mult = c(0, 0)),
                                limits     = NULL,
                                ...) {
  stopifnot(is.character(cols), length(cols) > 0)
  
  breaks <- seq_along(cols) * col_gap
  
  if (is.null(col_labels)) {
    labels <- cols
  } else if (is.null(names(col_labels))) {
    stopifnot(length(col_labels) == length(cols))
    labels <- col_labels
  } else {
    missing_nm <- setdiff(cols, names(col_labels))
    if (length(missing_nm))
      col_labels <- c(col_labels, stats::setNames(missing_nm, missing_nm))
    labels <- unname(col_labels[cols])
  }
  
  if (is.null(limits))
    limits <- c(min(breaks) - 0.5 * col_gap, max(breaks) + 0.5 * col_gap)
  
  ggplot2::scale_x_continuous(
    breaks   = breaks,
    labels   = labels,
    position = position,
    expand   = expand,
    limits   = limits,
    ...
  )
}



#' Draw alternating row stripes for forest plots and forest tables
#'
#' `geom_foreststripe()` adds alternating horizontal background bands (“zebra
#' striping”) to forest plots or forest tables. The stripes are drawn using the
#' y-axis scale (rather than the data itself), which makes the geom robust to
#' missing rows, summary rows, and faceting.
#'
#' This geom is designed to work with both the forest figure panel (estimates and
#' confidence intervals) and the forest table panel produced by
#' `geom_foresttable()`. Because stripes are computed from the panel scales, the
#' two panels remain visually aligned when composed with \pkg{patchwork}.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()]. Usually
#'   left empty; the geom does not require data-driven aesthetics.
#' @param data The data to be displayed in this layer. If `NULL`, the layer
#'   inherits the plot data; only the y scale is used.
#' @param position Position adjustment. Defaults to `"identity"`.
#' @param ... Additional arguments passed to the underlying `GeomForestStripe`.
#' @param n_cols Integer. Number of table columns. Used to determine the horizontal
#'   extent of the stripes.
#' @param col_gap Numeric. Spacing between table columns on the x axis.
#' @param start Integer. Index of the first row to stripe (counting from the top
#'   of the plot). Defaults to `2L`, so striping begins on the second row.
#' @param fill Fill colour for the stripes.
#' @param colour Border colour for the stripes. Defaults to `NA` (no border).
#' @param na.rm Logical. If `TRUE`, silently ignores missing values.
#' @param show.legend Logical. Whether this layer should be included in the legend.
#'   Defaults to `FALSE`.
#' @param inherit.aes Logical. If `FALSE`, the layer does not inherit aesthetic
#'   mappings from the parent plot.
#'
#' @details
#' The stripes are computed from the panel’s y range rather than the data rows.
#' For discrete y scales, this corresponds to the integer row positions used by
#' ggplot2. Every other row is selected starting at `start`.
#'
#' Horizontal extents are derived from the panel x range, allowing the stripes to
#' span both table columns and forest-plot panels without requiring explicit
#' xmin/xmax aesthetics.
#'
#' @return A ggplot2 layer object that can be added to a plot.
#'
#' @examples
#' # geom_foreststripe() is designed to be used as part of gg_forest().
#' # See ?gg_forest for a complete runnable example.
#'
#' @seealso
#' geom_foresttable(), geom_forestpoint(), gg_forest()
#'
#' @importFrom ggplot2 layer ggproto Geom aes draw_key_rect alpha
#' @importFrom grid rectGrob gTree gList gpar viewport nullGrob
#'
#' @export
geom_foreststripe <- function(mapping     = NULL,
                              data        = NULL,
                              position    = "identity",
                              ...,
                              n_cols,
                              col_gap     = 1,
                              start       = 2L,
                              fill        = "grey92",
                              colour      = NA,
                              na.rm       = FALSE,
                              show.legend = FALSE,
                              inherit.aes = TRUE) {
  
  if (missing(n_cols))
    stop("geom_foreststripe(): `n_cols` must be supplied (number of table columns).")
  
  ggplot2::layer(
    stat        = "identity",
    geom        = GeomForestStripe,
    mapping     = mapping,
    data        = data,
    position    = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params      = list(
      n_cols  = n_cols,
      col_gap = col_gap,
      start   = as.integer(start),
      fill    = fill,
      colour  = colour,
      na.rm   = na.rm,
      ...
    )
  )
}

#' @rdname geom_foreststripe
#' @format NULL
#' @usage NULL
#' @export
GeomForestStripe <- ggplot2::ggproto(
  "GeomForestStripe", ggplot2::Geom,
  
  # No required_aes — we only need the y scale, read inside draw_panel.
  required_aes = character(0),
  
  default_aes = ggplot2::aes(
    fill      = "grey92",
    colour    = NA,
    alpha     = 1,
    linetype  = 1,
    linewidth = 0
  ),
  
  draw_key = ggplot2::draw_key_rect,
  
  draw_panel = function(data, panel_params, coord,
                        n_cols  = 4L,
                        col_gap = 1,
                        start   = 2L,
                        fill    = "grey92",
                        colour  = NA,
                        na.rm   = FALSE) {
    
    # All unique y positions come from the panel y scale, not the data.
    # For a discrete scale this is panel_params$y.range (the integer range).
    y_range <- panel_params$y.range          # e.g. c(0.4, 13.6) for 13 levels
    y_vals  <- seq(floor(y_range[1] + 0.5),  # integer y positions (1, 2, …, n)
                   floor(y_range[2] + 0.5))
    
    # Sort descending so row 1 = top of plot (highest y = top on screen)
    y_vals <- sort(y_vals, decreasing = TRUE)
    n      <- length(y_vals)
    
    # Pick every-other row starting at `start`
    stripe_idx <- seq(from = start, to = (n-start+1), by = 2L)
    if (length(stripe_idx) == 0L) return(grid::nullGrob())
    y_stripe <- y_vals[stripe_idx]
    
    xmin_val <- panel_params$x.range[1]
    xmax_val <- panel_params$x.range[2]
    # xmin_val <- 0.5  * col_gap
    # xmax_val <- (n_cols + 0.5) * col_gap
    
    # Build a minimal data frame and transform through coord
    stripe_df <- data.frame(
      x    = (xmin_val + xmax_val) / 2,
      y    = y_stripe,
      xmin = xmin_val,
      xmax = xmax_val,
      ymin = y_stripe - 0.5,
      ymax = y_stripe + 0.5
    )
    
    coords <- coord$transform(stripe_df, panel_params)
    
    rects <- grid::rectGrob(
      x      = coords$x,
      y      = coords$y,
      width  = coords$xmax - coords$xmin,
      height = coords$ymax - coords$ymin,
      default.units = "native",
      gp = grid::gpar(
        fill = ggplot2::alpha(fill,   data$alpha[1] %||% 1),
        col  = colour,
        lty  = data$linetype[1]  %||% 1,
        lwd  = data$linewidth[1] %||% 0
      )
    )
    
    grid::gTree(
      children = grid::gList(rects),
      vp = grid::viewport(clip = "off")
    )
  }
)

# Internal helper (mirrors rlang's %||% without the dependency)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

#' Build a paired forest plot + companion table (patchwork-ready)
#'
#' `gg_forest()` creates two aligned ggplot objects: (1) a forest plot with
#' confidence intervals and weighted points (including summary diamonds), and
#' (2) a “table” rendered as text in a ggplot panel with column headers on the
#' top axis. The returned objects share the same y values so they can be combined
#' with \pkg{patchwork}.  This uses lower level functions in the package like `geom_forestpooint`, 
#' `geom_foreststripe`, `geom_foresttable` and `scale_x_foresttable` that could
#' be used to customize the look of the forest plot further. 
#'
#' Column arguments are provided as strings and are evaluated safely using
#' `.data[[...]]`.
#'
#' @param data A data frame containing one row per forest row (study and optional
#'   summary rows) and all referenced columns.
#' @param y String. Column name used for the forest-row coordinate.
#' @param x String. Column name containing the point estimate.
#' @param data_cols A names character vector of column names to display in the table.  
#' The values identify the variable names and the names are the header labels that 
#' will be used in the table. 
#' @param xmin_std,xmax_std Strings. Column names for the standard CI bounds.
#' @param size_prop Numeric scalar or string column name giving point-size weights.
#' If a numeric scalar, the all points will recieve the same weight. 
#' @param max_size Numeric. Maximum size for points (passed to
#'   [ggplot2::scale_size_area()]).  This will control the overall size of the points and may 
#'   need some tweaking depending on the size of the figure and the number of rows.
#' @param is_summary Logical scalar or string column name indicating summary rows.
#' If a logical scalar, all rows will be get square points (if `FALSE`) or diamond points (if `TRUE`).
#' @param xmin_ici,xmax_ici Optional strings for the names of the inferential CI bounds.
#' @param vline Numeric. Value at which a dashed vertical line will be drawn.  If `NULL`, no line is drawn. 
#' @param ci_colors Named character vector with entries `"std"` and `"ici"`.
#' @param use_log_scale Logical. If `TRUE`, uses a log-scaled x axis.
#' @param stripe_table Logical. Draw alternating row stripes in the table.
#' @param stripe_figure Logical. Draw alternating row stripes in the forest plot.
#' @param start_stripe Integer. Row index at which striping begins.
#' @param col_nudge Numeric. Horizontal nudge applied to table text. This is particularly 
#' useful for left-aligned columns where the header text and cell text do not line up natively. 
#' @param table_format_list Optional named list of formatting functions for
#'   `data_cols`.  This defaults to `as.character` for factors and characters, 
#'   `sprintf("%d")` for integers, and `sprintf("%.2f")` for numerics.  
#'   The names of the list should match the values of `data_cols`.
#' @param col_align Optional character vector giving alignment for each table
#'   column (`"left"`, `"center"`, `"right"`).
#' @param table_header_size Numeric. Font size for table headers.
#' @param table_text_size Numeric. Font size for table body text.
#' @param fig_xlab Character. X-axis label for the forest plot.
#' @param ... Additional arguments passed to `geom_forestpoint()`.
#'
#' @return An object of class `"gg_forest"`: a list with two ggplot objects
#'   (`forest`, `table`) suitable for composition with \pkg{patchwork}.
#'
#' @examples
#' # Load Packages
#' library(emmeans)
#' library(VizTest)
#' library(dplyr)
#' 
#' # Use built-in Esophageal Cancer Data
#' data(esoph)
#' 
#' # Aggregate data by age group
#' ag_data <- aggregate(esoph[,c("ncases", "ncontrols")], list(age = esoph$age), sum)
#' 
#' # Turn counts into integerss (not required, but makes printing nicer)
#' ag_data$ncases <- as.integer(ag_data$ncases)
#' ag_data$ncontrols <- as.integer(ag_data$ncontrols)
#' 
#' # Make age into unordered factor
#' ag_data$age <- factor(as.character(ag_data$age), 
#'                       levels=levels(esoph$age))
#'                       
#' # Estimate model of prevalence by age and overall (the summary model)                       
#' model1 <- glm(cbind(ncases, ncontrols) ~ age,
#'               data = ag_data, family = binomial())
#' model_sum <- glm(cbind(ncases, ncontrols) ~ 1,
#'                  data = ag_data, family = binomial())
#' 
#' # Make data frame of results for plotting using emmeans
#' fit <- emmeans(model1, "age")
#' fit_ci <- confint(fit)
#' 
#' # add in original count data
#' ag_data <- cbind(fit, ag_data[,c("ncases", "ncontrols")])
#' 
#' # turn coefficients and confidence intervals into odds ratio scale 
#' ag_data$or <- exp(ag_data$emmean)
#' ag_data$lower <- exp(ag_data$asymp.LCL)
#' ag_data$upper <- exp(ag_data$asymp.UCL)

#' # Make summary data frame that we can use for plotting                  
#' fit_sum <- data.frame(age= "Summary", emmean = coef(model_sum), 
#'   SE = unname(sqrt(vcov(model_sum))), or = exp(coef(model_sum)))
#' sum_ci <- confint(model_sum)
#' fit_sum$lower <- exp(sum_ci[1])
#' fit_sum$upper <- exp(sum_ci[2])
#' fit_sum$ncases <- sum(ag_data$ncases)
#' fit_sum$ncontrols <- sum(ag_data$ncontrols)
#' rownames(fit_sum) <- NULL
#' 
#' # Find the optimal visual testing intervals
#' viztest(fit, include_zero=FALSE, make_plot=FALSE, test_level = .05)
#' 
#' # Add inferential CIs to data (not for summary, though)
#' fit_ici <- confint(fit, level = .75)
#' ag_data$lower_ici <- exp(fit_ici$asymp.LCL)
#' ag_data$upper_ici <- exp(fit_ici$asymp.UCL)
#' 
#' # bind together the age-specific and summary data frames for plotting
#' ag_data <- dplyr::bind_rows(ag_data, fit_sum)
#' 
#' # identify summary row
#' ag_data$is_sum <- ag_data$age == "Summary"
#' 
#' # add point-size weight
#' ag_data$pt_size <- 1/ag_data$SE^2
#' 
#' # make age_label such that ages plot smallest at top and summary at bottom
#' ag_data$age_label <- factor(ag_data$age, levels=rev(ag_data$age))
#' 
#' # Make gg forest plot
#' out <- gg_forest(ag_data, 
#'     y = "age_label", 
#'     x = "or", 
#'     xmin_std = "lower", 
#'     xmax_std = "upper", 
#'     xmin_ici = "lower_ici", 
#'     xmax_ici = "upper_ici",
#'     size_prop = "pt_size", 
#'     is_summary = "is_sum", 
#'     use_log_scale = TRUE, 
#'     data_cols = c("Age" = "age_label", 
#'                    "Controls" = "ncontrols", 
#'                    "Cases"="ncases", 
#'                    "OR" = "or"), 
#'     max_size=5,
#'     table_header_size = 16, 
#'     table_text_size = 5, 
#'     col_nudge=c(-.085, 0,0,0), 
#'     diamond_aspect=15, 
#'     diamond_row_frac = .9)
#'    
#' # print plot
#' plot(out, widths=1, 1)
#'
#' @seealso
#' geom_forestpoint(), geom_foresttable(), geom_foreststripe(),
#' scale_x_foresttable()
#'
#' @importFrom ggplot2 ggplot aes labs theme theme_classic theme_void element_text
#' @importFrom ggplot2 element_blank element_line element_rect guides
#' @importFrom ggplot2 geom_linerange scale_color_manual scale_size_area geom_vline
#' @importFrom ggplot2 scale_x_log10
#' @importFrom ggplot2 .data
#' @importFrom grid unit
#'
#' @export
gg_forest <- function(data, 
                      y, 
                      x,
                      data_cols, 
                      xmin_std, 
                      xmax_std, 
                      size_prop = 1, 
                      max_size=15, 
                      is_summary = FALSE, 
                      xmin_ici = NULL, 
                      xmax_ici = NULL, 
                      vline = NULL, 
                      ci_colors=c("ici" = "gray65", "std" = "black"), 
                      use_log_scale = FALSE, 
                      stripe_table = TRUE, 
                      stripe_figure = TRUE, 
                      start_stripe = 3, 
                      col_nudge = 0, 
                      table_format_list = NULL, 
                      col_align = NULL, 
                      table_header_size = 16, 
                      table_text_size = 4, 
                      fig_xlab = "Estimate",
                      ...){
  if(!requireNamespace("patchwork")){
    stop("The `patchwork` package is required to use `gg_forest()`. Please install it with `install.packages('patchwork')`.", call. = FALSE)
  }
  dots <- list(...)
  if(!is.character(is_summary)){
    data$is_smry <- is_summary
    is_summary <- "is_smry"
  }
  if(!is.character(size_prop)){
    data$size_prop <- size_prop
    size_prop <- "size_prop"
  }
  if(!is.null(xmin_ici) & !is.null(xmax_ici)){
    col_vals <- unname(ci_colors[c( "ici", "std")])
  }else{
    col_vals <- unname(ci_colors["std"])
  }
  fg <- ggplot(data, aes(x=.data[[x]], y =.data[[y]]))  
  if(stripe_figure){
    fg <- fg + geom_foreststripe(n_cols = length(data_cols), col_gap = 1, start = start_stripe) 
  }
  fg <- fg + geom_linerange(aes(xmin = .data[[xmin_std]], xmax = .data[[xmax_std]], colour="Standard CI"),  linewidth=.75, show.legend=!is.null(xmin_ici)) 
  if(!is.null(xmin_ici) & !is.null(xmax_ici)){
    fg <- fg + geom_linerange(aes(xmin=.data[[xmin_ici]], xmax=.data[[xmax_ici]], colour="Inferential CI"),  linewidth=5)
  }
  fg <- fg +   geom_forestpoint(aes(is_summary = .data[[is_summary]], 
                                    xmin = .data[[xmin_std]], 
                                    xmax = .data[[xmax_std]], 
                                    size = .data[[size_prop]]), 
                                fill="black", show.legend=FALSE) + guides(size = "none")
  
  if(use_log_scale) fg <- fg + scale_x_log10()
  if(!is.null(vline)) fg <- fg + geom_vline(xintercept = vline, linetype="dashed", color="black")
  fg <- fg + scale_size_area(max_size = max_size) + 
    theme_classic() + 
    theme(legend.position="bottom", 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.title.y = element_blank(), 
          panel.border = element_blank(), 
          axis.line.y = element_blank()) + 
    labs(x = fig_xlab, color="", size="") + 
    scale_color_manual(values=col_vals)
  
  if(is.null(table_format_list)){
    flist <- vector(mode = "list", length = length(data_cols))
    names(flist) <- data_cols
    for(i in seq_along(data_cols)){
      z <- data[[data_cols[i]]]
      if(inherits(z, "factor") | inherits(z, "character")) flist[[i]] <- as.character
      if(inherits(z, "integer"))  flist[[i]] <- \(x) ifelse(is.na(x), "", sprintf("%d", x))
      if(inherits(z, "numeric")) flist[[i]] <- \(x) ifelse(is.na(x), "", sprintf("%.2f", x))
    }
  }else{
    flist <- table_format_list
  }
  if(is.null(col_align)) col_align <- c("left", rep("center", length(data_cols)-1))
  
  ft <- ggplot(data, aes(y = .data[[y]])) 
  if(stripe_table){
    ft <- ft + geom_foreststripe(n_cols = length(data_cols), col_gap = 1, start = start_stripe) 
  }
  ft <- ft + 
    geom_foresttable(
      cols = data_cols,
      col_nudge = col_nudge, 
      fmt = flist,
      col_align = col_align,
      size = table_text_size
    ) +
    scale_x_foresttable(data_cols, names(data_cols),
                        col_gap = 1) +   # returns list(scale, theme) — ggplot2 handles lists fine
    theme_void() +
    theme(
      axis.text.x.top   = element_text(face = "bold", size = table_header_size, margin = margin(t = 8)),
      axis.line.x.top   = element_blank(),
      axis.ticks.x.top  = element_blank(),
      axis.ticks.length = unit(0, "pt")
    )
  
  out <- list(fg, ft)
  class(out) <- c("gg_forest", class(out))
  out
}

#' Plotting method for `gg_forest` objects. 
#' 
#' Uses patchwork to combine the forest plot and table components.  
#' @param x An object of class `"gg_foest"`, produces with the `gg_forest()` function.
#' @param widths Numeric vector of length 2 giving the relative widths of the table and forest plot components. 
#'  The default is `c(1.4, 1)`, which gives the table slightly more width than the forest plot.  
#' @param ... Other arguments passed to `theme()` to control the appearance of the combined plot.  
#' For example, you might want to use `legend.position = "none"` if you don't want a legend for the CI colors.
#' 
#' @return A ggplot object combining the forest plot and table components.
#' 
#' @importFrom ggplot2 theme
#' @method plot gg_forest
#' @export
plot.gg_forest <- function(x, ..., widths=c(1.4, 1)){
  if(!requireNamespace("patchwork")){
    stop("The `patchwork` package is required to plot `gg_forest` objects. Please install it with `install.packages('patchwork')`.", call. = FALSE)
  }
  other_args <- list(...)
  if(!"legend.position" %in% names(other_args)){
    other_args$legend.position <- "bottom"
  }
  thm <- do.call(theme, other_args)
  pwand <- function (e1, e2) {
    if (TRUE) {
      if (TRUE) {
        e1$patches$annotation$theme <- e1$patches$annotation$theme + 
          e2
      }
      e1$patches$plots <- lapply(e1$patches$plots, function(p) {
        if (FALSE) {
          p <- p & e2
        }
        else {
          p <- p + e2
        }
        p
      })
    }
    e1 + e2
  }
  e1 <- (x[[2]] + x[[1]] + patchwork::plot_layout(widths = widths, guides="collect"))
  pwand(e1, thm)
}


#' Extract representative difficulty levels from a VizTest result
#'
#' `get_viztest_levels()` selects one or more representative stimulus levels
#' from a VizTest result object based on empirical difficulty and agreement.
#' Levels are chosen from rows with the maximum agreement probability
#' (`psame`) and can be returned individually (e.g., lowest, highest) or
#' collectively.
#'
#' @param x A VizTest result object containing a component `tab`, a data frame
#'   with at least the columns `level`, `psame`, and `easy`.
#' @param method Character string specifying which level(s) to return.
#'   One of `"easiest"`, `"all"`, `"lowest"`, `"highest"`, or `"middle"`.
#' @param ... Reserved for future extensions; currently unused.
#'
#' @details
#' The function first subsets `x$tab` to rows achieving the maximum value of
#' `psame` (ignoring missing values). From this subset:
#' \itemize{
#'   \item `"lowest"` returns the smallest level.
#'   \item `"highest"` returns the largest level.
#'   \item `"middle"` returns the median level.
#'   \item `"easiest"` returns the level with the maximum `easy` value.
#'   \item `"all"` returns a named vector containing all four selections.
#' }
#'
#' @return
#' A numeric value when `method` is one of `"lowest"`, `"highest"`, `"middle"`,
#' or `"easiest"`. When `method = "all"`, a named numeric vector with elements
#' `lowest`, `highest`, `middle`, and `easiest`.
#'
#' @examples
#' data(mtcars)
#' mtcars$cyl <- as.factor(mtcars$cyl)
#' mtcars$hp <- scale(mtcars$hp)
#' mtcars$wt <- scale(mtcars$wt)
#' mod <- lm(qsec ~ hp + wt + cyl, data=mtcars)
#' v <- viztest(mod)
#' v
#' 
#' get_viztest_levels(v, "easiest")
#' 
#' @seealso
#' VizTest
#'
#' @importFrom stats quantile
#'
#' @export
get_viztest_levels <- function(x, method=c("easiest", "all", "lowest", "highest", "middle"), ...){
  meth <- match.arg(method)
  tab <- x$tab
  tab <- tab[which(tab$psame == max(tab$psame, na.rm=TRUE)), ]
  l <- tab$level[1]
  h <- tab$level[nrow(tab)]
  m <- quantile(tab$level, .5, na.rm=TRUE)
  e <- tab$level[which.max(tab$easy)]
  if(meth == "all"){
    out <- c("lowest" = l, "highest" = h, "middle" = m, "easiest" = e)
  }
  if(meth == "lowest") out <- l
  if(meth == "highest") out <- h
  if(meth == "middle") out <- m
  if(meth == "easiest") out <- e
  out
}


