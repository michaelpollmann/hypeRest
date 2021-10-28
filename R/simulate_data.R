#' Simulate a data set of hypothetical evaluations
#'
#' @details
#' This function simulates data that can be used to illustrate to illustrate the estimation methods included in hypeRest package.
#' The data generating process is **not** calibrated to resemble any real data set.
#' The data are simulated such that the simple difference in means estimator is (severely) biased.
#'
#' @param n integer, sample size, must be large enough to allow regression on hypotheticals and fixed characteristics.
#' @param num_h integer, number of hypothetical evaluation variables, must be at least 1.
#' @param num_x integer, number of fixed characteristics, can be 0.
#' @param frac_coeff_0 numeric, between 0 and 1, probability of each coefficient being 0.
#'
#' @returns list:
#'          \item{dat_wide}{data.frame with `n` rows and `4+3*num_h+num_x` columns.
#'                            Column Y holds observed outcomes.
#'                            Column W holds realized treatment state.
#'                            Columns H1 through H<num_h> hold hypothetical evaluations corresponding to the realized treatment state.
#'                            Columns H1_0 through H<num_h>_0 hold hypothetical evaluations corresponding to the control state.
#'                            Columns H1_1 through H<num_h>_1 hold hypothetical evaluations corresponding to the treated state.
#'                            Columns X1 through X<num_x> hold fixed characteristics.
#'                            Columns Y1 and Y0 hold potential outcomes (under treatment and control) that are not typically available in real data sets}
#'          \item{dat_long}{data.frame with `2*n` rows and `2+num_h+num_x` columns.
#'                            The first `n` rows contain data under the control state of each observation,
#'                            the last `n` rows contain data under the treated state of each observation;
#'                            observations are in the same order in each half.
#'                            Column Y holds the outcome (if observed under the treatment state, otherwise NA).
#'                            Column W indicates the treatment state for which the row contains data.
#'                            Columns H1 through H<num_h> hold hypothetical evaluations for the treatment state of the row.
#'                            Columns X1 through X<num_x> hold fixed characteristics.}
#'          \item{tau_i}{vector of length `n`. Individual-level treatment effects (in the same order as the observations in `dat_wide` and `dat_long`.}
#'          \item{tau}{numeric. The true in-sample average treatment effect.}
#'
#' @examples
#' # simulate a data set with 500 observations, 3 hypothetical variables, 2 fixed characteristics
#' sims <- simulate_data(n=500, num_h=3, num_x=2)
#' # extract data
#' dat_wide <- sims$dat_wide
#' # estimate effects using hypothetical evaluations; match column names generated inside simulate_data
#' est_hyp <- hyperest(sims$dat_wide,
#'                     Y_name="Y",
#'                     W_name="W",
#'                     H_names=paste0("H",1:3),
#'                     H0_names=paste0("H",1:3,"_0"),
#'                     H1_names=paste0("H",1:3,"_1"),
#'                     X_names=paste0("X",1:2))
#' # compare true (in-sample) effect, estimates using hypothetical evaluations, and difference in means
#' round(c(sims$tau,
#'         est_hyp$tau,
#'         mean(dat_wide$Y[dat_wide$Y==TRUE]) - mean(dat_wide$Y[dat_wide$Y==FALSE])),
#'       2)
#'
#' @export
simulate_data <- function(n=300, num_h = 5, num_x = 5, frac_coeff_0=0) {

  if (!is.numeric(n)) { stop("n must be an integer") }
  if (!(n==as.integer(n))) { stop("n must be an integer") }

  if (!is.numeric(num_h)) { stop("num_h must be an integer >= 1") }
  if (!(num_h==as.integer(num_h))) { stop("num_h must be an integer >= 1") }
  if (num_h < 1) { stop("num_h must be an integer >= 1") }

  if (!is.numeric(num_x)) { stop("num_x must be an integer >= 0") }
  if (!(num_x==as.integer(num_x))) { stop("num_x must be an integer >= 0") }
  if (num_x < 0) { stop("num_x must be an integer >= 0") }

  if (n < num_h + num_x + 2) { stop("n must be an integer >= num_h + num_x + 2") }

  if (!is.numeric(frac_coeff_0)) { stop("frac_coeff_0 must be a numeric >= 0 and <= 1") }
  if ((frac_coeff_0 < 0) | (frac_coeff_0>1)) { stop("frac_coeff_0 must be a numeric >= 0 and <= 1") }

  # simulate hypothetical evaluations
  Sigma_h <- matrix(0.25,nrow=num_h,ncol=num_h)
  diag(Sigma_h) <- 1
  H0 <- MASS::mvrnorm(n=n, mu=rep.int(0,num_h), Sigma=Sigma_h)
  H1 <- 0.8*H0 + 0.25 + sqrt(1-1*0.8^2) * MASS::mvrnorm(n=n, mu=rep.int(0,num_h), Sigma=Sigma_h)
  # simulate fixed characteristics
  Sigma_x <- matrix(0.25,nrow=num_x,ncol=num_x)
  diag(Sigma_x) <- 1
  X <- (H0 %*% matrix(stats::runif(n=num_h*num_x),nrow=num_h,ncol=num_x)
        + MASS::mvrnorm(n=n, mu=rep.int(0,num_x), Sigma=Sigma_x))

  # treatment assignment
  gamma_h <- stats::runif(n=num_h, min=0.1) * (stats::runif(n=num_h) > frac_coeff_0)
  gamma_x <- stats::runif(n=num_x, min=0.1) * (stats::runif(n=num_x) > frac_coeff_0)
  W <- scale((H1 - H0) %*% gamma_h + X %*% gamma_x) + stats::rnorm(n=n, sd=2)
  W <- W >= stats::median(W)

  # hypothetical evaluation of observed treatment state
  H <- apply(H1,2,function(h) h*W) + apply(H0,2,function(h) h*(1-W))

  # outcome
  beta_h <- 0.2 * gamma_h + stats::runif(n=num_h) * (stats::runif(n=num_h) > frac_coeff_0)
  beta_x <- 0.2 * gamma_x + stats::runif(n=num_x) * (stats::runif(n=num_x) > frac_coeff_0)
  Y0 <- H0 %*% beta_h + X %*% beta_x + stats::rnorm(n=n,sd=1)
  Y1 <- H1 %*% beta_h + X %*% beta_x + stats::rnorm(n=n,sd=1)
  Y <- W*Y1 + (1-W)*Y0

  dat_wide <- as.data.frame(cbind(Y,W,H,H0,H1,X,Y0,Y1))
  col_names <- c("Y","W",paste0("H",1:num_h),paste0("H",1:num_h,"_0"),paste0("H",1:num_h,"_1"))
  if (num_x > 0) {
    col_names <- c(col_names,paste0("X",1:num_x))
  }
  col_names <- c(col_names, "Y0", "Y1")
  colnames(dat_wide) <- col_names

  dat_long <- as.data.frame(rbind(cbind(ifelse(W==FALSE,Y,NA),FALSE,H0,X),
                                  cbind(ifelse(W==TRUE,Y,NA),TRUE,H1,X)))
  col_names <- c("Y","W",paste0("H",1:num_h))
  if (num_x > 0) {
    col_names <- c(col_names,paste0("X",1:num_x))
  }
  colnames(dat_long) <- col_names

  return(list(dat_wide=dat_wide,
              dat_long=dat_long,
              tau_i = Y1 - Y0,
              tau = mean(Y1 - Y0)))
}
