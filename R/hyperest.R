#' Estimate treatment effects using hypothetical evaluations
#'
#' This function estimates treatment effects and standard errors using hypothetical evaluations.
#'
#' @details
#' Step 1 estimates a linear regression:
#'   \deqn{Y ~ H + X}
#' where the outcome variable Y is given in the input string Y_name,
#' the hypothetical evaluations for the realized treatment state are H given in the input string vector of column names H_names,
#' and optional fixed characteristics (not varying by treatment state) are given by the input string vector of column names X_names.
#'
#' Step 2 estimates treatment effects as the difference in predictions for the treated and control state:
#'   \deqn{\tau = (H(1) - H(0)) \beta}
#' where H(1) are the hypothetical evaluations for the treated state given by the input string column names H1_names,
#' H(0) are the hypothetical evaluations for the control state given by the input column names H0_names,
#' and \eqn{\beta} are the estimated coefficients on H in the regression of step 1.
#'
#' The function estimates the ATE, ATT, ATC, and their standard errors.
#' It also predicts individual level treatment effects, but no standard errors are provided
#' because theoretical properties are only available for average effects within groups (that grow with sample size).
#'
#' @references
#'
#' B. Douglas Bernheim, Daniel Björkegren, Jeffrey Naecker, and Michael Pollmann. "Causal Inference from Hypothetical Evaluations." 2021.
#'
#'
#' @param dat_wide dataframe or equivalent, holding data on separate (assumed to be independent) observations as rows,
#'           with column corresponding to different variables.
#' @param Y_name string, the name of the column of dat_wide that holds the outcome data.
#' @param W_name string, the name of the column of dat_wide that holds the treatment indicator. It is only used to estimate the ATT and ATC.
#' @param H_names string vector, the names of columns holding hypothetical evaluations for the realized treatment state.
#' @param H0_names string vector, the names of columns holding hypothetical evaluations in the absence of treatment.
#'           IMPORTANT: order of variables must be the same as in H_names.
#' @param H1_names string vector, the names of columns holding hypothetical evaluations in the presence of treatment.
#'           IMPORTANT: order of variables must be the same as in H_names.
#' @param X_names (optional) string vector, the names of columns holding fixed characteristics that do not vary by treatment state.
#' @param ... (optional) parameters passed to vcovHC of the sandwich package,
#'           determining how to calculate heteroskedasticity-robust standard errors for the linear regression in step 1.
#'           If no parameters are supplied, the defaults of sandwich::vcovHC are used.
#'
#' @return list:
#'         \item{tau}{estimate of the average treatment effect}
#'         \item{tau_att}{estimate of the average treatment effect on the treated (ATT)}
#'         \item{tau_atc}{estimate of the average treatment effect on the control (ATC)}
#'         \item{se}{standard error of the ATE estimator}
#'         \item{se_att}{standard error of the ATT estimator}
#'         \item{se_atc}{standard error of the ATC estimator}
#'         \item{se_naive_only_s1}{standard error of the ATE estimator taking into account only the step 1 regression variance (holding fixed hypothetical evaluations in step 2)}
#'         \item{se_naive_only_s2}{standard error of the ATE estimator taking into account only the step 2 variance due to sampling of hypothetical evaluations (holding fixed the step 1 coefficient estimates)}
#'         \item{reg}{`lm` regression object of step 1}
#'         \item{tau_i}{predictions of individual-level treatment effects}
#' @examples
#' # simulate a data set with 500 observations, 3 hypothetical variables, 2 fixed characteristics
#' sims <- simulate_data(n=500, num_h=3, num_x=2)
#' # extract data
#' dat_wide <- sims$dat_wide
#' # estimate effects using hypothetical evaluations; match column names generated inside simulate_data
#' est_hyp <- hyperest(dat_wide,
#'                     Y_name="Y",
#'                     W_name="W",
#'                     H_names=paste0("H",1:3),
#'                     H0_names=paste0("H",1:3,"_0"),
#'                     H1_names=paste0("H",1:3,"_1"),
#'                     X_names=paste0("X",1:2))
#' # compare true (in-sample) effect, estimates using hypothetical evaluations, and difference in means
#' round(c(sims$tau,
#'         est_hyp$tau,
#'         mean(sims$dat_wide$Y[dat_wide$W==TRUE]) - mean(sims$dat_wide$Y[dat_wide$W==FALSE])),
#'       2)
#'
#' # for additional verification, plot predicted treatment effects against individual-level effects
#' graphics::plot(est_hyp$tau_i,sims$tau_i)
#' reg <- stats::lm(tau ~ tau_hat, data = data.frame(tau = sims$tau_i, tau_hat = est_hyp$tau_i))
#' reg$coefficients
#' graphics::abline(reg)
#'
#' @export
hyperest <- function(dat_wide
                     , Y_name
                     , W_name
                     , H_names
                     , H0_names
                     , H1_names
                     , X_names=NULL
                     , ...
) {

  # number of observations
  N <- nrow(dat_wide)

  # create formula for linear regression
  # regress Y on S
  formula_string <- paste(Y_name,
                          paste(H_names, collapse=' + '),
                          sep=' ~ ')
  # control for X if specified
  if (!is.null(X_names)) {
    formula_string <- paste(formula_string,
                            paste(X_names, collapse=' + '),
                            sep=' + ')
  }

  formula <- stats::as.formula(formula_string)

  reg <- stats::lm(formula, data=dat_wide)
  # summary(reg)


  # get the regression coefficients of S
  gamma <- as.vector(reg$coefficients[H_names])
  # get the difference in subjective responses
  S1mS0 <- as.matrix(dat_wide[,H1_names] - dat_wide[,H0_names],ncol=length(H_names))
  # individual predictions
  tau_i <- S1mS0 %*% gamma
  # difference in conditional expectations
  tau <- mean(tau_i)
  tau_att <- mean(tau_i[dat_wide[[W_name]]==TRUE])
  tau_atc <- mean(tau_i[dat_wide[[W_name]]==FALSE])


  # estimate standard error

  # expectation of squared first moment
  g <- as.vector(tau - S1mS0 %*% gamma)
  Eg2 <-  mean(g^2)

  # expectation of derivative of first moment w.r.t. gamma (coefficients of S)
  Edgdgamma <- as.vector(colMeans(S1mS0))

  # robust standard errors (covariance matrix) of gamma hat
  Vgamma <- N*sandwich::vcovHC(reg,...)[2:(length(H_names)+1),2:(length(H_names)+1)]

  # expectation of derivative of first moment w.r.t. all regression coefficients
  Edgdtheta <- c(Edgdgamma, 0, rep(0,times=length(X_names)))  # zeros correspond to constant term and X

  # expectation of derivative of second moments (OLS) w.r.t. all regression coefficients
  if (!is.null(X_names)) {
    x <- as.matrix(cbind(dat_wide[,H_names],rep(1,N), dat_wide[,X_names]))
  }
  else {
    x <- as.matrix(cbind(dat_wide[,H_names],rep(1,N)))
  }
  Edmdtheta <- - t(x)%*%x / N

  # expectation of cross-term of moments
  e <- reg$residuals
  m <- x * e
  Egm <- colMeans(g * m)

  # putting it all together
  V_tau <- Eg2 + Edgdgamma %*% Vgamma %*% Edgdgamma - 2 * Edgdtheta %*% solve(Edmdtheta) %*% Egm

  # s.e. for ATT
  p <- mean(dat_wide[[W_name]])
  Eg2W <-  mean((dat_wide[[W_name]]*g)^2)
  EgmW <- colMeans(dat_wide[[W_name]] * g * m)
  EdgdgammaW <- as.vector(colSums(as.matrix(S1mS0[dat_wide[[W_name]]==TRUE,],ncol=length(H_names)))/N)
  EdgdthetaW <- c(EdgdgammaW, 0, rep(0,times=length(X_names)))  # zeros correspond to constant term and X
  V_tau_att <- (Eg2W + EdgdgammaW %*% Vgamma %*% EdgdgammaW
                   - 2 * EdgdthetaW %*% solve(Edmdtheta) %*% EgmW)/(p^2)
  # s.e. for ATC
  p <- mean(dat_wide[[W_name]])
  Eg2W0 <-  mean(((1-dat_wide[[W_name]])*g)^2)
  EgmW0 <- colMeans((1-dat_wide[[W_name]]) * g * m)
  EdgdgammaW0 <- as.vector(colSums(as.matrix(S1mS0[dat_wide[[W_name]]==FALSE,],ncol=length(H_names)))/N)
  EdgdthetaW0 <- c(EdgdgammaW0, 0, rep(0,times=length(X_names)))  # zeros correspond to constant term and X
  V_tau_atc <- (Eg2W0 + EdgdgammaW0 %*% Vgamma %*% EdgdgammaW0
                   - 2 * EdgdthetaW0 %*% solve(Edmdtheta) %*% EgmW0)/((1-p)^2)


  # variance of naive estimator
  tau_se_naive_only_s1 <- as.vector(sqrt(as.vector(colMeans(S1mS0)) %*% Vgamma %*% as.vector(colMeans(S1mS0)) / N))
  tau_se_naive_only_s2 <- as.vector(sqrt(t(gamma)%*%stats::var(S1mS0)%*%gamma / N))

  out <- list(tau = as.numeric(tau)
              , tau_att = as.numeric(tau_att)
              , tau_atc = as.numeric(tau_atc)
              , se = as.numeric(sqrt(diag((V_tau/N))))
              , se_att = as.numeric(sqrt(diag((V_tau_att/N))))
              , se_atc = as.numeric(sqrt(diag((V_tau_atc/N))))
              , se_naive_only_s1 = as.numeric(tau_se_naive_only_s1)
              , se_naive_only_s2 = as.numeric(tau_se_naive_only_s2)
              , reg = reg
              , tau_i = as.numeric(tau_i))
  return(out)
}
