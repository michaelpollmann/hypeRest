#' Estimate treatment effects using hypothetical evaluations in high-dimensional specifications
#'
#' This function estimates treatment effects and standard errors
#' using hypothetical evaluations in high-dimensional specifications.
#' It uses cross-validated LASSO (glmnet::cv.glmnet) to estimate high-dimensional step 1 regression
#' and approximate residual balancing based on Athey, Imbens, and Wager (2018, implemented in balanceHD).
#'
#' @details
#' Step 1 estimates a LASSO regression based on the formula input, with penalty parameter chosen using cross-validation.
#'
#' Step 1b calculates an approximate residual balancing correction using balanceHD::approx.balance (Athey, Imbens, and Wager, 2018).
#'
#' Step 2 estimates treatment effects as the difference in predictions for the treated and control state:
#'   \deqn{\tau = (\mu(H(1),X) - \mu(H(0),X)) + r}
#' where \eqn{\mu} is the prediction function estimated in step 1 and \eqn{r} is the residual balancing term estimated in step 1b.
#'
#' The function currently estimates the ATE the standard errors of the estimator.
#' It also predicts individual level treatment effects, but no standard errors are provided
#' because theoretical properties are only available for average effects within groups (that grow with sample size).
#'
#'
#' @references
#'
#' B. Douglas Bernheim, Daniel Björkegren, Jeffrey Naecker, and Michael Pollmann. "Causal Inference from Hypothetical Evaluations." 2021.
#'
#' Susan Athey, Guido W. Imbens, and Stefan Wager. "Approximate Residual Balancing: De-Biased Inference of Average Treatment Effects in High Dimensions." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 80, no. 4 (2018): 597-623.
#'
#'
#' @param formula formula (or string interpretable as formula),
#'                 specifies the regression of outcome on hypothetical evaluations and fixed characteristics.
#' @param W_name string, the name of the column of dat_long that holds the treatment indicator.
#' @param dat_long data.frame or equivalent. **Should have 2 rows for each observation**:
#'                 one for data corresponding to the treated state (`[[W_name]] == TRUE`)
#'                 and one for data corresponding to the control state (`[[W_name]] == FALSE`).
#'                 Each row contains outcome data (`NA` if the observations is unobserved under that treatment state),
#'                 the column named `W_name` indicating whether the row corresponds to treatment or control,
#'                 hypothetical evaluations for the treatment state of the row,
#'                 fixed characteristics (should appear in the rows for both treatment states).
#'                 **NOTE**: to interpret individual-level effects, the order of observations
#'                 within `[[W_name]] == TRUE` should be the same as the order of observations within `[[W_name]] == FALSE`.
#'                 For the average treatment effect estimate and standard error, the order of observations is irrelevant.
#' @param dont_penalize (optional) string vector, *variable names* (from the formula, see example) that should not be penalized in the step 1 LASSO regression.
#' @param args.approx.balance (optional) list, arguments passed to balanceHD::approximate.balance, see example for setting the optimization routine
#' @param ... (optional) parameters passed to `glmnet::cv.glmnet` and to `predict` such as number of folds, and criterion for selecting the penalty parameter.
#'
#' @return list of three elements:
#'         \item{tau}{estimate of the average treatment effect (ATE)}
#'         \item{se}{standard error of the ATE estimator}
#'         \item{tau_i}{predictions of individual-level treatment effects (see comment on input argument dat_long to ensure correct order of rows)}
#'
#' @examples
#' # simulate a data set with 500 observations, 20 hypothetical variables, 10 fixed characteristics,
#' #   where approximately 3/4 of coefficients are 0
#' num_h <- 20
#' num_x <- 10
#' sims <- simulate_data(n=500, num_h=num_h, num_x=num_x, frac_coeff_0 = 3/4)
#' # extract data
#' dat_long <- sims$dat_long
#' # create variable names and formula using all interactions
#' all_vars <- c(paste0("H",1:num_h),paste0("X",1:num_x))
#' formula <- paste("Y","~",paste0("(",paste(all_vars,collapse=" + "),")^2"))
#' W_name <- "W"
#'
#' est_hyp_arb <- hyperest_arb(formula=formula,
#'                             W_name="W",
#'                             dat_long=dat_long,
#'                             # don't penalize base hypotheticals
#'                             dont_penalize = paste0("H",1:num_h),
#'                             # how to pick lambda in cross-validation
#'                             s="lambda.min",
#'                             # which optimizer to use for balanceHD::approx.balance
#'                             args.approx.balance=list(optimizer = "quadprog"))
#' # compare true (in-sample) effect, estimates using hypothetical evaluations, and difference in means
#' round(c(sims$tau,
#'         est_hyp_arb$tau,
#'         mean(sims$dat_wide$Y[sims$dat_wide$W==TRUE])
#'           - mean(sims$dat_wide$Y[sims$dat_wide$W==FALSE])),
#'       2)
#'
#' # for additional verification, plot predicted treatment effects against individual-level effects
#' graphics::plot(est_hyp_arb$tau_i,sims$tau_i)
#' reg <- stats::lm(tau ~ tau_hat, data = data.frame(tau = sims$tau_i, tau_hat = est_hyp_arb$tau_i))
#' reg$coefficients
#' graphics::abline(reg)
#'
#' @export
hyperest_arb <- function(formula
                         , W_name
                         , dat_long
                         , dont_penalize=NULL
                         , args.approx.balance=list()
                         , ...
) {
  # regression model to get outcome and covariates for observed treatment state
  reg_obs <- stats::lm(formula=formula, data=dat_long, y=TRUE, x=TRUE, model=FALSE, qr=FALSE)

  # regression model to get (transformations of) covariates for both treatment states
  # use a dependent variable that exists irrespective of observed treatment
  f_reg_all <- paste(W_name,"~",paste(labels(stats::terms(stats::as.formula(formula))),collapse=" + "))
  reg_all <- stats::lm(formula=f_reg_all, data=dat_long, x=TRUE, model=FALSE, qr=FALSE)

  # residual balancing weights
  # balance target
  mean_x1 <- colMeans(reg_all$x[dat_long[[W_name]]==TRUE,])
  mean_x0 <- colMeans(reg_all$x[dat_long[[W_name]]==FALSE,])
  # balance
  gamma1 <- do.call(balanceHD::approx.balance,
                    c(list(M = reg_obs$x,
                           balance.target = mean_x1),
                      args.approx.balance))
  gamma0 <- do.call(balanceHD::approx.balance,
                    c(list(M = reg_obs$x,
                           balance.target = mean_x0),
                      args.approx.balance))

  # cross-validated LASSO
  reg_lasso <- glmnet::cv.glmnet(x=reg_obs$x, y=reg_obs$y,
                                 penalty.factor = !(colnames(reg_obs$x) %in% dont_penalize))
  # predict outcomes
  Y1_hat <- stats::predict(reg_lasso, newx = reg_all$x[dat_long[[W_name]]==TRUE,], ...)
  Y0_hat <- stats::predict(reg_lasso, newx = reg_all$x[dat_long[[W_name]]==FALSE,], ...)

  # then re-weight the residuals
  Y_hat <- stats::predict(reg_lasso, newx = reg_obs$x, ...)
  w_resid <- sum((gamma1-gamma0)*(reg_obs$y - Y_hat))

  # estimate effects and standard error
  tau_i <- Y1_hat - Y0_hat + w_resid
  tau <- mean(tau_i)
  se <- sqrt(sum((gamma1-gamma0)^2*(reg_obs$y - Y_hat)^2) *
               # degrees of freedom correction
               length(gamma1) / max(1, (length(gamma1) - sum(stats::coef(reg_lasso,...) != 0))))


  out <- list(tau = as.numeric(tau),
              se = as.numeric(se),
              tau_i = as.numeric(tau_i))

  return(out)
}
