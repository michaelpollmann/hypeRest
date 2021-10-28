# hypeRest

**Est**imation of treatment effects using **hyp**othetical **e**valuations (in **R**), as proposed by Bernheim, Björkegren, Naecker, and Pollmann (2021).

To install this package in `R`, run the following commands:  

```R
library(devtools) 
install_github("michaelpollmann/hypeRest")
```



## Usage

**Note: functions, arguments, and returns may change in future versions. Comments are welcome.**

The package contains a function `simulate_data` to quickly simulate data sets with user-specified dimensions and the data.frame shapes used by the functions in this package.
For low-dimensional specifications, use the function `hyperest`, which requires each variable in the regression to be part of the data.frame.
For high-dimensional specifications, use the function `hypererst_arb`, which uses a formula object to specify the relationship between the outcome and hypothetical evaluations and fixed characteristics.

### Low-dimensional specifications

```r
# simulate a data set with 500 observations, 3 hypothetical variables, 2 fixed characteristics
sims <- simulate_data(n=500, num_h=3, num_x=2)
# extract data
dat_wide <- sims$dat_wide
# estimate effects using hypothetical evaluations; match column names generated inside simulate_data
est_hyp <- hyperest(dat_wide,
                    Y_name="Y",
                    W_name="W",
                    H_names=paste0("H",1:3),
                    H0_names=paste0("H",1:3,"_0"),
                    H1_names=paste0("H",1:3,"_1"),
                    X_names=paste0("X",1:2))
# compare true (in-sample) effect, estimates using hypothetical evaluations, and difference in means
round(c(sims$tau,
        est_hyp$tau,
        mean(sims$dat_wide$Y[dat_wide$W==TRUE]) - mean(sims$dat_wide$Y[dat_wide$W==FALSE])),
      2)

# for additional verification, plot predicted treatment effects against individual-level effects
graphics::plot(est_hyp$tau_i,sims$tau_i)
reg <- stats::lm(tau ~ tau_hat, data = data.frame(tau = sims$tau_i, tau_hat = est_hyp$tau_i))
reg$coefficients
graphics::abline(reg)
```


### High-dimensional specifications

```r
# simulate a data set with 500 observations, 20 hypothetical variables, 10 fixed characteristics,
#   where approximately 3/4 of coefficients are 0
num_h <- 20
num_x <- 10
sims <- simulate_data(n=500, num_h=num_h, num_x=num_x, frac_coeff_0 = 3/4)
# extract data
dat_long <- sims$dat_long
# create variable names and formula using all interactions
all_vars <- c(paste0("H",1:num_h),paste0("X",1:num_x))
formula <- paste("Y","~",paste0("(",paste(all_vars,collapse=" + "),")^2"))
W_name <- "W"

est_hyp_arb <- hyperest_arb(formula=formula,
                            W_name="W",
                            dat_long=dat_long,
                            # don't penalize base hypotheticals
                            dont_penalize = paste0("H",1:num_h),
                            # how to pick lambda in cross-validation
                            s="lambda.min",
                            # which optimizer to use for balanceHD::approx.balance
                            args.approx.balance=list(optimizer = "quadprog"))
# compare true (in-sample) effect, estimates using hypothetical evaluations, and difference in means
round(c(sims$tau,
        est_hyp_arb$tau,
        mean(sims$dat_wide$Y[sims$dat_wide$W==TRUE])
          - mean(sims$dat_wide$Y[sims$dat_wide$W==FALSE])),
      2)

# for additional verification, plot predicted treatment effects against individual-level effects
graphics::plot(est_hyp_arb$tau_i,sims$tau_i)
reg <- stats::lm(tau ~ tau_hat, data = data.frame(tau = sims$tau_i, tau_hat = est_hyp_arb$tau_i))
reg$coefficients
graphics::abline(reg)
```



## Brief Description

Our method proceeds in 2 steps.

1. Estimate the relationship between outcome `Y` and hypothetical evaluations and fixed characteristics `H,X`:
        `Y ~ mu(H, X, theta)`
    where `H` are hypothetical evaluations corresponding to the realized treatment state, and `mu` is the specification, for instance linear or including interaction terms.
2. Predict treatment effects for both treatment states:
        `tau = mu(H(1), X, theta) - mu(H(0), X, theta)`

For high-dimensional specifications with LASSO penalty, `hyperest_arb`, an *approximate residual balancing* (Athey et al., 2018) term is added. Approximate residual balancing is implemented in the `R` package [balanceHD](https://github.com/swager/balanceHD).


## Reference

B. Douglas Bernheim, Daniel Björkegren, Jeffrey Naecker, and Michael Pollmann. **Causal Inference from Hypothetical Evaluations**. 2021.

Susan Athey, Guido W. Imbens, and Stefan Wager. **Approximate Residual Balancing: De-Biased Inference of Average Treatment Effects in High Dimensions**. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)* 80, no. 4 (2018): 597-623. [[arXiv](https://arxiv.org/pdf/1604.07125.pdf)]
