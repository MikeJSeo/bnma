# bnma 1.5.1

* Modified contrast.network.run to avoid error

# bnma 1.5.0

* When calculating treatment effect in the natural scale (i.e. absolute treatment effect coded T in JAGS code), added option to specify covariate/baseline risk that treatment effect should be calculated at.

# bnma 1.4.0

* Added different treatment comparisons for RR, RD, and NNT
* Added Treat.order variable in ume model
* Revised the function network.leverage.plot to include per arm contributions
* Added a plot that compares posterior mean deviance of inconsistency model and consistency model
* When fitting meta regression on baseline risk, made adjustments when the baseline treatment is not the control treatment in each trial (i.e. by adding a fictitious row).

# bnma 1.3.0

* Added finding absolute treatment effect in binomial and normal outcomes
* Added RNG.inits in the network.run function for reproducibility

# bnma 1.2.0

* Added modeling baseline risk and covariate effects when fitting a fixed effects model
* Modified code for network.forest.plot that allows plotting only the comparison with the reference treatment and fixed labeling issue

# bnma 1.1.2

* Added number needed to treat

# bnma 1.1.1

* Updated vignette and added description section

# bnma 1.1.0

* Added an extension to calculate relative risk and risk difference (using odds ratio and estimated placebo event rate) when the outcome is binomial

# bnma 1.0.0

* Added a `NEWS.md` file to track changes to the package.

