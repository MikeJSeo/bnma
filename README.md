
<!-- README.md is generated from README.Rmd. Please edit that file -->
bnma (Bayesian network meta analysis)
=====================================

This Package is quite similar and has been inspired by the popular Bayesian NMA package gemtc. Some additional features of bnma include:

-   bnma allows binomial, normal, and multinomial outcomes.
-   bnma adds modelling baseline risk
-   bnma automatically checks for convergence using gelman-rubin diagnostics before sampling full iteration amount
-   bnma generates reasonable and dispersed initial values if left unspecified
-   bnma adds an option to report relative risk, risk difference, and number needed to treat using external placebo event rate when the outcome is binomial
-   bnma adds new summary plots: rank plots, forest plot, etc

Please see the vignette for detailed examples of using this package.
