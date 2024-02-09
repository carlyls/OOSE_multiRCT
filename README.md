# OOSE_multiRCT

This repo contains code for methods that estimate prediction intervals for the conditional average treatment effect (CATE) in a target sample based on a model developed using multiple randomized controlled trials (RCTs). The purpose of this work is to be able to effectively assess uncertainty in a target sample that does not come from a particular trial, when we allow for between-study heterogeneity in the trial model. Code in the R folder does the following:

-   MDD_Generation_OOSEst.R: Generates simulated datasets based off of real trials comparing treatments for major depressive disorder (MDD). Contains the gen_mdd() function that creates datasets with varying covariate distributions according to trial membership.

-   MA_OOSEst.R: Constructs prediction intervals for the treatment effect based on a mixed effects meta-analysis.

-   CF_OOSEst.R: Constructs prediction intervals for the treatment effect based on a causal forest with pooling with trial indicator.

-   BART_OOSEst.R: Constructs prediction intervals for the treatment effect based on Bayesian Additive Regression Trees (BART) with pooling with trial indicator.

-   Comparing_OOSEst.R: Runs all modeling approaches on simulated data according to different data generation parameters and outputs performance results.

-   Simulations_OOSEst.R: Applies the functions from all files using different data generation parameters and scenarios, repeats multiple iterations, and saves the results.

-   OOSEst_cluster.sh: Applies the simulations in a computing cluster.
