# OOSE_multiRCT

This repo contains code for methods that estimate confidence intervals for the conditional average treatment effect (CATE) based on a mixed effects meta-analysis of multiple RCTs. The purpose of this work is to be able to effectively assess uncertainty in a target sample that does not come from a particular trial, when we allow for between-study heterogeneity in the meta-analysis model. Code in the R folder does the following:

-   MDD_Generation_OOSEst.R: Generates simulated datasets based somewhat off of trials comparing treatments for major depressive disorder (MDD). Contains the gen_mdd() function that creates datasets with varying covariate distributions according to trial membership.

-   MA_OOSEst.R: Fits meta-analysis models to the simulated data and then creates confidence intervals for every individual's CATE according to three different approaches:

    1.  Considering the variance of fixed effect coefficients only using the glht() function;

    2.  Manually computing the variance of each CATE according to the fixed and random effect distributions;

    3.  Bootstrapping the fixed and random effects to estimate the CATE and then derive a 95% confidence interval

-   Simulations_OOSEst.R: Applies the functions from the previous two files using different data generation parameters and scenarios, repeats multiple iterations, and saves the results.
