# intreg2
Generalized interval regression

This is a Stata package that fits a model of depvar on indepvars using maximum likelihood where the dependent variable can be point data, interval data, right-censored data, or left-censored data.This is a generalization of the built in STATA command intreg and will yield identical estimates if the normal distribution option is used.  Unlike intreg, intreg2 allows the underlying variable of interest to be distributed according to a more general distribution including all distributions in the Skewed Generalized T family and Generalized Beta of the Second Kind tree.


Installation

The development version from github:

net install gb2reg, from(https://raw.githubusercontent.com/jakesurf/intreg2/master/) force
