Function Regression in Ecology and Evolution (FREE)
===========================================================================
This README is for the software FREE (Function Regression in Ecology and Evolution).

Copyright (C) 2014, Jian Yen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
===========================================================================


FREE is a collection of helper functions for fitting regression models where response and/or predictor variables are functions, rather than scalar quantities. All functions in FREE are written in R 3.0.1.

Several additional packages are required (see Installation, below); and some of these packages are not available through CRAN.

Currently only functional response variables are considered but future updates will introduce functional predictors and a more complete formula interface.


Created 13 February 2014
Updated 21 February 2014


Installation
============
The distribution is 16 files that must be placed in a local directory. The R working directory must be set to this local directory and the file "FREEmain.R" must be sourced within R.

Several additional packages must be installed within R:
1. rstan <http://mc-stan.org/rstan.html>
2. fda <http://cran.r-project.org/web/packages/fda/index.html>
3. mboost <http://cran.r-project.org/web/packages/mboost/index.html>
4. INLA <http://www.r-inla.org/>
5. R2WinBUGS* <http://cran.r-project.org/web/packages/R2WinBUGS/index.html>
6. All dependencies for these packages
See the links above for relevant details. Note that rstan and INLA are not currently available through the CRAN and must be installed according to the instructions on their websites.

*R2WinBUGS is a package for calling the program WinBUGS from R. WinBUGS <http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/contents.shtml> or OpenBUGS still need to be installed locally on your computer. This only works on computers running Windows operating systems. Installing WinBUGS is not necessary to use other methods within package FREE.

Usage
=====
Once "FREEmain.R" is sourced in R there are two main functions, six methods for objects of class FREEfit, six helper functions and six internal functions for model fitting.

Exact details for each implementation are discussed in:
Yen JDL, et al. (in preparation) Functional regression models for ecology and evolution.

The two main functions are:
1. FREEfit       - fit a single function regression model
2. FREEfitCV     - perform k-fold cross validation for function regression                   
                   models

Both FREEfit and FREEfitCV take a formula and a named list of all required data. The argument "method" selects the method used for fitting the function regression model. Options are "gamboost" (default), "fda", "INLA", and "stan". Two further options ("BUGS" and "FREE") are not implemented yet but will introduce alternative MCMC approaches.

Other arguments are passed directly to the model fitting procedure and are specific to "method". See the internal functions (below) for further details.

Outputs:
FREEfit outputs an object of class FREEfit, which contains fitted values, observed values, coefficient estimates, coefficient standard deviations (where available), r2, the model family, and the bin data.

FREEfitCV outputs an object of class FREEfitCV, which contains fitted values, observed values and the cross-validation r2 value. No methods are implemented for FREEfitCV objects yet.

FREEfit objects are best explored using the following methods:
1. coef       - extract model coefficients and upper and lower bounds
2. plot       - plot fitted vs observed as well as all coefficients
3. print      - print the model call, coefficients and r2
4. summary    - summary the model coefficients (mean +/- sd) and print the
                resulting summary, including r2 and model call
5. predict    - predict the response variable for new predictor data
6. fitted     - returns the fitted values of the model
7. residuals

Inputs for all FREEfit methods are the same as those for the generic functions they are mimicking.

Diagnostic functions:
1. FREEfitTest    - runs a set of diagnostic tests for time and
                    cross-validation performance using different methods
2. FREEfitTime    - records computing time for a method or vector of
                    methods fitted repeatedly to the same dataset

The helper functions are:
1. ConvertB2A         - convert data from a matrix format to a vector
                        format
2. MakeInlaFormula    - create a formula for a specific INLA model
3. MakeMboostFormula  - create a formula for a specific gamboost model
4. MakeTridiag        - create a tridiagonal from three vectors or numerics
5. MakeMassiveTridiag - create a large pseudo-diagonal matrix from a
                        smaller matrix
6. stanCodeDefault    - code default for stan MCMC model
7. MakeCovMat

The internal model fitting functions are:
1. FREEgamboost     - fit a FREE model using gamboost from package Mboost
  Optional parameters and data:
    bins
    coord.data
    model.int
    model.pred
    model.site
    rand.eff
    family
    spatial
    cvm.set
    weights
    deg.m.int
    df.m.int
    diff.m.int
    deg.m.pred
    df.m.pred
    diff.m.pred
    df.spat
    df.mrf
    nu.m
    mstop
    trace
    offset
2. FREEfda          - fit a FREE model using the fda package
  Optional parameters and data:
    bins
    y.basis
    nbasis.y
    norder.y
    beta.basis
    nbasis.beta
    norder.beta
    loglam
3. FREEinla         - fit a FREE model using the Rinla package
  Optional parameters and data:
    bins
    model.int
    model.pred
    model.site
    model.eij
    family.resp
    order
    diag.int
    diag.pred
    prec.prior
    control.predictor.set
    control.compute.set
    group.mean
    n.groups
    group.vars
    n.groups.var
    verbose
4. FREEstan         - fit a FREE model using the rstan package
  Optional parameters and data:
    bins
    stan.file
    stan.model
    Kt
    iid.er
    n.chains
    n.iters
    n.burnin
    n.thin
    verbose
5. FREEbugs         - *NOT YET IMPLEMENTED* - fit using WinBUGS
6. FREEfree         - *NOT YET IMPLEMENTED* - fit using compiled MCMC code


Feedback
========
Please send comments and bug reports to
jdl.yen@gmail.com

