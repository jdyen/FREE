# Function Regression in Ecology and Evolution (FREE)

This README is for the R package FREE (Function Regression in Ecology and Evolution) (see also
Yen JDL, et al. (2015) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution, 6:17-26).

Copyright &copy; 2015, Jian Yen

*****

## Licence details
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

*****

## Overview
FREE is a collection of R functions for fitting regression models where response and/or predictor variables are functions, rather than scalar quantities. The emphasis is on easy model fitting and straightforward interfaces, with a focus on ecological and evolutionary applications. All functions in FREE are written in R 3.1.2.

The latest updates to the FREE package have introduced several new models. The FREE package now supports both function-valued response models and function-valued predictor models with multiple function-valued predictor variables. All new models are implemented using a purpose-built Gibbs sampler. All models now also accommodate clustering variables (similar to random intercepts in a mixed effects model), which can account for clustering of subjects in space or time (e.g., several responses measured on one individual or in one site). Models with function-valued response variables now support unequally spaced observations, with correlated errors inversely proportional to the distance between any two observations.

The latest version of the FREE package no longer includes the fda, gamboost, INLA, stan or BUGS methods. The old version of FREE, which includes these methods, is still available, but is no longer being tested or updated.  The latest verison of FREE has many fewer dependencies than the original FREE package, and all models are fitted using the same Gibbs sampling approach. Updated plotting methods and outputs provide greater model flexibility, as well as easier model checking (e.g., using Gelman-Rubin statistics and log-likelihood traces).


Created 13 February 2014

Updated 11 May 2017

*****

## Installation
FREE is distributed as an R package in source form. FREE is not currently available through the CRAN. The FREE package has been tested on Windows 7, Windows 10, Windows XP, and OSX 10.6, 10.7 and 10.10.

Installing FREE requires an appropriate C/C++ compiler to be installed. Easily installed options are [gcc](https://github.com/kennethreitz/osx-gcc-installer/) (OSX users) and [Rtools](https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows) (Windows users). If you're unsure of whether you need to install a C/C++ compiler, you can try installing the FREE package anyway; if you do not get any errors then no compiler is needed.

- the easiest (and recommended) method for installation is to use the devtools package to install directly from the GitHub source:
```
if (!require(devtools)) {
  install.packages("devtools")
}
# install the current version of FREE
devtools::install_github("jdyen/FREE/FREE")
# alternative; use if you want to use the old version of the FREE package
#devtools::install_github("jdyen/FREE/FREE", ref = 'original')
```

NOTE: If you use the old version of FREE, you will need to install the INLA package, which is not available through the CRAN. See the [INLA](http://www.r-inla.org) website for details. If you wish to use the WinBUGS methods in the old version of FREE you will need to install WinBUGS 1.4 and the jump add-in. See the [WinBUGS](http://www2.mrc-bsu.cam.ac.uk/bugs/) and [rjMCMC](http://www.winbugs-development.org.uk/rjmcmc.html) websites for details. WinBUGS does not install easily on non-Windows operating systems. Note that installing WinBUGS is not necessary to install the FREE package or to use other methods within the FREE package.

WARNING: We are no longer testing or updating the old version of FREE; changes to the fda, mboost, INLA, stan or WinBUGS software may cause errors. Please download and install the current version of the software for the most up-to-date version of the FREE package.

*****

## Usage
Once FREE has been installed there are two main functions to use: `FREEfit` and `FREEfitCV`. Both of these functions have a formula interface (`FREEfit.formula` and `FREEfitCV.formula`) and information about their use can be found by typing `?FREEfit` and `?FREEfitCV` in the R console.

The `FREEfit` and `FREEfitCV` functions should determine automatically whether the model to be fitted is a function-valued response or function-valued predictor model. At this stage there might be some errors/warnings in this process; see the relevant help files for details.

Mathematical and statistical details for function regression models are discussed in:
Yen JDL, et al. (2015) Function regression in ecology and evolution. Methods in Ecology and Evolution, 6:17-26.

Several applications of FREE to ecological data are:   
Yen JDL, et al. (2017) Balancing generality and specificity in ecological gradient analysis with species abundance distributions and individual size distributions. Global Ecology and Biogeography, 26:318–332.   
Yen JDL, et al. (2017) How do different aspects of biodiversity change through time? A case study on an Australian bird community. Ecography, 40:642–650.

Models fitted using `FREEfit` are of class FREEfit and have several S3 methods available: `print`, `summary`, `plot`, `coef`, `predict`, `fitted` and `residuals`. Models fitted using `FREEfitCV` are of class FREEfitCV and have several S3 methods available: `print`, `summary`, `plot` and `residuals`. There is also a `plotPretty` function, which displays cleaner plots of fitted model coefficients.

*****

## Feedback
Please send comments and bug reports to
<jdl.yen@gmail.com>
