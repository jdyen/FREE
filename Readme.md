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

Several additional packages are required (see Installation, below); some of these packages are not available through CRAN.

The FREE package has been updated to include several new methods. The FREE package now supports both function-valued response models and function-valued predictor models with one function-valued predictor variable. Function-valued predictor models are implemented using a purpose-built Gibbs sampler. Both types of models now also accommodate clustering variables (similar to random intercepts in a mixed effects model), which can account for clustering of subjects in space or time (e.g., several responses measures on one individual or in one site).

The new methods are implemented with a purpose-built Gibbs sampler, which is intended to replace the other methods. The Gibbs sampler is still being tested, so the other methods have not yet been removed. We are releasing FREE in its current form because functional data analysis has the potential to provide new insight into ecological and evolutionary questions and we hope to receive as much feedback as possible on the FREE package.


Created 13 February 2014

Updated 14 October 2015

*****

## Installation
FREE is distributed as an R package in source and binary form. FREE is not currently available through the CRAN. The FREE package has been tested on Windows 7, Windows XP, and OSX 10.6, 10.7 and 10.10.

There are four installation options. The first three require an appropriate C and C++ compiler installed. Easily installed options are [gcc](https://github.com/kennethreitz/osx-gcc-installer/) (OSX users) and [Rtools](https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows) (Windows users). If you're unsure of whether you need to install a C/C++ compiler, you can try installing the FREE package anyway; if you do not get any errors then no compiler is needed.

- install/load the devtools package and install directly from the source package on GitHub:
```
devtools::install_github("jdyen/FREE/FREE")
```
*****

- Download the source package FREE_x.x.tar.gz into a local directory on your computer (replace x.x with the appropriate version number, currently 2.0). Download and source the installFREE.R script:
```
source("installFREE.R")
installFREE()
```
*****

- Install all dependencies manually. First, install the INLA package (not available through CRAN; details at [INLA](http://www.r-inla.org/)). FREE imports functions from several packages (see Depends and Imports in the DESCRIPTION file) and these packages must be installed for FREE to install and load correctly. With the exception of INLA, all other packages are available through the CRAN and should be easy to install. Once the appropriate packages have been installed, you simply need to place the file FREE_x.x.tar.gz in the current working directory (where x.x is replaced with the current version number) and use
```
install.packages("FREE_x.x.tar.gz", repos=NULL, type="source")
```
to install the FREE R package. This package then can be loaded in R using `library(FREE)`.

*****

- install all dependencies manually and install the FREE package from binary. For OSX systems, use
```
install.packages("FREE_x.x.tgz", repos=NULL)
```
and, for Windows systems, use
```
install.packages("FREE.zip", repos=NULL)
```

*****

Errors during installation often are related to the installation of required packages. Restarting R and making sure all required packages can be loaded using `library(pkgName)` will identify missing or incorrectly installed packages.

An error might occur if rstan has been installed only for 64-bit architecture. There are two possible solutions to this problem:

- install FREE for 64-bit architecture only:
```
install.packages("FREE_x.x.tar.gz", repos = NULL, type = "source", INSTALL_opts=c("--no-multiarch"))
```
- install rstan for 32- and 64-bit architectures:
```
install.packages('rstan', type = 'source', INSTALL_opts = "--merge-multiarch")
```

Either one of these options should work.

NOTE: The BUGS methods require WinBUGS 1.4 and its jump add-in to be installed locally. See the [WinBUGS](http://www2.mrc-bsu.cam.ac.uk/bugs/) and [rjMCMC](http://www.winbugs-development.org.uk/rjmcmc.html) websites for details. WinBUGS does not install easily on non-Windows operating systems. Note that installing WinBUGS is not necessary to use other methods within package FREE.

*****

## Usage
Once FREE has been installed there are two main functions to use: `FREEfit` and `FREEfitCV`. Both of these functions have a formula interface (`FREEfit.formula` and `FREEfitCV.formula`) and information about their use can be found by typing `?FREEfit` and `?FREEfitCV` in the R console.

The `FREEfit` and `FREEfitCV` functions should determine automatically whether the model to be fitted is a function-valued response or function-valued predictor model. At this stage we expect some errors/warnings in this process; see the relevant help files for details.

Mathematical and statistical details for each implementation, as well as a comparison of their performance, are discussed in:
Yen JDL, et al. (2015) Function regression in ecology and evolution. Methods in Ecology and Evolution, 6:17-26.

Models fitted using `FREEfit` are of class FREEfit and have several S3 methods available: `print`, `summary`, `plot`, `coef`, `predict`, `fitted` and `residuals`. Models fitted using `FREEfitCV` are of class FREEfitCV and have several S3 methods available: `print`, `summary`, `plot` and `residuals`.

*****

## Feedback
Please send comments and bug reports to
<jdl.yen@gmail.com>