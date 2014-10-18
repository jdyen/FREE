# Function Regression in Ecology and Evolution (FREE)

This README is for the R package FREE (Function Regression in Ecology and Evolution).
Yen JDL, et al. (in press) Function regression in ecology and evolution. Methods in Ecology and Evolution.

Copyright &copy; 2014, Jian Yen

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
FREE is a collection of R helper functions for fitting regression models where response and/or predictor variables are functions, rather than scalar quantities. The emphasis is on easy model fitting and straightforward interfaces, with a focus on ecological and evolutionary applications. All functions in FREE are written in R 3.0.2.

Several additional packages are required (see Installation, below); and some of these packages are not available through CRAN.

The FREE package focuses on function-valued response models, but also supports basic function-valued predictor models (single intercept, several variables). Function-valued predictor models are implemented using a purpose-built Gibbs sampler, and this method is still in development, so has not been tested extensively.

The purpose-built Gibbs sampler is intended to replace the other supported methods, but is in development, so other methods have not yet been removed. We are releasing FREE in its current form because functional data analysis has the potential to provide new insight into ecological and evolutionary questions and we hope to receive as much feedback as possible on the FREE package.


Created 13 February 2014

Updated 11 September 2014

*****

## Installation
FREE is distributed as an R package (in source form) but is not currently available through the CRAN. There are two versions of the FREE package, one for Windows systems (tested on Windows 7 and Windows XP) and one for OSX/UNIX systems (tested only on OSX 10.6).

For Windows users: simply download the file FREE_x.x.tar.gz into a local directory on your computer (replace x.x with the appropriate version number).

For OSX users: navigate to the folder OSX/ and download the file FREE_x.x.tar.gz into a local directory on your computer. If you install this version, you do not need to install the rstan, R2WinBUGS, coda or boot packages

Because FREE is installed from source you will also need an appropriate C and C++ compiler installed. Easily installed options are [gcc](https://github.com/kennethreitz/osx-gcc-installer/) (OSX users) and [Rtools](https://github.com/stan-dev/rstan/wiki/Install-Rtools-for-Windows) (Windows users). If you're unsure of whether you need to install a C/C++ compiler, you can try installing the FREE package anyway; if you do not get any errors then no compiler is needed.

#### Experimental
A simple install procedure has been developed, but is currently experimental and may not work on all systems. You will still need to download the appropriate FREE_x.x.tar.gz file (OSX or Windows) to your working directory.

Once you've installed a C/C++ compiler, download the installFREE.R script to your working directory and enter the following into your R console (setting OSX.install=TRUE for OSX users).
```
source("installFREE.R")
installFREE(OSX.install=FALSE)
```
*****

#### Standard approach
Before installing this R package several additional packages must be installed within R:

1. [rstan](http://mc-stan.org/rstan.html) (Windows users only)
2. [INLA](http://www.r-inla.org/)

See the links above for relevant details. Note that rstan and INLA are not currently available through the CRAN and must be installed according to the instructions on their websites.

FREE imports functions from several packages (see Depends and Imports in the DESCRIPTION file) and these packages must be installed for FREE to install and load correctly. There is a separate DESCRIPTION file for OSX users in the OSX/ folder. With the exception of rstan and INLA, all other packages are available through the CRAN and should be easy to install.

If package `maptools` is not installed correctly from the CRAN, try
```
install.packages("maptools", repos="http://R-Forge.R-project.org")
```

Once the appropriate packages have been installed, you simply need to place the file FREE_x.x.tar.gz in the current working directory (where x.x is replaced with the current version number) and use
```
install.packages("FREE_x.x.tar.gz", repos=NULL, type="source")
```
to install the FREE R package. This package then can be loaded in R using `library(FREE)`.

Errors during installation often are related to the installation of required packages. Restarting R and making sure all required packages can be loaded using `library(pkgName)` will identify missing or incorrectly installed packages.

An error might occur if rstan has been installed only for 64-bit architecture. There are two possible solutions to this problem:

- install FREE for 64-bit architecture only:
```
install.packages("FREE_1.0.tar.gz", repos = NULL, type = "source", INSTALL_opts=c("--no-multiarch"))
```
- install rstan for 32- and 64-bit architectures:
```
install.packages('rstan', type = 'source', INSTALL_opts = "--merge-multiarch")
```

Either one of these options should work.

NOTE: The BUGS methods require WinBUGS 1.4 and its jump add-in to be installed locally. See the [WinBUGS](http://www2.mrc-bsu.cam.ac.uk/bugs/) and [rjMCMC](http://www.winbugs-development.org.uk/rjmcmc.html) websites for details. WinBUGS does not install easily on non-Windows operating systems. Note that installing WinBUGS is not necessary to use other methods within package FREE.

NOTE: The BUGS and stan methods are not available in the OSX version. The required packages either are not available or do not install correctly on OSX systems. If you want to get the stan method working and can install and load the [rstan](http://mc-stan.org/rstan.html) package then it is possible to install the Windows version of FREE on OSX systems. The `stan` function within the `rstan` package causes R to crash on all OSX systems we have used.

*****

## Usage
Once FREE has been installed there are two main functions to use: `FREEfit` and `FREEfitCV`. Both of these functions have a formula interface (`FREEfit.formula` and `FREEfitCV.formula`) and information about their use can be found by typing `?FREEfit` and `?FREEfitCV` in the R console.

The `FREEfit` and `FREEfitCV` functions should determine automatically whether the model to be fitted is a function-valued response or function-valued predictor model. At this stage we expect some errors/warnings in this process; see the relevant help files for details.

Mathematical and statistical details for each implementation, as well as a comparison of their performance, are discussed in:
Yen JDL, et al. (in press) Function regression in ecology and evolution. Methods in Ecology and Evolution.

Models fitted using `FREEfit` are of class FREEfit and have several S3 methods available: `print`, `summary`, `plot`, `coef`, `predict`, `fitted` and `residuals`. Models fitted using `FREEfitCV` are of class FREEfitCV and have several S3 methods available: `print`, `summary`, `plot` and `residuals`.

*****

## Feedback
Please send comments and bug reports to
<jdl.yen@gmail.com>