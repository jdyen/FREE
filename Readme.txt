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


FREE is a collection of helper functions for fitting regression models where response and/or predictor variables are functions, rather than scalar quantities. The emphasis is on easy model fitting and straightforward interfaces, with a focus on ecological and evolutionary applications. All functions in FREE are written in R 3.0.1.

Several additional packages are required (see Installation, below); and some of these packages are not available through CRAN.

Currently only functional response variables are considered but future updates will introduce functional predictors.

A purpose-built MCMC algorithm is in development to replace the need for other packages. We intend to incorporate both functional predictors and responses into this package. We are releasing FREE in its current form because functional data analysis has the potential to provide new insight into ecological and evolutionary questions and the testing and development of our purpose-built Bayesian approach could take several years.


Created 13 February 2014
Updated 27 February 2014


Installation
============
FREE is distributed as an R package but is not currently available through the CRAN. Simply download the file FREE_x.x.tar.gz into a local directory on your computer (replace x.x with the appropriate version number).

Before installing this R package several additional packages must be installed within R:
1. rstan <http://mc-stan.org/rstan.html>
2. INLA <http://www.r-inla.org/>
See the links above for relevant details. Note that rstan and INLA are not currently available through the CRAN and must be installed according to the instructions on their websites.

FREE imports functions from several packages (see DESCRIPTION file) and these packages must be installed for FREE to install and load correctly. With the exception of rstan and INLA, all other packages are available through the CRAN and should be easy to install.

If package 'maptools' is not installed correctly from the CRAN, try
> install.packages("maptools", repos="http://R-Forge.R-project.org")

Once the appropriate packages have been installed, you simply need to install the FREE package in R, using
>> install.packages("FREE_x.x.tar.gz", repos=NULL, type="source")
when FREE_x.x.tar.gz is in the current working directory (and x.x is replaced with the appropriate version number).

NOTE: The BUGS methods require WinBUGS 1.4 and the jump add-in to be installed locally. See <http://www.mrc-bsu.cam.ac.uk/bugs/winbugs/contents.shtml> for details. WinBUGS does not install easily onto non-Windows operating systems. Note that installing WinBUGS is not necessary to use other methods within package FREE.


Usage
=====
Once FREE has been installed there are two main functions to use: FREEfit and FREEfitCV. Both of these functions have a formula interface (FREEfit.formula and FREEfitCV.formula) and information about their use can be found by typing ?FREEfit and ?FREEfitCV in the R console.

Mathematical and statistical details for each implementation, as well as a comparison of their performance, are discussed in:
Yen JDL, et al. (in preparation) Functional regression in ecology and evolution.

Models fitted using FREEfit are of class FREEfit and have several S3 methods available: print, summary, plot, coef, predict, fitted and residuals. Models fitted using FREEfitCV are of class FREEfitCV and have several S3 methods available: print, summary, plot and residuals.


Feedback
========
Please send comments and bug reports to
jdl.yen@gmail.com

