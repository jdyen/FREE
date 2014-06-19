# Function Regression in Ecology and Evolution (FREE)

This is a to-do list for the R package FREE (Function Regression in Ecology and Evolution).

Maintainer: Jian Yen

*****

### Tasks to complete:
#### Major tasks
- Add scalar-valued responses with function-valued predictors using the FREE method
- Add subject-level errors with the FREE method
- Add random effects with the FREE method
- Add function-valued response and predictor models using the FREE method
- Add model families for the FREE method
- Calculate and return DIC for the FREE method

#### Minor tasks
- Check rstan install on UNIX systems. Currently crashes on stan() call on OSX
- Implement model family specification for stan and BUGS models
- Implement error families for fda methods
- Check and fix errors in different basis functions for fda method
- Add identity basis in stan method
- Test memory use for MakeMassiveTridiag and gamboost with AR1 errors. How many sites are too many?
- Add a preferred citation to the readme file
- Add citations for other software packages
- Check errors in the number of bins used by the gamboost method
