# Load data from raw, assemble into data set, place in interim
source("./data/make_dataset.R")

# Subset data into interim data sets, place in interim
# source("./data/make_subsets.R")

## Build features from subsets, place in processed
source("./features/build_features.R")
# source('./features/build_features_quittozero.R')


## Analysis 1 ------
# fit/predict analysis 1 (LCMM)
source("./models/smoking_trajectories.R")
source("./models/tableone_lca.R")

## Analysis 2 ---- 
# fit/predict class membership
source("./models/predict_class_ordinal.R")
source("./models/predict_class_results_ord.R")
# source("./models/predict_class_site.R")

## Analysis 3 ---- 
# fit/predict smoking cessation at 1-year follow-up
source("./models/smoking_cessation_logistic.R")
source("./models/smoking_cessation_results_logistic.R")

## Listed Protocol 11ppm -------
# source("./models/smoking_cessation_logistic_11ppm.R")
# source("./models/smoking_cessation_results_11ppm.R")