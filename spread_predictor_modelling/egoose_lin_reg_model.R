## RUNNING BASIC LINEAR REGRESSION MODEL TRAINED ON FLORIDA TRAINING DATA AND 
## TESTED ON FLORIDA TEST AND TEXAS DATA FOR COMPARISON WITH RF ALL PARAMS MODEL

# Author = John Gray
# Email = greyjohn15@gmail.com
# Last edit = 31/08/2023

## Load data ----

fl <- read.csv("data/training_data_overall_rf_model_florida_v5.csv")

## create additive linear regression model, fitted to florida training data ----
## that uses the same input parameters as the all parameters rf model

all_params_model <- lm(formula = fl$establishment_score ~ fl$tmax_ring_1 + 
                         fl$prcp_ring_1 + fl$water_ring_1 + fl$forest_ring_1 + 
                         fl$grass_ring_1 + fl$wetland_ring_1 + fl$farming_ring_1 +
                         fl$urban_ring_1 + fl$barren_ring_1 + fl$previous_establishment + 
                         fl$tmax_ring_2 + fl$prcp_ring_2 + fl$water_ring_2 + 
                         fl$forest_ring_2 + fl$grass_ring_2 + fl$wetland_ring_2 + 
                         fl$wetland_ring_2 + fl$farming_ring_2 + fl$urban_ring_2 + 
                         fl$barren_ring_2 + fl$bs_ring_1 + fl$bs_ring_2)


# apply to florida test data ----

# Load in Florida test data
fl_test <- read.csv("data/test_data_overall_rf_model_flordia_v5.csv")

# remove extraneous variables from florida test data
fl_test_all_params_data <- fl_test[,c(-1,-2,-24)]

# generate linear regression model predictions for florida test data
fl_test_all_params_data$predictions <- all_params_model$coefficients[1] +
  all_params_model$coefficients[2]*fl_test_all_params_data$tmax_ring_1 +
  all_params_model$coefficients[3]*fl_test_all_params_data$prcp_ring_1 + 
  all_params_model$coefficients[4]*fl_test_all_params_data$water_ring_1 +
  all_params_model$coefficients[5]*fl_test_all_params_data$forest_ring_1 + 
  all_params_model$coefficients[6]*fl_test_all_params_data$grass_ring_1 +
  all_params_model$coefficients[7]*fl_test_all_params_data$wetland_ring_1 +
  all_params_model$coefficients[8]*fl_test_all_params_data$farming_ring_1 +
  all_params_model$coefficients[9]*fl_test_all_params_data$urban_ring_1 +
  all_params_model$coefficients[10]*fl_test_all_params_data$barren_ring_1 +
  all_params_model$coefficients[11]*fl_test_all_params_data$previous_establishment +
  all_params_model$coefficients[12]*fl_test_all_params_data$tmax_ring_2 +
  all_params_model$coefficients[13]*fl_test_all_params_data$prcp_ring_2 + 
  all_params_model$coefficients[14]*fl_test_all_params_data$water_ring_2 +
  all_params_model$coefficients[15]*fl_test_all_params_data$forest_ring_2 + 
  all_params_model$coefficients[16]*fl_test_all_params_data$grass_ring_2 +
  all_params_model$coefficients[17]*fl_test_all_params_data$wetland_ring_2 +
  all_params_model$coefficients[18]*fl_test_all_params_data$farming_ring_2 +
  all_params_model$coefficients[19]*fl_test_all_params_data$urban_ring_2 +
  all_params_model$coefficients[20]*fl_test_all_params_data$barren_ring_2 +
  all_params_model$coefficients[21]*fl_test_all_params_data$bs_ring_1 +
  all_params_model$coefficients[22]*fl_test_all_params_data$bs_ring_2

## calculate MSE and variance explained for the linear model

# mse
mse_florida_all_params <- mean((fl_test_all_params_data$predictions-fl_test_all_params_data$establishment_score)^2)

# variance explained

mean_test_response <- mean(fl_test$establishment_score)

TSS <- sum((fl_test$establishment_score - mean_test_response)^2)

RSS_all_params <- sum((fl_test_all_params_data$establishment_score - fl_test_all_params_data$predictions)^2)

var_explained_all_params <- (1 - RSS_all_params / TSS)


# save performance statistic results in a csv file

model <- "all_parameters"
var_explained <- var_explained_all_params
mse <- mse_florida_all_params

lm_performance_florida <- data.frame(model, mse, var_explained)

write.csv(lm_performance_florida, "data/lm_performance_florida.csv")

#### REPEAT FOR TEXAS ----

# Load in texas data
tx <- read.csv("data/rf_model_v5_all_parameters_data_texas.csv")

# remove extraneous variables
tx_all_params_data <- tx[,c(-1,-24)]

# make predictions for texas data using trained linear regression model
tx_all_params_data$predictions <- all_params_model$coefficients[1] +
  all_params_model$coefficients[2]*tx_all_params_data$tmax_ring_1 +
  all_params_model$coefficients[3]*tx_all_params_data$prcp_ring_1 + 
  all_params_model$coefficients[4]*tx_all_params_data$water_ring_1 +
  all_params_model$coefficients[5]*tx_all_params_data$forest_ring_1 + 
  all_params_model$coefficients[6]*tx_all_params_data$grass_ring_1 +
  all_params_model$coefficients[7]*tx_all_params_data$wetland_ring_1 +
  all_params_model$coefficients[8]*tx_all_params_data$farming_ring_1 +
  all_params_model$coefficients[9]*tx_all_params_data$urban_ring_1 +
  all_params_model$coefficients[10]*tx_all_params_data$barren_ring_1 +
  all_params_model$coefficients[11]*tx_all_params_data$previous_establishment +
  all_params_model$coefficients[12]*tx_all_params_data$tmax_ring_2 +
  all_params_model$coefficients[13]*tx_all_params_data$prcp_ring_2 + 
  all_params_model$coefficients[14]*tx_all_params_data$water_ring_2 +
  all_params_model$coefficients[15]*tx_all_params_data$forest_ring_2 + 
  all_params_model$coefficients[16]*tx_all_params_data$grass_ring_2 +
  all_params_model$coefficients[17]*tx_all_params_data$wetland_ring_2 +
  all_params_model$coefficients[18]*tx_all_params_data$farming_ring_2 +
  all_params_model$coefficients[19]*tx_all_params_data$urban_ring_2 +
  all_params_model$coefficients[20]*tx_all_params_data$barren_ring_2 +
  all_params_model$coefficients[21]*tx_all_params_data$bs_ring_1 +
  all_params_model$coefficients[22]*tx_all_params_data$bs_ring_2


## calculate MSE and variance explained for linear models

# mse
mse_texas_all_params <- mean((tx_all_params_data$predictions-tx_all_params_data$establishment_score)^2)


# variance explained

mean_test_response <- mean(tx$establishment_score)

TSS <- sum((tx$establishment_score - mean_test_response)^2)

RSS_all_params <- sum((tx_all_params_data$establishment_score - tx_all_params_data$predictions)^2)

var_explained_tx_all_params <- (1 - RSS_all_params / TSS)

# save results in performance comparison table

tx_model <- "all_parameters"
tx_var_explained <- var_explained_tx_all_params
tx_mse <- mse_texas_all_params

lm_performance_texas <- data.frame(tx_model, tx_mse, tx_var_explained)

write.csv(lm_performance_texas, "data/lm_performance_texas.csv")

