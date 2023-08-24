## Running basic linear regression models for comparison with data

## Load data ----

fl <- read.csv("data/training_data_overall_rf_model_florida_v5.csv")

## create models

all_params_model <- lm(formula = fl$establishment_score ~ fl$tmax_ring_1 + 
                         fl$prcp_ring_1 + fl$water_ring_1 + fl$forest_ring_1 + 
                         fl$grass_ring_1 + fl$wetland_ring_1 + fl$farming_ring_1 +
                         fl$urban_ring_1 + fl$barren_ring_1 + fl$previous_establishment + 
                         fl$tmax_ring_2 + fl$prcp_ring_2 + fl$water_ring_2 + 
                         fl$forest_ring_2 + fl$grass_ring_2 + fl$wetland_ring_2 + 
                         fl$wetland_ring_2 + fl$farming_ring_2 + fl$urban_ring_2 + 
                         fl$barren_ring_2 + fl$bs_ring_1 + fl$bs_ring_2)

summary(all_params_model)

abiotic_pe_model <- lm(formula = fl$establishment_score ~ fl$tmax_ring_1 + 
                         fl$prcp_ring_1 + fl$water_ring_1 + fl$forest_ring_1 + 
                         fl$grass_ring_1 + fl$wetland_ring_1 + fl$farming_ring_1 +
                         fl$urban_ring_1 + fl$barren_ring_1 + fl$previous_establishment +
                         fl$tmax_ring_2 + fl$prcp_ring_2 + fl$water_ring_2 + 
                         fl$forest_ring_2 + fl$grass_ring_2 + fl$wetland_ring_2 + 
                         fl$wetland_ring_2 + fl$farming_ring_2 + fl$urban_ring_2 + 
                         fl$barren_ring_2)

abiotic_model <- lm(formula = fl$establishment_score ~ fl$tmax_ring_1 + 
                      fl$prcp_ring_1 + fl$water_ring_1 + fl$forest_ring_1 + 
                      fl$grass_ring_1 + fl$wetland_ring_1 + fl$farming_ring_1 +
                      fl$urban_ring_1 + fl$barren_ring_1 +
                      fl$tmax_ring_2 + fl$prcp_ring_2 + fl$water_ring_2 + 
                      fl$forest_ring_2 + fl$grass_ring_2 + fl$wetland_ring_2 + 
                      fl$wetland_ring_2 + fl$farming_ring_2 + fl$urban_ring_2 + 
                      fl$barren_ring_2)

biotic_pe_model <- lm(formula = fl$establishment_score ~ fl$previous_establishment + 
                       fl$bs_ring_1 + fl$bs_ring_2)

biotic_model <- lm(formula = fl$establishment_score ~ fl$bs_ring_1 + fl$bs_ring_2)

# apply to florida training data

fl_test <- read.csv("data/test_data_overall_rf_model_flordia_v5.csv")

fl_test_all_params_data <- fl_test[,c(-1,-2,-24)]

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

fl_test_abiotic_pe_data <- fl_test[,c(-1,-22,-23,-24)]

fl_test_abiotic_pe_data$predictions <- abiotic_pe_model$coefficients[1] +
  abiotic_pe_model$coefficients[2]*fl_test_abiotic_pe_data$tmax_ring_1 +
  abiotic_pe_model$coefficients[3]*fl_test_abiotic_pe_data$prcp_ring_1 + 
  abiotic_pe_model$coefficients[4]*fl_test_abiotic_pe_data$water_ring_1 +
  abiotic_pe_model$coefficients[5]*fl_test_abiotic_pe_data$forest_ring_1 + 
  abiotic_pe_model$coefficients[6]*fl_test_abiotic_pe_data$grass_ring_1 +
  abiotic_pe_model$coefficients[7]*fl_test_abiotic_pe_data$wetland_ring_1 +
  abiotic_pe_model$coefficients[8]*fl_test_abiotic_pe_data$farming_ring_1 +
  abiotic_pe_model$coefficients[9]*fl_test_abiotic_pe_data$urban_ring_1 +
  abiotic_pe_model$coefficients[10]*fl_test_abiotic_pe_data$barren_ring_1 +
  abiotic_pe_model$coefficients[11]*fl_test_abiotic_pe_data$previous_establishment +
  abiotic_pe_model$coefficients[12]*fl_test_abiotic_pe_data$tmax_ring_2 +
  abiotic_pe_model$coefficients[13]*fl_test_abiotic_pe_data$prcp_ring_2 + 
  abiotic_pe_model$coefficients[14]*fl_test_abiotic_pe_data$water_ring_2 +
  abiotic_pe_model$coefficients[15]*fl_test_abiotic_pe_data$forest_ring_2 + 
  abiotic_pe_model$coefficients[16]*fl_test_abiotic_pe_data$grass_ring_2 +
  abiotic_pe_model$coefficients[17]*fl_test_abiotic_pe_data$wetland_ring_2 +
  abiotic_pe_model$coefficients[18]*fl_test_abiotic_pe_data$farming_ring_2 +
  abiotic_pe_model$coefficients[19]*fl_test_abiotic_pe_data$urban_ring_2 +
  abiotic_pe_model$coefficients[20]*fl_test_abiotic_pe_data$barren_ring_2

fl_test_abiotic_data <- fl_test[,c(-1,-12,-22,-23,-24)]

fl_test_abiotic_data$predictions <- abiotic_model$coefficients[1] +
  abiotic_model$coefficients[2]*fl_test_abiotic_data$tmax_ring_1 +
  abiotic_model$coefficients[3]*fl_test_abiotic_data$prcp_ring_1 + 
  abiotic_model$coefficients[4]*fl_test_abiotic_data$water_ring_1 +
  abiotic_model$coefficients[5]*fl_test_abiotic_data$forest_ring_1 + 
  abiotic_model$coefficients[6]*fl_test_abiotic_data$grass_ring_1 +
  abiotic_model$coefficients[7]*fl_test_abiotic_data$wetland_ring_1 +
  abiotic_model$coefficients[8]*fl_test_abiotic_data$farming_ring_1 +
  abiotic_model$coefficients[9]*fl_test_abiotic_data$urban_ring_1 +
  abiotic_model$coefficients[10]*fl_test_abiotic_data$barren_ring_1 +
  abiotic_model$coefficients[11]*fl_test_abiotic_data$tmax_ring_2 +
  abiotic_model$coefficients[12]*fl_test_abiotic_data$prcp_ring_2 + 
  abiotic_model$coefficients[13]*fl_test_abiotic_data$water_ring_2 +
  abiotic_model$coefficients[14]*fl_test_abiotic_data$forest_ring_2 + 
  abiotic_model$coefficients[15]*fl_test_abiotic_data$grass_ring_2 +
  abiotic_model$coefficients[16]*fl_test_abiotic_data$wetland_ring_2 +
  abiotic_model$coefficients[17]*fl_test_abiotic_data$farming_ring_2 +
  abiotic_model$coefficients[18]*fl_test_abiotic_data$urban_ring_2 +
  abiotic_model$coefficients[19]*fl_test_abiotic_data$barren_ring_2

fl_test_biotic_pe_data <- fl_test[,c(2, 12, 22, 23)]

fl_test_biotic_pe_data$predictions <- biotic_pe_model$coefficients[1] +
  biotic_pe_model$coefficients[2]*fl_test_biotic_pe_data$previous_establishment + 
  biotic_pe_model$coefficients[3]*fl_test_biotic_pe_data$bs_ring_1 +
  biotic_pe_model$coefficients[4]*fl_test_biotic_pe_data$bs_ring_2

fl_test_biotic_data <- fl_test[,c(2,22,23)]

fl_test_biotic_data$predictions <- biotic_model$coefficients[1] +
  biotic_model$coefficients[2]*fl_test_biotic_data$bs_ring_1 +
  biotic_model$coefficients[3]*fl_test_biotic_data$bs_ring_2

fl_test_all_params_data$establishment_score <- fl_test$establishment_score

## calculate MSE and variance explained for linear models

# mse
mse_florida_all_params <- mean((fl_test_all_params_data$predictions-fl_test_all_params_data$establishment_score)^2)
mse_florida_abiotic_pe <- mean((fl_test_abiotic_pe_data$predictions-fl_test_abiotic_pe_data$establishment_score)^2)
mse_florida_abiotic <- mean((fl_test_abiotic_data$predictions-fl_test_abiotic_data$establishment_score)^2)
mse_florida_biotic_pe <- mean((fl_test_biotic_pe_data$predictions-fl_test_biotic_pe_data$establishment_score)^2)
mse_florida_biotic <- mean((fl_test_biotic_data$predictions-fl_test_biotic_data$establishment_score)^2)

# variance explained

mean_test_response <- mean(fl_test$establishment_score)

TSS <- sum((fl_test$establishment_score - mean_test_response)^2)

RSS_all_params <- sum((fl_test_all_params_data$establishment_score - fl_test_all_params_data$predictions)^2)
RSS_abiotic_pe <- sum((fl_test_abiotic_pe_data$establishment_score - fl_test_abiotic_pe_data$predictions)^2)
RSS_abiotic <- sum((fl_test_abiotic_data$establishment_score - fl_test_abiotic_data$predictions)^2)
RSS_biotic_pe <- sum((fl_test_biotic_pe_data$establishment_score - fl_test_biotic_pe_data$predictions)^2)
RSS_biotic <- sum((fl_test_biotic_data$establishment_score - fl_test_biotic_data$predictions)^2)


var_explained_all_params <- (1 - RSS_all_params / TSS)
var_explained_abiotic_pe <- (1 - RSS_abiotic_pe / TSS)
var_explained_abiotic <- (1 - RSS_abiotic / TSS)
var_explained_biotic_pe <- (1 - RSS_biotic_pe / TSS)
var_explained_biotic <- (1 - RSS_biotic / TSS)

# save results

model <- c("all_parameters", "abiotic_pe", "abiotic", "biotic_pe", "biotic")
var_explained <- c(var_explained_all_params, var_explained_abiotic_pe, 
                   var_explained_abiotic, var_explained_biotic_pe, 
                   var_explained_biotic)
mse <- c(mse_florida_all_params, mse_florida_abiotic_pe, mse_florida_abiotic,
         mse_florida_biotic_pe, mse_florida_biotic)

lm_performance_florida <- data.frame(model, mse, var_explained)

write.csv(lm_performance_florida, "data/lm_performance_florida.csv")

#### REPEAT FOR TEXAS ----

tx <- read.csv("data/rf_model_v5_all_parameters_data_texas.csv")

tx_all_params_data <- tx[,c(-1,-24)]

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

tx_abiotic_pe_data <- tx[,c(-1,-22,-23,-24)]

tx_abiotic_pe_data$predictions <- abiotic_pe_model$coefficients[1] +
  abiotic_pe_model$coefficients[2]*tx_abiotic_pe_data$tmax_ring_1 +
  abiotic_pe_model$coefficients[3]*tx_abiotic_pe_data$prcp_ring_1 + 
  abiotic_pe_model$coefficients[4]*tx_abiotic_pe_data$water_ring_1 +
  abiotic_pe_model$coefficients[5]*tx_abiotic_pe_data$forest_ring_1 + 
  abiotic_pe_model$coefficients[6]*tx_abiotic_pe_data$grass_ring_1 +
  abiotic_pe_model$coefficients[7]*tx_abiotic_pe_data$wetland_ring_1 +
  abiotic_pe_model$coefficients[8]*tx_abiotic_pe_data$farming_ring_1 +
  abiotic_pe_model$coefficients[9]*tx_abiotic_pe_data$urban_ring_1 +
  abiotic_pe_model$coefficients[10]*tx_abiotic_pe_data$barren_ring_1 +
  abiotic_pe_model$coefficients[11]*tx_abiotic_pe_data$previous_establishment +
  abiotic_pe_model$coefficients[12]*tx_abiotic_pe_data$tmax_ring_2 +
  abiotic_pe_model$coefficients[13]*tx_abiotic_pe_data$prcp_ring_2 + 
  abiotic_pe_model$coefficients[14]*tx_abiotic_pe_data$water_ring_2 +
  abiotic_pe_model$coefficients[15]*tx_abiotic_pe_data$forest_ring_2 + 
  abiotic_pe_model$coefficients[16]*tx_abiotic_pe_data$grass_ring_2 +
  abiotic_pe_model$coefficients[17]*tx_abiotic_pe_data$wetland_ring_2 +
  abiotic_pe_model$coefficients[18]*tx_abiotic_pe_data$farming_ring_2 +
  abiotic_pe_model$coefficients[19]*tx_abiotic_pe_data$urban_ring_2 +
  abiotic_pe_model$coefficients[20]*tx_abiotic_pe_data$barren_ring_2

tx_abiotic_data <- tx[,c(-1,-12,-22,-23,-24)]

tx_abiotic_data$predictions <- abiotic_model$coefficients[1] +
  abiotic_model$coefficients[2]*tx_abiotic_data$tmax_ring_1 +
  abiotic_model$coefficients[3]*tx_abiotic_data$prcp_ring_1 + 
  abiotic_model$coefficients[4]*tx_abiotic_data$water_ring_1 +
  abiotic_model$coefficients[5]*tx_abiotic_data$forest_ring_1 + 
  abiotic_model$coefficients[6]*tx_abiotic_data$grass_ring_1 +
  abiotic_model$coefficients[7]*tx_abiotic_data$wetland_ring_1 +
  abiotic_model$coefficients[8]*tx_abiotic_data$farming_ring_1 +
  abiotic_model$coefficients[9]*tx_abiotic_data$urban_ring_1 +
  abiotic_model$coefficients[10]*tx_abiotic_data$barren_ring_1 +
  abiotic_model$coefficients[11]*tx_abiotic_data$tmax_ring_2 +
  abiotic_model$coefficients[12]*tx_abiotic_data$prcp_ring_2 + 
  abiotic_model$coefficients[13]*tx_abiotic_data$water_ring_2 +
  abiotic_model$coefficients[14]*tx_abiotic_data$forest_ring_2 + 
  abiotic_model$coefficients[15]*tx_abiotic_data$grass_ring_2 +
  abiotic_model$coefficients[16]*tx_abiotic_data$wetland_ring_2 +
  abiotic_model$coefficients[17]*tx_abiotic_data$farming_ring_2 +
  abiotic_model$coefficients[18]*tx_abiotic_data$urban_ring_2 +
  abiotic_model$coefficients[19]*tx_abiotic_data$barren_ring_2

tx_biotic_pe_data <- tx[,c(2, 12, 22, 23)]

tx_biotic_pe_data$predictions <- biotic_pe_model$coefficients[1] +
  biotic_pe_model$coefficients[2]*tx_biotic_pe_data$previous_establishment + 
  biotic_pe_model$coefficients[3]*tx_biotic_pe_data$bs_ring_1 +
  biotic_pe_model$coefficients[4]*tx_biotic_pe_data$bs_ring_2

tx_biotic_data <- tx[,c(2,22,23)]

tx_biotic_data$predictions <- biotic_model$coefficients[1] +
  biotic_model$coefficients[2]*tx_biotic_data$bs_ring_1 +
  biotic_model$coefficients[3]*tx_biotic_data$bs_ring_2

## calculate MSE and variance explained for linear models

# mse
mse_texas_all_params <- mean((tx_all_params_data$predictions-tx_all_params_data$establishment_score)^2)
mse_texas_abiotic_pe <- mean((tx_abiotic_pe_data$predictions-tx_abiotic_pe_data$establishment_score)^2)
mse_texas_abiotic <- mean((tx_abiotic_data$predictions-tx_abiotic_data$establishment_score)^2)
mse_texas_biotic_pe <- mean((tx_biotic_pe_data$predictions-tx_biotic_pe_data$establishment_score)^2)
mse_texas_biotic <- mean((tx_biotic_data$predictions-tx_biotic_data$establishment_score)^2)

# variance explained

mean_test_response <- mean(tx$establishment_score)

TSS <- sum((tx$establishment_score - mean_test_response)^2)

RSS_all_params <- sum((tx_all_params_data$establishment_score - tx_all_params_data$predictions)^2)
RSS_abiotic_pe <- sum((tx_abiotic_pe_data$establishment_score - tx_abiotic_pe_data$predictions)^2)
RSS_abiotic <- sum((tx_abiotic_data$establishment_score - tx_abiotic_data$predictions)^2)
RSS_biotic_pe <- sum((tx_biotic_pe_data$establishment_score - tx_biotic_pe_data$predictions)^2)
RSS_biotic <- sum((tx_biotic_data$establishment_score - tx_biotic_data$predictions)^2)


var_explained_tx_all_params <- (1 - RSS_all_params / TSS)
var_explained_tx_abiotic_pe <- (1 - RSS_abiotic_pe / TSS)
var_explained_tx_abiotic <- (1 - RSS_abiotic / TSS)
var_explained_tx_biotic_pe <- (1 - RSS_biotic_pe / TSS)
var_explained_tx_biotic <- (1 - RSS_biotic / TSS)

# save results

tx_model <- c("all_parameters", "abiotic_pe", "abiotic", "biotic_pe", "biotic")
tx_var_explained <- c(var_explained_tx_all_params, var_explained_tx_abiotic_pe, 
                   var_explained_tx_abiotic, var_explained_tx_biotic_pe, 
                   var_explained_tx_biotic)
tx_mse <- c(mse_texas_all_params, mse_texas_abiotic_pe, mse_texas_abiotic,
         mse_texas_biotic_pe, mse_texas_biotic)

lm_performance_texas <- data.frame(tx_model, tx_mse, tx_var_explained)

write.csv(lm_performance_texas, "data/lm_performance_texas.csv")

