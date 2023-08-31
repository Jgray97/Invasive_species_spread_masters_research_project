### TESTING + ANALYSIS OF RF MODEL ON TEXAS DATASET

# Author = John Gray
# Email = greyjohn15@gmail.com
# Last edit = 31/08/23

# Load packages ----

library(randomForest)
library(ggplot2)
library(dplyr)
library(pROC)

# Load in RF MODEL INPUT AND TARGET VARIABLE DATASETS FOR TEXAS ----

# Load in abiotic input data for texas (this dataset also contains target
# target variable establishment score and previous establishment input values)
abiotic_data <- read.csv("data/rf_model_abiotic_data_texas.csv")

# Load in biotic variable input data
biotic_data <- read.csv("data/texas_rf_biotic_data.csv")

# bind abiotic and biotic input datasets together
new_model_data <- cbind(abiotic_data, biotic_data)

# get rid of extraneous and repeated columns
new_model_data <- new_model_data[,-c(1,25,26,27,28,29)]

# filter dataset to only consider predictions for cells where there is enough
# data as these are the only predictions which can be validated
refined_data <- filter(new_model_data, enough_checklists == 1)

# filter dataframe down to only cases where it is possible that the invasive
# species could become established due to spread (i.e. within 2 cells of an 
# established cell) in order to remove class imbalances
refined_data$interesting_cases <- refined_data$tmax_ring_1 + 
  refined_data$previous_establishment + refined_data$tmax_ring_2 + 
  refined_data$establishment_score

refined_data <- filter(refined_data, interesting_cases != 0)

# get rid of 2010 from data because this is when the goose established itself 
# first and we are not trying to predict cases of initial establishment
refined_data <- filter(refined_data, year != 2010)

#### Now loaded in rf models that have been pre-trained on Florida training ----
#### data ----

load(file = "all_parameters_rf_v5.RData")

load(file = "ring1_pe_rf_v5.RData")

load(file = "ring1_rf_v5.RData")

load(file = "ring2_pe_rf_v5.RData")

load(file = "ring2_rf_v5.RData")

load(file = "biotic_previous_establishment_rf_v5.RData")

load(file = "biotic_rf_v5.RData")

load(file = "abiotic_previous_establishment_rf_v5.RData")

load(file = "abiotic_rf_v5.RData")

# Apply all parameters version of the model to texas data and plot ----
# predicted establishment score against true establishment score ----

rf_data <- refined_data[,-c(1,2,4,26)]

rf_data$predicted_establishment <- predict(rf.fit, rf_data)

texas_pred <- ggplot(rf_data, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  theme_classic() +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12)) +
  scale_y_continuous(breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12)) +
  labs(title = "Texas Dataset",
       x = "Predicted establishment score (SP model output)",
       y = "True establishment score (ET model output)") +
  coord_cartesian(xlim = c(0, 0.12), ylim = c(0, 0.12))

# save plot
ggsave("plots_texas_rf_predictions_v2/all_parameters_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate_mse
mse_test_all_parameters <- mean((rf_data$predicted_establishment-rf_data$establishment_score)^2)

# calculate percentage Variance (equivalent of r2)
mean_test_all_parameters_response <- mean(rf_data$establishment_score)

TSS <- sum((rf_data$establishment_score - mean_test_all_parameters_response)^2)

RSS_all_parameters <- sum((rf_data$establishment_score - rf_data$predicted_establishment)^2)

var_explained_all_parameters <- (1 - RSS_all_parameters / TSS)

# Apply ring 1 + previous establishment version of the model to texas data ----
# and plot predicted establishment score against true establishment score ----

rf_data_ring_1_pe <- refined_data[,-c(1,2,4,15,16,17,18,19,20,21,22,23,25,26)]

rf_data_ring_1_pe$predicted_establishment <- predict(rf.fit_ring1_pe, rf_data)

ggplot(rf_data_ring_1_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "RF model V5 - 1st ring and previous establishment", 
       subtitle = "Predicted establishment vs actual establishment for Egyptian goose in texas")

# save plot
ggsave("plots_texas_rf_predictions_v2/ring_1_pe_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate_mse
mse_test_ring_1_pe <- mean((rf_data_ring_1_pe$predicted_establishment-rf_data_ring_1_pe$establishment_score)^2)

# calculate percentage Variance (equivalent of r2)
mean_test_response <- mean(rf_data$establishment_score)

TSS <- sum((rf_data$establishment_score - mean_test_response)^2)

RSS_ring_1_pe <- sum((rf_data_ring_1_pe$establishment_score - rf_data_ring_1_pe$predicted_establishment)^2)

var_explained_ring_1_pe <- (1 - RSS_ring_1_pe / TSS)

# Apply ring 1 version of the model to texas data ----
# and plot predicted establishment score against true establishment score ----

rf_data_ring1 <- refined_data[,-c(1,2,4,14,15,16,17,18,19,20,21,22,23,25,26)]

rf_data_ring1$predicted_establishment <- predict(rf.fit_ring1, rf_data_ring1)

ggplot(rf_data_ring1, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "RF model V5 - 1st ring", 
       subtitle = "Predicted establishment vs actual establishment for Egyptian goose in texas")

# save plot
ggsave("plots_texas_rf_predictions_v2/ring_1_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate_mse
mse_test_ring_1 <- mean((rf_data_ring1$predicted_establishment-rf_data_ring1$establishment_score)^2)

# calculate percentage Variance (equivalent of r2)

mean_ring_1_test_response <- mean(rf_data_ring1$establishment_score)

TSS_ring_1 <- sum((rf_data_ring1$establishment_score - mean_ring_1_test_response)^2)

RSS_ring_1 <- sum((rf_data_ring1$establishment_score - rf_data_ring1$predicted_establishment)^2)

var_explained_ring_1 <- (1 - RSS_ring_1 / TSS_ring_1)

# Apply ring 2 + previous establishment version of the model to texas data ----
# and plot predicted establishment score against true establishment score ----

rf_data_ring_2_pe <- refined_data[,c(3,14,15,16,17,18,19,20,21,22,23,25)]

rf_data_ring_2_pe$predicted_establishment <- predict(rf.fit_ring2_pe, rf_data)

ggplot(rf_data_ring_2_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "RF model V5 - 2nd ring and previous establishment", 
       subtitle = "Predicted establishment vs actual establishment for Egyptian goose in texas")

# save plot
ggsave("plots_texas_rf_predictions_v2/ring_2_pe_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate_mse

mse_test_ring_2_pe <- mean((rf_data_ring_2_pe$predicted_establishment-rf_data_ring_2_pe$establishment_score)^2)

# calculate percentage Variance (equivalent of r2)

mean_test_response <- mean(rf_data$establishment_score)

TSS <- sum((rf_data$establishment_score - mean_test_response)^2)

RSS_ring_2_pe <- sum((rf_data_ring_2_pe$establishment_score - rf_data_ring_2_pe$predicted_establishment)^2)

var_explained_ring_2_pe <- (1 - RSS_ring_2_pe / TSS)

# Apply ring 2 version of the model to texas data ----
# and plot predicted establishment score against true establishment score ----

rf_data_ring2 <- refined_data[,c(3,15,16,17,18,19,20,21,22,23,25)]

rf_data_ring2$predicted_establishment <- predict(rf.fit_ring2, rf_data_ring2)

ggplot(rf_data_ring2, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "RF model V5 - 2nd ring", 
       subtitle = "Predicted establishment vs actual establishment for Egyptian goose in texas")

# save plot
ggsave("plots_texas_rf_predictions_v2/ring_2_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate_mse

mse_test_ring_2 <- mean((rf_data_ring2$predicted_establishment-rf_data_ring2$establishment_score)^2)

# calculate percentage Variance (equivalent of r2)

mean_ring_2_test_response <- mean(rf_data_ring2$establishment_score)

TSS_ring_2 <- sum((rf_data_ring2$establishment_score - mean_ring_2_test_response)^2)

RSS_ring_2 <- sum((rf_data_ring2$establishment_score - rf_data_ring2$predicted_establishment)^2)

var_explained_ring_2 <- (1 - RSS_ring_2 / TSS_ring_2)

# Apply biotic + previous establishment version of the model to texas data ----
# and plot predicted establishment score against true establishment score ----

rf_data_biotic_pe <- refined_data[,c(3,14,24,25)]

rf_data_biotic_pe$predicted_establishment <- predict(rf.fit_biotic_pe, rf_data_biotic_pe)

ggplot(rf_data_biotic_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "RF model V5 - Biotic parameters + previous establishment", 
       subtitle = "Predicted establishment vs actual establishment for Egyptian goose in texas")

# save plot
ggsave("plots_texas_rf_predictions_v2/biotic_pe_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate_mse
mse_test_biotic_pe <- mean((rf_data_biotic_pe$predicted_establishment-rf_data_biotic_pe$establishment_score)^2)

# calculate percentage Variance (equivalent of r2)

mean_test_response <- mean(rf_data$establishment_score)

TSS <- sum((rf_data$establishment_score - mean_test_response)^2)

RSS_biotic_pe <- sum((rf_data_biotic_pe$establishment_score - rf_data_biotic_pe$predicted_establishment)^2)

var_explained_biotic_pe <- (1 - RSS_biotic_pe / TSS)


# Apply biotic version of the model to texas data ----
# and plot predicted establishment score against true establishment score ----

rf_data_biotic <- refined_data[,c(3,24,25)]

rf_data_biotic$predicted_establishment <- predict(rf.fit_biotic, rf_data_biotic)

ggplot(rf_data_biotic, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "RF model V5 - Biotic parameters", 
       subtitle = "Predicted establishment vs actual establishment for Egyptian goose in texas")

# save plot
ggsave("plots_texas_rf_predictions_v2/biotic_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate_mse

mse_test_biotic <- mean((rf_data_biotic$predicted_establishment-rf_data_biotic$establishment_score)^2)

# calculate percentage Variance (equivalent of r2)

mean_test_response <- mean(rf_data$establishment_score)

TSS <- sum((rf_data$establishment_score - mean_test_response)^2)

RSS_biotic <- sum((rf_data_biotic$establishment_score - rf_data_biotic$predicted_establishment)^2)

var_explained_biotic <- (1 - RSS_biotic / TSS)

# Apply abiotic + previous establishment version of the model to texas data ----
# and plot predicted establishment score against true establishment score ----

rf_data_abiotic_pe <- refined_data[,c(3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)]

rf_data_abiotic_pe$predicted_establishment <- predict(rf.fit_abiotic_pe, rf_data_abiotic_pe)

ggplot(rf_data_abiotic_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "RF model V5 - Abiotic parameters + previous establishment", 
       subtitle = "Predicted establishment vs actual establishment for Egyptian goose in texas")

# save plot
ggsave("plots_texas_rf_predictions_v2/abiotic_pe_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate_mse

mse_test_abiotic_pe <- mean((rf_data_abiotic_pe$predicted_establishment-rf_data_abiotic_pe$establishment_score)^2)

# calculate percentage Variance (equivalent of r2) for test data

mean_test_response <- mean(rf_data$establishment_score)

TSS <- sum((rf_data$establishment_score - mean_test_response)^2)

RSS_abiotic_pe <- sum((rf_data_abiotic_pe$establishment_score - rf_data_abiotic_pe$predicted_establishment)^2)

var_explained_abiotic_pe <- (1 - RSS_abiotic_pe / TSS)


# Apply abiotic version of the model to texas data ----
# and plot predicted establishment score against true establishment score ----

rf_data_abiotic <- refined_data[,c(3,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23)]

rf_data_abiotic$predicted_establishment <- predict(rf.fit_abiotic, rf_data_abiotic)

ggplot(rf_data_abiotic, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "RF model V5 - Abiotic parameters", 
       subtitle = "Predicted establishment vs actual establishment for Egyptian goose in texas")

# save plot
ggsave("plots_texas_rf_predictions_v2/abiotic_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate_mse

mse_test_abiotic <- mean((rf_data_abiotic$predicted_establishment-rf_data_abiotic$establishment_score)^2)

# calculate percentage Variance (equivalent of r2)

mean_test_response <- mean(rf_data$establishment_score)

TSS <- sum((rf_data$establishment_score - mean_test_response)^2)

RSS_abiotic <- sum((rf_data_abiotic$establishment_score - rf_data_abiotic$predicted_establishment)^2)

var_explained_abiotic <- (1 - RSS_abiotic / TSS)


#### Straightforward performance comparison ----

texas_rf_stats <- data.frame("Model" = c("all_parameters", 
                                   "ring_1_previous_establishment", "ring_1",
                                   "ring_2_previous_establishment", "ring_2",
                                   "biotic_previous_establishment", "biotic",
                                   "abiotic_previous_establishment", "abiotic"),
                       "MSE" = c(mse_test_all_parameters,
                                                  mse_test_ring_1_pe,
                                                  mse_test_ring_1,
                                 mse_test_ring_2_pe, mse_test_ring_2,
                                                  mse_test_biotic_pe,
                                                  mse_test_biotic,
                                                  mse_test_abiotic_pe,
                                                  mse_test_abiotic),
                       "Percentage_var_explained" = c(var_explained_all_parameters,
                                                           var_explained_ring_1_pe,
                                                           var_explained_ring_1,
                                                      var_explained_ring_2_pe,
                                                      var_explained_ring_2,
                                                           var_explained_biotic_pe,
                                                           var_explained_biotic,
                                                           var_explained_abiotic_pe,
                                                           var_explained_abiotic))

# save straightforward performance comparison table
write.csv(texas_rf_stats, "data/rf_texas_v2_performance_comparison.csv")

#### Analysis of capabilities of spread prediction ----

# We want to know: 
# percentage/number of spread events correctly predicted
# percentage/number of spread events incorrectly predicted
# accuracy when predicting spread events 

# creating table to analyse spread prediction capabilities

model = c("all_parameters", 
          "ring_1_previous_establishment", "ring_1",
          "ring_2_previous_establishment", "ring_2",
          "biotic_previous_establishment", "biotic",
          "abiotic_previous_establishment", "abiotic")

correct_spread_predictions = numeric(9)

correct_spread_prediction_percentage <- numeric(9)

inaccurate_spread_predictions <- numeric(9)

inaccurate_spread_prediction_percentage <- numeric(9)

detailed_analysis <- data.frame(model, correct_spread_predictions, 
                                correct_spread_prediction_percentage,
                                inaccurate_spread_predictions,
                                inaccurate_spread_prediction_percentage)

# create a dataframe which contains actual establishment scores and predicted
# establishment for every version of the rf model for texas
test_data_predictions <- cbind(rf_data[,c(1,11,23)],
                               rf_data_ring_1_pe$predicted_establishment,
                               rf_data_ring1$predicted_establishment,
                               rf_data_ring_2_pe$predicted_establishment,
                               rf_data_ring2$predicted_establishment,
                               rf_data_biotic_pe$predicted_establishment,
                               rf_data_biotic$predicted_establishment,
                               rf_data_abiotic_pe$predicted_establishment,
                               rf_data_abiotic$predicted_establishment)

# calculate number of correct spread predictions for each version of the model
detailed_analysis$correct_spread_predictions[1] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[1] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))


detailed_analysis$correct_spread_predictions[2] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (rf_data_ring_1_pe$predicted_establishment > 0.005)))


detailed_analysis$correct_spread_prediction_percentage[2] <- 100 * nrow(filter(test_data_predictions, 
                                                                               (establishment_score > 0.005) &
                                                                                 (previous_establishment < 0.005) &
                                                                                 (rf_data_ring_1_pe$predicted_establishment > 0.005))) /
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[3] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (rf_data_ring1$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[3] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (rf_data_ring1$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[4] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (rf_data_ring_2_pe$predicted_establishment > 0.005)))


detailed_analysis$correct_spread_prediction_percentage[4] <- 100 * nrow(filter(test_data_predictions, 
                                                                               (establishment_score > 0.005) &
                                                                                 (previous_establishment < 0.005) &
                                                                                 (rf_data_ring_2_pe$predicted_establishment > 0.005))) /
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[5] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (rf_data_ring2$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[5] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (rf_data_ring2$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[6] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (rf_data_biotic_pe$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[6] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (rf_data_biotic_pe$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[7] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (rf_data_biotic$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[7] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (rf_data_biotic$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[8] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (rf_data_abiotic_pe$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[8] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (rf_data_abiotic_pe$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[9] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (rf_data_abiotic$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[9] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (rf_data_abiotic$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

# rename detailed analysis columns
colnames(detailed_analysis)[3] <- "correct_spread_prediction_percentage"

colnames(detailed_analysis)[4] <- "correct_no_spread_predictions"

colnames(detailed_analysis)[5] <- "correct_no_spread_prediction_percentage"

# calculate number of correct predictions of no spread
detailed_analysis$correct_no_spread_predictions[1] <- nrow(filter(test_data_predictions, 
                                                                  (establishment_score < 0.005) &
                                                                    (previous_establishment < 0.005) &
                                                                    (predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[1] <- 100*nrow(filter(test_data_predictions, 
                                                                                (establishment_score < 0.005) &
                                                                                  (previous_establishment < 0.005) &
                                                                                  (predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))


detailed_analysis$correct_no_spread_predictions[2] <- nrow(filter(test_data_predictions, 
                                                                  (establishment_score < 0.005) &
                                                                    (previous_establishment < 0.005) &
                                                                    (rf_data_ring_1_pe$predicted_establishment < 0.005)))


detailed_analysis$correct_no_spread_prediction_percentage[2] <- 100 * nrow(filter(test_data_predictions, 
                                                                                  (establishment_score < 0.005) &
                                                                                    (previous_establishment < 0.005) &
                                                                                    (rf_data_ring_1_pe$predicted_establishment < 0.005))) /
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[3] <- nrow(filter(test_data_predictions, 
                                                                  (establishment_score < 0.005) &
                                                                    (previous_establishment < 0.005) &
                                                                    (rf_data_ring1$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[3] <- 100*nrow(filter(test_data_predictions, 
                                                                                (establishment_score < 0.005) &
                                                                                  (previous_establishment < 0.005) &
                                                                                  (rf_data_ring1$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[4] <- nrow(filter(test_data_predictions, 
                                                                  (establishment_score < 0.005) &
                                                                    (previous_establishment < 0.005) &
                                                                    (rf_data_ring_2_pe$predicted_establishment < 0.005)))


detailed_analysis$correct_no_spread_prediction_percentage[4] <- 100 * nrow(filter(test_data_predictions, 
                                                                                  (establishment_score < 0.005) &
                                                                                    (previous_establishment < 0.005) &
                                                                                    (rf_data_ring_2_pe$predicted_establishment < 0.005))) /
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[5] <- nrow(filter(test_data_predictions, 
                                                                  (establishment_score < 0.005) &
                                                                    (previous_establishment < 0.005) &
                                                                    (rf_data_ring2$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[5] <- 100*nrow(filter(test_data_predictions, 
                                                                                (establishment_score < 0.005) &
                                                                                  (previous_establishment < 0.005) &
                                                                                  (rf_data_ring2$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[6] <- nrow(filter(test_data_predictions, 
                                                                  (establishment_score < 0.005) &
                                                                    (previous_establishment < 0.005) &
                                                                    (rf_data_biotic_pe$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[6] <- 100*nrow(filter(test_data_predictions, 
                                                                                (establishment_score < 0.005) &
                                                                                  (previous_establishment < 0.005) &
                                                                                  (rf_data_biotic_pe$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[7] <- nrow(filter(test_data_predictions, 
                                                                  (establishment_score < 0.005) &
                                                                    (previous_establishment < 0.005) &
                                                                    (rf_data_biotic$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[7] <- 100*nrow(filter(test_data_predictions, 
                                                                                (establishment_score < 0.005) &
                                                                                  (previous_establishment < 0.005) &
                                                                                  (rf_data_biotic$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[8] <- nrow(filter(test_data_predictions, 
                                                                  (establishment_score < 0.005) &
                                                                    (previous_establishment < 0.005) &
                                                                    (rf_data_abiotic_pe$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[8] <- 100*nrow(filter(test_data_predictions, 
                                                                                (establishment_score < 0.005) &
                                                                                  (previous_establishment < 0.005) &
                                                                                  (rf_data_abiotic_pe$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[9] <- nrow(filter(test_data_predictions, 
                                                                  (establishment_score < 0.005) &
                                                                    (previous_establishment < 0.005) &
                                                                    (rf_data_abiotic$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[9] <- 100*nrow(filter(test_data_predictions, 
                                                                                (establishment_score < 0.005) &
                                                                                  (previous_establishment < 0.005) &
                                                                                  (rf_data_abiotic$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

# turn these values into numeric values to enable mathematical manipulation with
# them
detailed_analysis$correct_spread_predictions <- as.numeric(detailed_analysis$correct_spread_predictions)
detailed_analysis$correct_spread_prediction_percentage <- as.numeric(detailed_analysis$correct_spread_prediction_percentage)
detailed_analysis$correct_no_spread_predictions <- as.numeric(detailed_analysis$correct_no_spread_predictions)
detailed_analysis$correct_no_spread_prediction_percentage <- as.numeric(detailed_analysis$correct_no_spread_prediction_percentage)

# calculate total number of spread cases
detailed_analysis$no_spread <- 100*detailed_analysis$correct_no_spread_predictions/detailed_analysis$correct_no_spread_prediction_percentage

# calculate spread prediction precision
detailed_analysis$precision <- detailed_analysis$correct_spread_predictions / (detailed_analysis$correct_spread_predictions+detailed_analysis$no_spread-detailed_analysis$correct_no_spread_predictions)

# rename correct prediction percentage to spread recall
colnames(detailed_analysis)[3] <- "spread_recall"

# insert column for total number of spread events
detailed_analysis$spread <- 11 

# calculate overall spread prediction accuracy for all models
detailed_analysis$accuracy <- 100*(detailed_analysis$correct_spread_predictions + detailed_analysis$correct_no_spread_predictions) / (detailed_analysis$no_spread + detailed_analysis$spread)

# save detailed analysis table as a csv
write.csv(detailed_analysis, "data/rf_model_v5_texas_predictive_analysis.csv")

# save texas data predictions table as a csv
write.csv(test_data_predictions, "data/rf_model_v5_texas_predictions.csv")

### calculating AUC SCORE ----

# create auc dataframe for texas predictions to be populated with scores
auc <- data.frame(model, "auc_all_data" = numeric(9), "auc_spread_data" = numeric(9))

## To calculate AUC score we need to assign binary establishment values for each
## establishment score (we will take any non 0 actual establishment score to be
## indicative of establishment)

# First create new column for binary establishment value
test_data_predictions$established <- 0

# then assign binary establishment values to each entry accordingly
for (i in 1:nrow(test_data_predictions)) {
  if (test_data_predictions$establishment_score[i] < 0.005) {
    test_data_predictions$established[i] <- 0
  } else {
    test_data_predictions$established[i] <- 1
  }
}

# we calculate AUC scores both for all Texas data, and for only the cells
# where there is a possibility of spread, which provides more insight into the
# model's spread prediction capabilities

# filter test data to just cells where there is a possibility of spread
test_spread_data_predictions <- filter(test_data_predictions, previous_establishment < 0.005)

# calculate auc scores for each version of the model for both overall and just
# possible spread datasets
auc$auc_all_data[1] <- auc(test_data_predictions$established, test_data_predictions$predicted_establishment)
auc$auc_spread_data[1] <- auc(test_spread_data_predictions$established, test_spread_data_predictions$predicted_establishment)

auc$auc_all_data[2] <- auc(test_data_predictions$established, test_data_predictions[,4])
auc$auc_spread_data[2] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,4])

auc$auc_all_data[3] <- auc(test_data_predictions$established, test_data_predictions[,5])
auc$auc_spread_data[3] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,5])

auc$auc_all_data[4] <- auc(test_data_predictions$established, test_data_predictions[,6])
auc$auc_spread_data[4] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,6])

auc$auc_all_data[5] <- auc(test_data_predictions$established, test_data_predictions[,7])
auc$auc_spread_data[5] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,7])

auc$auc_all_data[6] <- auc(test_data_predictions$established, test_data_predictions[,8])
auc$auc_spread_data[6] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,8])

auc$auc_all_data[7] <- auc(test_data_predictions$established, test_data_predictions[,9])
auc$auc_spread_data[7] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,9])

auc$auc_alldata[8] <- auc(test_data_predictions$established, test_data_predictions[,10])
auc$auc_spread_data[8] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,10])

auc$auc_all_data[9] <- auc(test_data_predictions$established, test_data_predictions[,11])
auc$auc_spread_data[9] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,11])

# save auc scores in a csv file
write.csv(auc, "data/texas_data_auc_scores.csv")

# save overall rf input + target variable data in a csv file
write.csv(rf_data, "data/rf_model_v5_all_parameters_data.csv")

# create the roc curves for best performing models on spread and plot them ----

# comment out stuff below if not picking up from here
predictions <- read.csv("data/rf_model_v5_texas_predictions.csv")
predictions <- predictions[,-1]
test_spread_data_predictions <- filter(predictions, previous_establishment < 0.005)


# create roc curve plots for the 2 versions of the model with the highest AUC
# scores for the texas spread possibility data
roc.biotic <- roc(test_spread_data_predictions$established, test_spread_data_predictions$rf_data_biotic.predicted_establishment)
biotic_roc_test <- data.frame(roc.biotic$specificities, roc.biotic$sensitivities)
biotic_roc <- ggplot(biotic_roc_test, aes(x =roc.biotic.specificities, y = roc.biotic.sensitivities)) +
  geom_path() +
  scale_x_reverse() +
  theme_classic() +
  geom_abline(aes(intercept = 1, slope = 1, color = "red")) +
  theme(legend.position = "none") + 
  labs(title = "Biotic variables model",
       x = "Specificity (1 - false positive rate)",
       y = "Sensitivity (true positive rate)")


roc.abiotic <- roc(test_spread_data_predictions$established, test_spread_data_predictions$rf_data_abiotic.predicted_establishment)
abiotic_roc_test <- data.frame(roc.abiotic$specificities, roc.abiotic$sensitivities)
abiotic_roc <- ggplot(abiotic_roc_test, aes(x =roc.abiotic.specificities, y = roc.abiotic.sensitivities)) +
  geom_path() +
  scale_x_reverse() +
  theme_classic() +
  geom_abline(aes(intercept = 1, slope = 1, color = "red")) +
  theme(legend.position = "none") + 
  labs(title = "Abiotic variables model",
       x = "Specificity (1 - false positive rate)",
       y = "Sensitivity (true positive rate)")


# save these plots in a grid image
grid.arrange(biotic_roc, abiotic_roc, nrow = 1)
g <- arrangeGrob(biotic_roc, abiotic_roc, nrow = 1)
ggsave(file = "final_establishment_plots/combined_plots/roc_curves.png", g, 
       width = 10, height = 5)

