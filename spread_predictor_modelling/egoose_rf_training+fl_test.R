### TRAINING OF RF MODEL ON FLORIDA TRAINING DATASET + TESTING OF RF MODEL ON
### FLORIDA TEST DATASET

# Author = John Gray
# Email = greyjohn15@gmail.com
# Last edit = 31/08/23

# Load packages ----

library(randomForest)
library(ggplot2)
library(dplyr)
library(pROC)

# Load in RF MODEL INPUT AND TARGET VARIABLE DATASETS FOR FLORIDA ----

# Load in abiotic input data for florida (this dataset also contains target
# target variable establishment score and previous establishment input values)
abiotic_data <- read.csv("data/rf_model_abiotic_data.csv")

# Load in biotic variable input data
biotic_data <- read.csv("data/rf_biotic_data.csv")

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
refined_data$interesting_cases <- refined_data$tmax_ring_1 + refined_data$tmax_ring_2 + refined_data$previous_establishment + refined_data$establishment_score
refined_data <- filter(refined_data, interesting_cases != 0)

# get rid of 2003 from data because this is when the goose established itself 
# first and we are not trying to predict cases of initial establishment
refined_data <- filter(refined_data, year != 2003)

### load in train_test_split ----

# this is a csv file containing a table with 70% 1s and 30% 2s, which will be
# used to assign input data to a training or test dataset
train_test_split <- read.csv("data/train_test_split_standard.csv")
train_test_split <- train_test_split$x

#### Now we train models on these data ----

## All parameters version of rf model ----

# Choose a random see to keep results consistent across different models
set.seed(100)

# Remove extraneous variables
rf_data <- refined_data[,-c(1,2,4,26)]

# split data into training andd test datasets according to train/test split
training_data <- rf_data[train_test_split==1,]
test_data <- rf_data[train_test_split==2,]

# fit all parameter model to training data
rf.fit <- randomForest(formula = establishment_score ~ ., data=training_data, 
                       ntree=1000, importance = TRUE, proximity = TRUE)

# have a look at the results
rf.fit
importance(rf.fit)

# create objects for important stats
mse_training_all_parameters <- rf.fit$mse[1000]
var_explained_training_all_parameters <- rf.fit$rsq[1000]

# examine predictions on training data with plot comparing predicted
# establishment according to rf model and actual establishment score
training_data$predicted_establishment <- predict(rf.fit, training_data)

ggplot(training_data, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - all parameters", 
       subtitle = "Training data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/all_parameters_training_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# examine predictions on test data with plot comparing predicted establishment
# according to rf model and actual establishment score
test_data$predicted_establishment <- predict(rf.fit, test_data)

florida_test_pred <- ggplot(test_data, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  theme_classic() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "Florida Test Dataset ", 
       x = "Predicted establishment score (SP model output)",
       y = "True establishment score (ET model output)") +
  coord_cartesian(xlim = c(0,0.25), ylim = c(0, 0.25))

# save this plot
ggsave("rf_model_v5_plots/all_parameters_test_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate mse for test data
mse_test_all_parameters <- mean((test_data$predicted_establishment-test_data$establishment_score)^2)

# Calculate percentage Variance (equivalent of r2) for test data
mean_test_response <- mean(test_data$establishment_score)

TSS <- sum((test_data$establishment_score - mean_test_response)^2)

RSS <- sum((test_data$establishment_score - test_data$predicted_establishment)^2)

var_explained_test_all_parameters <- (1 - RSS / TSS) 

## just ring 1 and previous establishment version of rf model ----

# reduce input data down to appropriate variables
rf_data_ring_1_pe <- refined_data[,-c(1,2,4,15,16,17,18,19,20,21,22,23,25,26)]

# split data into training and test datasets according to train/test split
training_data_ring_1_pe <- rf_data_ring_1_pe[train_test_split==1,]
test_data_ring_1_pe <- rf_data_ring_1_pe[train_test_split==2,]

# fit ring 1 and pe model to training data
rf.fit_ring1_pe <- randomForest(formula = establishment_score ~ ., data=training_data_ring_1_pe, 
                       ntree=1000, importance = TRUE, proximity = TRUE)

# have a look at the results
rf.fit_ring1_pe

importance(rf.fit_ring1_pe)

# create objects for important stats
mse_training_ring1_pe <- rf.fit_ring1_pe$mse[1000]

var_explained_training_ring1_pe <- rf.fit_ring1_pe$rsq[1000]

# examine predictions on training data with plot comparing predicted
# establishment according to rf model and actual establishment score
training_data_ring_1_pe$predicted_establishment <- predict(rf.fit_ring1_pe, training_data_ring_1_pe)

ggplot(training_data_ring_1_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - 1st ring and previous establishment", 
       subtitle = "Training data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/ring1_pe_training_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# examine predictions on test data with plot comparing predicted establishment
# according to rf model and actual establishment score
test_data_ring_1_pe$predicted_establishment <- predict(rf.fit_ring1_pe, test_data_ring_1_pe)

ggplot(test_data_ring_1_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - 1st ring and previous establishment", 
       subtitle = "Test data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/ring1_pe_test_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate mse for test data
mse_test_ring1_pe <- mean((test_data_ring_1_pe$predicted_establishment-test_data_ring_1_pe$establishment_score)^2)

# Calculate percentage Variance (equivalent of r2) for test data
mean_test_ring1_pe_response <- mean(test_data_ring_1_pe$establishment_score)

TSS <- sum((test_data_ring_1_pe$establishment_score - mean_test_ring1_pe_response)^2)

RSS_ring_1_pe <- sum((test_data_ring_1_pe$establishment_score - test_data_ring_1_pe$predicted_establishment)^2)

var_explained_test_ring_1_pe <- (1 - RSS_ring_1_pe / TSS)

## Just ring 1 version of rf model ----

# reduce input data down to appropriate variables
rf_data_ring1 <- refined_data[,-c(1,2,4,14,15,16,17,18,19,20,21,22,23,25,26)]

# split data into training and test datasets according to train/test split
training_data_ring1 <- rf_data_ring1[train_test_split==1,]
test_data_ring1 <- rf_data_ring1[train_test_split==2,]

# fit ring 1 model to training data
rf.fit_ring1 <- randomForest(formula = establishment_score ~ ., data=training_data_ring1, 
                                ntree=1000, importance = TRUE, proximity = TRUE)

# have a look at the results
rf.fit_ring1

importance(rf.fit_ring1)

# create objects for important stats
mse_training_ring1 <- rf.fit_ring1$mse[1000]

var_explained_training_ring1 <- rf.fit_ring1$rsq[1000]

# examine predictions on training data with plot comparing predicted
# establishment according to rf model and actual establishment score
training_data_ring1$predicted_establishment <- predict(rf.fit_ring1, training_data_ring1)

ggplot(training_data_ring1, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - 1st ring", 
       subtitle = "Training data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/ring1_training_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# examine predictions on test data with plot comparing predicted establishment
# according to rf model and actual establishment score
test_data_ring1$predicted_establishment <- predict(rf.fit_ring1, test_data_ring1)

ggplot(test_data_ring1, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - 1st ring", 
       subtitle = "Test data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/ring1_test_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate mse for test data
mse_test_ring1 <- mean((test_data_ring1$predicted_establishment-test_data_ring1$establishment_score)^2)

# Calculate percentage Variance (equivalent of r2) for test data
mean_test_ring1_response <- mean(test_data_ring1$establishment_score)

TSS <- sum((test_data_ring1$establishment_score - mean_test_ring1_response)^2)

RSS_ring1 <- sum((test_data_ring1$establishment_score - test_data_ring1$predicted_establishment)^2)

var_explained_test_ring1 <- (1 - RSS_ring1 / TSS)

## Just ring 2 and previous establishment version of rf model ----

# reduce input data down to appropriate variables
rf_data_ring2_pe <- refined_data[,c(3,14,15,16,17,18,19,20,21,22,23,25)]

# split data into training and test datasets according to train/test split
training_data_ring2_pe <- rf_data_ring2_pe[train_test_split==1,]
test_data_ring2_pe <- rf_data_ring2_pe[train_test_split==2,]

# fit ring 2 and pe model to training data
rf.fit_ring2_pe <- randomForest(formula = establishment_score ~ ., data=training_data_ring2_pe, 
                             ntree=1000, importance = TRUE, proximity = TRUE)

# have a look at the results
rf.fit_ring2_pe

importance(rf.fit_ring2_pe)

# create objects for important stats
mse_training_ring2_pe <- rf.fit_ring2_pe$mse[1000]

var_explained_training_ring2_pe <- rf.fit_ring2_pe$rsq[1000]

# examine predictions on training data with plot comparing predicted
# establishment according to rf model and actual establishment score
training_data_ring2_pe$predicted_establishment <- predict(rf.fit_ring2_pe, training_data_ring2_pe)

ggplot(training_data_ring2_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - 2nd ring and previous establishment", 
       subtitle = "Training data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/ring2_pe_training_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# examine predictions on test data with plot comparing predicted establishment
# according to rf model and actual establishment score
test_data_ring2_pe$predicted_establishment <- predict(rf.fit_ring2_pe, test_data_ring2_pe)

ggplot(test_data_ring2_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - 2nd ring and previous establishment", 
       subtitle = "Test data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/ring2_pe_test_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate mse for test data
mse_test_ring2_pe <- mean((test_data_ring2_pe$predicted_establishment-test_data_ring2_pe$establishment_score)^2)

# Calculate percentage Variance (equivalent of r2) for test data
mean_test_ring2_pe_response <- mean(test_data_ring2_pe$establishment_score)

TSS <- sum((test_data_ring2_pe$establishment_score - mean_test_ring2_pe_response)^2)

RSS_ring2_pe <- sum((test_data_ring2_pe$establishment_score - test_data_ring2_pe$predicted_establishment)^2)

var_explained_test_ring2_pe <- (1 - RSS_ring2_pe / TSS)

## Just ring 2 version of rf model ----

# reduce input data down to appropriate variables
rf_data_ring2 <- refined_data[,c(3,15,16,17,18,19,20,21,22,23,25)]

# split data into training and test datasets according to train/test split
training_data_ring2 <- rf_data_ring2[train_test_split==1,]
test_data_ring2 <- rf_data_ring2[train_test_split==2,]

# fit ring 2 model to training data
rf.fit_ring2 <- randomForest(formula = establishment_score ~ ., data=training_data_ring2, 
                                ntree=1000, importance = TRUE, proximity = TRUE)

# have a look at the results
rf.fit_ring2

importance(rf.fit_ring2)

# create objects for important stats
mse_training_ring2 <- rf.fit_ring2$mse[1000]

var_explained_training_ring2 <- rf.fit_ring2$rsq[1000]

# examine predictions on training data with plot comparing predicted
# establishment according to rf model and actual establishment score
training_data_ring2$predicted_establishment <- predict(rf.fit_ring2, training_data_ring2)

ggplot(training_data_ring2, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - 2nd ring", 
       subtitle = "Training data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/ring2_training_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# examine predictions on test data with plot comparing predicted establishment
# according to rf model and actual establishment score
test_data_ring2$predicted_establishment <- predict(rf.fit_ring2, test_data_ring2)

ggplot(test_data_ring2, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - 2nd ring", 
       subtitle = "Test data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/ring2_test_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate mse for test data
mse_test_ring2 <- mean((test_data_ring2$predicted_establishment-test_data_ring2$establishment_score)^2)

# Calculate percentage Variance (equivalent of r2) for test data
mean_test_ring2_response <- mean(test_data_ring2$establishment_score)

TSS <- sum((test_data_ring2$establishment_score - mean_test_ring2_response)^2)

RSS_ring2 <- sum((test_data_ring2$establishment_score - test_data_ring2$predicted_establishment)^2)

var_explained_test_ring2 <- (1 - RSS_ring2 / TSS)

## JUST BIOTIC and Previous establishment version of rf model ----

# reduce input data down to appropriate variables
rf_data_biotic_pe <- refined_data[,c(3,14,24,25)]

# split data into training and test datasets according to train/test split
training_data_biotic_pe <- rf_data_biotic_pe[train_test_split==1,]
test_data_biotic_pe <- rf_data_biotic_pe[train_test_split==2,]

# fit biotic and previous establishment model to training data
rf.fit_biotic_pe <- randomForest(formula = establishment_score ~ ., data=training_data_biotic_pe, 
                             ntree=1000, importance = TRUE, proximity = TRUE)

# have a look at the results
rf.fit_biotic_pe

importance(rf.fit_biotic_pe)

# create objects for important stats
mse_training_biotic_pe <- rf.fit_biotic_pe$mse[1000]

var_explained_training_biotic_pe <- rf.fit_biotic_pe$rsq[1000]

# examine predictions on training data with plot comparing predicted
# establishment according to rf model and actual establishment score
training_data_biotic_pe$predicted_establishment <- predict(rf.fit_biotic_pe, training_data_biotic_pe)

ggplot(training_data_biotic_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - Biotic and previous establishment parameters", 
       subtitle = "Training data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/biotic_pe_training_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# examine predictions on test data with plot comparing predicted establishment
# according to rf model and actual establishment score
test_data_biotic_pe$predicted_establishment <- predict(rf.fit_biotic_pe, test_data_biotic_pe)

ggplot(test_data_biotic_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - Biotic and previous establishment parameters", 
       subtitle = "Test data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/biotic_pe_test_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate mse for test data
mse_test_biotic_pe <- mean((test_data_biotic_pe$predicted_establishment-test_data_biotic_pe$establishment_score)^2)

# Calculate percentage Variance (equivalent of r2) for test data
mean_test_biotic_pe_response <- mean(test_data_biotic_pe$establishment_score)

TSS <- sum((test_data_biotic_pe$establishment_score - mean_test_biotic_pe_response)^2)

RSS_biotic_pe <- sum((test_data_biotic_pe$establishment_score - test_data_biotic_pe$predicted_establishment)^2)

var_explained_test_biotic_pe <- (1 - RSS_biotic_pe / TSS)

## JUST BIOTIC version of rf model ----

# reduce input data down to appropriate variables
rf_data_biotic <- refined_data[,c(3,24,25)]

# split data into training and test datasets according to train/test split
training_data_biotic <- rf_data_biotic[train_test_split==1,]
test_data_biotic <- rf_data_biotic[train_test_split==2,]

# fit biotic model to training data
rf.fit_biotic <- randomForest(formula = establishment_score ~ ., data=training_data_biotic, 
                                 ntree=1000, importance = TRUE, proximity = TRUE)

# have a look at the results
rf.fit_biotic

importance(rf.fit_biotic)

# create objects for important stats
mse_training_biotic <- rf.fit_biotic$mse[1000]

var_explained_training_biotic <- rf.fit_biotic$rsq[1000]

# examine predictions on training data with plot comparing predicted
# establishment according to rf model and actual establishment score
training_data_biotic$predicted_establishment <- predict(rf.fit_biotic, training_data_biotic)

ggplot(training_data_biotic, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - Biotic parameters", 
       subtitle = "Training data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/biotic_training_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# examine predictions on test data with plot comparing predicted establishment
# according to rf model and actual establishment score
test_data_biotic$predicted_establishment <- predict(rf.fit_biotic, test_data_biotic)

ggplot(test_data_biotic, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - Biotic parameters", 
       subtitle = "Test data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/biotic_test_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate mse for test data
mse_test_biotic <- mean((test_data_biotic$predicted_establishment-test_data_biotic$establishment_score)^2)

# Calculate percentage Variance (equivalent of r2) for test data
mean_test_biotic_response <- mean(test_data_biotic$establishment_score)

TSS <- sum((test_data_biotic$establishment_score - mean_test_biotic_response)^2)

RSS_biotic <- sum((test_data_biotic$establishment_score - test_data_biotic$predicted_establishment)^2)

var_explained_test_biotic <- (1 - RSS_biotic / TSS)

## JUST ABIOTIC and previous establishment version of rf model ----

# reduce input data down to appropriate variables
rf_data_abiotic_pe <- refined_data[,c(3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)]

# split data into training and test datasets according to train/test split
training_data_abiotic_pe <- rf_data_abiotic_pe[train_test_split==1,]
test_data_abiotic_pe <- rf_data_abiotic_pe[train_test_split==2,]

# fit abiotic and previous establishment model to training data
rf.fit_abiotic_pe <- randomForest(formula = establishment_score ~ ., data=training_data_abiotic_pe, 
                              ntree=1000, importance = TRUE, proximity = TRUE)

# have a look at the results
rf.fit_abiotic_pe

importance(rf.fit_abiotic_pe)

# create objects for important stats
mse_training_abiotic_pe <- rf.fit_abiotic_pe$mse[1000]

var_explained_training_abiotic_pe <- rf.fit_abiotic_pe$rsq[1000]

# examine predictions on training data with plot comparing predicted
# establishment according to rf model and actual establishment score
training_data_abiotic_pe$predicted_establishment <- predict(rf.fit_abiotic_pe, training_data_abiotic_pe)

ggplot(training_data_abiotic_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - abiotic and previous establishment parameters", 
       subtitle = "Training data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/abiotic_pe_training_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# examine predictions on test data with plot comparing predicted establishment
# according to rf model and actual establishment score
test_data_abiotic_pe$predicted_establishment <- predict(rf.fit_abiotic_pe, test_data_abiotic_pe)

ggplot(test_data_abiotic_pe, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - abiotic and previous establishment parameters", 
       subtitle = "Test data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/abiotic_pe_test_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate mse for test data
mse_test_abiotic_pe <- mean((test_data_abiotic_pe$predicted_establishment-test_data_abiotic_pe$establishment_score)^2)

# Calculate percentage Variance (equivalent of r2) for test data
mean_test_abiotic_pe_response <- mean(test_data_abiotic_pe$establishment_score)

TSS <- sum((test_data_abiotic_pe$establishment_score - mean_test_abiotic_pe_response)^2)

RSS_abiotic_pe <- sum((test_data_abiotic_pe$establishment_score - test_data_abiotic_pe$predicted_establishment)^2)

var_explained_test_abiotic_pe <- (1 - RSS_abiotic_pe / TSS)

## JUST ABIOTIC version of rf model ----

# reduce input data down to appropriate variables
rf_data_abiotic <- refined_data[,c(3,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23)]

# split data into training and test datasets according to train/test split
training_data_abiotic <- rf_data_abiotic[train_test_split==1,]
test_data_abiotic <- rf_data_abiotic[train_test_split==2,]

# fit abiotic model to training data
rf.fit_abiotic <- randomForest(formula = establishment_score ~ ., data=training_data_abiotic, 
                                  ntree=1000, importance = TRUE, proximity = TRUE)

# have a look at the results
rf.fit_abiotic

importance(rf.fit_abiotic)

# create objects for important stats
mse_training_abiotic <- rf.fit_abiotic$mse[1000]

var_explained_training_abiotic <- rf.fit_abiotic$rsq[1000]

# examine predictions on training data with plot comparing predicted
# establishment according to rf model and actual establishment score
training_data_abiotic$predicted_establishment <- predict(rf.fit_abiotic, training_data_abiotic)

ggplot(training_data_abiotic, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - abiotic parameters", 
       subtitle = "Training data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/abiotic_training_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# examine predictions on test data with plot comparing predicted establishment
# according to rf model and actual establishment score
test_data_abiotic$predicted_establishment <- predict(rf.fit_abiotic, test_data_abiotic)

ggplot(test_data_abiotic, aes(x = predicted_establishment, y = establishment_score)) +
  geom_point() + 
  geom_abline(aes(intercept = 0, slope = 1)) +
  labs(title = "RF model v5 - abiotic establishment parameters", 
       subtitle = "Test data: predicted establishment vs establishment")

# save this plot
ggsave("rf_model_v5_plots/abiotic_test_predicted_establishment_vs_establishment.png", 
       width = 10, height = 10)

# calculate mse for test data
mse_test_abiotic <- mean((test_data_abiotic$predicted_establishment-test_data_abiotic$establishment_score)^2)

# Calculate percentage Variance (equivalent of r2) for test data
mean_test_abiotic_response <- mean(test_data_abiotic$establishment_score)

TSS <- sum((test_data_abiotic$establishment_score - mean_test_abiotic_response)^2)

RSS_abiotic <- sum((test_data_abiotic$establishment_score - test_data_abiotic$predicted_establishment)^2)

var_explained_test_abiotic <- (1 - RSS_abiotic / TSS)

#### SAVE all RF MODEL versions ----

save(rf.fit, file = "all_parameters_rf_v5.RData")

save(rf.fit_ring1_pe, file = "ring1_pe_rf_v5.RData")

save(rf.fit_ring1, file = "ring1_rf_v5.RData")

save(rf.fit_ring2_pe, file = "ring2_pe_rf_v5.RData")

save(rf.fit_ring2, file = "ring2_rf_v5.RData")

save(rf.fit_biotic_pe, file = "biotic_previous_establishment_rf_v5.RData")

save(rf.fit_biotic, file = "biotic_rf_v5.RData")

save(rf.fit_abiotic_pe, file = "abiotic_previous_establishment_rf_v5.RData")

save(rf.fit_abiotic, file = "abiotic_rf_v5.RData")

#### create a dataframe which provides a straightforward performance ----
#### comparison for the different versions of the rf model ----

rf_stats <- data.frame("Model" = c("all_parameters", 
                                   "ring_1_previous_establishment", "ring_1",
                                   "ring_2_previous_establishment", "ring_2",
                                   "biotic_previous_establishment", "biotic",
                                   "abiotic_previous_establishment", "abiotic"),
                       "MSE_on_training_data" = c(mse_training_all_parameters,
                                                  mse_training_ring1_pe,
                                                  mse_training_ring1,
                                                  mse_training_ring2_pe,
                                                  mse_training_ring2,
                                                  mse_training_biotic_pe,
                                                  mse_training_biotic,
                                                  mse_training_abiotic_pe,
                                                  mse_training_abiotic),
                       "%_var_explained_training_data" = c(var_explained_training_all_parameters,
                                                           var_explained_training_ring1_pe,
                                                           var_explained_training_ring1,
                                                           var_explained_training_ring2_pe,
                                                           var_explained_training_ring2,
                                                           var_explained_training_biotic_pe,
                                                           var_explained_training_biotic,
                                                           var_explained_training_abiotic_pe,
                                                           var_explained_training_abiotic),
                       "MSE_on_test_data" = c(mse_test_all_parameters,
                                              mse_test_ring1_pe, mse_test_ring1,
                                              mse_test_ring2_pe, mse_test_ring2,
                                              mse_test_biotic_pe, mse_test_biotic, 
                                              mse_test_abiotic_pe, mse_test_abiotic),
                       "%_var_explained_test_data" = c(var_explained_test_all_parameters,
                                                       var_explained_test_ring_1_pe,
                                                       var_explained_test_ring1,
                                                       var_explained_test_ring2_pe,
                                                       var_explained_test_ring2,
                                                       var_explained_test_biotic_pe,
                                                       var_explained_test_biotic,
                                                       var_explained_test_abiotic_pe,
                                                       var_explained_test_abiotic))

# save this performance comparison df in a csv
write.csv(rf_stats, "data/rf_models_v5_comparison.csv")

#### Analysis of capabilities of spread prediction----

# We want to know: 
# percentage/number of spread events correctly predicted
# percentage/number of spread events incorrectly predicted
# accuracy when predicting spread events 

# ignore ring 2 models from now on cuz they're just a lot worse

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
# establishment for every version of the rf model
test_data_predictions <- cbind(test_data[,c(1,11,23)],
                               test_data_ring_1_pe$predicted_establishment,
                               test_data_ring1$predicted_establishment,
                               test_data_ring2_pe$predicted_establishment,
                               test_data_ring2$predicted_establishment,
                               test_data_biotic_pe$predicted_establishment,
                               test_data_biotic$predicted_establishment,
                               test_data_abiotic_pe$predicted_establishment,
                               test_data_abiotic$predicted_establishment)

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
                                                                 (test_data_ring_1_pe$predicted_establishment > 0.005)))


detailed_analysis$correct_spread_prediction_percentage[2] <- 100 * nrow(filter(test_data_predictions, 
                                                                               (establishment_score > 0.005) &
                                                                                 (previous_establishment < 0.005) &
                                                                                 (test_data_ring_1_pe$predicted_establishment > 0.005))) /
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[3] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_ring1$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[3] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_ring1$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[4] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_ring2_pe$predicted_establishment > 0.005)))


detailed_analysis$correct_spread_prediction_percentage[4] <- 100 * nrow(filter(test_data_predictions, 
                                                                               (establishment_score > 0.005) &
                                                                                 (previous_establishment < 0.005) &
                                                                                 (test_data_ring2_pe$predicted_establishment > 0.005))) /
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[5] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_ring2$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[5] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_ring2$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[6] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_biotic_pe$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[6] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_biotic_pe$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[7] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_biotic$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[7] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_biotic$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[8] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_abiotic_pe$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[8] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_abiotic_pe$predicted_establishment > 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score > 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_spread_predictions[9] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score > 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_abiotic$predicted_establishment > 0.005)))

detailed_analysis$correct_spread_prediction_percentage[9] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score > 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_abiotic$predicted_establishment > 0.005)))/
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
                                                                 (test_data_ring_1_pe$predicted_establishment < 0.005)))


detailed_analysis$correct_no_spread_prediction_percentage[2] <- 100 * nrow(filter(test_data_predictions, 
                                                                               (establishment_score < 0.005) &
                                                                                 (previous_establishment < 0.005) &
                                                                                 (test_data_ring_1_pe$predicted_establishment < 0.005))) /
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[3] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score < 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_ring1$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[3] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score < 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_ring1$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[4] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score < 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_ring2_pe$predicted_establishment < 0.005)))


detailed_analysis$correct_no_spread_prediction_percentage[4] <- 100 * nrow(filter(test_data_predictions, 
                                                                               (establishment_score < 0.005) &
                                                                                 (previous_establishment < 0.005) &
                                                                                 (test_data_ring2_pe$predicted_establishment < 0.005))) /
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[5] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score < 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_ring2$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[5] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score < 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_ring2$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[6] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score < 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_biotic_pe$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[6] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score < 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_biotic_pe$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[7] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score < 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_biotic$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[7] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score < 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_biotic$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[8] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score < 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_abiotic_pe$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[8] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score < 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_abiotic_pe$predicted_establishment < 0.005)))/
  nrow(filter(test_data_predictions, (establishment_score < 0.005) & (previous_establishment < 0.005)))

detailed_analysis$correct_no_spread_predictions[9] <- nrow(filter(test_data_predictions, 
                                                               (establishment_score < 0.005) &
                                                                 (previous_establishment < 0.005) &
                                                                 (test_data_abiotic$predicted_establishment < 0.005)))

detailed_analysis$correct_no_spread_prediction_percentage[9] <- 100*nrow(filter(test_data_predictions, 
                                                                             (establishment_score < 0.005) &
                                                                               (previous_establishment < 0.005) &
                                                                               (test_data_abiotic$predicted_establishment < 0.005)))/
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
detailed_analysis$spread <- 3 

# calculate overall spread prediction accuracy for all models
detailed_analysis$accuracy <- 100*(detailed_analysis$correct_spread_predictions + detailed_analysis$correct_no_spread_predictions) / (detailed_analysis$no_spread + detailed_analysis$spread)

# save detailed analysis table as a csv
write.csv(detailed_analysis, "data/rf_models_v5_predictive_analysis.csv")

# save test data predictions table as a csv
write.csv(test_data_predictions, "data/rf_models_v5_test_data_predictions.csv")

### calculating AUC SCORE ----

# create auc dataframe for Florida test predictions to be populated with scores
auc <- data.frame(model, "auc_all_test_data" = numeric(9), "auc_test_spread_data" = numeric(9))

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

# we calculate AUC scores both for all florida test data, and for only the cells
# where there is a possibility of spread, which provides more insight into the
# model's spread prediction capabilities

# filter test data to just cells where there is a possibility of spread
test_spread_data_predictions <- filter(test_data_predictions, previous_establishment < 0.005)

# calculate auc scores for each version of the model for both overall and just
# possible spread datasets
auc$auc_all_test_data[1] <- auc(test_data_predictions$established, test_data_predictions$predicted_establishment)
auc$auc_test_spread_data[1] <- auc(test_spread_data_predictions$established, test_spread_data_predictions$predicted_establishment)

auc$auc_all_test_data[2] <- auc(test_data_predictions$established, test_data_predictions[,4])
auc$auc_test_spread_data[2] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,4])

auc$auc_all_test_data[3] <- auc(test_data_predictions$established, test_data_predictions[,5])
auc$auc_test_spread_data[3] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,5])

auc$auc_all_test_data[4] <- auc(test_data_predictions$established, test_data_predictions[,6])
auc$auc_test_spread_data[4] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,6])

auc$auc_all_test_data[5] <- auc(test_data_predictions$established, test_data_predictions[,7])
auc$auc_test_spread_data[5] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,7])

auc$auc_all_test_data[6] <- auc(test_data_predictions$established, test_data_predictions[,8])
auc$auc_test_spread_data[6] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,8])

auc$auc_all_test_data[7] <- auc(test_data_predictions$established, test_data_predictions[,9])
auc$auc_test_spread_data[7] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,9])

auc$auc_all_test_data[8] <- auc(test_data_predictions$established, test_data_predictions[,10])
auc$auc_test_spread_data[8] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,10])

auc$auc_all_test_data[9] <- auc(test_data_predictions$established, test_data_predictions[,11])
auc$auc_test_spread_data[9] <- auc(test_spread_data_predictions$established, test_spread_data_predictions[,11])

# save auc scores in a csv file
write.csv(auc, "data/florida_test_data_auc_scores.csv")

# save overall rf input + target variable data in a csv file
write.csv(rf_data, "data/rf_data_overall_rf_model_florida_v5.csv")

# save florida training data in a csv file
write.csv(training_data, "data/training_data_overall_rf_model_florida_v5.csv")

# save florida test data in a csv file
write.csv(test_data, "data/test_data_overall_rf_model_flordia_v5.csv")

# save importance values for all parameters version of the model in a csv file
imp <- as.data.frame(importance(rf.fit))
write.csv(imp, "data/rf_model_v5_variable_importance.csv")

# create and save partial dependence plots for every all parameters model ----
# input variable ----

# bs ring 1
png(filename = "pdp_plot_bs_ring_1_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = bs_ring_1,
            rug = TRUE, 
            xlab = "Establishment-weighted biotic similarity in adjacent cells",
            ylab = "Establishment score",
            main = "Partial dependence of establishment score on biotic similarity in adjacent cells for RF Invasive spread model"))
dev.off()

# previous establishment
png(filename = "pdp_plot_previous_establishment_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = previous_establishment,
                  rug = TRUE, 
                  xlab = "Establishment score in previous year in target cell",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on establishment score in previous year"))
dev.off()

# bs ring 2
png(filename = "pdp_plot_bs_ring_2_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = bs_ring_2,
                  rug = TRUE, 
                  xlab = "Establishment-weighted biotic similarity in 1-degree_separation cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on biotic similarity in 1-degree-separation cells for RF Invasive spread model"))
dev.off()

# tmax ring 1
png(filename = "pdp_plot_tmax_ring_1_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = tmax_ring_1,
                  rug = TRUE, 
                  xlab = "Establishment-weighted temperature similarity in adjacent cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on temperature similarity in adjacent cells for RF Invasive spread model"))
dev.off()

# prcp ring 1
png(filename = "pdp_plot_prcp_ring_1_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = prcp_ring_1,
                  rug = TRUE, 
                  xlab = "Establishment-weighted precipitation similarity in adjacent cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on precipitation similarity in adjacent cells for RF Invasive spread model"))
dev.off()

# water ring 1
png(filename = "pdp_plot_water_ring_1_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = water_ring_1,
                  rug = TRUE, 
                  xlab = "Establishment-weighted water area similarity in adjacent cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on water area similarity in adjacent cells for RF Invasive spread model"))
dev.off()

# forest ring 1
png(filename = "pdp_plot_forest_ring_1_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = forest_ring_1,
                  rug = TRUE, 
                  xlab = "Establishment-weighted forest area similarity in adjacent cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on forest area similarity in adjacent cells for RF Invasive spread model"))
dev.off()

# grass ring 1
png(filename = "pdp_plot_grass_ring_1_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = grass_ring_1,
                  rug = TRUE, 
                  xlab = "Establishment-weighted grass area similarity in adjacent cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on grass area similarity in adjacent cells for RF Invasive spread model"))
dev.off()

# wetland ring 1
png(filename = "pdp_plot_wetland_ring_1_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = wetland_ring_1,
                  rug = TRUE, 
                  xlab = "Establishment-weighted wetland area similarity in adjacent cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on wetland area similarity in adjacent cells for RF Invasive spread model"))
dev.off()

# farming ring 1
png(filename = "pdp_plot_farming_ring_1_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = farming_ring_1,
                  rug = TRUE, 
                  xlab = "Establishment-weighted farming area similarity in adjacent cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on farming area similarity in adjacent cells for RF Invasive spread model"))
dev.off()

# urban ring 1
png(filename = "pdp_plot_urban_ring_1_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = urban_ring_1,
                  rug = TRUE, 
                  xlab = "Establishment-weighted urban area similarity in adjacent cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on urban area similarity in adjacent cells for RF Invasive spread model"))
dev.off()

# barren ring 1
png(filename = "pdp_plot_barren_ring_1_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = barren_ring_1,
                  rug = TRUE, 
                  xlab = "Establishment-weighted barren area similarity in adjacent cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on barren area similarity in adjacent cells for RF Invasive spread model"))
dev.off()

## 2nd RING PLOTS

# tmax ring 2
png(filename = "pdp_plot_tmax_ring_2_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = tmax_ring_2,
                  rug = TRUE, 
                  xlab = "Establishment-weighted temperature similarity in 1-degree-separation cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on temperature similarity in 1-degree-separation cells for RF Invasive spread model"))
dev.off()

# prcp ring 2
png(filename = "pdp_plot_prcp_ring_2_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = prcp_ring_2,
                  rug = TRUE, 
                  xlab = "Establishment-weighted precipitation similarity in 1-degree-separation cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on precipitation similarity in 1-degree-separation cells for RF Invasive spread model"))
dev.off()

# water ring 2
png(filename = "pdp_plot_water_ring_2_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = water_ring_2,
                  rug = TRUE, 
                  xlab = "Establishment-weighted water area similarity in 1-degree-separation cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on water area similarity in 1-degree-separation cells for RF Invasive spread model"))
dev.off()

# forest ring 2
png(filename = "pdp_plot_forest_ring_2_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = forest_ring_2,
                  rug = TRUE, 
                  xlab = "Establishment-weighted forest area similarity in 1-degree-separation cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on forest area similarity in 1-degree-separation cells for RF Invasive spread model"))
dev.off()

# grass ring 2
png(filename = "pdp_plot_grass_ring_2_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = grass_ring_2,
                  rug = TRUE, 
                  xlab = "Establishment-weighted grass area similarity in 1-degree-separation cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on grass area similarity in 1-degree-separation cells for RF Invasive spread model"))
dev.off()

# wetland ring 2
png(filename = "pdp_plot_wetland_ring_2_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = wetland_ring_2,
                  rug = TRUE, 
                  xlab = "Establishment-weighted wetland area similarity in 1-degree-separation cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on wetland area similarity in 1-degree-separation cells for RF Invasive spread model"))
dev.off()

# farming ring 2
png(filename = "pdp_plot_farming_ring_2_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = farming_ring_2,
                  rug = TRUE, 
                  xlab = "Establishment-weighted farming area similarity in 1-degree-separation cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on farming area similarity in 1-degree-separation cells for RF Invasive spread model"))
dev.off()

# urban ring 2
png(filename = "pdp_plot_urban_ring_2_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = urban_ring_2,
                  rug = TRUE, 
                  xlab = "Establishment-weighted urban area similarity in 1-degree-separation cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on urban area similarity in 1-degree-separation cells for RF Invasive spread model"))
dev.off()

# barren ring 2
png(filename = "pdp_plot_barren_ring_2_overall_rf.png", width = 1000, height = 1000)
print(partialPlot(x = rf.fit, pred.data = training_data, x.var = barren_ring_2,
                  rug = TRUE, 
                  xlab = "Establishment-weighted barren area similarity in 1-degree-separation cells",
                  ylab = "Establishment score",
                  main = "Partial dependence of establishment score on barren area similarity in 1-degree-separation cells for RF Invasive spread model"))
dev.off()

## create grid of noteworthy partial dependence results ---

# first need to turn the pdp's into dataframes so I can plot using ggplot which
# lets me save them in a grid
pe_df <- as.data.frame(previous_establishment_pdp)
prcp_df <- as.data.frame(prcp_pdp)
farm_cover_df <- as.data.frame(farm_cover_pdp)
forest_cover_df <- as.data.frame(forest_cover_pdp)

# create ggplots for interesting partial dependence
pe_plot <- ggplot(pe_df, aes(x=x, y=y)) +
  geom_path() +
  theme_classic() + 
  labs(title = "Previous establishment in target cell",
       x = "Establishment score in previous year in target cell",
       y = "Establishment score")

prcp_plot <- ggplot(prcp_df, aes(x=x, y=y)) +
  geom_path() +
  theme_classic() + 
  labs(title = "Precipitation in adjacent cells",
       x = "Establishment-weighted precipitation similarity in adjacent cells",
       y = "Establishment score")

farm_cover_plot <- ggplot(farm_cover_df, aes(x=x, y=y)) +
  geom_path() +
  theme_classic() + 
  labs(title = "Farmland cover in adjacent cells",
       x = "Establishment-weighted farming area similarity in adjacent cells",
       y = "Establishment score")

forest_cover_plot <- ggplot(forest_cover_df, aes(x=x, y=y)) +
  geom_path() +
  theme_classic() + 
  labs(title = "Forest cover in adjacent cells",
       x = "Establishment-weighted forest area similarity in adjacent cells",
       y = "Establishment score")

# create and save grid of plots for interesting pdps
grid.arrange(pe_plot, prcp_plot, farm_cover_plot, forest_cover_plot, nrow = 2)
g <- arrangeGrob(pe_plot, prcp_plot, farm_cover_plot, forest_cover_plot, nrow = 2)
ggsave(file = "final_establishment_plots/combined_plots/PDPs.png", g, 
       width = 10, height = 10)
