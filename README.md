# Invasive_species_spread_masters_research_project
## Author: John Gray
## Email: johnpatrickgray97@btinternet.com

### Summary:
The R scripts contained within this repository are the final versions of those which have been used to conduct analysis for my UCL Ecology & Data Science masters research project, entitled: 'Spread prediction unlocked: how citizen science and machine learning can advance
predicitons of biological invasions'. In this project, I have developed data-driven random forest models that use processed eBird data to make predictions about the future spread of invasive species. The models perform very well at making predictions of future spread for
invasions which they have been trained upon, and adequately for invasions where the model has been trained on a different invasion of the same species.

Scripts have been grouped into the folders provided here according to their function. These groups are outlined below:

### EGoose_establishment_calc_FL
These scripts use eBird observation data to calculate annual local establishment scores for the Egyptian goose in distinct hexagonal cells across Florida from 2002-2019 

### EGoose_establishment_calc_TX
These scripts use eBird observation data to calculate annual local establishment scores for the Egyptian goose in distinct hexagonal cells across Texas from 2002-2019

### abiotic_inputs_processing_FL
These scripts use various sources of environmental data, coupled with Florida establishment scores, to generate input data related to abiotic variables for a random forest model predicting spread of the Egyptian goose in Florida

### abiotic_inputs_processing_TX
These scripts use various sources of environmental data, coupled with Texas establishment scores, to generate input data related to abiotic variables for a random forest model predicting spread of the Egyptian goose in Texas

### biotic_inputs_processing_FL
These scripts use eBird observation data to determine native species assemblages in distinct hexagonal cells across Florida, which when coupled with Florida establishment scores, are used to generate input data related to biotic variables for a random forest model
predicting spread of the Egyptian goose in Florida

### biotic_inputs_processing_TX
These scripts use eBird observation data to determine native species assemblages in distinct hexagonal cells across Texas, which when coupled with Texas establishment scores, are used to generate input data related to biotic variables for a random forest model predicting 
spread of the Egyptian goose in Texas

### spread_predictor_modelling
These scripts train random forst models on a Florida training dataset, test their predictions on a Florida test dataset and Texas dataset, and compare their performance with an additive linear regression model for Florida test and Texas datasets
