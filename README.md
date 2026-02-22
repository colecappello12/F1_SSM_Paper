# A State-Space Approach to Modeling Tire Degradation in Formula 1 Racing

This repository contains code and data used to create the paper "A State-Space Approach to Modeling Tire Degradation in Formula 1 Racing. There are four folders containing different aspects of the paper creation process.

# Paper_Generation_Files

The file Tire_Deg_Paper_Review1.qmd generates the final paper, though we manually perform some final cleaning of the output tex file to ensure proper formatting. 

The folder itself shoud contain all the necessary results and files to generate the paper from Tire_Deg_Paper_Review1.qmd. 

## Data_Generation_Scripts

This folder contains the python script used to pull all lap time data for Lewis Hamilton in the 2025 season. The output is also contained in the folder.

## Cross_Validation_Scripts_and_Stan_Code

This folder contains two Quarto Markdown documents (.qmd) that generate the main results of the paper. The necessary stan code and dataset to fit the models is also included.

- Full_Race_CV.qmd contains code that runs the initial model selection procedure on the Austrian Grand Prix and outputs those results as .csv files.

- 2025_CV.qmd contains code that runs the cross-session analysis. It is currently set to run Cross Validation on one race due to the fact that the full CV pipeline takes many hours. The full pipeline can easily be run however by uncommenting a couple lines of code. There are two main codes chunks in this file, the first of which outputs race-level results, and the second aggregates the race-level results into a summary table for the full season. 

## Cross_Validation_Results

This folder contains the results of Full_Race_CV.qmd and 2025_CV.qmd.

- Model_Selection_Results.csv and CRPS_Results.csv contain the results of running Full_Race_CV.qmd

- The remaining files titled "racename"_results1.csv and All_CV_results1.csv contain race-level and aggregated results respectively, obtained by running 2025_CV.qmd


