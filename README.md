# mb4-analysis
ManyBabies4 analysis

# Description of Project Structure

## *final_analysis* folder

- MB4_part1_data_clean.qmd: This script aggregates the three main data files and prepares the dataset for analysis
- MB4_part2_Bayesian_analysis.qmd: This script performs the main analysis
- MB4_part3_Frequentist_analysis.qmd: This script performs the frequentist analysis as a supplement

## *pilot_simulations* folder

- analysis_structure.Rmd: This document details the pre-registered analysis plan for the MB4 study, and illustrates it with simulated data
- Simulations.Rmd: This document provides evidence from simulated data for choices made in the ManyBabies4 study
- MB4_pilot_analysis.Rmd: This script analyzes data from the initial pilot
- Other R scripts: Dependencies and/or functions for simuliations and plotting during pilot analyses

## *pilot_data* folder

- MB4_compiledData_June2019.csv: Dataset from the initial pilot

## *simulations_results* and *plots* folder

These folders contains results of the pilot analysis and simulations

# Instructions for Reproducing Analyses

Data from the main study is stored on the project's OSF repository. To fully reproduce the analyses performed in the qmd files in the *final_analysis* folder, follow these steps:

1. Download the entire RProject from this repository (in browser window, you can do so by clicking code -> Download ZIP, then unzipping the file with the software of your choice)

2. In the main *mb4-analysis* folder, create two folders: *main_data* and *intermediates*.

3. Download *clean_data.csv*, *cb_orders.csv*, *looking_time_reliability_data.csv*, and *contributing_lab_list.csv* from the OSF repository; place it in the *main_data* folder.

4. Run the codes in the *final_analysis* folder in the following order: *MB4_part1_data_clean.qmd* -> *MB4_part2_Bayesian_analysis.qmd* -> *MB4_part3_Frequentist_analysis.qmd*


