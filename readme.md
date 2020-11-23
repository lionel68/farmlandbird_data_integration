# Repository to reproduce analyze and main figures from Data Integration manuscript

## Description

In this repository can be found R-scripts to reproduce the analysis and the main figures shown in the manuscript Hertzog et al "XX".

Note that because of restriction in data sharing, only artificial data are provided here, to reproduce these analysis on other datasets make sure that these have the same structure (column names ...) as the artificial dataset .


## Structure of the repo

There are 4 folders in this repo:

	* data: folder to place the .csv data files
	* figures: folder where the figures will be saved
	* model_output: folder where the files resulting from the model fitting will be saved
	* script: folder that contains the following Rscripts:
		1. 01_function_scripts.r: R file with the helper function to set and fit the 6 different models (see Table 2 in the manuscript)
		2. 02_fitting_script.r: R script to load the data and fit to each of the species present in the data file the 6 different models
		3. 03_plotting_script.r: R script to gather the results from the model fitting and reproduce the main figures shown in the manuscript.

