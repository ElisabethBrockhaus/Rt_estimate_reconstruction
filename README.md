# Rt_estimate_reconstruction #

To plot the estimates as in the thesis run *plot_final_estimates.R*.

To run the estimations based on the incidence data run *parameter_influence_EpiEstim.R* and *compare_parameters_combinations.R*. The Script *parameter_influence_EpiEstim.R* runs the estimations using EpiEstim with different parameters and data. The Script *compare_parameters_combinations.R* runs the estimations using different estimation methods with different levels of adjustment. Some estimations have to be performed outside of this script, either in Python (globalrt, rtlive) or in a more time-intensive script (epiforecasts). Those estimates are contained in this repo and do not have to be calculated again. For all but rtlive, it is possible, though.
- **epiforecasts**: Run *epiforecasts/epinow_estimation.R*. Takes multiple hours depending on the length of input data.
- **globalrt**: Run *ArroyoMarioli/input_output_dataset/estimate_R_filter_and_smoother.py* (Might want to change the import between *dataset.csv* and *dataset_final.csv*. The former is the data used originally for the globalrt estimates. The latter is the correctly formatted aggregated version of the RKI line list data.)
- **rtlive**: The data on the number of tests performed is missing and not publicly available. Thus, the estimation cannot be performed on the basis of the repository.

### How to do all from scratch without using prepared data sets from repo? ###
- **ETH**: Before loading ETH data run the script *ETH/otherScripts/format_linelist_data.R* to calculate the delays. Then load the data with *new_deconvolution = TRUE* in the function *load_incidence_data()* with *method = 'ETHZ_sliding_window'*
- **globalrt**: Run the scripts from *ArroyoMarioli/input_output/dataset/* in the following order: *format_data_from_various_sources.R* -> *construct_dataset_adj.py* -> *construct_STAN_models.py* -> *estimate_R_filter_and_smoother.py*
- **rtlive**: The incidence data used for the main comparison results from running the script *rtlive/important_files_from_repo/notebooks/data_aggregation.ipynb*.
- **SDSC**: The smoothed data can be constructed using the script *SDSC/smoothing_example.ipynb*.

Note that for the execution the RKI's line list data might be necessary. These are to large too be saved to this repository. Thus, they have to be sources from https://hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0/explore. Older versions of these data can be obtained when filtering by "Meldedatum".

### Requirements ###
- *KITmetricslab/reproductive_numbers* repository must be stored in the same directory as this repository.
