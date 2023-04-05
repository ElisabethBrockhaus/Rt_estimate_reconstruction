# Rt_estimate_reconstruction #

## Within-model temporal coherence of real-time estimates ##

Figure 4 from the paper (real-time estimates from different groups) can be generated using *compare_real_time_estimates.R*. Figure 5 (within-method consistence) results from running *within_method_agreement.R*. First the consistence measures are calculated. The resulting .csv files are already contained in this repository (for days_until_final = 70) and can be used directly for plotting (see bottom of the script).

## Between-method agreement of retrospective estimates ##

To plot the estimates as in Section 4 of the paper run *plots_between_method_agreement.R*.

To run the estimations based on the incidence data run *parameter_influence_in_consensus_model.R* and *align_parameters_step_by_step.R*. The Script *parameter_influence_in_consensus_model.R* runs the estimations using EpiEstim varying parameters and data between the choices of the research groups. The Script *align_parameters_step_by_step.R* runs the estimations using the different estimation methods with different levels of parameter adjustment. Some estimations have to be performed outside of this script, either in Python (globalrt, rtlive) or in a more time-intensive script (epiforecasts). Those estimates are contained in this repo and do not have to be calculated again. For all but rtlive and HZI, it is possible, though.
- **epiforecasts**: Run *epiforecasts/epinow_estimation.R*. Takes multiple hours depending on the resources and the length of input data.
- **globalrt**: Run *ArroyoMarioli/input_output_dataset/estimate_R_filter_and_smoother.py* (Might want to change the import between *dataset.csv* and *dataset_final.csv*. The former is the data used originally for the globalrt estimates. The latter is the correctly formatted aggregated version of the RKI line list data.)
- **rtlive**: The data on the number of tests performed is missing and not publicly available. Thus, the estimation cannot be performed on the basis of the repository.

### How to do all from scratch without using prepared data sets from repo? ###
- **ETH**: Before loading ETH data run the script *ETH/otherScripts/format_linelist_data.R* to calculate the delays. Then load the data with *new_deconvolution = TRUE* in the function *load_incidence_data()* with *method = 'ETHZ_sliding_window'*
- **globalrt**: Run the scripts from *ArroyoMarioli/input_output/dataset/* in the following order: *format_data_from_various_sources.R* -> *construct_dataset_adj.py* -> *construct_STAN_models.py* -> *estimate_R_filter_and_smoother.py*
- **rtlive**: The incidence data used for the main comparison results from running the script *rtlive/important_files_from_repo/notebooks/data_aggregation.ipynb*.
- **SDSC**: The smoothed data can be constructed using the script *SDSC/smoothing_example.ipynb*.

Note that for the execution the RKI's line list data might be necessary. These are to large too be saved to this repository. Thus, they have to be sourced from https://hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0/explore. Older versions of these data can be obtained through filtering by "Meldedatum".

### Requirements ###
- *KITmetricslab/reproductive_numbers* repository must be stored in the same directory as this repository. It contains the realtime estimates as published by the resepective teams until early 2022.
