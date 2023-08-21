# Symptoms that differentiate patients with long covid from those with PVS

### Adith Arun, Fred Warner, Chenxi Huang, Harlan Krumholz

To identify the most important symptoms for differentiating participants experiencing Long Covid compared to those experiencing PVS, we trained a prediction model using information on presence or absence of all of the symptoms for each participants. We implemented a gradient boosted tree machine learning model for predicting whether a participant experienced internal tremor or not with 5-fold 5-repeat cross-validation. We selected the hyperparameters with highest area under the curve (AUC) from the internal cross-validation. Then, we computed the importance of each variable in differentiating participants using a permutation-based approach. We sorted the variables based on their importance and progressively excluded those with least importance from the model and evaluated the change in the AUC. We selected the best model and corresponding number of variables when the AUC first decreased by at least 1.5%. If a drop of this magnitude did not occur, we selected the model with the largest drop in AUC. 
To assess the robustness of this process to choice of modeling methods and variable importance metrics, this process was repeated with an XGBoost model with the gain in accuracy metric used to assess variable importance, as well as a XGBoost model with the Shapley value used to assess variable importance. We compare each of the three methods’ results for variable importance values using Pearson correlation coefficients. 

## How to use
Clone this github repository and navigate to the lc-vi-ml directory. Then, in your terminal, type the command `Rscript predict-lc-vi.R lc-vi/merged.csv LC lc-vi FALSE` to produce the output generated in the subfolder `lc-vi/`. The patient level binary symptom matrix can be found at `lc-vi/merged.csv`. 

If you have the `lc-vi/results.rds` object created and would like to use the outputs from that object, then, run `Rscript predict-lc-vi.R lc-vi/merged.csv LC lc-vi TRUE` to create the plots. 

The variable `LC` is the class 1 variable name used for computing AUC and other binary classification metrics
The variable `lc-vi` is the name of the folder to store outputs in

The variable `FALSE` performs the ML analysis from scratch. If set to `TRUE` then the script reads in `lc-vi/results.rds` object with pre-computed results that are plotted to form outputs in `lc-vi/`. 







