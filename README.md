# Model Variable Augmentation (MVA) for Diagnostic Assessment of Sensitivity Analysis Results
*by Juliane Mai and Bryan A Tolson (University of Waterloo, Canada)*

## Abstract
The method of Model Variable Augmentation (MVA) was introduced to assess the quality of SA results 
without performing any additional model runs or requiring bootstrapping. MVA is proven to perform 
well when only a small number of model runs was used to obtain the sensitivity indexes. 
MVA augments the original model input variables with additional variables of known properties. 
The sensitivities of the augmented model variables are used to draw conclusions on the reliability 
of the other "original" model parameters' sensitivities. The MVA method is already successfully 
tested with two global SA methods: the variance-based Sobol' method and the moment-independent PAWN method. The full paper can be found [here](https://agupubs.onlinelibrary.wiley.com/journal/19447973).

## Step-by-Step Tutorial
The step-by-step tutorial describes all the steps to estimate sensitivity indexes for (original) model variables and the augmented parameters. It also explains how to analyse these results and how to draw conclusions on the reliablility of the sensitivity indexes of the original model variables. Details can be found [here](https://github.com/julemai/MVA/wiki/Step-by-Step-Tutorial).

## Examples
We provide some case studies to show how MVA can help:
- to check the implementation of the sensitivity analysis method (see [here](https://github.com/julemai/MVA/wiki/Examples#sensitivity-analysis-method-implementation-check))
- to obtain a robust ranking of the model variables (see [here](https://github.com/julemai/MVA/wiki/Examples#check-for-convergence-of-input-variable-importance-ranking))
- to estimate the uncertainty of the sensitivity indexes without the necessity of bootstrapping (see [here](https://github.com/julemai/MVA/wiki/Examples#check-for-convergence-of-sensitivity-indexes))

## Citation
J Mai and BA Tolson (2019): <br>
Model Variable Augmentation (MVA) for Diagnostic Assessment of Sensitivity Analysis Results
Water Resources Research, ??, ??–??.
