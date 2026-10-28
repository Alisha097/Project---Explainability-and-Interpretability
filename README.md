Evaluating and Validating Causal Treatment Effect Estimation Method (Variant 5)
MSc Data Science Dissertation
King's College London

Student: Alisha
Student Number: K25118600
Supervisor: Dr. Jeroen Berrevoets

CONTENTS


This folder contains the source code used to reproduce the experiments
presented in the report.

Files included:

1. causal_forest.ipynb
   Python implementation of the complete Causal Forest experimental
   pipeline. The notebook performs:
   
   • Loading and preprocessing the cleaned LBIDD covariates
   • Semi-synthetic data generation following the benchmark framework
   • Causal Forest model training
   • SHAP explanation generation
   • Built-in feature importance extraction
   • Calculation of Attrpred and Attrprog
   • Calculation of PEHE and Normalised PEHE
   • Export of experimental results

2. bcf.R
   R implementation of the Bayesian Causal Forest (BCF) experiments.
   This script trains the Bayesian Causal Forest model, estimates
   heterogeneous treatment effects, computes permutation feature
   importance and evaluates explanation fidelity and estimation
   performance.

DATASET


The experiments are based on the Linked Births and Infant Deaths
Database (LBIDD).

The original dataset was obtained from the IBM Causal Inference
Benchmarking Framework repository:

https://github.com/IBM-HRL-MLHLS/IBM-Causal-Inference-Benchmarking-Framework/tree/master/data/LBIDD

The original LBIDD data were cleaned prior to experimentation by
removing redundant, near-constant and highly correlated variables, as
described in Chapter 4 of the report. The resulting cleaned
dataset was then used to generate the semi-synthetic benchmark datasets
used throughout the experiments.

SOFTWARE REQUIREMENTS


Python:
• Python 3.12.0
• NumPy 2.4.6
• pandas 3.0.3
• SciPy 1.17.1
• scikit-learn 1.6.1
• EconML 0.16.0
• SHAP 0.48.0
• Matplotlib 3.11.0

R:
• R 4.5.1
• bcf 2.0.2
• ranger 0.18.0
• ggplot2 4.0.3


OUTPUTS

Running the provided scripts generates the experimental results reported in the report, including:

Causal Forest
• causal_forest_results.csv
• causal_forest_summary.csv
• confounding_strength_sweep.csv
• confounding_strength_summary.csv

Bayesian Causal Forest
• bcf_estimation_results.csv
• bcf_estimation_summary_robust.csv
• bcf_explanation_results.csv
• bcf_explanation_summary.csv

NOTE:
The Python notebook should be executed before the R script, as it
performs data preparation and generates the semi-synthetic benchmark
datasets required by the Bayesian Causal Forest experiments

