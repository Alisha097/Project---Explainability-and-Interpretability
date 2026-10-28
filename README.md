# Evaluating and Validating Causal Treatment Effect Estimation Methods

How can we determine whether a machine learning model's estimates of causal effects are reliable? This project investigates the evaluation and validation of causal treatment effect estimators, with a focus on whether feature attribution methods correctly identify the covariates responsible for heterogeneous treatment effects.

**MSc Data Science Dissertation | 2025–2026**

## Overview

This project investigates the reliability of feature attribution methods for explaining heterogeneous treatment effect (HTE) predictions produced by tree-based causal models.

Although causal machine learning models can estimate conditional average treatment effects (CATE), accurate treatment effect estimation does not necessarily mean that a model correctly identifies the covariates responsible for treatment effect heterogeneity.

This dissertation extends the predictive-versus-prognostic attribution benchmarking framework proposed by [Crabbe et al.](https://arxiv.org/abs/2206.08363) to two tree-based causal estimators: **Causal Forest (CF)** and **Bayesian Causal Forest (BCF)**.

Using a semi-synthetic experimental framework based on real covariates from the **Linked Births and Infant Deaths Database (LBIDD)**, the project evaluates explanation fidelity against known ground truth while independently assessing treatment effect estimation accuracy.

## Research Question

How reliable are feature attribution methods for explaining heterogeneous treatment effect estimates produced by tree-based causal estimators, and under what conditions does this reliability break down?

### Research Objectives

- Evaluate whether feature attribution methods recover the true covariates responsible for treatment effect heterogeneity.
- Investigate how explanation fidelity changes with predictive signal strength.
- Examine the effects of predictive and prognostic confounding on explanation reliability.
- Compare built-in feature importance and SHAP for Causal Forest.
- Evaluate Bayesian Causal Forest using permutation feature importance.
- Assess the relationship between explanation fidelity and treatment effect estimation accuracy.

## Methodology

### Dataset and Experimental Framework

The study uses real covariates from LBIDD to construct a semi-synthetic dataset. Treatment assignments and potential outcomes are generated synthetically, allowing the true predictive and prognostic covariates to be identified.

The experimental design varies:

- **Predictive signal strength:** Weak, moderate, and strong treatment effect heterogeneity.
- **Confounding type:** No confounding, prognostic confounding, and predictive confounding.
- **Confounding strength:** Additional experiments investigate the effects of increasing confounding strength.

This framework allows explanation methods to be evaluated against known ground truth under controlled experimental conditions.

### Models and Explanation Methods

| Model | Implementation | Explanation method |
|---|---|---|
| Causal Forest (CF) | Python, EconML | Built-in feature importance and SHAP |
| Bayesian Causal Forest (BCF) | R, `bcf` | Permutation feature importance |

### Evaluation Metrics

**Precision in Estimation of Heterogeneous Effect (PEHE)**

Measures the discrepancy between estimated and true individual-level treatment effects. Lower PEHE indicates more accurate treatment effect estimation.

**Predictive Attribution (Attrpred)**

Measures the proportion of total feature attribution assigned to known predictive covariates that determine treatment effect heterogeneity.

**Prognostic Attribution (Attrprog)**

Measures the proportion of total feature attribution assigned to prognostic covariates that affect outcomes but do not directly determine treatment effect heterogeneity.

Together, these metrics assess treatment effect estimation accuracy and explanation fidelity as distinct dimensions of model performance.

## Key Findings

The experiments produced the following main findings:

- Explanation fidelity generally improves as predictive signal strength increases.
- Predictive confounding substantially reduces explanation reliability, even when treatment effect estimation remains relatively stable.
- Causal Forest's built-in feature importance generally identifies the true predictive covariates more accurately than SHAP across the evaluated experimental settings.
- Bayesian Causal Forest exhibits similar overall patterns and greater robustness in several confounding settings.
- Accurate treatment effect estimation does not necessarily imply trustworthy explanations of treatment effect heterogeneity.

These findings highlight the importance of evaluating explanation fidelity independently of conventional estimation accuracy.

## Technologies

- **Python** — Causal Forest implementation and data preprocessing.
- **EconML** — Causal Forest estimation using `CausalForestDML`.
- **SHAP** — Post-hoc feature attribution.
- **R** — Bayesian Causal Forest implementation.
- **bcf** — Bayesian Causal Forest estimation.
- **ranger** — Propensity score estimation.
- **ggplot2** — Statistical visualisation.
- **CSV** — Sharing the processed covariate matrix between Python and R.

## Repository Structure

```text
dissertation-project/
├── data/               # Processed covariates or data access instructions
├── notebooks/          # Exploratory analysis, if included
├── src/                # Python source code
├── R/                  # Bayesian Causal Forest implementation
├── results/             # Experimental outputs and figures
├── requirements.txt    # Python dependencies
└── README.md            # Project documentation
```


## Reproducibility

The experiments use a shared, cleaned covariate matrix for the Python and R implementations. Experimental runs vary the selected predictive and prognostic covariates, treatment assignments, outcome noise, and train-test splits across random seeds.

For reproducibility, the repository should document the required software versions, dependencies, data preparation steps, and experiment-specific configurations.

Where data cannot be redistributed, include instructions explaining how to obtain or prepare the required covariates.


## References

- Crabbe et al. — [Predictive-versus-prognostic attribution benchmarking framework (arXiv)](https://arxiv.org/abs/2206.08363)
- [EconML documentation](https://www.pywhy.org/EconML/)
- [Bayesian Causal Forest (`bcf`) package](https://cran.r-project.org/package=bcf)
- [SHAP documentation](https://shap.readthedocs.io/)

