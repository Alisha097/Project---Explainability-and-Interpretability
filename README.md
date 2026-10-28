# Evaluating and Validating Causal Treatment Effect Estimation Methods
How can we know whether a machine learning models estimate of a causal effect is actually correct? This project explores how to evaluate and validate causal treatment effect estimators such as those used in personalised medicine, economics and policy design. 

**MSc Data Science Dissertation | 2025–2026**

# Explainability and Interpretability - How can we explain heterogeneous treatment effect predictions?

## Overview
This project investigates the reliability of feature attribution methods for explaining heterogeneous treatment effect (HTE) predictions produced by tree-based causal models.

Although causal machine learning models can estimate conditional average treatment effects (CATE), accurate treatment effect estimation does not necessarily mean that the model correctly identifies the covariates responsible for treatment effect heterogeneity.

This dissertation extends the predictive-versus-prognostic attribution benchmarking framework proposed by **[Crabbe et al.](https://arxiv.org/abs/2206.08363)** to two tree-based causal estimators: **Causal Forest (CF)** and **Bayesian Causal Forest (BCF)**.
**Reference:** [Read the paper on arXiv](https://arxiv.org/abs/2206.08363) · [Download PDF](https://arxiv.org/pdf/2206.08363)

Using a semi-synthetic experimental framework based on real covariates from the Linked Births and Infant Deaths Database (LBIDD), the project evaluates explanation fidelity against known ground truth while independently assessing treatment effect estimation accuracy.

## Research Question

How reliable are feature attribution methods for explaining heterogeneous treatment effect estimates produced by tree-based causal estimators and under what conditions does this reliability break down?

### Research Objectives

- Evaluate whether feature attribution methods recover the true covariates responsible for treatment effect heterogeneity.
- Investigate how explanation fidelity changes with predictive signal strength.
- Examine the effects of predictive and prognostic confounding on explanation reliability.
- Compare explanation methods for Causal Forest and evaluate Bayesian Causal Forest using permutation feature importance.
- Assess the relationship between explanation fidelity and treatment effect estimation accuracy.

## Methodology

### Dataset and Experimental Framework

The study uses real covariates from LBIDD to construct a semi-synthetic dataset. Treatment assignments and potential outcomes are generated synthetically, allowing the true predictive and prognostic covariates to be identified.

The experimental design varies:

- **Predictive signal strength:** Weak, moderate and strong treatment effect heterogeneity.
- **Confounding type:** No confounding, prognostic confounding and predictive confounding.
- **Confounding strength:** Additional experiments investigate the effect of increasing confounding strength.

### Models and Explanation Methods

| Model | Implementation | Explanation method |
|---|---|---|
| Causal Forest | Python, EconML | Built-in feature importance and SHAP |
| Bayesian Causal Forest | R, `bcf` | Permutation feature importance |

### Evaluation Metrics

**Precision in Estimation of Heterogeneous Effect (PEHE)** measures the discrepancy between estimated and true treatment effects.

**Predictive attribution (Attrpred)** measures the proportion of total feature attribution assigned to the known predictive covariates.

**Prognostic attribution (Attrprog)** measures the proportion of total attribution assigned to prognostic covariates that affect outcomes but do not directly determine treatment effect heterogeneity.

Together, these metrics assess estimation accuracy and explanation fidelity as distinct dimensions of model performance.

## Key Findings

The experiments produced the following main findings:

- Explanation fidelity generally improves as predictive signal strength increases.
- Predictive confounding substantially reduces explanation reliability, even when treatment effect estimation remains relatively stable.
- Causal Forest's built-in feature importance generally identifies the true predictive covariates more accurately than SHAP across the evaluated experimental settings.
- Bayesian Causal Forest exhibits similar overall patterns and greater robustness in several confounding settings.
- Accurate treatment effect estimation does not necessarily imply trustworthy explanations of treatment effect heterogeneity.

These findings highlight the importance of evaluating explanation fidelity independently of conventional estimation accuracy.

## Technologies

- **Python:** Causal Forest implementation and data preprocessing.
- **EconML:** Causal Forest estimation using `CausalForestDML`.
- **SHAP:** Post-hoc feature attribution.
- **R:** Bayesian Causal Forest implementation.
- **bcf:** Bayesian Causal Forest estimation.
- **ranger:** Propensity score estimation.
- **ggplot2:** Statistical visualisation.
- **CSV:** Sharing the processed covariate matrix between Python and R.

## Repository Structure

```text
dissertation-project/
├── data/              # Processed covariates or data access instructions
├── notebooks/         # Exploratory analysis, if applicable
├── src/               # Python source code
├── R/                 # Bayesian Causal Forest implementation, if applicable
├── results/            # Experimental outputs and figures
├── requirements.txt   # Python dependencies, if applicable
└── README.md           # Project documentation
```

## Reproducibility

The experiments use a shared, cleaned covariate matrix for the Python and R implementations. Experimental runs vary the selected predictive and prognostic covariates, treatment assignments, outcome noise and train-test splits across random seeds.

To reproduce the experiments, consult the project-specific installation instructions, software versions and experimental configuration files.

## Dissertation

The full dissertation and supplementary materials can be linked here when they are available for public distribution.

- **Full dissertation:** [Add link]
- **Experimental results:** [Add link]
- **Supplementary materials:** [Add link]


## References

Please consult the dissertation's reference list for the complete academic sources underpinning the methodology and experimental framework.



