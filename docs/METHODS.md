# Methods — Detailed Descriptions

This document describes the statistical and machine learning methods implemented in the Genomic Causal Inference Pipeline.

---

## 1. Prediction Models

### 1.1 Ordinary Least Squares (OLS)

Standard linear regression minimizing the sum of squared residuals:

$$\hat{\beta} = (X^\top X)^{-1} X^\top Y$$

Used as a baseline for comparison with nonlinear models.

### 1.2 Random Forest (RF)

An ensemble of decision trees trained on bootstrap samples with random feature subsets. Predictions are averaged across all trees. Controlled by `rf_n_estimators` (default: 100).

### 1.3 K-Nearest Neighbors (KNN)

Predicts the outcome as the average of the K nearest training samples in feature space (Euclidean distance). Controlled by `knn_n_neighbors` (default: 5).

---

## 2. Two-Stage Instrumental Variable Models

These methods address **endogeneity** — the concern that the treatment SNP may be correlated with unobserved confounders.

### Causal Structure

```
Instruments (Z) ──→ Treatment (T) ──→ Outcome (Y)
                                        ↑
                          Confounders (U, observed & unobserved)
```

### Assignment

| Role | Variable |
|---|---|
| Treatment (T) | First (most-correlated) SNP: `X[:, 0]` |
| Instruments (Z) | Remaining top SNPs: `X[:, 1:]` |
| Confounders | Sex, Age |
| Outcome (Y) | Height |

### 2.1 Two-Stage Least Squares (2SLS)

- **Stage 1:** Regress T on Z using OLS → obtain predicted treatment T̂
- **Stage 2:** Regress Y on T̂ using OLS → obtain causal estimate

Both stages are fit on **training data only**; predictions are made on held-out test data.

### 2.2 Two-Stage Random Forest (2SRF)

Same structure as 2SLS but uses Random Forest in both stages.

### 2.3 Two-Stage KNN (2SKNN)

Same structure as 2SLS but uses KNN in both stages.

---

## 3. Mendelian Randomization (MR)

MR uses genetic variants as **instruments** to estimate causal effects, exploiting the fact that genotypes are assigned at conception (natural randomization).

### 3.1 Wald Ratio

For a single instrument Z:

$$\hat{\beta}_{causal} = \frac{\hat{\beta}_{ZY}}{\hat{\beta}_{ZX}}$$

where β_ZY is the instrument–outcome association and β_ZX is the instrument–treatment association. We select the strongest instrument (highest |β_ZX|).

Standard error via the delta method:

$$SE = |\hat{\beta}_{causal}| \sqrt{\left(\frac{SE_{ZY}}{\hat{\beta}_{ZY}}\right)^2 + \left(\frac{SE_{ZX}}{\hat{\beta}_{ZX}}\right)^2}$$

### 3.2 Inverse-Variance Weighted (IVW)

Fixed-effect meta-analysis across multiple instruments:

$$\hat{\beta}_{IVW} = \frac{\sum_j w_j \cdot (\hat{\beta}_{ZY,j} / \hat{\beta}_{ZX,j})}{\sum_j w_j}, \quad w_j = \frac{\hat{\beta}_{ZX,j}^2}{SE_{ZY,j}^2}$$

Assumes **all instruments are valid** (no horizontal pleiotropy).

### 3.3 MR-Egger

Allows for **directional pleiotropy** by including an intercept in the weighted regression of β_ZY on β_ZX:

$$\hat{\beta}_{ZY,j} = \alpha + \beta_{causal} \cdot \hat{\beta}_{ZX,j} + \epsilon_j$$

A non-zero intercept (α ≠ 0) indicates the presence of directional pleiotropy. The InSIDE assumption (Instrument Strength Independent of Direct Effect) must hold.

### 3.4 Weighted Median

Robust if up to **50% of instruments are invalid**. Computes per-instrument Wald ratios, weights them by precision, and takes the weighted median. Bootstrap standard errors (1000 iterations).

### Confounder Adjustment

All MR methods **residualise** the treatment and outcome on observed confounders (Sex, Age) before computing instrument–exposure and instrument–outcome associations.

---

## 4. Double Machine Learning (DML)

Implements the partially-linear model from **Chernozhukov et al. (2018)**:

$$Y = \theta \cdot T + g(W) + \varepsilon$$
$$T = m(W) + \eta$$

where T is the treatment SNP, W = [remaining SNPs, Sex, Age], and g(·), m(·) are unknown nuisance functions.

### Algorithm (Cross-Fitting)

1. Split data into K folds.
2. For each fold k:
   - Train nuisance models ĝ(W) and m̂(W) on all data **except** fold k.
   - Compute residuals on fold k:
     - Ỹ_k = Y_k − ĝ(W_k)
     - T̃_k = T_k − m̂(W_k)
3. Estimate θ by regressing pooled Ỹ on pooled T̃:

$$\hat{\theta} = \frac{\sum_i \tilde{T}_i \tilde{Y}_i}{\sum_i \tilde{T}_i^2}$$

4. Heteroskedasticity-robust (HC0) standard error:

$$\widehat{Var}(\hat{\theta}) = \frac{\sum_i \tilde{T}_i^2 \hat{\varepsilon}_i^2}{\left(\sum_i \tilde{T}_i^2\right)^2}$$

where ε̂_i = Ỹ_i − θ̂ · T̃_i.

Nuisance models use **Random Forest** regressors.

---

## 5. DoWhy Causal Model

Uses Microsoft's [DoWhy](https://github.com/py-why/dowhy) library for:

1. **Model specification** — Define treatment, outcome, instruments, and confounders as a causal graph.
2. **Identification** — Automatically find an estimand (IV or backdoor).
3. **Estimation** — Compute the causal effect.
4. **Refutation** — Test robustness via:
   - **Placebo treatment** — Replace the treatment with random noise. Estimate should → 0.
   - **Random common cause** — Add a random confounder. Estimate should remain stable.
   - **Data subset** — Re-estimate on a random subset. Estimate should remain stable.

---

## References

- Bowden, J. et al. (2015). Mendelian randomization with invalid instruments. *Statistical Methods in Medical Research*.
- Burgess, S. et al. (2013). Mendelian randomization analysis with multiple genetic variants. *International Journal of Epidemiology*.
- Chernozhukov, V. et al. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*.
- Sharma, A. & Kiciman, E. (2020). DoWhy: An end-to-end library for causal inference. *arXiv:2011.04216*.
