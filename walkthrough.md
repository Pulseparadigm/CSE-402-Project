# Comparative Report: Weerakoon–Fernando Method & Modern Upgrades

**Paper Title:** *A Variant of Newton's Method with Accelerated Third-Order Convergence*  
**Authors:** S. Weerakoon and T. G. I. Fernando (*Applied Mathematics Letters*, 13(8), pp. 87–93, 2000)  
**Implementation & Benchmark Suite:** [`src/weerakoon_fernando.py`](file:///f:/4-1/402/project/src/weerakoon_fernando.py)

---

## 📋 Executive Summary of Comparison

This benchmark suite evaluates the **Weerakoon–Fernando Method (VNM)** against the classical **Newton–Raphson** method and two high-level extensions:
1. **Simpson–VNM Upgrade:** Uses Simpson’s 1/3 quadrature rule to boost convergence to **Order 4 (Quartic)**.
2. **Steffensen–VNM Upgrade:** Replaces analytical derivatives with finite difference approximations for **Derivative-Free Order 3** convergence.
3. **Multivariate VNM System:** Extends VNM to solve non-linear 2D systems $F(X) = 0$.

---

## 📈 Method Comparison & Performance Summary

### 1. Theoretical vs. Empirical Metrics

| Method | Order ($p$) | Derivative Needed? | Formula | Iterations to $10^{-15}$ | NFEV Efficiency |
|--------|------------|--------------------|---------|--------------------------|-----------------|
| **Newton–Raphson** | $2.0$ | $f'(x)$ | $x - \frac{f(x)}{f'(x)}$ | 5 – 10 | Baseline (2 eval/iter) |
| **Weerakoon–Fernando (VNM)** | **$3.0$** | $f'(x)$ | $x - \frac{2f(x)}{f'(x) + f'(y_n)}$ | 3 – 5 | **High** (3 eval/iter) |
| **Simpson–VNM Upgrade** | **$4.0$** | $f'(x)$ | $x - \frac{6f(x)}{f'(x) + 4f'(z_n) + f'(y_n)}$ | **2 – 4** | **Highest** (4 eval/iter) |
| **Steffensen–VNM Upgrade** | **$3.0$** | ❌ **No** | Finite difference derivative approx. | 3 – 5 | High (3-4 eval/iter) |
| **Multivariate VNM System** | **$3.0$** | Jacobian $J(X)$ | $(J(X) + J(Y))\Delta X = 2F(X)$ | **3** | High (2 Jacobian evaluations) |

---

## 📊 Comparison Visualizations

### 1. Error Trajectories (Log Scale)
Demonstrates the rate of error decay per iteration. Higher-order methods (VNM, Simpson-VNM) drop to machine precision ($10^{-15}$) in significantly fewer iterations.

![Convergence Trajectories](plots/convergence_comparison.png)

---

### 2. Efficiency Comparison (Iterations & Function Evaluations)
Compares total iterations and total function evaluations (NFEVs) needed to reach $10^{-15}$ tolerance across all 6 test functions.

![Efficiency Summary](plots/efficiency_summary.png)

---

### 3. Computational Order of Convergence (COC)
Plots the empirical order of convergence $p \approx \frac{\ln|e_k/e_{k-1}|}{\ln|e_{k-1}/e_{k-2}|}$ across iterations.

![Computational Order of Convergence](plots/coc_comparison.png)

---

### 4. Multivariate System Trajectory Phase Plane
Visualizes the convergence trajectory of Multivariate VNM on the non-linear system $F(x,y) = [x^2+y^2-4, xy-1]^T = 0$ starting from $X_0 = [1.5, 0.5]^T$.

![Multivariate System Trajectory](plots/multivariate_trajectory.png)

---

## 📁 Repository Organization

```text
CSE-402-Project/
├── README.md                           # Main repository documentation
├── walkthrough.md                      # Complete comparative report
├── implementation_plans.md             # Specification & code details
├── .gitignore                          # Excluded files configuration
├── src/
│   └── weerakoon_fernando.py           # Core module & benchmark generator
├── paper/
│   ├── A Variant of Newton’s Method...pdf  # Original research paper PDF
│   ├── extract_pdf.py                  # PDF extraction utility
│   └── paper_text.txt                  # Extracted paper text
└── plots/
    ├── convergence_comparison.png      # Error log-scale plots
    ├── efficiency_summary.png          # Iteration & NFEV bar charts
    ├── coc_comparison.png              # Order of convergence plots
    ├── multivariate_trajectory.png     # 2D system phase plane
    └── root_paths.png                  # Function trajectory plots
```
