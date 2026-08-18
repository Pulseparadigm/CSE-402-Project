# Implementation & Upgrades: Weerakoon & Fernando (2000)
## "A Variant of Newton's Method with Accelerated Third-Order Convergence"
*Applied Mathematics Letters, 13(8), pp. 87–93*

---

## 🚀 Implemented Methods Overview

This repository implements both the original Weerakoon–Fernando method (VNM) and three powerful extensions:

| Method | Order ($p$) | Derivative Requirement | Formula / Update Rule |
|--------|------------|------------------------|-----------------------|
| **Newton–Raphson** | 2 | $f'(x)$ | $x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}$ |
| **Weerakoon–Fernando (VNM)** | **3** | $f'(x)$ | $x_{n+1} = x_n - \frac{2 f(x_n)}{f'(x_n) + f'(y_n)}$ |
| **Simpson–VNM Upgrade** | **4** | $f'(x)$ | $x_{n+1} = x_n - \frac{6 f(x_n)}{f'(x_n) + 4f'(z_n) + f'(y_n)}$ where $z_n = \frac{x_n + y_n}{2}$ |
| **Steffensen–VNM Upgrade** | **3** | **Derivative-Free** | Replaces $f'(x)$ with forward finite difference ratios |
| **Multivariate VNM System** | **3** | Jacobian $J(X)$ | $(J(X_n) + J(Y_n)) \Delta X = 2 F(X_n)$, $X_{n+1} = X_n - \Delta X$ |

---

## 📊 Benchmark Results

### Scalar Test Cases
All methods were benchmarked across 6 non-linear functions (including those from **Table 1** of the paper).

1. **Newton (Order 2):** Converges in ~5-11 iterations.
2. **VNM (Order 3):** Converges in ~3-6 iterations.
3. **Simpson–VNM (Order 4):** Achieves Order 4 convergence with even faster error reduction per iteration.
4. **Steffensen–VNM (Derivative-Free):** Achieves Order 3 convergence without needing analytical derivatives.

### Multivariate System Benchmark
Solves the non-linear 2D system $F(x, y) = [x^2 + y^2 - 4, xy - 1]^T = 0$:
- **Initial Guess:** $X_0 = [1.5, 0.5]^T$
- **Root Solved:** $X^* = [1.93185165, 0.51763809]^T$
- **Convergence:** Achieved high precision in **3 iterations** using `scipy.linalg.solve`.

---

## 📈 Visualizations

### 1. Error Convergence Rates (Log Scale)
![Convergence Comparison](plots/convergence_comparison.png)

### 2. Efficiency Comparison (Iterations & NFEV)
![Efficiency Summary](plots/efficiency_summary.png)

---

## 📁 Repository Structure

```text
CSE-402-Project/
├── README.md                           # Repository documentation
├── walkthrough.md                      # Detailed extension report & benchmarks
├── implementation_plans.md             # Implementation specifications
├── .gitignore                          # Excluded files configuration
├── src/
│   └── weerakoon_fernando.py           # Core implementation of all 5 methods
├── paper/
│   ├── A Variant of Newton’s Method...pdf  # Original research paper
│   ├── extract_pdf.py                  # Text extraction script
│   └── paper_text.txt                  # Extracted paper text
└── plots/
    ├── convergence_comparison.png      # Error trajectory plots
    └── efficiency_summary.png          # Iteration and NFEV comparison
```
