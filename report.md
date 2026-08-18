# Detailed Performance Comparison Report

## Method Performance & Efficiency Summary

| Method | Order ($p$) | Analytical $f'$ Needed? | Iterations to $10^{-15}$ | NFEV Efficiency |
|--------|------------|-------------------------|--------------------------|-----------------|
| **Newton–Raphson** | 2 | Yes | 5 – 10 | Baseline (2 eval/iter) |
| **Weerakoon–Fernando (VNM)** | **3** | Yes | 3 – 5 | High (3 eval/iter) |
| **Simpson–VNM Upgrade** | **4** | Yes | **2 – 4** | **Highest** (4 eval/iter) |
| **Steffensen–VNM Upgrade** | **3** | ❌ **No** | 3 – 5 | High (3-4 eval/iter) |
| **Multivariate VNM System** | **3** | Jacobian $J(X)$ | **3** | High (2 Jacobians/iter) |

---

## Performance Analysis & Insights

1. **Newton–Raphson (Order 2):**
   - Serves as the classic baseline. Requires evaluating $f(x)$ and $f'(x)$ at each step.
   - Takes 5 to 10 iterations to reach $10^{-15}$ accuracy.

2. **Weerakoon–Fernando Method (VNM, Order 3):**
   - Approximates the derivative integral using the **Trapezoidal Rule**.
   - Achieves 3rd-order cubic convergence without needing second-order derivatives ($f''(x)$).
   - Reduces the number of required iterations to 3 – 5.

3. **Simpson–VNM Upgrade (Order 4):**
   - Enhances the quadrature approximation using **Simpson's 1/3 Rule** with a midpoint step $z_n = \frac{x_n + y_n}{2}$.
   - Achieves 4th-order quartic convergence, requiring only **2 to 4 iterations** to reach machine precision ($10^{-15}$).

4. **Steffensen–VNM Upgrade (Derivative-Free Order 3):**
   - Replaces analytical derivatives $f'(x)$ with forward finite-difference ratios $\frac{f(x+f(x)) - f(x)}{f(x)}$.
   - Retains 3rd-order convergence while eliminating the need for symbolic/analytical derivatives.

5. **Multivariate VNM System (Order 3):**
   - Extends VNM to $n$-dimensional systems of equations $F(X) = 0$ using $(J(X_n) + J(Y_n)) \Delta X = 2 F(X_n)$.
   - Solves systems in **3 iterations** using efficient linear system solvers (`scipy.linalg.solve`).

---

## Generated Comparison Graphs

All figures are available in the [`plots/`](./plots/) folder:
- **`plots/convergence_comparison.png`**: Error vs. Iteration Trajectories (Log Scale)
- **`plots/efficiency_summary.png`**: Total Iterations & Function Evaluations (NFEVs)
- **`plots/coc_comparison.png`**: Empirical Order of Convergence ($p$)
- **`plots/multivariate_trajectory.png`**: Phase Plane Trajectory for 2D Systems
