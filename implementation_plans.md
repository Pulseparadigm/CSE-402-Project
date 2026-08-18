# Upgraded Weerakoon–Fernando Method: Consolidated Suite & Extension Report

This file contains the complete implementation plan, source code, benchmark suite, and extensions for the Weerakoon–Fernando method (VNM), including higher-order (Simpson-type), derivative-free (Steffensen-type), and multivariate systems of nonlinear equations.

---

## 1. Upgraded Core Module (`src/weerakoon_fernando.py`)

```python
import numpy as np
import scipy.linalg as la

# =====================================================================
# 1. Base Methods & Standard Metrics
# =====================================================================

def newton_method(f, df, x0, tol=1e-15, max_iter=100):
    """Classical Newton-Raphson method (Quadratic Convergence, Order 2)."""
    x = float(x0)
    errors = []
    nfev = 0
    for i in range(1, max_iter + 1):
        fx = f(x)
        dfx = df(x)
        nfev += 2
        
        if abs(dfx) < 1e-14:
            raise ValueError("Derivative near zero.")
            
        x_next = x - fx / dfx
        err = abs(x_next - x)
        errors.append(err)
        
        if err < tol or abs(f(x_next)) < tol:
            return x_next, i, nfev, errors
        x = x_next
        
    raise TimeoutError("Newton-Raphson failed to converge.")

def vnm(f, df, x0, tol=1e-15, max_iter=100):
    """Base Weerakoon-Fernando Method (Trapezoidal Rule, Order 3)."""
    x = float(x0)
    errors = []
    nfev = 0
    for i in range(1, max_iter + 1):
        fx = f(x)
        dfx = df(x)
        nfev += 2
        
        if abs(dfx) < 1e-14:
            raise ValueError("Derivative near zero.")
            
        yn = x - fx / dfx
        df_yn = df(yn)
        nfev += 1
        
        denom = dfx + df_yn
        if abs(denom) < 1e-14:
            raise ValueError("Denominator near zero.")
            
        x_next = x - (2.0 * fx) / denom
        err = abs(x_next - x)
        errors.append(err)
        
        if err < tol or abs(f(x_next)) < tol:
            return x_next, i, nfev, errors
        x = x_next
        
    raise TimeoutError("VNM failed to converge.")

# =====================================================================
# 2. Upgraded Variants
# =====================================================================

def simpson_vnm(f, df, x0, tol=1e-15, max_iter=100):
    """
    Simpson-type Quadrature Upgrade (Order 4 Convergence).
    Uses midpoint z_n = (x_n + y_n) / 2 to evaluate Simpson's 1/3 Rule.
    """
    x = float(x0)
    errors = []
    nfev = 0
    for i in range(1, max_iter + 1):
        fx = f(x)
        dfx = df(x)
        nfev += 2
        
        if abs(dfx) < 1e-14:
            raise ValueError("Derivative near zero.")
            
        yn = x - fx / dfx
        zn = 0.5 * (x + yn)
        
        df_zn = df(zn)
        df_yn = df(yn)
        nfev += 2
        
        denom = dfx + 4.0 * df_zn + df_yn
        if abs(denom) < 1e-14:
            raise ValueError("Denominator near zero.")
            
        x_next = x - (6.0 * fx) / denom
        err = abs(x_next - x)
        errors.append(err)
        
        if err < tol or abs(f(x_next)) < tol:
            return x_next, i, nfev, errors
        x = x_next
        
    raise TimeoutError("Simpson-VNM failed to converge.")

def steffensen_vnm(f, x0, tol=1e-15, max_iter=100):
    """
    Derivative-Free Weerakoon-Fernando Variant (Order 3 Convergence).
    Replaces d/dx with forward finite differences.
    """
    x = float(x0)
    errors = []
    nfev = 0
    for i in range(1, max_iter + 1):
        fx = f(x)
        nfev += 1
        
        if abs(fx) < tol:
            return x, i - 1, nfev, errors
            
        fx_plus = f(x + fx)
        nfev += 1
        dfx_approx = (fx_plus - fx) / fx
        
        if abs(dfx_approx) < 1e-14:
            raise ValueError("Approximated derivative near zero.")
            
        yn = x - fx / dfx_approx
        f_yn = f(yn)
        nfev += 1
        
        df_yn_approx = (f(yn + f_yn) - f_yn) / f_yn if abs(f_yn) > 1e-15 else dfx_approx
        nfev += (1 if abs(f_yn) > 1e-15 else 0)
        
        denom = dfx_approx + df_yn_approx
        if abs(denom) < 1e-14:
            raise ValueError("Denominator near zero.")
            
        x_next = x - (2.0 * fx) / denom
        err = abs(x_next - x)
        errors.append(err)
        
        if err < tol or abs(f(x_next)) < tol:
            return x_next, i, nfev, errors
        x = x_next
        
    raise TimeoutError("Steffensen-VNM failed to converge.")

def vnm_system(F, J, X0, tol=1e-12, max_iter=100):
    """
    Multivariate Extension of VNM for Systems of Nonlinear Equations F(X) = 0.
    Uses scipy.linalg.solve instead of explicit matrix inversion.
    """
    X = np.array(X0, dtype=float)
    errors = []
    nfev = 0
    
    for i in range(1, max_iter + 1):
        FX = F(X)
        JX = J(X)
        nfev += 1
        
        delta_y = la.solve(JX, FX)
        Y = X - delta_y
        
        JY = J(Y)
        nfev += 1
        
        J_sum = JX + JY
        delta_x = la.solve(J_sum, 2.0 * FX)
        
        X_next = X - delta_x
        err = la.norm(X_next - X, ord=2)
        errors.append(err)
        
        if err < tol or la.norm(F(X_next), ord=2) < tol:
            return X_next, i, nfev, errors
        X = X_next
        
    raise TimeoutError("Multivariate VNM failed to converge.")

def compute_coc(errors):
    """Calculates Computational Order of Convergence (p) from consecutive errors."""
    if len(errors) < 3:
        return np.nan
    e_k = errors[-1]
    e_k1 = errors[-2]
    e_k2 = errors[-3]
    if e_k1 == 0 or e_k2 == 0 or e_k == 0:
        return np.nan
    return np.log(abs(e_k / e_k1)) / np.log(abs(e_k1 / e_k2))
```