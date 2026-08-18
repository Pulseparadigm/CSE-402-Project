# Implementation: Weerakoon & Fernando (2000)
## "A Variant of Newton's Method with Accelerated Third-Order Convergence"
*Applied Mathematics Letters, 13(8), pp. 87–93*

---

## The Method

The **Weerakoon-Fernando VNM** solves `f(x) = 0` with cubic convergence — without needing `f''(x)`.

Each iteration has two steps:

| Step | Formula | Role |
|------|---------|------|
| Predictor | `y_n = x_n - f(x_n) / f'(x_n)` | Standard Newton step |
| Corrector | `x_{n+1} = x_n - 2·f(x_n) / (f'(x_n) + f'(y_n))` | Weerakoon-Fernando update |

**Why it works:** The corrector uses the trapezoidal rule to approximate the integral in Newton's method derivation, which upgrades convergence from order 2 → order **3**.

---

## Implementation: `src/weerakoon_fernando.py`

### Functions
- `newton_method(f, df, x0)` — Classical Newton-Raphson (order 2)
- `vnm(f, df, x0)` — Weerakoon-Fernando VNM (order 3)
- `computational_order_of_convergence(errors)` — Estimates COC from error sequence

### Test Functions (Table 1 from paper + extras)

| # | Function | x₀ | True Root |
|---|----------|----|-----------|
| 1 | x³ + 4x² − 10 | −0.5 | 1.36523001341448 |
| 2 | sin²(x) − x² + 1 | 1.0 | 1.40449164821621 |
| 3 | eˣ − 3x | 0.5 | 0.61906128673594 |
| 4 | cos(x) − x | 0.5 | 0.73908513321516 |
| 5 | x³ − 2x − 5 | 2.0 | 2.09455148154233 |
| 6 | ln(x) + x − 2 | 1.5 | 1.55714852908987 |

---

## Results & Visualizations

### Convergence Rate (Log Scale)
![Convergence Comparison](plots/convergence_comparison.png)

### Efficiency Comparison (Iterations & NFEV)
![Efficiency Summary](plots/efficiency_summary.png)

### Root-Finding Iteration Paths
![Root Finding Paths](plots/root_paths.png)

### Computational Order of Convergence (COC)
![Computational Order of Convergence](plots/coc_comparison.png)

---

## Key Takeaways

> [!IMPORTANT]
> VNM needs **1 extra function evaluation per iteration** vs Newton (evaluates `f'` at both `x_n` and `y_n`), but converges in **fewer total iterations**, so **total NFEVs are lower** overall.

- **NM COC ≈ 2.0** (quadratic) — confirmed by plots
- **VNM COC ≈ 3.0** (cubic) — confirmed by plots
- VNM requires **no second derivatives** — unlike other third-order methods (e.g., Halley's method)

---

## Repository Structure

```
CSE-402-Project/
├── README.md                           # Main repository overview
├── walkthrough.md                      # Detailed technical report
├── .gitignore                          # Excluded files configuration
├── src/
│   └── weerakoon_fernando.py           # Core implementation of NM & VNM
├── paper/
│   ├── A Variant of Newton’s Method...pdf  # Original research paper
│   ├── extract_pdf.py                  # Text extraction script
│   └── paper_text.txt                  # Extracted paper text
└── plots/
    ├── convergence_comparison.png      # Iteration error plots
    ├── efficiency_summary.png          # Bar chart comparison
    ├── root_paths.png                  # Function trajectory plots
    └── coc_comparison.png              # Convergence order analysis
```
