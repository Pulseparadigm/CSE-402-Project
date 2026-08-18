"""
=============================================================================
Upgraded Weerakoon-Fernando Method (VNM) & Extensions
-----------------------------------------------------------------------------
Implementation of the Weerakoon-Fernando method (2000) and its modern
upgrades:
1. Classical Newton-Raphson (Order 2)
2. Weerakoon-Fernando VNM (Trapezoidal, Order 3)
3. Simpson-type VNM (Simpson 1/3 Rule, Order 4)
4. Steffensen-type VNM (Derivative-Free, Order 3)
5. Multivariate System VNM (Systems of Nonlinear Equations F(X) = 0)

Authors of Paper: S. Weerakoon and T. G. I. Fernando
Refined Implementation: Python / NumPy / SciPy
=============================================================================
"""

import os
import math
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import warnings

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# 1. CORE & UPGRADED METHODS
# ─────────────────────────────────────────────────────────────────────────────

def newton_method(f, df, x0, tol=1e-15, max_iter=100):
    """Classical Newton-Raphson method (Quadratic Convergence, Order 2)."""
    x = float(x0)
    history = [x]
    errors = []
    nfev = 0
    for i in range(1, max_iter + 1):
        fx = f(x)
        dfx = df(x)
        nfev += 2
        
        if abs(dfx) < 1e-14:
            raise ValueError("Derivative near zero.")
            
        x_next = x - fx / dfx
        history.append(x_next)
        err = abs(x_next - x)
        errors.append(err)
        
        if err < tol or abs(f(x_next)) < tol:
            return x_next, history, errors, nfev, i
        x = x_next
        
    return x, history, errors, nfev, max_iter


def vnm(f, df, x0, tol=1e-15, max_iter=100):
    """Base Weerakoon-Fernando Method (Trapezoidal Rule, Order 3)."""
    x = float(x0)
    history = [x]
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
        history.append(x_next)
        err = abs(x_next - x)
        errors.append(err)
        
        if err < tol or abs(f(x_next)) < tol:
            return x_next, history, errors, nfev, i
        x = x_next
        
    return x, history, errors, nfev, max_iter


def simpson_vnm(f, df, x0, tol=1e-15, max_iter=100):
    """
    Simpson-type Quadrature Upgrade (Order 4 Convergence).
    Uses midpoint z_n = (x_n + y_n) / 2 to evaluate Simpson's 1/3 Rule.
    """
    x = float(x0)
    history = [x]
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
        history.append(x_next)
        err = abs(x_next - x)
        errors.append(err)
        
        if err < tol or abs(f(x_next)) < tol:
            return x_next, history, errors, nfev, i
        x = x_next
        
    return x, history, errors, nfev, max_iter


def steffensen_vnm(f, x0, tol=1e-15, max_iter=100):
    """
    Derivative-Free Weerakoon-Fernando Variant (Order 3 Convergence).
    Replaces d/dx with forward finite differences.
    """
    x = float(x0)
    history = [x]
    errors = []
    nfev = 0
    for i in range(1, max_iter + 1):
        fx = f(x)
        nfev += 1
        
        if abs(fx) < tol:
            return x, history, errors, nfev, i - 1
            
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
        history.append(x_next)
        err = abs(x_next - x)
        errors.append(err)
        
        if err < tol or abs(f(x_next)) < tol:
            return x_next, history, errors, nfev, i
        x = x_next
        
    return x, history, errors, nfev, max_iter


def vnm_system(F, J, X0, tol=1e-12, max_iter=100):
    """
    Multivariate Extension of VNM for Systems of Nonlinear Equations F(X) = 0.
    Uses scipy.linalg.solve instead of explicit matrix inversion.
    """
    X = np.array(X0, dtype=float)
    history = [X.copy()]
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
        history.append(X_next.copy())
        err = la.norm(X_next - X, ord=2)
        errors.append(err)
        
        if err < tol or la.norm(F(X_next), ord=2) < tol:
            return X_next, history, errors, nfev, i
        X = X_next
        
    return X, history, errors, nfev, max_iter


def compute_coc(errors):
    """Calculates Computational Order of Convergence (p) from consecutive errors."""
    if len(errors) < 3:
        return float("nan")
    e_k  = errors[-1]
    e_k1 = errors[-2]
    e_k2 = errors[-3]
    if e_k1 == 0 or e_k2 == 0 or e_k == 0:
        return float("nan")
    try:
        return np.log(abs(e_k / e_k1)) / np.log(abs(e_k1 / e_k2))
    except Exception:
        return float("nan")


# ─────────────────────────────────────────────────────────────────────────────
# 2. TEST FUNCTIONS & MULTIVARIATE SYSTEM BENCHMARKS
# ─────────────────────────────────────────────────────────────────────────────

TEST_CASES = [
    {
        "name"   : r"$f(x) = x^3 + 4x^2 - 10$",
        "label"  : "f1: x³+4x²−10",
        "f"      : lambda x: x**3 + 4*x**2 - 10,
        "df"     : lambda x: 3*x**2 + 8*x,
        "x0"     : -0.5,
        "root"   : 1.36523001341448,
    },
    {
        "name"   : r"$f(x) = \sin^2(x) - x^2 + 1$",
        "label"  : "f2: sin²(x)−x²+1",
        "f"      : lambda x: math.sin(x)**2 - x**2 + 1,
        "df"     : lambda x: 2*math.sin(x)*math.cos(x) - 2*x,
        "x0"     : 1.0,
        "root"   : 1.40449164821621,
    },
    {
        "name"   : r"$f(x) = e^x - 3x$",
        "label"  : "f3: eˣ−3x",
        "f"      : lambda x: math.exp(x) - 3*x,
        "df"     : lambda x: math.exp(x) - 3,
        "x0"     : 0.5,
        "root"   : 0.6190612867359450,
    },
    {
        "name"   : r"$f(x) = \cos(x) - x$",
        "label"  : "f4: cos(x)−x",
        "f"      : lambda x: math.cos(x) - x,
        "df"     : lambda x: -math.sin(x) - 1,
        "x0"     : 0.5,
        "root"   : 0.7390851332151607,
    },
    {
        "name"   : r"$f(x) = x^3 - 2x - 5$",
        "label"  : "f5: x³−2x−5",
        "f"      : lambda x: x**3 - 2*x - 5,
        "df"     : lambda x: 3*x**2 - 2,
        "x0"     : 2.0,
        "root"   : 2.09455148154233,
    },
    {
        "name"   : r"$f(x) = \ln(x) + x - 2$",
        "label"  : "f6: ln(x)+x−2",
        "f"      : lambda x: math.log(x) + x - 2,
        "df"     : lambda x: 1.0/x + 1,
        "x0"     : 1.5,
        "root"   : 1.55714852908987,
    },
]

# 2D System:
# F1(x, y) = x^2 + y^2 - 4 = 0
# F2(x, y) = x*y - 1 = 0
def system_F(X):
    x, y = X[0], X[1]
    return np.array([x**2 + y**2 - 4.0, x*y - 1.0])

def system_J(X):
    x, y = X[0], X[1]
    return np.array([
        [2.0*x, 2.0*y],
        [y,     x]
    ])


# ─────────────────────────────────────────────────────────────────────────────
# 3. BENCHMARK RUNNER
# ─────────────────────────────────────────────────────────────────────────────

def run_all_benchmarks():
    print("=" * 100)
    print(f"{'Function':<25} {'Method':<14} {'Iter':>5} {'NFEV':>6}  {'COC':>7}  {'Root'}")
    print("=" * 100)

    results = []
    methods = [
        ("Newton (Order 2)", lambda f, df, x0: newton_method(f, df, x0)),
        ("VNM (Order 3)",    lambda f, df, x0: vnm(f, df, x0)),
        ("Simpson (Order 4)",lambda f, df, x0: simpson_vnm(f, df, x0)),
        ("Steffensen (No-df)",lambda f, df, x0: steffensen_vnm(f, x0)),
    ]

    for tc in TEST_CASES:
        f, df, x0 = tc["f"], tc["df"], tc["x0"]
        tc_res = {"tc": tc, "runs": {}}

        for m_name, m_func in methods:
            try:
                root, history, errors, nfev, ni = m_func(f, df, x0)
                coc = compute_coc(errors)
                tc_res["runs"][m_name] = (root, history, errors, nfev, ni, coc)
                print(f"{tc['label']:<25} {m_name:<14} {ni:>5} {nfev:>6}  {coc:>7.4f}  {root:.15f}")
            except Exception as e:
                print(f"{tc['label']:<25} {m_name:<14} ERROR: {e}")

        print("-" * 100)
        results.append(tc_res)

    # System benchmark
    print("\n" + "=" * 100)
    print("  MULTIVARIATE SYSTEM VNM BENCHMARK: F(x, y) = [x² + y² - 4, xy - 1]^T")
    print("=" * 100)
    X0 = [1.5, 0.5]
    root_sys, hist_sys, errs_sys, nfev_sys, ni_sys = vnm_system(system_F, system_J, X0)
    coc_sys = compute_coc(errs_sys)
    print(f"  Initial Guess X0: {X0}")
    print(f"  Solved Root X*  : [{root_sys[0]:.15f}, {root_sys[1]:.15f}]")
    print(f"  Iterations      : {ni_sys},  NFEV: {nfev_sys},  COC: {coc_sys:.4f}")
    print("=" * 100)

    sys_res = {
        "X0": X0, "root": root_sys, "history": hist_sys, "errors": errs_sys,
        "nfev": nfev_sys, "ni": ni_sys, "coc": coc_sys
    }

    return results, sys_res


# ─────────────────────────────────────────────────────────────────────────────
# 4. PLOTTING & VISUALIZATION
# ─────────────────────────────────────────────────────────────────────────────

COLORS = {
    "nm_line"   : "#4FC3F7",   # light blue
    "vnm_line"  : "#FF8A65",   # orange
    "simp_line" : "#A5D6A7",   # green
    "steff_line": "#CE93D8",  # purple
    "bg"        : "#0D1117",
    "panel"     : "#161B22",
    "grid"      : "#21262D",
    "text"      : "#E6EDF3",
    "sub"       : "#8B949E",
    "accent"    : "#58A6FF",
}

plt.rcParams.update({
    "figure.facecolor"  : COLORS["bg"],
    "axes.facecolor"    : COLORS["panel"],
    "axes.edgecolor"    : COLORS["grid"],
    "axes.labelcolor"   : COLORS["text"],
    "axes.titlecolor"   : COLORS["text"],
    "xtick.color"       : COLORS["sub"],
    "ytick.color"       : COLORS["sub"],
    "grid.color"        : COLORS["grid"],
    "text.color"        : COLORS["text"],
    "legend.facecolor"  : COLORS["panel"],
    "legend.edgecolor"  : COLORS["grid"],
    "font.family"       : "DejaVu Sans",
    "font.size"         : 10,
})

def get_plot_dir():
    plot_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "plots"))
    os.makedirs(plot_dir, exist_ok=True)
    return plot_dir


def generate_plots(results, sys_res):
    plot_dir = get_plot_dir()
    n = len(results)
    cols = 3
    rows = math.ceil(n / cols)

    color_map = {
        "Newton (Order 2)": COLORS["nm_line"],
        "VNM (Order 3)": COLORS["vnm_line"],
        "Simpson (Order 4)": COLORS["simp_line"],
        "Steffensen (No-df)": COLORS["steff_line"],
    }

    # 1. Convergence Comparison
    fig, axes = plt.subplots(rows, cols, figsize=(18, rows * 4))
    fig.suptitle(
        "Convergence Comparison Across All Newton Variants\n"
        "NM (Ord 2)  ·  VNM (Ord 3)  ·  Simpson-VNM (Ord 4)  ·  Steffensen-VNM (Ord 3)",
        fontsize=14, fontweight="bold", color=COLORS["text"], y=1.01
    )

    for idx, res in enumerate(results):
        ax = axes.flat[idx]
        tc = res["tc"]

        for m_name, run_data in res["runs"].items():
            errors = run_data[2]
            xs = list(range(1, len(errors) + 1))
            safe_errs = [max(e, 1e-17) for e in errors]
            ax.semilogy(xs, safe_errs, "o-", color=color_map[m_name],
                        label=f"{m_name} ({run_data[4]} it)", linewidth=2, markersize=4)

        ax.set_title(tc["name"], fontsize=11, pad=8)
        ax.set_xlabel("Iteration", fontsize=9)
        ax.set_ylabel("Absolute Error |xₙ₊₁ - xₙ|", fontsize=9)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(fontsize=7, loc="upper right")

    for idx in range(len(results), rows * cols):
        axes.flat[idx].set_visible(False)

    plt.tight_layout()
    out1 = os.path.join(plot_dir, "convergence_comparison.png")
    plt.savefig(out1, dpi=150, bbox_inches="tight", facecolor=COLORS["bg"])
    print(f"[✓] Saved: {out1}")
    plt.close()

    # 2. Efficiency Bar Chart
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))
    fig.patch.set_facecolor(COLORS["bg"])
    labels = [r["tc"]["label"] for r in results]
    x = np.arange(len(labels))
    w = 0.2

    m_keys = ["Newton (Order 2)", "VNM (Order 3)", "Simpson (Order 4)", "Steffensen (No-df)"]
    offsets = [-1.5*w, -0.5*w, 0.5*w, 1.5*w]

    for ax, metric_idx, ylabel, title in [
        (ax1, 4, "Iterations to Converge", "Iteration Count Comparison"),
        (ax2, 3, "Total Function Evaluations (NFEV)", "Function Evaluations (NFEV)"),
    ]:
        ax.set_facecolor(COLORS["panel"])
        for k_idx, m_name in enumerate(m_keys):
            vals = [r["runs"][m_name][metric_idx] for r in results]
            bars = ax.bar(x + offsets[k_idx], vals, w, label=m_name,
                          color=color_map[m_name], alpha=0.85)
            ax.bar_label(bars, padding=2, fontsize=7, color=COLORS["text"])

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontweight="bold")
        ax.legend(fontsize=8)
        ax.grid(axis="y", linestyle="--", alpha=0.4)

    fig.suptitle("Efficiency Summary Across Upgraded Variants", fontsize=13, fontweight="bold")
    plt.tight_layout()
    out2 = os.path.join(plot_dir, "efficiency_summary.png")
    plt.savefig(out2, dpi=150, bbox_inches="tight", facecolor=COLORS["bg"])
    print(f"[✓] Saved: {out2}")
    plt.close()

    # 3. Computational Order of Convergence (COC) per Iteration
    fig, axes = plt.subplots(rows, cols, figsize=(18, rows * 4))
    fig.suptitle(
        "Computational Order of Convergence (COC) per Iteration\n"
        "NM ≈ 2  ·  VNM ≈ 3  ·  Simpson-VNM ≈ 4  ·  Steffensen-VNM ≈ 3",
        fontsize=14, fontweight="bold", color=COLORS["text"], y=1.01
    )

    for idx, res in enumerate(results):
        ax = axes.flat[idx]
        tc = res["tc"]

        for m_name, run_data in res["runs"].items():
            errors = run_data[2]
            coc_vals = []
            for k in range(2, len(errors)):
                e2, e1, e0 = errors[k-2], errors[k-1], errors[k]
                if e2 > 0 and e1 > 0 and e0 > 0:
                    try:
                        p = math.log(e0 / e1) / math.log(e1 / e2)
                        coc_vals.append(p)
                    except Exception:
                        pass
            if coc_vals:
                ax.plot(range(1, len(coc_vals) + 1), coc_vals, "o-",
                        color=color_map[m_name], label=f"{m_name}", linewidth=2, markersize=4)

        ax.axhline(2, color=COLORS["nm_line"],  linestyle="--", alpha=0.4, label="Order 2")
        ax.axhline(3, color=COLORS["vnm_line"], linestyle="--", alpha=0.4, label="Order 3")
        ax.axhline(4, color=COLORS["simp_line"],linestyle="--", alpha=0.4, label="Order 4")

        ax.set_title(tc["name"], fontsize=11, pad=8)
        ax.set_xlabel("Iteration", fontsize=9)
        ax.set_ylabel("COC (p)", fontsize=9)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_ylim(0, 5.5)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(fontsize=7, loc="upper right")

    for idx in range(len(results), rows * cols):
        axes.flat[idx].set_visible(False)

    plt.tight_layout()
    out3 = os.path.join(plot_dir, "coc_comparison.png")
    plt.savefig(out3, dpi=150, bbox_inches="tight", facecolor=COLORS["bg"])
    print(f"[✓] Saved: {out3}")
    plt.close()

    # 4. Multivariate System Trajectory Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor(COLORS["bg"])
    ax.set_facecolor(COLORS["panel"])

    x_vals = np.linspace(0.5, 2.5, 300)
    y_vals = np.linspace(0.1, 2.0, 300)
    X_grid, Y_grid = np.meshgrid(x_vals, y_vals)

    Z1 = X_grid**2 + Y_grid**2 - 4.0
    Z2 = X_grid * Y_grid - 1.0

    ax.contour(X_grid, Y_grid, Z1, levels=[0], colors=COLORS["nm_line"], linewidths=2)
    ax.contour(X_grid, Y_grid, Z2, levels=[0], colors=COLORS["vnm_line"], linewidths=2)

    hist_arr = np.array(sys_res["history"])
    ax.plot(hist_arr[:, 0], hist_arr[:, 1], "ro-", linewidth=2.5, markersize=8,
            label="Multivariate VNM Trajectory", zorder=5)

    ax.scatter([sys_res["X0"][0]], [sys_res["X0"][1]], color=COLORS["steff_line"], s=100, zorder=6, label="Start X0")
    ax.scatter([sys_res["root"][0]], [sys_res["root"][1]], color="#56d364", s=120, marker="*", zorder=6, label=f"Root X*")

    ax.set_title("Multivariate VNM Trajectory on F(x, y) = [x²+y²-4, xy-1]^T", fontweight="bold")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(loc="upper right")
    ax.grid(True, linestyle="--", alpha=0.4)

    plt.tight_layout()
    out4 = os.path.join(plot_dir, "multivariate_trajectory.png")
    plt.savefig(out4, dpi=150, bbox_inches="tight", facecolor=COLORS["bg"])
    print(f"[✓] Saved: {out4}")
    plt.close()


def main():
    print("\n" + "=" * 100)
    print("  UPGRADED WEERAKOON-FERNANDO METHOD (VNM) & EXTENSIONS BENCHMARK SUITE")
    print("  Weerakoon & Fernando (2000) with Order 4 Simpson & Steffensen Upgrades")
    print("=" * 100 + "\n")

    results, sys_res = run_all_benchmarks()
    print("\n[→] Generating updated performance plots...")
    generate_plots(results, sys_res)

    print("\n" + "=" * 100)
    print("  All benchmark suites completed and 4 plots saved in /plots directory!")
    print("  1. convergence_comparison.png  - Semi-log error trajectories")
    print("  2. efficiency_summary.png      - Iterations & NFEV bar charts")
    print("  3. coc_comparison.png          - Computational Order of Convergence")
    print("  4. multivariate_trajectory.png - 2D System trajectory phase plane")
    print("=" * 100 + "\n")


if __name__ == "__main__":
    main()
