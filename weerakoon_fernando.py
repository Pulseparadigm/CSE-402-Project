"""
=============================================================================
A Variant of Newton's Method with Accelerated Third-Order Convergence
-----------------------------------------------------------------------------
Implementation of the Weerakoon-Fernando method (2000), Published in
Applied Mathematics Letters, 13(8), pp. 87-93.

Authors of Paper: S. Weerakoon and T. G. I. Fernando
Implementation:   Python

Method (VNM):
  Given x_n, compute:
    y_n     = x_n - f(x_n) / f'(x_n)         [Newton predictor step]
    x_{n+1} = x_n - 2*f(x_n) / (f'(x_n) + f'(y_n))   [corrector step]

This achieves third-order (cubic) convergence without needing f''(x).
=============================================================================
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
import warnings

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
# CORE METHODS
# ─────────────────────────────────────────────────────────────────────────────

def newton_method(f, df, x0, tol=1e-15, max_iter=200):
    """
    Classical Newton-Raphson method (second-order convergence).
    x_{n+1} = x_n - f(x_n) / f'(x_n)
    """
    x = x0
    history = [x]
    errors  = []
    nfev    = 0   # number of function evaluations

    for i in range(max_iter):
        fx  = f(x);  nfev += 1
        dfx = df(x); nfev += 1

        if abs(dfx) < 1e-300:
            return None, history, errors, nfev, i + 1  # derivative vanished

        x_new = x - fx / dfx
        history.append(x_new)
        errors.append(abs(x_new - x))

        if abs(x_new - x) < tol and abs(fx) < tol:
            return x_new, history, errors, nfev, i + 1

        x = x_new

    return x, history, errors, nfev, max_iter


def vnm(f, df, x0, tol=1e-15, max_iter=200):
    """
    Weerakoon-Fernando Variant of Newton's Method (VNM) — third-order convergence.

    Predictor:  y_n     = x_n - f(x_n) / f'(x_n)
    Corrector:  x_{n+1} = x_n - 2*f(x_n) / (f'(x_n) + f'(y_n))
    """
    x = x0
    history = [x]
    errors  = []
    nfev    = 0

    for i in range(max_iter):
        fx  = f(x);  nfev += 1
        dfx = df(x); nfev += 1

        if abs(dfx) < 1e-300:
            return None, history, errors, nfev, i + 1

        # Predictor step (Newton)
        y    = x - fx / dfx

        # Corrector step (VNM)
        dfy  = df(y);  nfev += 1
        denom = dfx + dfy
        if abs(denom) < 1e-300:
            return None, history, errors, nfev, i + 1

        x_new = x - 2.0 * fx / denom
        history.append(x_new)
        errors.append(abs(x_new - x))

        if abs(x_new - x) < tol and abs(f(x_new)) < tol:
            nfev += 1
            return x_new, history, errors, nfev, i + 1

        x = x_new

    return x, history, errors, nfev, max_iter


def computational_order_of_convergence(errors):
    """
    Estimate the computational order of convergence (COC):
      COC ≈ log(|e_{n+1}| / |e_n|) / log(|e_n| / |e_{n-1}|)
    """
    coc_vals = []
    for k in range(2, len(errors)):
        e_prev2 = errors[k - 2]
        e_prev1 = errors[k - 1]
        e_cur   = errors[k]
        if e_prev2 > 0 and e_prev1 > 0 and e_cur > 0:
            num   = math.log(abs(e_cur   / e_prev1) + 1e-300)
            denom = math.log(abs(e_prev1 / e_prev2) + 1e-300)
            if abs(denom) > 1e-10:
                coc_vals.append(num / denom)
    return coc_vals[-1] if coc_vals else float("nan")


# ─────────────────────────────────────────────────────────────────────────────
# TEST FUNCTIONS   (matching Table 1 in the paper)
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


# ─────────────────────────────────────────────────────────────────────────────
# RUN ALL TESTS
# ─────────────────────────────────────────────────────────────────────────────

def run_all_tests():
    results = []
    print("=" * 90)
    print(f"{'Function':<30} {'Method':<6} {'Iter':>5} {'NFEV':>6}  {'COC':>7}  {'Root'}")
    print("=" * 90)

    for tc in TEST_CASES:
        f, df, x0, true_root = tc["f"], tc["df"], tc["x0"], tc["root"]

        # Newton's Method
        r_nm, h_nm, e_nm, nfev_nm, ni_nm = newton_method(f, df, x0)
        coc_nm = computational_order_of_convergence(e_nm)

        # Variant of Newton's Method
        r_vn, h_vn, e_vn, nfev_vn, ni_vn = vnm(f, df, x0)
        coc_vn = computational_order_of_convergence(e_vn)

        print(f"{tc['label']:<30} {'NM':<6} {ni_nm:>5} {nfev_nm:>6}  {coc_nm:>7.4f}  {r_nm:.15f}")
        print(f"{'':<30} {'VNM':<6} {ni_vn:>5} {nfev_vn:>6}  {coc_vn:>7.4f}  {r_vn:.15f}")
        print("-" * 90)

        results.append({
            "tc"    : tc,
            "nm"    : (r_nm, h_nm, e_nm, nfev_nm, ni_nm, coc_nm),
            "vnm"   : (r_vn, h_vn, e_vn, nfev_vn, ni_vn, coc_vn),
        })

    return results


# ─────────────────────────────────────────────────────────────────────────────
# VISUALISATION
# ─────────────────────────────────────────────────────────────────────────────

COLORS = {
    "nm_line" : "#4FC3F7",   # light blue
    "vnm_line": "#FF8A65",   # orange
    "nm_mark" : "#0288D1",
    "vnm_mark": "#E64A19",
    "bg"      : "#0D1117",
    "panel"   : "#161B22",
    "grid"    : "#21262D",
    "text"    : "#E6EDF3",
    "sub"     : "#8B949E",
    "accent"  : "#58A6FF",
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


def plot_convergence(results):
    """Plot error vs iteration for each test function."""
    n = len(results)
    cols = 3
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(18, rows * 4))
    fig.suptitle(
        "Convergence Comparison: Newton's Method vs Weerakoon-Fernando VNM\n"
        "S. Weerakoon & T.G.I. Fernando — Applied Mathematics Letters (2000)",
        fontsize=14, fontweight="bold", color=COLORS["text"], y=1.01
    )

    for idx, res in enumerate(results):
        ax  = axes.flat[idx]
        tc  = res["tc"]
        e_nm  = res["nm"][2]
        e_vn  = res["vnm"][2]

        xs_nm = list(range(1, len(e_nm) + 1))
        xs_vn = list(range(1, len(e_vn) + 1))

        # Filter out zeros for log scale
        def safe_log(lst):
            return [max(e, 1e-17) for e in lst]

        ax.semilogy(xs_nm, safe_log(e_nm), "o-", color=COLORS["nm_line"],
                    label=f"Newton (NM) — {res['nm'][4]} iter",
                    linewidth=2, markersize=5, markerfacecolor=COLORS["nm_mark"])
        ax.semilogy(xs_vn, safe_log(e_vn), "s-", color=COLORS["vnm_line"],
                    label=f"VNM — {res['vnm'][4]} iter",
                    linewidth=2, markersize=5, markerfacecolor=COLORS["vnm_mark"])

        ax.set_title(tc["name"], fontsize=11, pad=8)
        ax.set_xlabel("Iteration", fontsize=9)
        ax.set_ylabel("Absolute Error  |xₙ₊₁ − xₙ|", fontsize=9)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(fontsize=8, loc="upper right")
        ax.set_facecolor(COLORS["panel"])

    # Hide unused axes
    for idx in range(len(results), rows * cols):
        axes.flat[idx].set_visible(False)

    plt.tight_layout()
    plt.savefig("convergence_comparison.png", dpi=150, bbox_inches="tight",
                facecolor=COLORS["bg"])
    print("\n[✓] Saved: convergence_comparison.png")
    plt.close()


def plot_summary_bar(results):
    """Bar chart: iterations and NFEV comparison."""
    labels    = [r["tc"]["label"] for r in results]
    nm_iters  = [r["nm"][4]  for r in results]
    vn_iters  = [r["vnm"][4] for r in results]
    nm_nfev   = [r["nm"][3]  for r in results]
    vn_nfev   = [r["vnm"][3] for r in results]

    x    = np.arange(len(labels))
    w    = 0.35
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))
    fig.patch.set_facecolor(COLORS["bg"])

    for ax, nm_vals, vn_vals, ylabel, title in [
        (ax1, nm_iters, vn_iters, "Iterations", "Number of Iterations to Converge"),
        (ax2, nm_nfev,  vn_nfev,  "Function Evaluations", "Total Function Evaluations (NFEV)"),
    ]:
        ax.set_facecolor(COLORS["panel"])
        b1 = ax.bar(x - w/2, nm_vals, w, label="Newton (NM)",
                    color=COLORS["nm_line"], alpha=0.85, edgecolor=COLORS["nm_mark"])
        b2 = ax.bar(x + w/2, vn_vals, w, label="VNM",
                    color=COLORS["vnm_line"], alpha=0.85, edgecolor=COLORS["vnm_mark"])

        ax.bar_label(b1, padding=3, fontsize=9, color=COLORS["text"])
        ax.bar_label(b2, padding=3, fontsize=9, color=COLORS["text"])

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontweight="bold")
        ax.legend()
        ax.grid(axis="y", linestyle="--", alpha=0.4)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "Efficiency Summary: Newton's Method vs Weerakoon-Fernando VNM",
        fontsize=13, fontweight="bold", color=COLORS["text"]
    )
    plt.tight_layout()
    plt.savefig("efficiency_summary.png", dpi=150, bbox_inches="tight",
                facecolor=COLORS["bg"])
    print("[✓] Saved: efficiency_summary.png")
    plt.close()


def plot_function_roots(results):
    """Visualise each function and show Newton vs VNM iteration paths."""
    n = len(results)
    cols = 3
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(18, rows * 4))
    fig.suptitle("Root-Finding Paths: Newton's Method vs VNM",
                 fontsize=14, fontweight="bold", color=COLORS["text"], y=1.01)

    for idx, res in enumerate(results):
        ax  = axes.flat[idx]
        tc  = res["tc"]
        f   = tc["f"]
        root = tc["root"]

        h_nm = res["nm"][1]
        h_vn = res["vnm"][1]

        # x range for plotting
        all_x = h_nm + h_vn
        xmin  = min(all_x)
        xmax  = max(all_x)
        margin = max(abs(xmax - xmin) * 0.3, 0.5)
        xlo, xhi = xmin - margin, xmax + margin

        xs = np.linspace(xlo, xhi, 600)
        try:
            ys = [f(xi) for xi in xs]
        except Exception:
            ys = [0] * len(xs)

        ax.plot(xs, ys, color=COLORS["accent"], linewidth=2, label=tc["name"], zorder=2)
        ax.axhline(0, color=COLORS["sub"], linewidth=0.8, linestyle="--")
        ax.axvline(root, color="#56d364", linewidth=1.2, linestyle=":", label=f"Root ≈ {root:.6f}")

        # Plot iteration x-points on the x-axis
        ys_nm = np.zeros(len(h_nm))
        ys_vn = np.zeros(len(h_vn))
        ax.scatter(h_nm, ys_nm, color=COLORS["nm_mark"], s=40, zorder=5, label="NM iterates")
        ax.scatter(h_vn, ys_vn, color=COLORS["vnm_mark"], s=40, marker="D", zorder=5, label="VNM iterates")

        ax.set_title(tc["name"], fontsize=11, pad=8)
        ax.set_xlabel("x", fontsize=9)
        ax.set_ylabel("f(x)", fontsize=9)
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.set_facecolor(COLORS["panel"])

        # Clip y-axis for readability
        valid_ys = [y for y in ys if not math.isnan(y) and not math.isinf(y)]
        if valid_ys:
            ylo = min(valid_ys)
            yhi = max(valid_ys)
            pad = (yhi - ylo) * 0.1 + 0.5
            ax.set_ylim(ylo - pad, yhi + pad)

    for idx in range(len(results), rows * cols):
        axes.flat[idx].set_visible(False)

    plt.tight_layout()
    plt.savefig("root_paths.png", dpi=150, bbox_inches="tight",
                facecolor=COLORS["bg"])
    print("[✓] Saved: root_paths.png")
    plt.close()


def plot_coc_comparison(results):
    """Plot computational order of convergence per iteration."""
    n = len(results)
    cols = 3
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(18, rows * 4))
    fig.suptitle(
        "Computational Order of Convergence (COC) per Iteration\n"
        "NM ≈ 2 (quadratic)  ·  VNM ≈ 3 (cubic)",
        fontsize=14, fontweight="bold", color=COLORS["text"], y=1.01
    )

    for idx, res in enumerate(results):
        ax   = axes.flat[idx]
        tc   = res["tc"]
        e_nm = res["nm"][2]
        e_vn = res["vnm"][2]

        def coc_series(errors):
            series = []
            for k in range(2, len(errors)):
                e2, e1, e0 = errors[k-2], errors[k-1], errors[k]
                if e2 > 0 and e1 > 0 and e0 > 0:
                    num   = math.log(max(e0 / e1, 1e-300))
                    denom = math.log(max(e1 / e2, 1e-300))
                    if abs(denom) > 1e-10:
                        series.append(num / denom)
            return series

        coc_nm = coc_series(e_nm)
        coc_vn = coc_series(e_vn)

        if coc_nm:
            ax.plot(range(1, len(coc_nm) + 1), coc_nm, "o-", color=COLORS["nm_line"],
                    label="NM COC", linewidth=2, markersize=5)
        if coc_vn:
            ax.plot(range(1, len(coc_vn) + 1), coc_vn, "s-", color=COLORS["vnm_line"],
                    label="VNM COC", linewidth=2, markersize=5)

        ax.axhline(2, color=COLORS["nm_mark"],  linestyle="--", alpha=0.5, label="Order 2")
        ax.axhline(3, color=COLORS["vnm_mark"], linestyle="--", alpha=0.5, label="Order 3")

        ax.set_title(tc["name"], fontsize=11, pad=8)
        ax.set_xlabel("Iteration", fontsize=9)
        ax.set_ylabel("COC", fontsize=9)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_ylim(0, 5)
        ax.legend(fontsize=8)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.set_facecolor(COLORS["panel"])

    for idx in range(len(results), rows * cols):
        axes.flat[idx].set_visible(False)

    plt.tight_layout()
    plt.savefig("coc_comparison.png", dpi=150, bbox_inches="tight",
                facecolor=COLORS["bg"])
    print("[✓] Saved: coc_comparison.png")
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# =============================================================================
# INTERACTIVE DEMO
# =============================================================================

def demo_single(f_str, f, df, x0, true_root=None):
    """Run both methods on a single function and print a detailed trace."""
    print(f"\n{'='*70}")
    print(f"  Function : {f_str}")
    print(f"  x0       : {x0}")
    if true_root is not None:
        print(f"  True root: {true_root:.15f}")
    print(f"{'='*70}")

    for name, method in [("Newton's Method (NM)", newton_method),
                          ("Variant of Newton's Method (VNM)", vnm)]:
        root, history, errors, nfev, ni = method(f, df, x0)
        print(f"\n  >> {name}")
        print(f"    {'n':>3}  {'x_n':>22}  {'|x_{n+1}-x_n|':>18}")
        print(f"    {'-'*3}  {'-'*22}  {'-'*18}")
        for k, (xk, ek) in enumerate(zip(history[1:], errors)):
            print(f"    {k+1:>3}  {xk:>22.15f}  {ek:>18.6e}")
        coc = computational_order_of_convergence(errors)
        print(f"\n    Root  ~ {root:.15f}")
        print(f"    Iterations : {ni},  NFEV : {nfev},  COC ~ {coc:.4f}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "=" * 90)
    print("  A Variant of Newton's Method with Accelerated Third-Order Convergence")
    print("  Weerakoon & Fernando (2000) - Applied Mathematics Letters 13(8):87-93")
    print("=" * 90)

    # ── 1. Run all paper test cases ──────────────────────────────────────────
    results = run_all_tests()

    # ── 2. Detailed trace for first two cases (like Table 1 in paper) ───────
    for tc in TEST_CASES[:2]:
        demo_single(tc["label"], tc["f"], tc["df"], tc["x0"], tc["root"])

    # ── 3. Generate plots ────────────────────────────────────────────────────
    print("\n[→] Generating plots...")
    plot_convergence(results)
    plot_summary_bar(results)
    plot_function_roots(results)
    plot_coc_comparison(results)

    print("\n" + "=" * 90)
    print("  All done! Four plots saved in the project directory.")
    print("  convergence_comparison.png  - Error vs Iterations (log scale)")
    print("  efficiency_summary.png      - Iterations & NFEV bar chart")
    print("  root_paths.png              - Function plots with iteration paths")
    print("  coc_comparison.png          - Computational Order of Convergence")
    print("=" * 90 + "\n")


if __name__ == "__main__":
    main()
