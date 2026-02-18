#!/usr/bin/env python3
"""
Two-Horizon Entropy Sharing: SPARC Radial Acceleration Relation

Reproduces all numerical results and figures from:
  "The Radial Acceleration Relation from Two-Horizon Entropy Sharing
   with Zero Free Parameters" — K. Brodie (2026)

The model: f(a) = a/(a + cH₀/6), where the factor 1/6 arises from the
cos²θ mode overlap between planar (Unruh) and spherical (Hubble)
entanglement structures, integrated over the backward hemisphere.

Requires: numpy, scipy, matplotlib
Data: SPARC RAR (McGaugh, Lelli & Schombert 2016, PRL 117, 201101)
      Bundled as data/RAR.mrt
"""

import numpy as np
from scipy import integrate, optimize
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================
# Physical constants
# =============================================
c = 2.998e8          # m/s
H0 = 67.4e3 / 3.086e22  # 67.4 km/s/Mpc in 1/s
cH0 = c * H0        # ≈ 6.548e-10 m/s²
a0_mond = 1.2e-10   # Milgrom's empirical value

print("=" * 70)
print("TWO-HORIZON ENTROPY SHARING — SPARC RAR ANALYSIS")
print("=" * 70)
print(f"\ncH₀ = {cH0:.4e} m/s²")
print(f"a₀_MOND = {a0_mond:.1e} m/s²")

# =============================================
# The geometric integral
# =============================================
print("\n--- Geometric Factor ---")

# cos²θ over backward hemisphere:
# (1/4π) ∫_{backward} cos²θ dΩ = (1/2) ∫_{-1}^{0} u² du = 1/6
geometric_factor = 1/6

# Numerical verification
def cos2_integrand(theta):
    return np.cos(theta)**2 * np.sin(theta)

back_num, _ = integrate.quad(cos2_integrand, np.pi/2, np.pi)
back_num *= 2 * np.pi / (4 * np.pi)
print(f"Analytical: 1/6 = {geometric_factor:.6f}")
print(f"Numerical:        {back_num:.6f}")

# The predicted a₀
a0_predicted = cH0 * geometric_factor  # = cH₀/6
print(f"\nPredicted a₀ = cH₀/6 = {a0_predicted:.4e} m/s²")
print(f"Milgrom's a₀ =         {a0_mond:.4e} m/s²")
print(f"Ratio:                  {a0_predicted/a0_mond:.3f}")
print(f"Discrepancy:            {abs(a0_predicted/a0_mond - 1)*100:.1f}%")

# =============================================
# Model functions
# =============================================

def f_this_work(a):
    """f(a) = a/(a + cH₀/6) — zero free parameters."""
    return a / (a + a0_predicted)

def f_bare(a):
    """f(a) = a/(a + cH₀) — no geometric factor."""
    return a / (a + cH0)

def mond_simple(g_bar, a0=1.2e-10):
    """Standard MOND simple interpolation (force law)."""
    return 0.5 * (g_bar + np.sqrt(g_bar**2 + 4 * g_bar * a0))

def solve_gobs(g_bar, f_func):
    """Solve g_bar = f(g_obs) × g_obs for g_obs given g_bar."""
    if g_bar > 1e-7:
        return g_bar
    def residual(log_gobs):
        gobs = 10**log_gobs
        f = f_func(gobs)
        if f < 1e-15:
            return 1e10
        return np.log10(g_bar) - np.log10(f * gobs)
    try:
        return 10**optimize.brentq(residual, np.log10(g_bar) - 1,
                                   np.log10(g_bar) + 6, xtol=1e-8)
    except Exception:
        return g_bar

# =============================================
# Load SPARC data
# =============================================
data_path = os.path.join(os.path.dirname(__file__), 'data', 'RAR.mrt')
if not os.path.exists(data_path):
    # Try alternate location
    data_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'RAR.mrt')

data = np.loadtxt(data_path, skiprows=13)
log_gbar = data[:, 0]
log_gobs = data[:, 2]
g_bar_data = 10**log_gbar
g_obs_data = 10**log_gobs
print(f"\nLoaded SPARC: {len(g_bar_data)} data points, 153 galaxies")

# =============================================
# Compute predictions
# =============================================
print("\nComputing predictions...")

g_obs_this_work = np.array([solve_gobs(gb, f_this_work) for gb in g_bar_data])
g_obs_mond = mond_simple(g_bar_data)
g_obs_bare = np.array([solve_gobs(gb, f_bare) for gb in g_bar_data])

resid_this_work = log_gobs - np.log10(np.maximum(g_obs_this_work, 1e-20))
resid_mond = log_gobs - np.log10(g_obs_mond)
resid_bare = log_gobs - np.log10(np.maximum(g_obs_bare, 1e-20))

# =============================================
# Table 1: Model comparison
# =============================================
print("\n" + "=" * 70)
print("TABLE 1: Model comparison against SPARC RAR")
print("=" * 70)
print(f"\n{'Model':>45s} {'params':>6s} {'mean':>8s} {'σ':>8s} {'a₀':>12s}")
print("-" * 85)
for label, resid, a0, npar in [
    ('This work: f=a/(a+cH₀/6)', resid_this_work, f'{a0_predicted:.2e}', '0'),
    ('MOND fitted: a₀=1.2e-10', resid_mond, '1.20e-10', '1'),
    ('Bare: f=a/(a+cH₀)', resid_bare, f'{cH0:.2e}', '0')]:
    m = np.mean(resid)
    s = np.std(resid)
    print(f"{label:>45s} {npar:>6s} {m:+8.4f} {s:8.4f} {a0:>12s}")

# =============================================
# Table 2: Regime analysis
# =============================================
print("\n" + "=" * 70)
print("TABLE 2: Regime analysis")
print("=" * 70)
a0_boundary = a0_mond  # use Milgrom's value for regime boundaries
regimes = [
    ('Deep MOND (g < 3a₀)', g_obs_data < 3 * a0_boundary),
    ('Transition (3a₀ < g < 30a₀)',
     (g_obs_data >= 3 * a0_boundary) & (g_obs_data < 30 * a0_boundary)),
    ('Newtonian (g > 30a₀)', g_obs_data >= 30 * a0_boundary),
]

print(f"\n{'Regime':>35s} {'N':>6s} {'σ_this':>8s} {'σ_MOND':>8s}")
print("-" * 65)
for name, mask in regimes:
    n = np.sum(mask)
    s_tw = np.std(resid_this_work[mask]) if n > 1 else 0
    s_m = np.std(resid_mond[mask]) if n > 1 else 0
    print(f"{name:>35s} {n:>6d} {s_tw:8.3f} {s_m:8.3f}")

# =============================================
# Figures
# =============================================
print("\nGenerating figures...")
outdir = os.path.dirname(__file__)

# --- Figure 1: f(a) curves ---
fig1, ax = plt.subplots(figsize=(8, 5))
a_range = np.logspace(-12, -7, 300)

ax.semilogx(a_range, a_range / (a_range + a0_predicted), 'b-', lw=2.5,
            label=f'This work: $a/(a+cH_0/6)$,  $a_0={a0_predicted:.2e}$')
ax.semilogx(a_range, a_range / (a_range + a0_mond), 'r--', lw=2,
            label=f'MOND fitted: $a/(a+a_0)$,  $a_0=1.2\\times10^{{-10}}$')
ax.semilogx(a_range, a_range / (a_range + cH0), 'gray', lw=1.5, ls=':',
            label=f'Bare: $a/(a+cH_0)$,  $a_0={cH0:.2e}$')
ax.axvline(x=a0_predicted, color='blue', ls=':', alpha=0.3)
ax.axvline(x=a0_mond, color='red', ls=':', alpha=0.3)
ax.set_xlabel('$a$ [m/s²]', fontsize=12)
ax.set_ylabel('$f(a) = m_i / m_g$', fontsize=12)
ax.set_title('Inertia modification function', fontsize=13)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(-0.05, 1.05)
fig1.tight_layout()
fig1.savefig(os.path.join(outdir, 'fig_f_curves.png'), dpi=150)
print("  Saved fig_f_curves.png")

# --- Figure 2: Mode overlap geometry ---
fig2, ax = plt.subplots(figsize=(8, 5))
theta = np.linspace(0, np.pi, 200)
ax.plot(np.degrees(theta), np.cos(theta)**2, 'b-', lw=2, label='cos²θ (mode overlap)')
ax.plot(np.degrees(theta), np.abs(np.cos(theta)), 'g--', lw=1.5, label='|cosθ|')
ax.axvline(x=90, color='k', ls=':', alpha=0.3, label='equator')
ax.axvspan(90, 180, alpha=0.1, color='blue', label='backward hemisphere')
ax.set_xlabel('θ (degrees from acceleration axis)', fontsize=12)
ax.set_ylabel('Competition strength', fontsize=12)
ax.set_title('Entanglement mode overlap: planar vs spherical', fontsize=13)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
fig2.tight_layout()
fig2.savefig(os.path.join(outdir, 'fig_geometry.png'), dpi=150)
print("  Saved fig_geometry.png")

# --- Figure 3: SPARC RAR ---
fig3, ax = plt.subplots(figsize=(8, 6))
ax.scatter(log_gbar, log_gobs, s=1, alpha=0.2, color='gray', label='SPARC data (2,693 points)')

g_range_plot = np.logspace(-13, -8, 200)
g_pred_tw = np.array([solve_gobs(gb, f_this_work) for gb in g_range_plot])
g_pred_mond = mond_simple(g_range_plot)

ax.plot(np.log10(g_range_plot), np.log10(np.maximum(g_pred_tw, 1e-20)),
        'b-', lw=2.5, label=f'This work, 0 params (σ={np.std(resid_this_work):.3f} dex)')
ax.plot(np.log10(g_range_plot), np.log10(g_pred_mond),
        'r--', lw=2, label=f'MOND fitted (σ={np.std(resid_mond):.3f} dex)')
ax.plot(np.log10(g_range_plot), np.log10(g_range_plot), 'k:', alpha=0.3,
        label='Newtonian ($g_{obs}=g_{bar}$)')
ax.set_xlabel('$\\log_{10}(g_{\\rm bar}$ / m s$^{-2})$', fontsize=12)
ax.set_ylabel('$\\log_{10}(g_{\\rm obs}$ / m s$^{-2})$', fontsize=12)
ax.set_title('Radial Acceleration Relation', fontsize=13)
ax.legend(fontsize=9, loc='upper left')
ax.grid(True, alpha=0.3)
ax.set_xlim(-12.5, -8.5)
ax.set_ylim(-12, -8.5)
fig3.tight_layout()
fig3.savefig(os.path.join(outdir, 'fig_sparc_rar.png'), dpi=150)
print("  Saved fig_sparc_rar.png")

# --- Figure 4: Residual histograms ---
fig4, ax = plt.subplots(figsize=(8, 5))
ax.hist(resid_mond, bins=50, alpha=0.4, color='red',
        label=f'MOND fitted: σ={np.std(resid_mond):.3f} dex')
ax.hist(resid_this_work, bins=50, alpha=0.4, color='blue',
        label=f'This work (0 params): σ={np.std(resid_this_work):.3f} dex')
ax.axvline(x=0, color='k', ls=':', alpha=0.3)
ax.set_xlabel('Residual [dex]', fontsize=12)
ax.set_ylabel('Count', fontsize=12)
ax.set_title('Residual distributions vs SPARC', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
fig4.tight_layout()
fig4.savefig(os.path.join(outdir, 'fig_residuals.png'), dpi=150)
print("  Saved fig_residuals.png")

plt.close('all')

# =============================================
# Summary
# =============================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
  Geometric factor:  ∫_backward cos²θ dΩ/(4π) = 1/6
  Predicted a₀:      cH₀/6 = {a0_predicted:.4e} m/s²
  Milgrom's a₀:              {a0_mond:.4e} m/s²
  Discrepancy:                {abs(a0_predicted/a0_mond-1)*100:.1f}%

  SPARC σ (this work):       {np.std(resid_this_work):.4f} dex  (0 free parameters)
  SPARC σ (MOND fitted):     {np.std(resid_mond):.4f} dex  (1 free parameter)
  Mean residual (this work):  {np.mean(resid_this_work):+.4f} dex

  The factor 1/6 = (1/2) × (1/3):
    1/2 from backward hemisphere (Unruh trace direction)
    1/3 from cos²θ mode overlap (planar vs spherical entanglement)
""")
