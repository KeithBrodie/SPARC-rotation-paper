# The Radial Acceleration Relation from Two-Horizon Entropy Sharing with Zero Free Parameters

**Author:** Keith Brodie (2026)

## Summary

We derive MOND's interpolation function from first principles using Jacobson's thermodynamic gravity, with zero free parameters. The Hubble horizon competes with local Rindler horizons for vacuum entanglement. Their different geometries — planar (Unruh) and spherical (Hubble) — create angle-dependent competition described by the cos²θ mode overlap, integrated over the backward hemisphere, yielding a geometric factor of 1/6. The result:

$$f(a) = \frac{a}{a + cH_0/6}$$

with $a_0 = cH_0/6 = 1.09 \times 10^{-10}$ m/s² — 9% from Milgrom's empirical value.

## Key Result

| Model | Free params | a₀ (m/s²) | σ (dex) | Mean residual (dex) |
|-------|:-:|:-:|:-:|:-:|
| **This work** | **0** | **1.09 × 10⁻¹⁰** | **0.133** | **+0.005** |
| MOND (fitted) | 1 | 1.20 × 10⁻¹⁰ | 0.133 | −0.009 |
| Bare (cH₀) | 0 | 6.55 × 10⁻¹⁰ | 0.148 | −0.296 |

Tested against the SPARC Radial Acceleration Relation (2,693 data points, 153 galaxies), the zero-parameter prediction matches fitted MOND at the intrinsic scatter floor of the data.

## Files

- `Draft-v1.md` — Full paper text (Markdown with LaTeX math)
- `two_horizon_sparc.py` — Complete numerical analysis script (reproduces all tables and figures)
- `data/RAR.mrt` — SPARC RAR dataset (McGaugh, Lelli & Schombert 2016)
- `fig_f_curves.png` — Inertia modification function f(a)
- `fig_sparc_rar.png` — SPARC RAR with predictions
- `fig_residuals.png` — Residual histograms

## Running the Code

```bash
pip install numpy scipy matplotlib
python two_horizon_sparc.py
```

Requires only `numpy`, `scipy`, and `matplotlib`. Generates all four figures and prints Tables 1 and 2.

## Related Papers

1. **Paper 1:** K. Brodie, "Quantized Inertia as a Boundary Correction to Jacobson's Thermodynamic Spacetime" (2026). [DOI: 10.5281/zenodo.18664801](https://doi.org/10.5281/zenodo.18664801)

2. **Paper 2:** K. Brodie, "Karlsson's Redshift Periodicity as an Efimov Spectrum" (2026). [DOI: 10.5281/zenodo.18664931](https://doi.org/10.5281/zenodo.18664931)

3. **Paper 3:** K. Brodie, "Accelerated Structure Formation from Horizon Thermodynamics" (2026). [DOI: 10.5281/zenodo.18665076](https://doi.org/10.5281/zenodo.18665076)

## License

This work is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
