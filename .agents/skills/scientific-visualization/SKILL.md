---
name: scientific-visualization
description: Publication standards for generating clean, beautiful, high-data-density vector figures for ACS journals (JPCB/JCTC) using Python/matplotlib.
---

# Scientific Visualization & Plotting Protocol

## 1. Visual Aesthetics & Standards
- **Font Selection**: Use standard serif fonts (`DejaVu Serif`, `Times New Roman`, or `Computer Modern`) matching LaTeX typesetting.
- **Figure Dimensions**:
  - Single-column width: 3.25–3.33 inches (8.5 cm)
  - Two-column / double width: 6.75–7.00 inches (17.5 cm)
  - Table of Contents (TOC) Graphic: 3.25 inches wide $\times$ 1.75 inches tall (8.25 cm $\times$ 4.45 cm)
- **Palette**: Use distinct, colorblind-safe, curated palette (e.g. Deep Blue `#1e3a8a`, Forest Green `#065f46`, Rust Orange `#c2410c`, Plum `#701a75`, Slate `#334155`). Avoid neon defaults and low-contrast pastel washes.
- **Line & Marker Quality**:
  - Line widths: 1.5–2.0 pt.
  - Scatter markers: sized appropriately with subtle alpha transparency (0.7–0.85) to resolve overlapping densities.
  - Axes spines: 0.8–1.0 pt, ticks facing outwards or subtle, light gridlines (`linestyle=":"`, `alpha=0.5`).

## 2. High Data-to-Ink Ratio
- Strip away decorative chartjunk, excess borders, and redundant legends.
- Label key data points or curve features directly where possible.
- Include units on every axis without exception (e.g., $r$ (Å), $W(r)$ (kcal/mol), Time (ns), $k_{\rm on}$ ($\text{M}^{-1}\text{s}^{-1}$)).
- Provide uncertainty bounds, standard error shading, or error bars for computed free energies and rate constants.
