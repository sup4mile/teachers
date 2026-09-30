# First-stage estimates of the ability block (T2b)

Generated 2026-09-29 18:34:33 by `julia/spatial_model/calibration/spatial_first_stage.jl`. First pass outside the equilibrium: Q_l and t_l fixed, selection into teaching ignored for both genders (spatial_calibration.md §1, Table 2, *Ability block*). η = 0.08, α = 1.0, b = 1/(1 − η) = 1.08696, χ = log 1.12 = 0.11333.

## Estimates

| parameter | estimate | target | source |
|---|---|---|---|
| ρz | 0.5766 | latent mother–child correlation 0.577 (SE 0.014) | NLSY79/CNLSY; closed form |
| s_z | 0.1956 (SE 0.0071) | latent log-wage slope 0.2127 (SE 0.0077) = b α s_z | NLSY79; closed form |
| σξ = s_z √(1 − ρz²) | 0.1598 | — | implied |
| σϵ (Nz = 5) | 0.7861 | pooled within-occupation 90/10 3.7374 | ACS 2009–13; model 3.7374 |
| σϵ (continuous z) | 0.7870 | same | — |
| c = b s_z/χ | 1.8765 | — | score loading of log Q_l |

## Checks

- **χ gap:** b α s_z − χ = 0.2127 − 0.1133 = 0.0993 (SE 0.0077): the NLSY slope is 12.9 SEs above χ.
- **Siblings:** model latent sibling correlation ρz² = 0.3325, against 0.6063 (SE 0.0165) in the data (check only; siblings share inputs outside the model).
- **Corner (caveat 1):** z alone gives a 90/10 of 1.725, well below 3.737, so σϵ is interior.
- **Components:** at σϵ̂ the pooled 90/10 is 3.313 with ϵ alone (z shut off); z alone gives 1.725 for continuous z and 1.530 on the Nz = 5 lattice, where q10 and q90 fall on nodes. The no-selection lognormal benchmark (`sigma_eps_placeholder`) gives σϵ = 0.4309, so Roy selection within cells raises σϵ̂ by 82%.
- **Shares block:** the inversion reproduces the share targets to 4.4e-13 (max abs. error, both genders, after quadrature).
- **Occupational block at σϵ̂:** log A_i/A_HP ∈ [-1.526, -0.327]; log r_i ∈ [-1.577, 0.436] (market occupations; A_HP = r_HP = 1).
- **Cell 90/10s (untargeted):** model by gender 3.646 (men), 3.851 (women) vs. data 3.821, 3.633; weighted correlation of model and data cell ratios -0.188; model range [3.21, 4.15] vs. data [2.80, 6.59].

## Ability grid (Nz)

σϵ̂ on each Rouwenhorst grid, and each grid's pooled 90/10 at the continuous-limit σϵ̂ = 0.7870 (target 3.7374).

| Nz | σϵ̂ | Δ vs. limit | pooled 90/10 at limit σϵ̂ |
|---|---|---|---|
| 3 | 0.7852 | -0.0018 | 3.7468 |
| 5 | 0.7861 | -0.0009 | 3.7420 |
| 7 | 0.7864 | -0.0006 | 3.7404 |
| 9 | 0.7866 | -0.0004 | 3.7397 |
| 15 | 0.7867 | -0.0002 | 3.7387 |
| 25 | 0.7868 | -0.0001 | 3.7382 |
| continuous | 0.7870 | +0.0000 | 3.7374 |

## Table 5 variants

Each row changes one input and refits σϵ on the Nz = 5 grid. Rouwenhorst nodes and weights depend on s_z but not ρz, so the ρz rows change only ρz and σξ.

| variant | ρz | s_z | σξ | σϵ̂ | c | b s_z − χ | ρz² |
|---|---|---|---|---|---|---|---|
| Baseline | 0.5766 | 0.1956 | 0.1598 | 0.7861 | 1.876 | 0.0993 | 0.333 |
| ρz: cross-sectional mothers (0.529) | 0.5293 | 0.1956 | 0.1660 | 0.7861 | 1.876 | 0.0993 | 0.280 |
| ρz: unweighted (0.630) | 0.6300 | 0.1956 | 0.1519 | 0.7861 | 1.876 | 0.0993 | 0.397 |
| b s_z = χ | 0.5766 | 0.1043 | 0.0852 | 0.8435 | 1.000 | 0.0000 | 0.333 |
| NLSY97 slope (0.199) | 0.5766 | 0.1828 | 0.1494 | 0.7966 | 1.754 | 0.0854 | 0.333 |
| Age-residualized 90/10 (3.635) | 0.5766 | 0.1956 | 0.1598 | 0.7661 | 1.876 | 0.0993 | 0.333 |
| η = 0.073 | 0.5766 | 0.1971 | 0.1611 | 0.7921 | 1.876 | 0.0993 | 0.333 |
| η = 0.103 | 0.5766 | 0.1908 | 0.1558 | 0.7665 | 1.876 | 0.0993 | 0.333 |

## Occupational block and cell 90/10s at σϵ̂

A_i/A_HP from men's shares, r_i = Θ_{i,f}/(ωf A_i) from women's (1 − τω_{i,f} = ωf r_i; ωf is internal). Shares are over the 20 non-teaching occupations. Cell 90/10s are untargeted checks on a common σϵ.

| occupation | share m | share f | A_i/A_HP | r_i | 90/10 m: model | data | 90/10 f: model | data |
|---|---|---|---|---|---|---|---|---|
| Executives | 0.0702 | 0.0629 | 0.5666 | 0.7355 | 3.72 | 4.15 | 3.84 | 3.59 |
| Management related | 0.0319 | 0.0401 | 0.4131 | 0.8390 | 3.52 | 3.38 | 3.71 | 2.97 |
| Architects, engineers, CS | 0.0462 | 0.0139 | 0.4767 | 0.4930 | 3.60 | 2.92 | 3.48 | 2.99 |
| Scientists, arts | 0.0314 | 0.0371 | 0.4107 | 0.8181 | 3.51 | 3.94 | 3.69 | 3.24 |
| Doctors and lawyers | 0.0140 | 0.0135 | 0.3084 | 0.7545 | 3.36 | 6.59 | 3.48 | 5.48 |
| Nurses, therapists | 0.0158 | 0.0939 | 0.3214 | 1.5470 | 3.38 | 5.62 | 3.98 | 4.53 |
| Postsecondary teachers | 0.0078 | 0.0078 | 0.2547 | 0.7643 | 3.27 | 4.23 | 3.39 | 3.97 |
| Other teachers, librarians | 0.0046 | 0.0155 | 0.2173 | 1.1220 | 3.21 | 3.31 | 3.50 | 3.58 |
| Technicians | 0.0395 | 0.0333 | 0.4481 | 0.7194 | 3.57 | 3.90 | 3.66 | 3.44 |
| Sales | 0.0791 | 0.0681 | 0.5967 | 0.7225 | 3.76 | 4.71 | 3.86 | 4.62 |
| Administrative support | 0.0626 | 0.1407 | 0.5398 | 1.1159 | 3.69 | 3.27 | 4.15 | 2.80 |
| Fire, police, guards | 0.0322 | 0.0079 | 0.4143 | 0.4713 | 3.52 | 3.51 | 3.39 | 3.50 |
| Food, cleaning, personal | 0.0579 | 0.0740 | 0.5226 | 0.8552 | 3.67 | 3.50 | 3.89 | 3.86 |
| Farm, extraction | 0.0295 | 0.0047 | 0.4013 | 0.4132 | 3.50 | 4.51 | 3.32 | 4.19 |
| Mechanics, construction | 0.1201 | 0.0032 | 0.7213 | 0.2067 | 3.91 | 3.66 | 3.28 | 4.00 |
| Precision manufacturing | 0.0165 | 0.0051 | 0.3259 | 0.5240 | 3.39 | 3.19 | 3.34 | 3.31 |
| Manufacturing operators | 0.0263 | 0.0077 | 0.3843 | 0.5038 | 3.48 | 3.52 | 3.39 | 3.23 |
| Fabricators, handlers | 0.0434 | 0.0112 | 0.4649 | 0.4699 | 3.59 | 3.39 | 3.45 | 3.24 |
| Vehicle operators | 0.0417 | 0.0037 | 0.4578 | 0.3388 | 3.58 | 3.47 | 3.30 | 3.56 |
| Home production | 0.2293 | 0.3556 | 1 | 1 | — | — | — | — |

## Limitations of the first pass

- Teaching selection is ignored in shares and wage distributions for both genders; the §4 consistency check measures it after the first internal fit.
- Cross-location variation in Q_l and t_l and the goods moving cost (Ξ) are left out of wages; both add dispersion, so σϵ̂ is an upper bound in that respect.
- ACS wages include measurement error and transitory shocks, which also load on σϵ (caveat 2).
- The ACS block is 2009–13, the spatial targets late 2010s (caveat 3).
- Mean wages by occupation (Table 2 check) are not in the workbook and are not compared here.
