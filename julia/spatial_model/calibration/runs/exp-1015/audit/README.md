# Teacher-supply audit, final script (2026-10-15)

`spatial_supply_audit.jl` (snapshot r1015w, `--two-sided`; job 5965513) at each fitted run: `audit.txt` (report) and `audit.toml`
(all margins, shocks and checks). Columns: workplace nesting ρw and σw; the staffing-ratio elasticity ε22 − ε12 for a
suburban pay rise (central difference) as margins are added, (a) workplace only, (b) + residence, (c) + Λ re-solved,
(d) + occupation and education, (e) + stationary distribution, (f) full GE with its ε22 and ε12; the small-district
elasticity μy/(σw C); commuting directions, teachers (T) and others (NT): P(work suburb | live city), P(work city |
live suburb). ACS: T 0.0862, 0.0832; NT 0.0919, 0.1164. Built by `experiments-2026-10d/audit_table.py`.

| case | ρw | σw | (a) | (b) | (c) | (d) | (e) | (f) GE ratio | GE ε22 | GE ε12 | district | T: c→s, s→c | NT: c→s, s→c |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| c0 | nan | nan | nan | 1.4 | 1.4 | 4.0 | 4.1 | 1.3 | 1.51 | 0.25 | 4.4 | —, — | —, — |
| c1 | 0.066 | 0.017 | 28.5 | 29.5 | 29.5 | 30.9 | 30.9 | 10.9 | 5.54 | -5.37 | 66.9 | 0.159, 0.090 | 0.127, 0.092 |
| c3 | 0.066 | 0.017 | 28.4 | 30.2 | 30.3 | 31.0 | 31.1 | 12.4 | 6.26 | -6.12 | 68.5 | 0.162, 0.080 | 0.128, 0.092 |
| c4p | 0.068 | 0.017 | 28.1 | 30.0 | 30.1 | 30.8 | 30.9 | 12.3 | 5.88 | -6.40 | 68.7 | 0.145, 0.090 | 0.133, 0.088 |
| c4p_closure | 0.068 | 0.017 | 28.2 | 30.1 | 30.2 | 30.8 | 30.9 | 12.4 | 5.79 | -6.61 | 68.7 | 0.132, 0.101 | 0.133, 0.088 |
| zeta_closure | 0.066 | 0.016 | 32.9 | 34.4 | 34.5 | 35.1 | 35.2 | 16.3 | 7.18 | -9.15 | 68.9 | 0.075, 0.168 | 0.125, 0.091 |
| s1 | 1.000 | 0.250 | 1.8 | 3.7 | 3.8 | 4.6 | 4.7 | 1.7 | 1.64 | -0.02 | 4.7 | 0.142, 0.084 | 0.130, 0.091 |
| s1_u | 1.000 | 0.250 | 1.7 | 3.6 | 3.6 | 4.1 | 4.2 | 5.1 | 3.04 | -2.09 | 4.4 | 0.138, 0.087 | 0.131, 0.091 |
| s1_e | 1.000 | 0.250 | 2.1 | 3.9 | 4.0 | 4.6 | 4.7 | 1.8 | 1.66 | -0.18 | 4.7 | 0.095, 0.127 | 0.129, 0.090 |
