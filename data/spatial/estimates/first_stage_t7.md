# First stage on the 2016–19 occupational block (T7)

Generated 2026-09-29 18:34:39 by `spatial_first_stage.jl --t7`. Data: pooled ACS 1-year PUMS 2016–19 with the workbook's definitions (`acs_occupations.py`; [report](acs_occupations.md)); the ability block (ρz, s_z) is unchanged. Same first pass as T2b: Q_l and t_l fixed, teaching selection ignored.

| | 2009–13 | 2016–19 |
|---|---|---|
| Pooled non-teacher 90/10 (target) | 3.7374 | 3.8214 |
| σϵ̂ (Nz = 5) | 0.7861 | 0.8063 |
| K–12 teaching share, men / women | 0.0194 / 0.0599 | 0.0181 / 0.0571 |
| Home production share, men / women | 0.2293 / 0.3556 | 0.2056 / 0.3123 |
| Teachers' 90/10, men / women / pooled | 2.597 / 2.830 / 2.773 | 2.587 / 2.850 / 2.786 |
| Years of schooling, non-teachers / K–12 teachers | 13.64 / 16.04 | 13.87 / 15.93 |
| log A_i/A_HP range | [-1.526, -0.327] | [-1.446, -0.307] |
| log r_i range | [-1.577, 0.436] | [-1.499, 0.441] |

Across the 19 market occupations, log A_i/A_HP moves by +0.054 on average (SD 0.038; correlation 0.992) and log r_i by +0.036 (SD 0.036; correlation 0.998). Market A_i rise against home production because home-production shares fell (men -2.4 pp, women -4.3 pp); the internal levels (κ̃, ωf) move with them.

| occupation | log A_i 09–13 | log A_i 16–19 | log r_i 09–13 | log r_i 16–19 | 90/10 m: model / data 16–19 | 90/10 f: model / data 16–19 |
|---|---|---|---|---|---|---|
| Executives | -0.568 | -0.463 | -0.307 | -0.247 | 3.86 / 4.26 | 3.99 / 3.68 |
| Management related | -0.884 | -0.805 | -0.175 | -0.190 | 3.62 / 3.45 | 3.78 / 3.17 |
| Architects, engineers, CS | -0.741 | -0.620 | -0.707 | -0.674 | 3.74 / 3.06 | 3.60 / 3.14 |
| Scientists, arts | -0.890 | -0.815 | -0.201 | -0.155 | 3.62 / 4.18 | 3.80 / 3.36 |
| Doctors and lawyers | -1.176 | -1.157 | -0.282 | -0.232 | 3.43 / 6.82 | 3.55 / 5.79 |
| Nurses, therapists | -1.135 | -1.017 | 0.436 | 0.441 | 3.50 / 5.52 | 4.11 / 4.48 |
| Postsecondary teachers | -1.368 | -1.365 | -0.269 | -0.239 | 3.33 / 4.34 | 3.44 / 4.11 |
| Other teachers, librarians | -1.526 | -1.446 | 0.115 | 0.163 | 3.30 / 3.26 | 3.61 / 3.69 |
| Technicians | -0.803 | -0.707 | -0.329 | -0.338 | 3.68 / 4.60 | 3.75 / 3.79 |
| Sales | -0.516 | -0.498 | -0.325 | -0.292 | 3.83 / 4.86 | 3.93 / 4.80 |
| Administrative support | -0.617 | -0.560 | 0.110 | 0.074 | 3.79 / 3.36 | 4.19 / 2.87 |
| Fire, police, guards | -0.881 | -0.863 | -0.752 | -0.731 | 3.59 / 3.55 | 3.45 / 3.59 |
| Food, cleaning, personal | -0.649 | -0.600 | -0.156 | -0.135 | 3.76 / 3.60 | 3.97 / 3.77 |
| Farm, extraction | -0.913 | -0.916 | -0.884 | -0.761 | 3.56 / 4.30 | 3.41 / 4.29 |
| Mechanics, construction | -0.327 | -0.307 | -1.577 | -1.499 | 3.99 / 3.53 | 3.36 / 3.92 |
| Precision manufacturing | -1.121 | -1.072 | -0.646 | -0.594 | 3.47 / 3.17 | 3.42 / 3.24 |
| Manufacturing operators | -0.956 | -0.929 | -0.686 | -0.663 | 3.55 / 3.43 | 3.45 / 3.22 |
| Fabricators, handlers | -0.766 | -0.710 | -0.755 | -0.704 | 3.68 / 3.36 | 3.54 / 3.25 |
| Vehicle operators | -0.781 | -0.734 | -1.082 | -1.008 | 3.67 / 3.45 | 3.38 / 3.70 |
