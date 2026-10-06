# PREREG 118 -- clusters and groups of the same size: does visible matter alone set their gravity? (Coalesce)
Committed BEFORE running iter118_pairs.py.
Data: clusters.csv (Vikhlinin+2006 10 clusters; Sun+2009 23 groups with direct M500; 31 systems after duplicates).
Boost B = M500 / (M_gas + M_stars); stars f* = 0.025 (M500/1e14)^-0.37 (memory, Gonzalez+2013-like; sensitivity f* = 0
and f* = 0.02). Pulls at r500: g_obs = G M500/r500^2, g_bar = g_obs / B.
NOTE (definitional): r500 encloses 500x the critical density, so at a given size the total pull is fixed by construction;
what can vary at the same size is the visible share, i.e. the boost.
Rule tested -- 'visible matter alone sets gravity' in its galaxy form g_obs = g_bar nu(g_bar/a), nu = 1/(1-exp(-sqrt x)):
 T1 best-fit acceleration scale a for these 31 systems vs galaxies' 1.2e-10. Expect a ~ 1e-9 (x5-20 galaxies).
 T2 same size: in bins of r500, the rule fixes the boost exactly (one g_obs -> one g_bar). Measure intrinsic scatter of
    log B beyond errors. Expect 0.05-0.12 dex; significance uncertain with these errors (report it either way).
 T3 trend: the rule (deep regime) predicts visible fraction rising ~linearly with size (slope ~1 in log f_vis vs log r500,
    ~0.7 in the transition). Expect measured slope 0.2-0.6 (shallower).
Free choices: stellar fractions (1, with sensitivity), simple nu (1). Hydrostatic masses ~10-20% low (known bias) -- noted.
