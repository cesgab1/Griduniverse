# PREREG 119 -- cluster points without a boundary convention (Coalesce: 'measure each point independently')
Committed BEFORE running iter119_points.py.
Each measured radius in each system is one data point: r, total enclosed mass (gravity / hot-gas balance), visible
enclosed mass (gas glow + stars). The radius labels (r2500, r500) only choose WHERE to sample; values at a radius do not
depend on the convention. Compare each point with the galaxy rule g_obs = g_bar nu(g_bar/1.2e-10) -- exactly as for
galaxy rotation-curve points (expedition leg 1).
Data: Vikhlinin+2006 (M2500, fg2500, M500, fg500; z from their Table 1), Sun+2009 (r2500, fg2500 for all groups with
r2500; r500, M500, fg500 where direct; z from their Table 1). Duplicates: Sun kept. Enclosed mass at r2500 for Sun groups
from r2500 and the critical density at their z (that step uses the definition, but only to recover the measured mass).
Stars: same three assumptions as iteration 118.
Expectations: points sit ABOVE the galaxy rule at both radii; offset larger in the inner points (r2500: x1.6-2.5) than
the outer (r500: x1.3-1.6) -- i.e. the shortfall is worst toward cluster centres. Galaxy rule predicts x1 everywhere.
