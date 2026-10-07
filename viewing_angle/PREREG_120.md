# PREREG 120 -- does the galaxy pattern depend on our viewing angle? (Coalesce: 2-D sky -> 3-D convention)
Committed BEFORE running iter120_angle.py.
Rotation speeds are measured along our line of sight and divided by sin(inclination) to get the true speed; a tilt error
changes speed^2 (the pull) by 1/sin^2. Physics cannot care how we view a galaxy; an artefact of the correction would.
Data: SPARC rotation curves + main-table inclinations (browser, checksummed: 175 rows, sum 10344 deg).
Split: low tilt 30-50 deg (most sensitive to tilt errors), middle 50-70, high 70-90 (edge-on: speeds need little
correction but disks are seen through themselves). Below 30 deg excluded as usual (huge correction).
Measure per group: median boost in bins of visible pull (default methods: mass-to-light 0.5, also gas-dominated points),
and the scatter around the all-galaxy curve.
Expectation: the pattern is the same in all three groups to within ~0.1 dex in boost; scatter larger in the low-tilt
group (bigger correction errors). A systematic difference > 0.15 dex between groups would flag the tilt convention.
