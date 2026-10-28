# Experimental higher-relation searches

These scripts explore the action of lasso twists on sliders in quotients with a
specified two-strand untwisting relation. They use specialized bases and Krylov
calculations to propose annihilating or minimal polynomials. They are research
sources outside the maintained installed package.

- `theta3_min_poly.py` records the earlier three-strand search.
- `theta_min_poly.py` extends the computational search to larger strand counts.
- `gap_scripts/theta3_exact.gap` in the repository root is an independent exact
  computational route retained with its historical evidence.

Agreement of specialized dimensions supplies evidence for candidate generic
patterns. Small finite-field coefficients do not establish exact rational
identities. An annihilating relation needs an exact witness;
minimality and dimension need additional certified-basis or representation
evidence. A successful run at finitely many parameters does not supply either
conclusion by itself.

The [mathematical specification](../../docs/mathematical-specification/main.tex)
defines what the maintained certificate checks establish.
