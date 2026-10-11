# Frozen PRN14 elevation-counterfactual PVT protocol

Status: post-hoc counterfactual geometry analysis. Freeze this file before
running the three requested elevation cases.

## Fixed inputs

Reuse `PVT_GEOMETRY_PROPAGATION_PROTOCOL.md` without changing:

- the nine-satellite SignalSim geometry and PRN14 azimuth;
- the existing 150 PRN14 delay-error trials;
- the three methods: gated, always single, and always two;
- the gate rule `P_cross <= 0.20`;
- the unaffected-satellite noise model;
- the ENU-plus-clock WLS formula and elevation weighting;
- primary zenith code sigma 1.0 m and sensitivity values 0.5/2.0 m.

## Counterfactual variable

Replace only PRN14 elevation by exactly 15, 30, and 45 degrees in three
separate runs. Keep its azimuth fixed at the SignalSim value 211.408198280931
degrees. The original 4.108709808829-degree result is the reference point.

These are synthetic counterfactual geometries, not claims that PRN14 occupied
those elevations at the frozen epoch.

## Primary endpoint

At zenith code sigma 1.0 m, use the gated-versus-always-single 3D position RMSE
reduction, reported in metres and percent. Horizontal and vertical reductions
are secondary endpoints.

## Frozen interpretation rules

1. **Structural elevation support:** the primary 3D relative improvement must
   strictly increase over the sequence 4.1087, 15, 30, and 45 degrees. The
   45-degree improvement must also be at least five times the original
   4.1087-degree relative improvement.
2. **Practical PVT support:** in addition to structural support, the 45-degree
   gated-versus-single 3D RMSE reduction must be at least 0.010 m and at least
   0.5%.
3. If rule 1 passes but rule 2 fails, conclude that elevation leverage is
   mechanistically supported but insufficient for a meaningful single-channel
   PVT claim.
4. If rule 1 fails, do not claim a monotonic elevation effect from this
   geometry. Report the component-wise nonmonotonicity without changing the
   tested elevations or endpoint.

Always-two results remain a safety comparison and are not used to decide
whether the elevation hypothesis passes.
