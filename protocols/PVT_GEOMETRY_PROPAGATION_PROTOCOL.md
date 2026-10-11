# Frozen PVT geometry-propagation calculation

Status: post-hoc engineering propagation using an already frozen 3x50
single-channel pilot. This calculation does not change the channel estimator,
gate, threshold, or sample set.

## Geometry

Use SignalSim's own `CNavData`, `GpsSatPosSpeedEph`, and `SatElAz` code with
the exact test configuration:

- UTC: 2021-06-19 10:05:40;
- receiver: latitude 37.352721 deg, longitude -121.915773 deg, altitude 20 m;
- navigation file: `BRDC00IGS_R_20211700000_01D_MN.rnx`;
- GPS elevation mask: 3 deg.

Use every healthy visible GPS satellite. PRN14 is the only multipath-affected
channel. Do not choose or remove satellites after viewing PVT outcomes.

## Channel errors

Use the existing 150 direct-delay errors converted by
`c / 1.023e6 m/chip`:

- always single: saved single-path tau estimate;
- always two: saved full two-path tau estimate;
- gated: two-path iff the already frozen `P_cross <= 0.20`, otherwise single.

No additional random error is added to PRN14 because its saved tau error is
treated as the complete affected-channel pseudorange error.

## Other-satellite baseline noise

For every unaffected satellite, use independent zero-mean Gaussian code noise
with standard deviation

`sigma_i = sigma_zenith / sin(elevation_i)`.

The primary result uses `sigma_zenith = 1.0 m`. Sensitivity values 0.5 m and
2.0 m are reported without selecting among them.

Do not draw extra Monte Carlo samples. Because linearized WLS is linear in
pseudorange errors, analytically add the unaffected-satellite noise covariance
to the squared position errors from each of the existing 150 PRN14 errors.

## WLS

In local ENU coordinates use one row per satellite:

`[-cos(el) sin(az), -cos(el) cos(az), -sin(el), 1]`.

Use `W_ii = 1 / sigma_i^2`, solve position plus receiver clock by
`(H'WH)^-1 H'W`, and use the same geometry and weights for all methods.

## Endpoints

Report horizontal, vertical, and 3D position RMSE over the existing 150
channel-error trials after analytic averaging over unaffected-channel noise.
Also report per-cell endpoints, satellite geometry, DOP diagnostics, PRN14
WLS influence coefficients, and gate-minus-single absolute and relative
changes. These are geometry-propagation proxies, not a complete nonlinear PVT
receiver experiment.
