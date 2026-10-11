# PRN28 medium-elevation v2 validation — final report

## Scope

This is an independent medium-elevation validation at the frozen SignalSim
geometry, using the current receiver binary, a newly calibrated PRN28 response
library, the frozen shared-ratio model, and the frozen 3-of-5 persistence rule.
Truth was not used for response calibration, existence testing, persistence, or
PVT method selection; it was used only for the final position score.

## Infrastructure audit and repair

The first PVT attempt with the original windows stopped before scoring because
the source RINEX stream contained 6, then 9, then 11 satellites as navigation
data became available.  The original clean window therefore could not satisfy
the frozen 11-satellite common-observable requirement.  The original run is
retained as a failed infrastructure batch.

The versioned v2 amendment moved the three 16-window scenes to relative
seconds 91–106, 108–123, and 125–140, where the source stream contains the
complete frozen set at every scored epoch.  No estimator, threshold, rule, or
PVT solver was changed.

## Frozen geometry and checks

- Time: `2021-06-19 12:25:30 UTC`; geometry audit at `12:25:40 UTC`;
- target: GPS L1 C/A PRN28;
- target elevation/azimuth: `45.932503638994° / 244.110613445362°`;
- frozen common set: `G05,G07,G08,G09,G13,G14,G17,G20,G21,G28,G30`;
- precheck: passed (`locked_fraction=0.97028`, C/N0 median `43.361 dB-Hz`,
  C/N0 5th percentile `41.624 dB-Hz`, longest bad run `63 ms`);
- current-binary response shape errors for banks 98/116/132:
  `3.39e-5`, `5.19e-5`, and `4.63e-5`;
- confirmation/PVT seeds: `2026101130`–`2026101134`;
- PVT methods: single, always shared-ratio, gated 3-of-5 shared-ratio.

All five PVT runs had 11 satellites at every scored epoch, native GNSS-SDR
PVT output, and only PRN28 C1C differed between the three common observation
streams.

## Gate result

The gate activated `0/80` windows in each of clean, 0.05, and 0.10 scenes.
Consequently, the gated method is exactly the single-path method in this
medium-elevation batch.  This is a safety/non-inferiority result, not a
positive correction-gain result.

## PVT results

Values are means over five independent noise seeds.  “Improvement” is defined
as single-path 3D RMSE minus the competing method's 3D RMSE.

| Scene | Method | Horizontal RMSE (m) | Vertical RMSE (m) | 3D RMSE (m) | 95% 3D (m) | Gross rate | PVT availability |
|---|---|---:|---:|---:|---:|---:|---:|
| clean | single | 1.8816 | 2.3512 | 3.0198 | 4.7505 | 0 | 1.000 |
| clean | always shared-ratio | 2.7134 | 2.7838 | 3.8937 | 5.8230 | 0 | 0.975 |
| clean | gated 3-of-5 | 1.8816 | 2.3512 | 3.0198 | 4.7505 | 0 | 1.000 |
| 0.05 | single | 1.7150 | 2.4160 | 2.9777 | 5.1211 | 0 | 1.000 |
| 0.05 | always shared-ratio | 2.9706 | 3.9505 | 4.9445 | 7.6149 | 0.0286 | 0.925 |
| 0.05 | gated 3-of-5 | 1.7150 | 2.4160 | 2.9777 | 5.1211 | 0 | 1.000 |
| 0.10 | single | 1.9414 | 2.6814 | 3.3198 | 5.2586 | 0.0125 | 1.000 |
| 0.10 | always shared-ratio | 2.4532 | 3.5289 | 4.3036 | 6.4617 | 0 | 0.975 |
| 0.10 | gated 3-of-5 | 1.9414 | 2.6814 | 3.3198 | 5.2586 | 0.0125 | 1.000 |

The paired 3D-RMSE improvement of always shared-ratio versus single was:

- clean: `-0.8738 m`, bootstrap 95% CI `[-1.1726, -0.5094]`;
- amplitude 0.05: `-1.9667 m`, bootstrap 95% CI `[-2.2042, -1.7383]`;
- amplitude 0.10: `-0.9837 m`, bootstrap 95% CI `[-1.3437, -0.5183]`.

The gated-versus-single difference was exactly `0.0000 m` in all three scenes
and all five seeds because the frozen gate never activated.

## Conclusion

PRN28 supplies the requested independent 40°–50° validation point at
`45.93°`.  Under this geometry, the current model and frozen 3-of-5 gate are
safe: they never introduce the degradation caused by always applying the
shared-ratio correction.  The batch does not show a medium-elevation PVT
improvement because the gate correctly rejected all 80 windows.  Therefore
the defensible paper claim is “cross-geometry protection against unsafe
correction,” while the positive PVT-improvement claim remains supported only
by the previously validated geometries/scenarios where the gate actually
activated.
