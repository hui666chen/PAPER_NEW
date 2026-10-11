# PRN28 independent medium-elevation validation protocol — v2

Status: frozen before v2 calibration and confirmation data generation.

## Reason for the versioned amendment

The original frozen PRN28 protocol used score windows beginning at relative
second 61.  The first common-observable source run showed that the receiver
had acquired all channels, but its RINEX stream contained only 6 satellites
through relative second 59, 9 satellites from second 60, and the complete
11-satellite set from second 90.  This is a source-observation readiness issue,
not an estimator result.  The v2 amendment moves the three scene windows into
the already observed 11-satellite interval.  No model, threshold, persistence
rule, correction formula, or PVT solver is changed.

The original protocol, source run, and failure log remain retained as an
infrastructure audit and are not silently replaced by this amendment.

## Frozen geometry

- Navigation file: `<LOCAL_USER>/chh/GNSS/SignalSim-main/EphData/BRDC00IGS_R_20211700000_01D_MN.rnx`
- Receiver LLA: `37.352721 deg, -121.915773 deg, 20 m`
- Signal start: `2021-06-19 12:25:30 UTC`
- Geometry audit epoch: `2021-06-19 12:25:40 UTC`
- Target: GPS L1 C/A PRN28
- Target elevation: `45.932503638994 deg`
- Target azimuth: `244.110613445362 deg`
- Visible GPS set at the audit epoch: `G05,G07,G08,G09,G13,G14,G17,G20,G21,G28,G30`
- Frozen common-observable set: the same 11 satellites
- Elevation mask: `3 deg`

## Precheck acceptance rule

The v1 PRN28 precheck already passed with the current receiver binary.  The
same thresholds remain in force for the v2 source stream: after a one-second
acquisition allowance, lock fraction at least 0.95, C/N0 median at least 35
dB-Hz, C/N0 fifth percentile at least 30 dB-Hz, maximum bad run at most 100
ms, and at least 10 seconds of post-acquisition state rows.

## Formal v2 experiment

- Current-binary PRN28 response calibration windows: `98`, `116`, `132`;
- the response-bank GPU shape error must be at most `0.01` for every bank;
- frozen shared-ratio model and null threshold;
- frozen persistence rule: `3-of-5`;
- scenes, each with 16 scoring windows:
  - clean: relative seconds `91–106`;
  - amplitude ratio 0.05: relative seconds `108–123`;
  - amplitude ratio 0.10: relative seconds `125–140`;
- delay `0.8` chip;
- relative phase `8.833159347728436` degrees;
- independent seeds `2026101130–2026101134`;
- common-observable PVT comparison: single, always shared-ratio, and gated
  persistence shared-ratio, changing only PRN28 C1C;
- the source observation stream and satellite set must be identical across
  all three PVT methods;
- truth is used only for final PVT scoring, never for response/gate decisions.

The claim remains limited to this frozen PRN28/time geometry.  Passing v2
establishes medium-elevation safety/non-inferiority under the frozen gate; it
does not by itself claim a positive correction gain if the gate remains
inactive.
