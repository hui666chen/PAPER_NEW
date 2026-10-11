# Claims and limitations for manuscript writing

## Safe claims

1. The continuous-shift response assumption produces measurable fractional-delay
   mismatch in the tested receiver chain.
2. A calibrated receiver-chain-consistent response substantially reduces this
   mismatch in the controlled 45-case experiment.
3. The classical dual-path correction improves PVT in the tested strong,
   identifiable multipath condition.
4. The 3-of-5 gate prevents unsafe correction in tested clean, weak, and
   medium-elevation cases.
5. The common-observable PVT protocol removes satellite-search timing as a
   comparison confounder.

## Claims to avoid

- “works for all satellite elevations”;
- “improves PVT for all multipath amplitudes”;
- “the gate is a new statistical theory”;
- “P_global is directly available in a receiver” when using the oracle truth-
  based analysis;
- “hardware validated” unless the exact reported result used the RF/USRP chain;
- “AI-based” or “intelligent” unless a real AI method is actually part of the
  final experiment.

## Main limitations

- response banks are calibrated for the tested receiver binary and states;
- medium-elevation batches demonstrate safety but not positive gain;
- historical PRN30 records contain a 3-of-4 protocol version;
- PRN9 batches with different protocol versions are not pooled;
- weak/short-delay reflection parameters may be statistically unidentifiable;
- absolute pseudorange/PVT benefit depends on satellite geometry and the
  corrected PRN's leverage in the full constellation.

