# Selected scripts

These are the analysis and experiment scripts retained for provenance and
future reproduction. They are not a turnkey package: local SignalSim,
GNSS-SDR, CUDA, navigation files, and response banks were intentionally not
uploaded. Paths were replaced with placeholders before publication.

The scripts most relevant to the paper are:

- `receiver_chain_response_control.py`: response-model control experiment;
- `test_receiver_chain_response.py`: receiver-response checks;
- `test_profile_interpolation.py`: interpolation/fractional-delay checks;
- `compare_tau_pseudorange_proxy.py`: direct-path delay and pseudorange-proxy comparison;
- `pvt_geometry_propagation.py`: geometry/PVT propagation utilities;
- `crossfit_basin_confidence.py`: exploratory observation-only confidence analysis;
- `fit_dual_component_delay.py` and `fit_dual_component_with_template.py`: dual-path fitting prototypes.

The frozen classical-baseline package is under `scripts/baselines/`:

- `correlator_baselines.py`: EML, normalized Double-Delta, and MEDLL ideal-vs-chain template fits on retained 401-point observations;
- `run_prn09_baseline_observation.py`: opt-in PRN9 rerun that retains large mixed response CSVs outside Git;
- `gate_roc.py`: single-window and K-of-N operating-point tables;
- `design_parameter_sweep.py`: deterministic delay/amplitude/phase/C/N0 design;
- `runtime_benchmark.py`: reference algorithm timing plus a separately entered CUDA real-time factor.
- `baseline_schedule.py`: deterministic median-per-window adapter into the existing `correction_m` PVT schedule format.
- `crb_identifiability_grid.py`: profiled M0/M1/M2 CRB grid with explicit reference-C/N0 scaling.

The protocol and exact settings are in
`protocols/BASELINE_AND_GENERALIZATION_PROTOCOL.md`.  These additions do not
claim numerical baseline results until a fresh retained observation batch has
been processed.

Before running any script, replace `<LOCAL_USER>`, `<TMP_DIR>`, and other
placeholders with local paths and supply the required response/data artifacts.
