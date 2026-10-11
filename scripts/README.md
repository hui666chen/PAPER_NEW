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

Before running any script, replace `<LOCAL_USER>`, `<TMP_DIR>`, and other
placeholders with local paths and supply the required response/data artifacts.

