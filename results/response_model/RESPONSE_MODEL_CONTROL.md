# Receiver-chain response control experiment

## Purpose

Test whether the simple assumption “the same clean template can be shifted
continuously” describes the response produced by the actual receiver chain.
The estimator, search range, and fitting procedure were held fixed; only the
response representation was changed.

## Result

| Metric | Simple continuous-shift template | Receiver-chain-consistent response |
|---|---:|---:|
| Delay RMSE | 0.006924 chip | 0.0001455 chip |
| Maximum absolute error | 0.012985 chip | 0.000255 chip |
| Cases within 0.01 chip | 40/45 | 45/45 |

## Interpretation

The result identifies response-model mismatch as a real source of systematic
delay error. It does not prove that the response bank is valid for every
receiver state, PRN, or dynamic condition. The paper should describe the bank
as calibrated for the tested receiver binary and validation geometry.

## Related code

- `scripts/receiver_chain_response_control.py`
- `scripts/test_receiver_chain_response.py`
- `scripts/test_profile_interpolation.py`

