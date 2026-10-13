# Short-delay pilot infrastructure amendment

Amended date: 2026-10-13 (Asia/Shanghai)  
Status: frozen before generating any amended validation IQ or reading any
pseudorange value from the failed run.

## Reason

Validation seed `2026101310` completed signal generation, injection QC,
correlator extraction, and receiver execution, but scoring stopped before any
pseudorange metric was calculated because G09 was absent at the originally
scheduled score seconds. A presence-only RINEX audit (no C1C values read)
showed that G09 was continuously available in both paired receiver runs only
from relative second 31 through 67. The failure is therefore classified as a
tracking/acquisition timing infrastructure failure under the frozen seed
policy. The run is retained and is not scored.

## Sole design amendment

The signal factors, estimators, zero-point rule, endpoint, and pass criterion
remain unchanged. Only recording duration and the absolute placement of the
four non-overlapping scenes are changed so all score windows occur after G09
acquisition:

| Cell | Segment (s) | Score windows (s) | Response-bank window |
|---|---:|---:|---:|
| d005_a025 | 34--46 | 36--43 | 39 |
| d005_a050 | 48--60 | 50--57 | 53 |
| d010_a025 | 62--74 | 64--71 | 67 |
| d010_a050 | 76--88 | 78--85 | 81 |

Recording duration is increased from 68 to 96 seconds. The amended isolated
clean calibration recording uses the same duration and the same four score
windows. Response banks and method-specific zero points are recalibrated from
that independent recording; no bank or zero point from the superseded early
windows is reused.

## Seeds and blindness

Seed `2026101310` is retained as an infrastructure-failed run. In accordance
with the already-frozen replacement rule, the five scored seeds are
`2026101311`, `2026101312`, `2026101313`, `2026101314`, and replacement seed
`2026101320`. No numerical pseudorange outcome from seed `2026101310` was
computed or inspected before this amendment.

This amendment does not permit any further estimator, baseline, parameter,
threshold, phase, delay, amplitude, metric, or pass-rule change.
