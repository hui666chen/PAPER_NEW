"""Observation-only adapter for the existing epoch-response CSV format."""
import argparse
import json
from pathlib import Path

from scripts.baselines.correlator_baselines import read_epoch_scans
from scripts.gss_multipath.gss import Config, SampledCorrelation, estimate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-taps", type=int, default=401)
    parser.add_argument("--axis", choices=("receiver-shift", "physical-delay"),
                        default="receiver-shift")
    args = parser.parse_args()
    scans = read_epoch_scans(args.input, expected_taps=args.expected_taps)
    config = Config()
    rows = []
    for scan in scans:
        x, z = scan.delay_chip, scan.mixed
        if args.axis == "receiver-shift":
            x, z = -x[::-1], z[::-1]
        correlation = SampledCorrelation(x, z)
        correlation.validate_support(config)
        result = estimate(correlation, config=config)
        result.update(window=scan.window, epoch=scan.epoch,
                      input_mode="sampled_linear_interpolation",
                      input_axis=args.axis, output_axis="physical-delay",
                      epoch_count=scan.metadata.get("epoch_count"),
                      reference="unit_ideal_triangle")
        rows.append(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive create protects previous diagnostic results.
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump({"input": str(args.input.resolve()), "results": rows}, stream,
                  indent=2, allow_nan=False)
    print(f"Wrote {len(rows)} observation estimates: {args.output}")


if __name__ == "__main__":
    main()
