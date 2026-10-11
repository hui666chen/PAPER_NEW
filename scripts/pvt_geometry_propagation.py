"""Propagate frozen PRN14 delay errors through SignalSim GPS geometry and WLS."""

import argparse
import csv
import hashlib
import io
import json
import subprocess
from pathlib import Path

import numpy as np

from blind_response_experiment import write_new
from run_blind_transfer import ROOT


HERE = Path(__file__).resolve().parent
SIGNALSIM = Path("<LOCAL_USER>/chh/GNSS/SignalSim-main")
RINEX = SIGNALSIM / "EphData/BRDC00IGS_R_20211700000_01D_MN.rnx"
SOURCE = ROOT / "crossfit_basin_pilot_3x50/results.json"
OUT = ROOT / "pvt_geometry_propagation_3x50"
CPP = HERE / "extract_signalsim_gps_geometry.cpp"
PROTOCOL = HERE / "PVT_GEOMETRY_PROPAGATION_PROTOCOL.md"
COUNTERFACTUAL_PROTOCOL = HERE / "PVT_ELEVATION_COUNTERFACTUAL_PROTOCOL.md"
METRES_PER_CHIP = 299_792_458.0 / 1_023_000.0
THRESHOLD = .20
SIGMA_ZENITH = (0.5, 1.0, 2.0)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extract_geometry():
    binary = Path("<TMP_DIR>/extract_signalsim_gps_geometry_pvt")
    sources = [
        "NavData.cpp", "Rinex.cpp", "GnssTime.cpp", "Coordinate.cpp",
        "Almanac.cpp", "FastMath.cpp", "MessageOutput.cpp",
    ]
    command = [
        "g++", "-std=c++11", "-O2", f"-I{SIGNALSIM / 'inc'}", str(CPP),
        *[str(SIGNALSIM / "src" / name) for name in sources],
        "-o", str(binary),
    ]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    extracted = subprocess.run(
        [str(binary), str(RINEX)], check=True, capture_output=True, text=True
    ).stdout
    rows = []
    for row in csv.DictReader(io.StringIO(extracted)):
        rows.append(dict(
            prn=int(row["prn"]), elevation_deg=float(row["elevation_deg"]),
            azimuth_deg=float(row["azimuth_deg"]),
            los_ecef=[float(row["los_x"]), float(row["los_y"]),
                      float(row["los_z"])],
        ))
    if 14 not in {row["prn"] for row in rows}:
        raise RuntimeError("PRN14 is not visible in extracted SignalSim geometry")
    return rows, completed.stderr


def design(geometry, sigma_zenith):
    elevation = np.radians([row["elevation_deg"] for row in geometry])
    azimuth = np.radians([row["azimuth_deg"] for row in geometry])
    east = np.cos(elevation) * np.sin(azimuth)
    north = np.cos(elevation) * np.cos(azimuth)
    up = np.sin(elevation)
    h = np.column_stack((-east, -north, -up, np.ones(len(geometry))))
    sigma = sigma_zenith / np.sin(elevation)
    weight = np.diag(1 / sigma ** 2)
    normal = h.T @ weight @ h
    mapping = np.linalg.solve(normal, h.T @ weight)
    unweighted_q = np.linalg.inv(h.T @ h)
    return dict(
        h=h, sigma=sigma, mapping=mapping,
        gdop=float(np.sqrt(np.trace(unweighted_q))),
        pdop=float(np.sqrt(np.trace(unweighted_q[:3, :3]))),
        hdop=float(np.sqrt(unweighted_q[0, 0] + unweighted_q[1, 1])),
        vdop=float(np.sqrt(unweighted_q[2, 2])),
        condition=float(np.linalg.cond(normal)),
    )


def expected_metrics(errors_m, mapping, sigma, affected_index):
    errors_m = np.asarray(errors_m, float)
    other = np.arange(mapping.shape[1]) != affected_index
    covariance_other = (
        mapping[:, other] @ np.diag(sigma[other] ** 2)
        @ mapping[:, other].T
    )
    affected_response = mapping[:, affected_index, None] * errors_m[None, :]
    component_mse = np.mean(affected_response ** 2, axis=1) + np.diag(covariance_other)
    return dict(
        east_rmse_m=float(np.sqrt(component_mse[0])),
        north_rmse_m=float(np.sqrt(component_mse[1])),
        up_rmse_m=float(np.sqrt(component_mse[2])),
        horizontal_rmse_m=float(np.sqrt(component_mse[0] + component_mse[1])),
        vertical_rmse_m=float(np.sqrt(component_mse[2])),
        position_3d_rmse_m=float(np.sqrt(component_mse[:3].sum())),
        clock_rmse_m=float(np.sqrt(component_mse[3])),
        prn14_error_rmse_m=float(np.sqrt(np.mean(errors_m ** 2))),
        analytic_other_noise_covariance=covariance_other.tolist(),
    )


def trial_errors(source):
    cells = []
    pooled = {name: [] for name in ("gated", "single", "two")}
    for cell in source["cells"]:
        two = np.asarray([row["full_tau"] for row in cell["records"]])
        single = np.asarray([row["single_tau"] for row in cell["records"]])
        risk = np.asarray([row["P_cross"] for row in cell["records"]])
        gated = np.where(risk <= THRESHOLD + 1e-12, two, single)
        values = dict(
            gated=gated * METRES_PER_CHIP,
            single=single * METRES_PER_CHIP,
            two=two * METRES_PER_CHIP,
        )
        for name in pooled:
            pooled[name].append(values[name])
        cells.append(dict(
            delay=cell["delay"], amp=cell["amp"], N_epoch=cell["N_epoch"],
            errors=values,
        ))
    return cells, {name: np.concatenate(parts) for name, parts in pooled.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument(
        "--prn14-elevation-deg", type=float,
        help="counterfactually replace only PRN14 elevation; keep azimuth fixed",
    )
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    geometry, compiler_stderr = extract_geometry()
    original_prn14_elevation = next(
        row["elevation_deg"] for row in geometry if row["prn"] == 14
    )
    if args.prn14_elevation_deg is not None:
        if not 3.0 <= args.prn14_elevation_deg < 90.0:
            raise ValueError("counterfactual PRN14 elevation must be in [3, 90) deg")
        for row in geometry:
            if row["prn"] == 14:
                row["original_elevation_deg"] = row["elevation_deg"]
                row["elevation_deg"] = float(args.prn14_elevation_deg)
                row["counterfactual_elevation"] = True
    affected_index = next(
        index for index, row in enumerate(geometry) if row["prn"] == 14
    )
    source = json.loads(SOURCE.read_text())
    cells, pooled = trial_errors(source)
    sensitivity = []
    for sigma_zenith in SIGMA_ZENITH:
        model = design(geometry, sigma_zenith)
        methods = {
            name: expected_metrics(errors, model["mapping"], model["sigma"],
                                   affected_index)
            for name, errors in pooled.items()
        }
        per_cell = []
        for cell in cells:
            per_cell.append(dict(
                delay=cell["delay"], amp=cell["amp"], N_epoch=cell["N_epoch"],
                methods={
                    name: expected_metrics(errors, model["mapping"], model["sigma"],
                                           affected_index)
                    for name, errors in cell["errors"].items()
                },
            ))
        comparisons = {}
        for endpoint in ("horizontal_rmse_m", "vertical_rmse_m", "position_3d_rmse_m"):
            gate = methods["gated"][endpoint]
            single = methods["single"][endpoint]
            comparisons[endpoint] = dict(
                gate_minus_single_m=gate - single,
                relative_improvement_percent=100 * (single - gate) / single,
            )
        sensitivity.append(dict(
            sigma_zenith_m=sigma_zenith,
            dop={key: model[key] for key in ("gdop", "pdop", "hdop", "vdop", "condition")},
            prn14_mapping_ENU_clock=model["mapping"][:, affected_index].tolist(),
            methods=methods, comparisons=comparisons, per_cell=per_cell,
        ))

    payload = dict(
        schema=1, status="post-hoc frozen-error geometry propagation",
        protocol=str(PROTOCOL), protocol_sha256=sha256(PROTOCOL),
        counterfactual_protocol=(
            str(COUNTERFACTUAL_PROTOCOL)
            if args.prn14_elevation_deg is not None else None
        ),
        counterfactual_protocol_sha256=(
            sha256(COUNTERFACTUAL_PROTOCOL)
            if args.prn14_elevation_deg is not None else None
        ),
        source=str(SOURCE), source_sha256=sha256(SOURCE),
        geometry_source=dict(
            extractor=str(CPP), extractor_sha256=sha256(CPP),
            signalsim_root=str(SIGNALSIM), rinex=str(RINEX), rinex_sha256=sha256(RINEX),
            utc="2021-06-19T10:05:40Z",
            receiver_lla_deg_m=[37.352721, -121.915773, 20.0],
            elevation_mask_deg=3.0,
            prn14_original_elevation_deg=original_prn14_elevation,
            prn14_counterfactual_elevation_deg=args.prn14_elevation_deg,
            prn14_azimuth_held_fixed=True,
            compiler_warnings_present=bool(compiler_stderr.strip()),
        ),
        geometry=geometry, affected_prn=14, trials=150,
        threshold=THRESHOLD, metres_per_chip=METRES_PER_CHIP,
        baseline_noise_model="independent sigma_zenith/sin(elevation); analytic expectation",
        sensitivity=sensitivity,
        limitations=[
            "linearized single-epoch ENU plus clock WLS, not nonlinear receiver PVT",
            "only PRN14 carries saved multipath-estimator error",
            "other-satellite errors are an ideal independent Gaussian baseline",
            "same static geometry is used for all 150 channel-error trials",
        ],
    )
    write_new(args.output / "summary.json", payload)
    print(json.dumps({
        "satellites": [(r["prn"], r["elevation_deg"], r["azimuth_deg"])
                       for r in geometry],
        "primary_sigma_zenith_1m": sensitivity[1],
    }, indent=2))
    print(f"COMPLETE={args.output}")


if __name__ == "__main__":
    main()
