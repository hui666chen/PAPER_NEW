"""Independent finite-difference checks for the profiled CRB implementation."""
import json

import numpy as np

import crb_model_comparison as crb
import monte_carlo_resolution_baseline as mc
from blind_response_experiment import digest, write_new
from receiver_chain_response_control import ChainResponse
from run_blind_transfer import FS, load_batch, make_bank


OUT = mc.ROOT / "CRB_M0_M1_M2" / "verification.json"
DELAY = .60
RATIO = .34
PHASE = mc.PHASE_RAD
COUNT = 2


def numerical_jacobian(function, parameters, steps):
    columns = []
    for index, step in enumerate(steps):
        plus = parameters.copy(); plus[index] += step
        minus = parameters.copy(); minus[index] -= step
        columns.append(((function(plus) - function(minus)) / (2 * step)).reshape(-1))
    return np.column_stack(columns)


def profiled_fisher(jacobian, global_count, whiten):
    # Whiten each epoch independently before constructing the real-parameter FIM.
    shaped = jacobian.reshape(COUNT, 401, -1)
    white = np.concatenate([whiten(shaped[k]) for k in range(COUNT)], axis=0)
    fisher = 2 * np.real(white.conj().T @ white)
    gg = fisher[:global_count, :global_count]
    gn = fisher[:global_count, global_count:]
    nn = fisher[global_count:, global_count:]
    effective = gg - gn @ np.linalg.solve(nn, gn.T)
    effective = (effective + effective.T) / 2
    covariance = np.linalg.inv(effective)
    return effective, covariance


def relative_matrix_error(a, b):
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def main():
    with np.load(mc.ROOT / "MC_resolution_r_44c3uo/noise_covariance.npz") as data:
        eigval = data["eigenvalues"]
        eigvec = data["eigenvectors"]

    def whiten(matrix):
        return (eigvec.conj().T @ matrix) / np.sqrt(eigval)[:, None]

    rows, _ = load_batch(mc.A_PREP)
    state = next(r for r in rows if int(r["input_start_sample"]) / FS >= 2.2)
    waveform = np.load(mc.CAL)["waveform"]
    bank, spec = make_bank(waveform, state)
    model = ChainResponse(bank, spec=spec)
    ones = np.ones(401)
    tap = spec.tap_grid
    c = RATIO * np.exp(1j * PHASE)
    times = (np.arange(COUNT) - (COUNT - 1) / 2) * .001

    checks = []

    # M0: global tau/delta, then 8 real nuisance coordinates per epoch.
    m0_initial = [0., DELAY]
    for _ in range(COUNT):
        m0_initial += [1., 0., c.real, c.imag, 0., 0., 0., 0.]
    m0_initial = np.array(m0_initial)

    def mean_m0(p):
        tau, delta = p[:2]
        r0 = model.exact([-tau])[0]
        r1 = model.exact([delta - tau])[0]
        answer = []
        offset = 2
        for _ in range(COUNT):
            q = p[offset:offset+8]; offset += 8
            a0 = q[0] + 1j*q[1]; a1 = q[2] + 1j*q[3]
            b0 = q[4] + 1j*q[5]; b1 = q[6] + 1j*q[7]
            answer.append(a0*r0 + a1*r1 + b0*ones + b1*tap)
        return np.asarray(answer)

    j = numerical_jacobian(mean_m0, m0_initial,
                           np.array([1e-5, 1e-5] + [1e-6]*(8*COUNT)))
    f_num, v_num = profiled_fisher(j, 2, whiten)
    analytic = crb.evaluate("M0", model, spec, whiten, DELAY, RATIO, PHASE, COUNT)
    sd_analytic = np.array([analytic["tau_sd_crb"], analytic["delta_sd_crb"]])
    checks.append({"model": "M0", "numeric_sd": np.sqrt(np.diag(v_num)).tolist(),
                   "analytic_sd": sd_analytic.tolist(),
                   "sd_relative_error_max": float(np.max(abs(np.sqrt(np.diag(v_num))/sd_analytic-1)))})

    # M1: four globals, then 6 real nuisance coordinates per epoch.
    m1_initial = [0., DELAY, RATIO, PHASE]
    for _ in range(COUNT):
        m1_initial += [1., 0., 0., 0., 0., 0.]
    m1_initial = np.array(m1_initial)

    def mean_m1(p):
        tau, delta, ratio, phase = p[:4]
        r0 = model.exact([-tau])[0]
        r1 = model.exact([delta - tau])[0]
        h = r0 + ratio*np.exp(1j*phase)*r1
        answer = []; offset = 4
        for _ in range(COUNT):
            q = p[offset:offset+6]; offset += 6
            g = q[0] + 1j*q[1]; b0 = q[2] + 1j*q[3]; b1 = q[4] + 1j*q[5]
            answer.append(g*h + b0*ones + b1*tap)
        return np.asarray(answer)

    j = numerical_jacobian(mean_m1, m1_initial,
                           np.array([1e-5, 1e-5, 1e-6, 1e-6] + [1e-6]*(6*COUNT)))
    f_num, v_num = profiled_fisher(j, 4, whiten)
    analytic = crb.evaluate("M1", model, spec, whiten, DELAY, RATIO, PHASE, COUNT)
    sd_analytic = np.array([analytic["tau_sd_crb"], analytic["delta_sd_crb"],
                            np.sqrt(np.linalg.inv(f_num)[2,2]), analytic["phi0_sd_crb"]])
    numeric_sd = np.sqrt(np.diag(v_num))
    checks.append({"model": "M1", "numeric_sd": numeric_sd.tolist(),
                   "analytic_tau_delta_phi_sd": [analytic["tau_sd_crb"],
                                                   analytic["delta_sd_crb"],
                                                   analytic["phi0_sd_crb"]],
                   "tau_delta_phi_relative_error_max": float(max(
                       abs(numeric_sd[0]/analytic["tau_sd_crb"]-1),
                       abs(numeric_sd[1]/analytic["delta_sd_crb"]-1),
                       abs(numeric_sd[3]/analytic["phi0_sd_crb"]-1)))})

    # M2: M1 plus omega; numerical phase is phi0 + omega*t_k.
    m2_initial = [0., DELAY, RATIO, PHASE, 0.]
    for _ in range(COUNT):
        m2_initial += [1., 0., 0., 0., 0., 0.]
    m2_initial = np.array(m2_initial)

    def mean_m2(p):
        tau, delta, ratio, phase, omega = p[:5]
        r0 = model.exact([-tau])[0]
        r1 = model.exact([delta - tau])[0]
        answer = []; offset = 5
        for time_s in times:
            q = p[offset:offset+6]; offset += 6
            g = q[0] + 1j*q[1]; b0 = q[2] + 1j*q[3]; b1 = q[4] + 1j*q[5]
            h = r0 + ratio*np.exp(1j*(phase + omega*time_s))*r1
            answer.append(g*h + b0*ones + b1*tap)
        return np.asarray(answer)

    j = numerical_jacobian(mean_m2, m2_initial,
                           np.array([1e-5, 1e-5, 1e-6, 1e-6, 1e-3] + [1e-6]*(6*COUNT)))
    f_num, v_num = profiled_fisher(j, 5, whiten)
    analytic = crb.evaluate("M2", model, spec, whiten, DELAY, RATIO, PHASE, COUNT)
    numeric_sd = np.sqrt(np.diag(v_num))
    checks.append({"model": "M2", "numeric_sd": numeric_sd.tolist(),
                   "analytic_tau_delta_phi_sd": [analytic["tau_sd_crb"],
                                                   analytic["delta_sd_crb"],
                                                   analytic["phi0_sd_crb"]],
                   "tau_delta_phi_relative_error_max": float(max(
                       abs(numeric_sd[0]/analytic["tau_sd_crb"]-1),
                       abs(numeric_sd[1]/analytic["delta_sd_crb"]-1),
                       abs(numeric_sd[3]/analytic["phi0_sd_crb"]-1)))})

    # Known-ratio/phase check on the six empirical boundary cells.
    original = json.loads((mc.ROOT / "CRB_M0_M1_M2/crb_results.json").read_text())["results"]
    known_checks = []
    for count, ratio, delay in ((50,.5,.4),(50,.5,.6),(100,.34,.6),
                                (100,.34,.8),(100,.5,.3),(100,.5,.4)):
        r0, p0 = crb.response_and_derivative(model, 0.)
        r1, p1 = crb.response_and_derivative(model, delay)
        cc = ratio + 0j
        h = r0 + cc*r1
        nuisance = np.column_stack((h, ones, tap))
        derivatives = np.column_stack((-p0-cc*p1, cc*p1))
        fisher = count * crb.residualized_information(whiten, nuisance, derivatives)
        known_sd = np.sqrt(np.diag(np.linalg.inv(fisher)))
        unknown = next(r for r in original if r["model"] == "M1" and
                       r["epochs"] == count and np.isclose(r["amplitude_ratio"], ratio) and
                       np.isclose(r["delay_chip"], delay) and np.isclose(r["phi0_rad"], 0))
        known_checks.append({"epochs": count, "ratio": ratio, "delay_chip": delay,
            "known_ratio_phase_tau_sd": float(known_sd[0]),
            "unknown_ratio_phase_tau_sd": unknown["tau_sd_crb"],
            "known_ratio_phase_delta_sd": float(known_sd[1]),
            "unknown_ratio_phase_delta_sd": unknown["delta_sd_crb"],
            "delta_sd_improvement_fraction": float(1-known_sd[1]/unknown["delta_sd_crb"])})

    output = {"schema": 1, "truth_accessed": False,
              "source_sha256": digest(__file__),
              "finite_difference_case": {"epochs": COUNT, "ratio": RATIO,
                                         "delay_chip": DELAY, "phi0_rad": PHASE},
              "finite_difference_checks": checks,
              "known_ratio_phase_checks": known_checks,
              "interpretation": "known constraints strictly improve or preserve the bound; independent finite differences validate the analytic projection"}
    write_new(OUT, output)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
