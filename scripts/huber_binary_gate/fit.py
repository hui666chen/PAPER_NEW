"""Complex radial Huber, shared-ratio grid and per-epoch profiled gains.

Joint IRLS uses a complete grid search at each weighted-LS step; this is
not a guarantee of the global Huber optimum. No simulator labels are read.
"""
import numpy as np


def gain_residual(y, template, pixel_weights=None):
    w = np.ones(y.shape) if pixel_weights is None else pixel_weights
    gain = np.sum(w * y * template.conj(), axis=1) / np.maximum(
        np.sum(w * abs(template)**2, axis=1), 1e-30)
    return y - gain[:, None] * template


def huber(r, delta):
    a = abs(r)
    return np.where(a <= delta, a*a, 2*delta*a-delta*delta)


class Grid:
    def __init__(self, dictionary, direct, reflected, ratios):
        self.d = np.asarray(dictionary, complex)
        self.direct = np.asarray(direct, int)
        self.reflected = np.asarray(reflected, int)
        self.ratios = np.asarray(ratios, complex)
        self.gram = self.d.conj().T @ self.d

    def template(self, index):
        ratio_index, row = divmod(index, len(self.direct))
        return (self.d[:, self.direct[row]]
                + self.ratios[ratio_index] * self.d[:, self.reflected[row]])

    def search(self, y, epoch_weights, pixel_weights):
        # G[e,i,j] = D_i^H W_e D_j; no squared-condition Gram inverse.
        gram = (np.broadcast_to(self.gram, (len(y),)+self.gram.shape)
                if np.all(pixel_weights == 1)
                else self.d.conj().T[None] @ (pixel_weights[:, :, None]*self.d[None]))
        corr = (pixel_weights*y) @ self.d.conj()
        energy = np.sum(pixel_weights*abs(y)**2, axis=1)
        d, r = self.direct, self.reflected
        best_score, best_index = np.inf, -1
        for k, rho in enumerate(self.ratios):
            norm = (gram[:, d, d].real + abs(rho)**2*gram[:, r, r].real
                    + 2*(rho*gram[:, d, r]).real)
            rhs = corr[:, d] + rho.conjugate()*corr[:, r]
            scores = np.sum(epoch_weights[:, None] *
                            (energy[:, None] - abs(rhs)**2/np.maximum(norm, 1e-30)), axis=0)
            index = int(np.argmin(scores))
            if scores[index] < best_score:
                best_score, best_index = float(scores[index]), k*len(d)+index
        return best_index

    def fit(self, y, weights, delta=None, max_steps=6):
        pixel = np.ones(y.shape)
        index = self.search(y, weights, pixel)
        residual = gain_residual(y, self.template(index))
        ls_index = index
        if delta is None:
            return dict(index=index, ls_index=index, residual=residual,
                        converged=True, objective=[float(np.sum(weights[:, None]*abs(residual)**2))])
        history = [float(np.sum(weights[:, None]*huber(residual, delta)))]
        converged = False
        for _ in range(max_steps):
            pixel = np.minimum(1., delta/np.maximum(abs(residual), 1e-30))
            index = self.search(y, weights, pixel)
            residual = gain_residual(y, self.template(index), pixel)
            cost = float(np.sum(weights[:, None]*huber(residual, delta)))
            if cost > history[-1] + 1e-8*max(1., history[-1]):
                raise RuntimeError("IRLS objective increased")
            converged = abs(cost-history[-1]) <= 1e-5*max(history[-1], 1e-30)
            history.append(cost)
            if converged:
                break
        return dict(index=index, ls_index=ls_index, residual=residual,
                    converged=bool(converged), objective=history)


def evaluate(y, weights, template, delta=None):
    residual = gain_residual(y, template)
    if delta is not None:
        for _ in range(30):
            updated = gain_residual(y, template,
                                   np.minimum(1., delta/np.maximum(abs(residual), 1e-30)))
            if np.max(abs(updated-residual)) < 1e-9:
                residual = updated
                break
            residual = updated
    return float(np.sum(weights[:, None]*(abs(residual)**2 if delta is None
                                         else huber(residual, delta))))
