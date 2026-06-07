#include "FrameCore/ModalDynamics.h"
#include "PreparedSystemImpl.h"

#include <cmath>

namespace frame {

ModalTimeHistory solveModalStepResponse(const PreparedSystem& prepared, const FrameModel& model,
                                        const ModalResult& modes, const ModalDynamicsOptions& opts) {
    ModalTimeHistory R;
    R.dt = opts.dt;
    const PreparedSystem::Impl& S = *prepared.impl;
    if (S.singular) { R.singular = true; R.diagnostic = S.diagnostic; return R; }
    const int N = S.N, nm = static_cast<int>(modes.modes.size());
    if (nm == 0 || opts.nSteps <= 0 || opts.dt <= 0) { R.singular = true; R.diagnostic = "bad modal/step input"; return R; }

    // global step force vector
    VecX F = VecX::Zero(N);
    for (const auto& nl : model.nodalLoads) {
        const int ni = model.nodeIndex(nl.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) F(gdof(ni, d)) += nl.comp[d];
    }

    // per-mode quantities (modal mass = 1; modal force f_n = phi_n^T F; c = 2 zeta w; k = w^2)
    std::vector<VecX> phi(nm);
    std::vector<real> w(nm), c(nm), k(nm), p(nm), q(nm, 0.0), qd(nm, 0.0), qdd(nm, 0.0);
    for (int i = 0; i < nm; ++i) {
        VecX ph = VecX::Zero(N);
        for (int g = 0; g < N; ++g) ph(g) = modes.modes[i].shape[static_cast<size_t>(g)];
        phi[i] = ph;
        w[i]   = modes.modes[i].omega;
        c[i]   = 2.0 * opts.zeta * w[i];
        k[i]   = w[i] * w[i];
        p[i]   = ph.dot(F);
        qdd[i] = p[i];                       // m*qdd0 = p - c*0 - k*0
    }

    const real dt = opts.dt, beta = 0.25, gamma = 0.5;
    R.u.reserve(static_cast<size_t>(opts.nSteps) + 1);
    auto record = [&]() {
        std::vector<real> u(static_cast<size_t>(N), 0.0);
        for (int i = 0; i < nm; ++i)
            for (int g = 0; g < N; ++g) u[static_cast<size_t>(g)] += phi[i](g) * q[i];
        R.u.push_back(std::move(u));
    };
    record();   // t = 0

    for (int step = 0; step < opts.nSteps; ++step) {
        for (int i = 0; i < nm; ++i) {
            // Newmark average-acceleration update for m*qdd + c*qd + k*q = p (m = 1, p constant)
            const real keff = 1.0 / (beta * dt * dt) + gamma * c[i] / (beta * dt) + k[i];
            const real rhs  = p[i]
                + (1.0 / (beta * dt * dt)) * q[i] + (1.0 / (beta * dt)) * qd[i] + (1.0 / (2.0 * beta) - 1.0) * qdd[i]
                + c[i] * ((gamma / (beta * dt)) * q[i] + (gamma / beta - 1.0) * qd[i] + dt * (gamma / (2.0 * beta) - 1.0) * qdd[i]);
            const real qn   = rhs / keff;
            const real qddn = (1.0 / (beta * dt * dt)) * (qn - q[i]) - (1.0 / (beta * dt)) * qd[i] - (1.0 / (2.0 * beta) - 1.0) * qdd[i];
            const real qdn  = qd[i] + dt * ((1.0 - gamma) * qdd[i] + gamma * qddn);
            q[i] = qn; qd[i] = qdn; qdd[i] = qddn;
        }
        record();
    }
    return R;
}

}  // namespace frame
