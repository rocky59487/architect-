#include "FrameCore/ResponseSpectrum.h"
#include "PreparedSystemImpl.h"

#include <cmath>
#include <algorithm>

namespace frame {

ResponseSpectrumResult solveResponseSpectrum(const PreparedSystem& prepared, const ModalResult& modes,
                                             const Spectrum& spectrum, int excDof,
                                             SpectrumCombo combo, real zeta) {
    ResponseSpectrumResult R;
    const PreparedSystem::Impl& S = *prepared.impl;
    if (S.singular) { R.singular = true; R.diagnostic = S.diagnostic; return R; }
    const int N = S.N;
    R.u.assign(static_cast<size_t>(std::max(0, N)), 0.0);
    if (modes.modes.empty()) { R.singular = true; R.diagnostic = "no modes for response spectrum"; return R; }

    // mass matrix + influence vector r (rigid-body unit motion along excDof)
    std::vector<Triplet> mtrips;
    for (const auto& el : S.elems) el->assembleMass(mtrips);
    SpMat M(N, N); M.setFromTriplets(mtrips.begin(), mtrips.end()); M.makeCompressed();

    VecX r = VecX::Zero(N);
    for (int node = 0; node * 6 < N; ++node) r(node * 6 + excDof) = 1.0;
    const VecX Mr = M * r;
    R.totalMass = r.dot(Mr);                              // total directional mass r^T M r

    const int nm = static_cast<int>(modes.modes.size());
    std::vector<VecX> umode(nm);                          // per-mode peak displacement
    std::vector<real> Vmode(nm);                          // per-mode base shear
    R.effMass.assign(static_cast<size_t>(nm), 0.0);

    for (int i = 0; i < nm; ++i) {
        VecX phi = VecX::Zero(N);
        for (int g = 0; g < N; ++g) phi(g) = modes.modes[i].shape[static_cast<size_t>(g)];
        const real Gamma = phi.dot(Mr);                  // phi^T M r (phi mass-normalized)
        R.effMass[i] = Gamma * Gamma;                    // effective modal mass
        const real w = modes.modes[i].omega;
        const real T = (w > 0) ? twoPi / w : 0;
        const real Sa = spectrum.at(T);
        const real Sd = (w > 0) ? Sa / (w * w) : 0;      // spectral displacement
        Vmode[i] = R.effMass[i] * Sa;                    // modal base shear = M* Sa
        umode[i] = (Gamma * Sd) * phi;                   // modal peak displacement contribution
    }

    auto rho = [&](int i, int j) -> real {
        if (combo == SpectrumCombo::SRSS) return (i == j) ? 1.0 : 0.0;
        const real wi = modes.modes[i].omega, wj = modes.modes[j].omega;
        if (wi <= 0 || wj <= 0) return (i == j) ? 1.0 : 0.0;
        const real b = wi / wj;
        const real num = 8.0 * zeta * zeta * (1.0 + b) * std::pow(b, 1.5);
        const real den = (1.0 - b * b) * (1.0 - b * b) + 4.0 * zeta * zeta * b * (1.0 + b) * (1.0 + b);
        return (den > 0) ? num / den : ((i == j) ? 1.0 : 0.0);
    };

    for (int g = 0; g < N; ++g) {
        real s = 0;
        for (int i = 0; i < nm; ++i)
            for (int j = 0; j < nm; ++j) s += rho(i, j) * umode[i](g) * umode[j](g);
        R.u[static_cast<size_t>(g)] = (s > 0) ? std::sqrt(s) : 0;
    }
    real sv = 0;
    for (int i = 0; i < nm; ++i)
        for (int j = 0; j < nm; ++j) sv += rho(i, j) * Vmode[i] * Vmode[j];
    R.baseShear = (sv > 0) ? std::sqrt(sv) : 0;
    return R;
}

}  // namespace frame
