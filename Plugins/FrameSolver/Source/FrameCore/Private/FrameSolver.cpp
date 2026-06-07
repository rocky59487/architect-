#include "FrameCore/FrameSolver.h"
#include "FrameEigen.h"
#include "IElement.h"
#include "BeamColumnElement.h"
#include "MITC4ShellElement.h"

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>

namespace frame {

SolveResult solve(const FrameModel& model, const SolveOptions& opts) {
    SolveResult R;

    std::string why;
    if (!model.validate(why)) {
        R.singular = true;
        R.diagnostic = "invalid model: " + why;
        const size_t n = (size_t)std::max(0, model.dofCount());
        R.u.assign(n, 0.0);
        R.reactions.assign(n, 0.0);
        return R;
    }

    const int N = model.dofCount();
    R.u.assign((size_t)N, 0.0);
    R.reactions.assign((size_t)N, 0.0);

    // ---- build the element list (seam #1). The solver no longer hard-codes the
    // 12x12 beam: each member becomes a BeamColumnElement, each shell facet a
    // MITC4ShellElement; both slot in behind the same IElement interface. ----
    std::vector<std::unique_ptr<IElement>> elems;
    elems.reserve(model.members.size() + model.shells.size());
    for (size_t e = 0; e < model.members.size(); ++e)
        elems.push_back(std::make_unique<BeamColumnElement>((int)e));
    for (size_t s = 0; s < model.shells.size(); ++s)
        elems.push_back(std::make_unique<MITC4ShellElement>((int)s));

    // ---- prepare each element: local stiffness + fixed-end forces + (optional)
    // member-end release condensation. An ill-posed release set (singular released
    // sub-block) flags the model singular with a clear diagnostic. ----
    bool prepSingular = false;
    for (auto& el : elems) {
        std::string ewhy;
        if (!el->prepare(model, opts, ewhy)) {
            prepSingular = true;
            R.singular   = true;
            R.diagnostic = ewhy;   // last failing element wins (matches prior behavior)
        }
    }
    if (prepSingular) {            // do not trust a solve with an ill-posed release set
        R.u.assign((size_t)N, 0.0);
        R.reactions.assign((size_t)N, 0.0);
        return R;
    }

    // ---- assemble global stiffness K from each element's (possibly condensed) kl ----
    std::vector<Triplet> trips;
    size_t tripCap = 0;
    for (const auto& el : elems) { const size_t d = (size_t)el->localDof(); tripCap += d * d; }
    trips.reserve(tripCap);
    for (const auto& el : elems) el->assemble(trips);

    // ---- global load vector F: nodal loads + equivalent nodal loads from each element ----
    VecX F = VecX::Zero(N);
    for (const auto& nl : model.nodalLoads) {
        const int ni = model.nodeIndex(nl.node);
        if (ni < 0) continue;
        for (int d = 0; d < 6; ++d) F(gdof(ni, d)) += nl.comp[d];
    }
    for (const auto& el : elems) el->addEquivalentNodalLoads(F);

    // ---- assemble sparse K ----
    SpMat K(N, N);
    K.setFromTriplets(trips.begin(), trips.end());
    K.makeCompressed();

    // ---- free-DOF map ----
    std::vector<int>  fmap(N, -1);
    std::vector<real> presc(N, 0.0);
    int nf = 0;
    for (size_t k = 0; k < model.nodes.size(); ++k)
        for (int d = 0; d < 6; ++d) {
            const int g = gdof((int)k, d);
            if (model.nodes[k].fixed[d]) presc[g] = model.nodes[k].prescribed[d];
            else fmap[g] = nf++;
        }
    if (nf == 0) { R.singular = true; R.diagnostic = "fully constrained (no free DOF)"; return R; }

    // ---- reduced system  K_ff u_f = F_f - K_fc u_c ----
    std::vector<Triplet> ftrips;
    VecX Ff = VecX::Zero(nf);
    for (int c = 0; c < N; ++c) {
        for (SpMat::InnerIterator it(K, c); it; ++it) {
            const int r = it.row();
            if (fmap[r] < 0) continue;                 // constrained row -> reactions later
            const real v = it.value();
            if (fmap[c] >= 0) ftrips.emplace_back(fmap[r], fmap[c], v);
            else if (presc[c] != 0.0) Ff(fmap[r]) -= v * presc[c];
        }
    }
    for (int g = 0; g < N; ++g) if (fmap[g] >= 0) Ff(fmap[g]) += F(g);

    SpMat Kff(nf, nf);
    Kff.setFromTriplets(ftrips.begin(), ftrips.end());
    Kff.makeCompressed();

    // ---- factor + mechanism detection (from the factorization, NOT connectivity) ----
    LDLTSolver ldlt;
    ldlt.compute(Kff);

    bool sing = false;
    std::string diag;
    if (ldlt.info() != Eigen::Success) {
        sing = true;
        diag = "LDLT factorization failed (info != Success): rank-deficient stiffness / mechanism";
    } else {
        const VecX D = ldlt.vectorD();
        real maxAbs = 0;
        for (int i = 0; i < D.size(); ++i) maxAbs = std::max(maxAbs, std::abs(D(i)));
        const real tol = opts.pivotTol * std::max(real(1), maxAbs);   // pivotTol defaults to 1e-12
        for (int i = 0; i < D.size(); ++i) {
            if (std::abs(D(i)) < tol) {
                sing = true;
                diag = "near-zero pivot in LDLT D (mechanism) at reduced dof #" + std::to_string(i);
                break;
            }
            if (D(i) < 0) {
                sing = true;
                diag = "negative pivot in LDLT D (indefinite stiffness / mechanism) at reduced dof #" + std::to_string(i);
                break;
            }
        }
    }

    VecX uf = VecX::Zero(nf);
    if (!sing) {
        uf = ldlt.solve(Ff);
        if (ldlt.info() != Eigen::Success || !uf.allFinite()) {
            sing = true;
            diag = "solve produced non-finite displacements (mechanism)";
        }
    }

    R.singular = sing;
    R.diagnostic = diag;

    // scatter to full u (constrained dofs = prescribed)
    VecX u = VecX::Zero(N);
    for (int g = 0; g < N; ++g) u(g) = (fmap[g] >= 0 && !sing) ? uf(fmap[g]) : presc[g];
    for (int g = 0; g < N; ++g) R.u[(size_t)g] = u(g);

    if (sing) return R;   // do not trust downstream forces on a mechanism

    // reactions  R = K u - F  (nonzero at constrained DOF)
    const VecX Rv = K * u - F;
    for (int g = 0; g < N; ++g) R.reactions[(size_t)g] = Rv(g);

    // member/shell internal forces (local): each element recovers its own results
    // (beams -> R.memberForces[memberIndex], shells -> R.shellForces[shellIndex]).
    R.memberForces.resize(model.members.size());
    R.shellForces.resize(model.shells.size());
    for (const auto& el : elems) el->recover(u, R);

    return R;
}

} // namespace frame
