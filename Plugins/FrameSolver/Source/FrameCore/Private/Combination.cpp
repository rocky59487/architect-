#include "FrameCore/Combination.h"
#include <algorithm>

namespace frame {

static void addScaledEnd(MemberEndForces& a, const MemberEndForces& b, real f) {
    a.N += f * b.N; a.Vy += f * b.Vy; a.Vz += f * b.Vz;
    a.T += f * b.T; a.My += f * b.My; a.Mz += f * b.Mz;
}

SolveResult combine(const std::vector<SolveResult>& cases, const std::vector<real>& factors) {
    SolveResult out;
    if (cases.empty()) return out;

    const SolveResult& base = cases[0];
    out.u.assign(base.u.size(), 0.0);
    out.reactions.assign(base.reactions.size(), 0.0);
    out.memberForces.resize(base.memberForces.size());
    out.shellForces.resize(base.shellForces.size());
    for (size_t e = 0; e < out.memberForces.size(); ++e) out.memberForces[e].member = base.memberForces[e].member;
    for (size_t s = 0; s < out.shellForces.size(); ++s)  out.shellForces[s].shell  = base.shellForces[s].shell;

    const size_t nC = std::min(cases.size(), factors.size());
    bool sing = false;
    for (size_t c = 0; c < nC; ++c) {
        const real f = factors[c];
        const SolveResult& R = cases[c];
        if (R.singular) sing = true;
        for (size_t k = 0; k < out.u.size() && k < R.u.size(); ++k) out.u[k] += f * R.u[k];
        for (size_t k = 0; k < out.reactions.size() && k < R.reactions.size(); ++k) out.reactions[k] += f * R.reactions[k];
        for (size_t e = 0; e < out.memberForces.size() && e < R.memberForces.size(); ++e) {
            addScaledEnd(out.memberForces[e].endI, R.memberForces[e].endI, f);
            addScaledEnd(out.memberForces[e].endJ, R.memberForces[e].endJ, f);
        }
        for (size_t s = 0; s < out.shellForces.size() && s < R.shellForces.size(); ++s) {
            ShellElementForces& o = out.shellForces[s];
            const ShellElementForces& r = R.shellForces[s];
            o.Mxx += f * r.Mxx; o.Myy += f * r.Myy; o.Mxy += f * r.Mxy;
            o.Qx  += f * r.Qx;  o.Qy  += f * r.Qy;
            o.Nxx += f * r.Nxx; o.Nyy += f * r.Nyy; o.Nxy += f * r.Nxy;
        }
    }
    out.singular = sing;
    return out;
}

}  // namespace frame
