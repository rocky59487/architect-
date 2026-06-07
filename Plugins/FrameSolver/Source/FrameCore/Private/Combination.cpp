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

static void envEnd(MemberEndForces& hi, MemberEndForces& lo, const MemberEndForces& v) {
    hi.N = std::max(hi.N, v.N); lo.N = std::min(lo.N, v.N);
    hi.Vy = std::max(hi.Vy, v.Vy); lo.Vy = std::min(lo.Vy, v.Vy);
    hi.Vz = std::max(hi.Vz, v.Vz); lo.Vz = std::min(lo.Vz, v.Vz);
    hi.T = std::max(hi.T, v.T); lo.T = std::min(lo.T, v.T);
    hi.My = std::max(hi.My, v.My); lo.My = std::min(lo.My, v.My);
    hi.Mz = std::max(hi.Mz, v.Mz); lo.Mz = std::min(lo.Mz, v.Mz);
}

static void envShell(ShellElementForces& hi, ShellElementForces& lo, const ShellElementForces& v) {
    hi.Mxx = std::max(hi.Mxx, v.Mxx); lo.Mxx = std::min(lo.Mxx, v.Mxx);
    hi.Myy = std::max(hi.Myy, v.Myy); lo.Myy = std::min(lo.Myy, v.Myy);
    hi.Mxy = std::max(hi.Mxy, v.Mxy); lo.Mxy = std::min(lo.Mxy, v.Mxy);
    hi.Qx  = std::max(hi.Qx,  v.Qx);  lo.Qx  = std::min(lo.Qx,  v.Qx);
    hi.Qy  = std::max(hi.Qy,  v.Qy);  lo.Qy  = std::min(lo.Qy,  v.Qy);
    hi.Nxx = std::max(hi.Nxx, v.Nxx); lo.Nxx = std::min(lo.Nxx, v.Nxx);
    hi.Nyy = std::max(hi.Nyy, v.Nyy); lo.Nyy = std::min(lo.Nyy, v.Nyy);
    hi.Nxy = std::max(hi.Nxy, v.Nxy); lo.Nxy = std::min(lo.Nxy, v.Nxy);
}

ResultEnvelope envelope(const std::vector<SolveResult>& cases) {
    ResultEnvelope E;
    if (cases.empty()) return E;
    const SolveResult& b = cases[0];
    E.uMax = b.u;        E.uMin = b.u;
    E.reactMax = b.reactions; E.reactMin = b.reactions;
    E.endIMax.resize(b.memberForces.size()); E.endIMin.resize(b.memberForces.size());
    E.endJMax.resize(b.memberForces.size()); E.endJMin.resize(b.memberForces.size());
    E.shellMax.resize(b.shellForces.size()); E.shellMin.resize(b.shellForces.size());
    for (size_t e = 0; e < b.memberForces.size(); ++e) {
        E.endIMax[e] = E.endIMin[e] = b.memberForces[e].endI;
        E.endJMax[e] = E.endJMin[e] = b.memberForces[e].endJ;
    }
    for (size_t s = 0; s < b.shellForces.size(); ++s) {
        E.shellMax[s] = E.shellMin[s] = b.shellForces[s];
    }
    for (size_t c = 1; c < cases.size(); ++c) {
        const SolveResult& R = cases[c];
        for (size_t k = 0; k < E.uMax.size() && k < R.u.size(); ++k) { E.uMax[k] = std::max(E.uMax[k], R.u[k]); E.uMin[k] = std::min(E.uMin[k], R.u[k]); }
        for (size_t k = 0; k < E.reactMax.size() && k < R.reactions.size(); ++k) { E.reactMax[k] = std::max(E.reactMax[k], R.reactions[k]); E.reactMin[k] = std::min(E.reactMin[k], R.reactions[k]); }
        for (size_t e = 0; e < E.endIMax.size() && e < R.memberForces.size(); ++e) {
            envEnd(E.endIMax[e], E.endIMin[e], R.memberForces[e].endI);
            envEnd(E.endJMax[e], E.endJMin[e], R.memberForces[e].endJ);
        }
        for (size_t s = 0; s < E.shellMax.size() && s < R.shellForces.size(); ++s) {
            envShell(E.shellMax[s], E.shellMin[s], R.shellForces[s]);
        }
    }
    return E;
}

}  // namespace frame
