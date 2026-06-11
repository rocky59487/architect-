#include "frame_cli_core.h"

#include "FrameCore/FrameSolver.h"
#include "FrameCore/ModalAnalysis.h"
#include "FrameCore/PDeltaAnalysis.h"
#include "FrameCore/TensionOnly.h"
#include "FrameCore/SizeOpt.h"
#include "FrameCore/DynamicCollapse.h"
#include "FrameCore/MemberGeometry.h"

#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdarg>
#include <cmath>
#include <algorithm>

using namespace frame;

#ifndef FRAMECORE_BUILD_SHA
#define FRAMECORE_BUILD_SHA "unknown"   // overridden by the build script via /D (git short SHA)
#endif

namespace {

// Append a printf-formatted line to a string (every protocol line fits comfortably in 512 chars:
// 12 doubles * ~20 chars). Output goes to a buffer, not stdout, so the DLL can return it.
void appendf(std::string& s, const char* fmt, ...) {
    char buf[512];
    va_list ap; va_start(ap, fmt);
    int n = std::vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    if (n > 0) s.append(buf, (size_t)std::min(n, (int)sizeof(buf) - 1));
}

struct RawMat { real E, G, rho, nu; bool hasCap; real capC, capT, capS; };
struct RawSec { real A, Iy, Iz, J, cy, cz, Asy, Asz; };
struct RawNode { int id; real x, y, z; int f[6]; real p[6]; };
struct RawMem { int id, i, j, mat, sec; real rx, ry, rz; int active; int tonly; };
struct RawShell { int id, n[4], mat; real t; int active; };
struct RawNL { int node; real c[6]; };
struct RawUDL { int member; real wx, wy, wz; };
struct RawSP { int shell; real p; };
struct RawHinge { int member, dof; real Mp; };

struct Block {
    std::vector<RawMat> mats; std::vector<RawSec> secs;
    std::vector<RawNode> nodes; std::vector<RawMem> mems; std::vector<RawShell> shes;
    std::vector<RawNL> nls; std::vector<RawUDL> udls; std::vector<RawSP> sps;
    std::vector<RawHinge> hins;
    SolveOptions opt;
    int nModes = 0;
    int pdelta = -1;
    std::string analysis;
    int  toMaxIter = 32, toAllowReact = 1;
    real soAmin = 0; int soMaxIter = 40; real soDcTol = 1e-8;
    real dcDt = 1e-3, dcMaxTime = 0.5;  std::vector<int> dcRemovals;
    bool empty = true;
};

void parseLine(Block& b, const std::string& tag, std::istringstream& ss) {
    b.empty = false;
    if (tag == "MAT") {
        RawMat m{}; m.nu = 0; m.hasCap = false; ss >> m.E >> m.G >> m.rho;
        real c; if (ss >> c) { m.hasCap = true; m.capC = c; m.capT = (ss >> c) ? c : m.capC; m.capS = (ss >> c) ? c : real(0.6) * m.capC; }
        b.mats.push_back(m);
    } else if (tag == "SMAT") {
        RawMat m{}; m.rho = 0; m.hasCap = false; ss >> m.E >> m.nu >> m.G;
        real c; if (ss >> c) { m.hasCap = true; m.capC = c; m.capT = (ss >> c) ? c : m.capC; m.capS = (ss >> c) ? c : real(0.6) * m.capC; }
        b.mats.push_back(m);
    } else if (tag == "SEC") { RawSec s{}; ss >> s.A >> s.Iy >> s.Iz >> s.J >> s.cy >> s.cz >> s.Asy >> s.Asz; b.secs.push_back(s); }
    else if (tag == "NODE") { RawNode n{}; ss >> n.id >> n.x >> n.y >> n.z; for (int k=0;k<6;++k) ss >> n.f[k]; for (int k=0;k<6;++k) ss >> n.p[k]; b.nodes.push_back(n); }
    else if (tag == "MEMBER") { RawMem mm{}; mm.active = 1; mm.tonly = 0; ss >> mm.id >> mm.i >> mm.j >> mm.mat >> mm.sec >> mm.rx >> mm.ry >> mm.rz;
                                int act; if (ss >> act) mm.active = act;
                                int to;  if (ss >> to)  mm.tonly  = to;
                                b.mems.push_back(mm); }
    else if (tag == "SHELL") { RawShell s{}; s.active = 1; ss >> s.id >> s.n[0] >> s.n[1] >> s.n[2] >> s.n[3] >> s.mat >> s.t;
                               int act; if (ss >> act) s.active = act;
                               b.shes.push_back(s); }
    else if (tag == "NLOAD") { RawNL l{}; ss >> l.node; for (int k=0;k<6;++k) ss >> l.c[k]; b.nls.push_back(l); }
    else if (tag == "UDL") { RawUDL u{}; ss >> u.member >> u.wx >> u.wy >> u.wz; b.udls.push_back(u); }
    else if (tag == "SPRESS") { RawSP s{}; ss >> s.shell >> s.p; b.sps.push_back(s); }
    else if (tag == "HINGE") { RawHinge h{}; ss >> h.member >> h.dof >> h.Mp; b.hins.push_back(h); }
    else if (tag == "OPT") { int er=0, ut=0; real pt=1e-12; ss >> er >> ut >> pt; b.opt.enableReleases=er!=0; b.opt.useTimoshenko=ut!=0; b.opt.pivotTol=pt;
                             int im; if (ss >> im) b.opt.useIncompatibleMembrane = im!=0;   // S8-8a (optional, back-compat)
                             int dk; if (ss >> dk) b.opt.useDKQPlate            = dk!=0; }  // S8-8b (optional, back-compat)
    else if (tag == "EIGEN") { ss >> b.nModes; }
    else if (tag == "PDELTA") { ss >> b.pdelta; }
    else if (tag == "TONLY")  { b.analysis = "TONLY";  int v; if (ss >> v) b.toMaxIter = v; if (ss >> v) b.toAllowReact = v; }
    else if (tag == "SIZEOPT"){ b.analysis = "SIZEOPT"; real a; int mi; real dt; if (ss >> a) b.soAmin = a; if (ss >> mi) b.soMaxIter = mi; if (ss >> dt) b.soDcTol = dt; }
    else if (tag == "DYNC")   { b.analysis = "DYNC";   real d; if (ss >> d) b.dcDt = d; if (ss >> d) b.dcMaxTime = d; int rid; while (ss >> rid) b.dcRemovals.push_back(rid); }
    // unknown tags are ignored (forward-compatible)
}

FrameModel buildModel(const Block& b) {
    FrameModel model;
    model.materials.reserve(b.mats.size());
    model.sections.reserve(b.secs.size());
    for (const auto& m : b.mats) {
        Material fm(m.E, m.G, m.rho); fm.nu = m.nu;
        fm.cap = m.hasCap ? Capacity::make(m.capC, m.capT, m.capS) : Capacity::make(300, 300, 180);
        model.materials.push_back(fm);
    }
    for (const auto& s : b.secs) {
        Section fs; fs.A=s.A; fs.Iy=s.Iy; fs.Iz=s.Iz; fs.J=s.J; fs.cy=s.cy; fs.cz=s.cz; fs.Asy=s.Asy; fs.Asz=s.Asz;
        model.sections.push_back(fs);
    }
    model.nodes.reserve(b.nodes.size());
    for (const auto& n : b.nodes) {
        Node fn(n.id, n.x, n.y, n.z);
        for (int k=0;k<6;++k) { fn.fixed[k] = (n.f[k]!=0); fn.prescribed[k] = n.p[k]; }
        model.nodes.push_back(fn);
    }
    model.members.reserve(b.mems.size());
    for (const auto& mm : b.mems) {
        Member fmem(mm.id, mm.i, mm.j, mm.mat, mm.sec);
        fmem.refVec = Vec3(mm.rx, mm.ry, mm.rz);
        fmem.active = (mm.active != 0);
        fmem.tensionOnly = (mm.tonly != 0);
        model.members.push_back(fmem);
    }
    model.shells.reserve(b.shes.size());
    for (const auto& s : b.shes) {
        ShellQuad sq(s.id, s.n[0], s.n[1], s.n[2], s.n[3], s.mat, s.t);
        sq.active = (s.active != 0);
        model.shells.push_back(sq);
    }
    for (const auto& l : b.nls) { NodalLoad nl; nl.node=l.node; for (int k=0;k<6;++k) nl.comp[k]=l.c[k]; model.nodalLoads.push_back(nl); }
    for (const auto& u : b.udls) { MemberUDL mu; mu.member=u.member; mu.w_local=Vec3(u.wx,u.wy,u.wz); model.memberUDLs.push_back(mu); }
    for (const auto& s : b.sps) { ShellPressure sp; sp.shell=s.shell; sp.p=s.p; model.shellPressures.push_back(sp); }
    for (const auto& h : b.hins) { model.hinges.push_back(PlasticHinge{ h.member, h.dof, h.Mp }); }
    return model;
}

void printState(std::string& out, const FrameModel& model, const SolveResult& r) {
    appendf(out, "SINGULAR %d\n", r.singular ? 1 : 0);
    for (size_t k = 0; k < model.nodes.size(); ++k) {
        appendf(out, "DISP %d %.12g %.12g %.12g %.12g %.12g %.12g\n", model.nodes[k].id,
                r.disp((int)k,Ux), r.disp((int)k,Uy), r.disp((int)k,Uz),
                r.disp((int)k,Rx), r.disp((int)k,Ry), r.disp((int)k,Rz));
        appendf(out, "RF %d %.12g %.12g %.12g %.12g %.12g %.12g\n", model.nodes[k].id,
                r.reaction((int)k,Ux), r.reaction((int)k,Uy), r.reaction((int)k,Uz),
                r.reaction((int)k,Rx), r.reaction((int)k,Ry), r.reaction((int)k,Rz));
    }
    for (size_t e = 0; e < r.memberForces.size(); ++e) {
        const MemberForcePair& mf = r.memberForces[e];
        appendf(out, "MF %d %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n", mf.member,
                mf.endI.N, mf.endI.Vy, mf.endI.Vz, mf.endI.T, mf.endI.My, mf.endI.Mz,
                mf.endJ.N, mf.endJ.Vy, mf.endJ.Vz, mf.endJ.T, mf.endJ.My, mf.endJ.Mz);
    }
    for (size_t e = 0; e < r.shellForces.size(); ++e) {
        const ShellElementForces& sf = r.shellForces[e];
        appendf(out, "SF %d %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n", sf.shell,
                sf.Mxx, sf.Myy, sf.Mxy, sf.Qx, sf.Qy, sf.Nxx, sf.Nyy, sf.Nxy);
    }
}

void runBlock(std::string& out, const Block& b) {
    appendf(out, "VERSION %s\n", FRAMECORE_BUILD_SHA);
    if (b.empty || b.nodes.empty()) return;   // bare END = handshake only

    FrameModel model = buildModel(b);

    if (b.analysis == "TONLY") {
        TensionOnlyOptions to; to.maxIter = b.toMaxIter; to.allowReactivation = (b.toAllowReact != 0); to.solve = b.opt;
        const TensionOnlyResult R = runTensionOnly(model, to);
        appendf(out, "TONLY %d %d %d\n", R.converged ? 1 : 0, R.cycled ? 1 : 0, R.iterations);
        out += "SLACK";
        for (MemberId id : R.slack) appendf(out, " %d", id);
        out += "\n";
        printState(out, model, R.finalState);
        return;
    }
    if (b.analysis == "SIZEOPT") {
        SizeOptOptions so; so.Amin = b.soAmin; so.maxIter = b.soMaxIter; so.dcTol = b.soDcTol; so.solve = b.opt;
        const SizeOptResult R = runSizeOptimization(model, so);
        appendf(out, "SIZEOPT %d %d %d\n", R.converged ? 1 : 0, R.iterations, R.singular ? 1 : 0);
        real vol = 0;
        for (size_t e = 0; e < model.members.size(); ++e) {
            const real A  = e < R.finalAreas.size() ? R.finalAreas[e] : real(0);
            const real dc = e < R.finalDC.size()    ? R.finalDC[e]    : real(0);
            appendf(out, "AREA %d %.12g %.12g\n", model.members[e].id, A, dc);
            const int ni = model.nodeIndex(model.members[e].i), nj = model.nodeIndex(model.members[e].j);
            if (ni >= 0 && nj >= 0)
                vol += A * norm(model.nodes[(size_t)nj].pos - model.nodes[(size_t)ni].pos);
        }
        appendf(out, "WEIGHTVOL %.12g\n", vol);
        return;
    }
    if (b.analysis == "DYNC") {
        DynCollapseOptions dco; dco.dt = b.dcDt; dco.maxTime = b.dcMaxTime; dco.solve = b.opt;
        for (int rid : b.dcRemovals) dco.initialRemovals.push_back((MemberId)rid);
        const DynCollapseHistory H = runDynamicCollapse(model, dco);
        const real tend = H.frames.empty() ? real(0) : H.frames.back().t;
        appendf(out, "DYNC %d %d %d %.12g\n", (int)H.outcome, (int)H.events.size(), (int)H.frames.size(), tend);
        for (const DynCollapseEvent& ev : H.events)
            appendf(out, "DEVENT %.12g %d %d %d\n", ev.t, (int)ev.mode,
                    (int)ev.removedMembers.size(), (int)ev.detached.size());
        for (const DynCollapseFrame& fr : H.frames) {
            real mx = 0; for (real x : fr.u) mx = std::max(mx, std::fabs(x));
            appendf(out, "DFRAME %.12g %.12g\n", fr.t, mx);
        }
        return;
    }

    SolveResult r;
    if (b.pdelta >= 0) {
        PDeltaOptions po; po.refactorPath = (b.pdelta != 0); po.maxIter = 5000; po.tolU = 1e-13; po.solve = b.opt;
        const PDeltaResult pr = runPDelta(model, po);
        appendf(out, "PDSTATUS %d %d %d\n", pr.converged ? 1 : 0, pr.diverged ? 1 : 0, pr.iterations);
        r = pr.finalState;
    } else {
        r = solve(model, b.opt);
    }
    printState(out, model, r);
    if (b.nModes > 0) {
        const PreparedSystem ps = assembleAndFactor(model, b.opt);
        ModalOptions mo; mo.numModes = b.nModes;
        const ModalResult mr = solveModal(ps, mo);
        appendf(out, "FREQ %d", (int)mr.modes.size());
        for (const auto& md : mr.modes) appendf(out, " %.12g", md.omega);
        out += "\n";
    }
}

}  // namespace

namespace frame_cli {

std::string processAll(const std::string& input) {
    std::string out;
    std::istringstream in(input);
    std::string line;
    Block b;
    bool inBlock = false;
    while (std::getline(in, line)) {
        std::istringstream ss(line);
        std::string tag; if (!(ss >> tag)) continue;
        if (tag == "END") {
            runBlock(out, b);
            out += "EOR\n";
            b = Block{};
            inBlock = false;
            continue;
        }
        inBlock = true;
        parseLine(b, tag, ss);
    }
    if (inBlock) { runBlock(out, b); out += "EOR\n"; }
    return out;
}

}  // namespace frame_cli
