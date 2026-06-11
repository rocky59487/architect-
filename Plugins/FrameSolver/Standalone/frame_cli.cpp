// frame_cli -- stdin/stdout driver around the FrameCore analyses (OpenSees harness #14 + the S6
// Grasshopper / external-client bridge). Reads a line-based model description and prints the solved
// state at full precision. Engine-only; never linked into the game. Units: N, mm, MPa.
//
// DAEMON MODE (S6 J1.5): the driver loops over model BLOCKS. Each block ends at `END`; the driver
// solves it, prints its output, then a lone `EOR` (end-of-response) line + flushes, and resets for
// the next block. A client that keeps the pipe open streams many blocks through ONE process (the
// PreparedSystem/ReSolveSession reuse is a future optimisation; results are already identical to an
// independent cli per block). A single-shot client (model + END + EOF) is the one-block case and is
// byte-for-byte unchanged; the trailing EOR is an unknown token the OpenSees harness ignores.
//
// INPUT (whitespace/line tokens; MAT/SMAT/SEC must precede the elements that index them):
//   MAT    E G rho  [capComp capTens capShear]         (beam material; nu=0. Optional allowable cap
//                                                        for SIZEOPT / D-C; default make(300,300,180).
//                                                        1 cap value -> tens=comp, shear=0.6*comp)
//   SMAT   E nu G    [capComp capTens capShear]         (shell material; carries nu)
//   SEC    A Iy Iz J cy cz Asy Asz
//   NODE   id x y z  fUx fUy fUz fRx fRy fRz  [pUx pUy pUz pRx pRy pRz]
//   MEMBER id i j matIdx secIdx  refx refy refz  [active [tonly]]  (active default 1; tonly default 0)
//   SHELL  id n0 n1 n2 n3 matIdx t  [active]
//   NLOAD  node Fx Fy Fz Mx My Mz
//   UDL    member wx wy wz                              (local)
//   SPRESS shellId p
//   HINGE  member dof Mp                                (dof 4/5/10/11, signed Mp)
//   OPT    enableReleases useTimoshenko pivotTol
//   EIGEN  nModes
//   PDELTA path                                         (0=frozen reuse, 1=K_T ref; absent/<0=linear)
//   TONLY  [maxIter [allowReact]]                       (tension-only eliminator on MEMBER ... tonly)
//   SIZEOPT Amin maxIter dcTol                          (fully-stressed sizing of every active member)
//   DYNC   dt maxTime [rid...]                          (dynamic collapse; trailing ids = removals)
//   END                                                 (solve this block; daemon resets after)
// MAT and SMAT append to ONE material pool in input order; matIdx indexes that pool.
//
// OUTPUT (per block):
//   VERSION <sha>                                       (always the first line of a block -- handshake)
//   SINGULAR <0|1>
//   DISP nodeId ux uy uz rx ry rz                       (one per node)
//   RF   nodeId Fx Fy Fz Mx My Mz                       (one per node)
//   MF   id Ni Vyi Vzi Ti Myi Mzi Nj Vyj Vzj Tj Myj Mzj (one per member)
//   SF   id Mxx Myy Mxy Qx Qy Nxx Nyy Nxy               (one per shell)
//   FREQ n omega1 ...                                   (when EIGEN given)
//   PDSTATUS conv div iters                             (when PDELTA given)
//   TONLY conv cycled iters / SLACK id... / <state>     (when TONLY given)
//   SIZEOPT conv iters singular / AREA id A DC / WEIGHTVOL v   (when SIZEOPT given)
//   DYNC outcome nEvents nFrames Tend / DEVENT t mode nRem nDet / DFRAME t maxAbsU   (when DYNC given)
//   EOR                                                 (end-of-response sentinel; flushes the block)
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ModalAnalysis.h"
#include "FrameCore/PDeltaAnalysis.h"
#include "FrameCore/TensionOnly.h"
#include "FrameCore/SizeOpt.h"
#include "FrameCore/DynamicCollapse.h"
#include "FrameCore/MemberGeometry.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cmath>

using namespace frame;

#ifndef FRAMECORE_BUILD_SHA
#define FRAMECORE_BUILD_SHA "unknown"   // overridden by the build script via /D (git short SHA)
#endif

namespace {
struct RawMat { real E, G, rho, nu; bool hasCap; real capC, capT, capS; };
struct RawSec { real A, Iy, Iz, J, cy, cz, Asy, Asz; };
struct RawNode { int id; real x, y, z; int f[6]; real p[6]; };
struct RawMem { int id, i, j, mat, sec; real rx, ry, rz; int active; int tonly; };
struct RawShell { int id, n[4], mat; real t; int active; };
struct RawNL { int node; real c[6]; };
struct RawUDL { int member; real wx, wy, wz; };
struct RawSP { int shell; real p; };
struct RawHinge { int member, dof; real Mp; };

// All the state of ONE model block (reset between daemon requests).
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
    bool empty = true;   // no model tokens seen -> a bare END is a handshake-only request
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
    else if (tag == "OPT") { int er=0, ut=0; real pt=1e-12; ss >> er >> ut >> pt; b.opt.enableReleases=er!=0; b.opt.useTimoshenko=ut!=0; b.opt.pivotTol=pt; }
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

// Standard solved-state output (shared by the linear / P-Delta / tension-only paths).
void printState(const FrameModel& model, const SolveResult& r) {
    std::printf("SINGULAR %d\n", r.singular ? 1 : 0);
    for (size_t k = 0; k < model.nodes.size(); ++k) {
        std::printf("DISP %d %.12g %.12g %.12g %.12g %.12g %.12g\n", model.nodes[k].id,
                    r.disp((int)k,Ux), r.disp((int)k,Uy), r.disp((int)k,Uz),
                    r.disp((int)k,Rx), r.disp((int)k,Ry), r.disp((int)k,Rz));
        std::printf("RF %d %.12g %.12g %.12g %.12g %.12g %.12g\n", model.nodes[k].id,
                    r.reaction((int)k,Ux), r.reaction((int)k,Uy), r.reaction((int)k,Uz),
                    r.reaction((int)k,Rx), r.reaction((int)k,Ry), r.reaction((int)k,Rz));
    }
    for (size_t e = 0; e < r.memberForces.size(); ++e) {
        const MemberForcePair& mf = r.memberForces[e];
        std::printf("MF %d %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n", mf.member,
                    mf.endI.N, mf.endI.Vy, mf.endI.Vz, mf.endI.T, mf.endI.My, mf.endI.Mz,
                    mf.endJ.N, mf.endJ.Vy, mf.endJ.Vz, mf.endJ.T, mf.endJ.My, mf.endJ.Mz);
    }
    for (size_t e = 0; e < r.shellForces.size(); ++e) {
        const ShellElementForces& sf = r.shellForces[e];
        std::printf("SF %d %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n", sf.shell,
                    sf.Mxx, sf.Myy, sf.Mxy, sf.Qx, sf.Qy, sf.Nxx, sf.Nyy, sf.Nxy);
    }
}

void runBlock(const Block& b) {
    // Version handshake first (a bare END = handshake-only request: VERSION then EOR, no solve).
    std::printf("VERSION %s\n", FRAMECORE_BUILD_SHA);
    if (b.empty || b.nodes.empty()) return;

    FrameModel model = buildModel(b);

    if (b.analysis == "TONLY") {
        TensionOnlyOptions to; to.maxIter = b.toMaxIter; to.allowReactivation = (b.toAllowReact != 0); to.solve = b.opt;
        const TensionOnlyResult R = runTensionOnly(model, to);
        std::printf("TONLY %d %d %d\n", R.converged ? 1 : 0, R.cycled ? 1 : 0, R.iterations);
        std::printf("SLACK");
        for (MemberId id : R.slack) std::printf(" %d", id);
        std::printf("\n");
        printState(model, R.finalState);
        return;
    }
    if (b.analysis == "SIZEOPT") {
        SizeOptOptions so; so.Amin = b.soAmin; so.maxIter = b.soMaxIter; so.dcTol = b.soDcTol; so.solve = b.opt;
        const SizeOptResult R = runSizeOptimization(model, so);
        std::printf("SIZEOPT %d %d %d\n", R.converged ? 1 : 0, R.iterations, R.singular ? 1 : 0);
        real vol = 0;
        for (size_t e = 0; e < model.members.size(); ++e) {
            const real A  = e < R.finalAreas.size() ? R.finalAreas[e] : real(0);
            const real dc = e < R.finalDC.size()    ? R.finalDC[e]    : real(0);
            std::printf("AREA %d %.12g %.12g\n", model.members[e].id, A, dc);
            const int ni = model.nodeIndex(model.members[e].i), nj = model.nodeIndex(model.members[e].j);
            if (ni >= 0 && nj >= 0)
                vol += A * norm(model.nodes[(size_t)nj].pos - model.nodes[(size_t)ni].pos);
        }
        std::printf("WEIGHTVOL %.12g\n", vol);
        return;
    }
    if (b.analysis == "DYNC") {
        DynCollapseOptions dco; dco.dt = b.dcDt; dco.maxTime = b.dcMaxTime; dco.solve = b.opt;
        for (int rid : b.dcRemovals) dco.initialRemovals.push_back((MemberId)rid);
        const DynCollapseHistory H = runDynamicCollapse(model, dco);
        const real tend = H.frames.empty() ? real(0) : H.frames.back().t;
        std::printf("DYNC %d %d %d %.12g\n", (int)H.outcome, (int)H.events.size(), (int)H.frames.size(), tend);
        for (const DynCollapseEvent& ev : H.events)
            std::printf("DEVENT %.12g %d %d %d\n", ev.t, (int)ev.mode,
                        (int)ev.removedMembers.size(), (int)ev.detached.size());
        // J1b: stream a compact per-frame timeline (peak |displacement|) for client-side playback.
        for (const DynCollapseFrame& fr : H.frames) {
            real mx = 0; for (real x : fr.u) mx = std::max(mx, std::fabs(x));
            std::printf("DFRAME %.12g %.12g\n", fr.t, mx);
        }
        return;
    }

    SolveResult r;
    if (b.pdelta >= 0) {
        PDeltaOptions po; po.refactorPath = (b.pdelta != 0); po.maxIter = 5000; po.tolU = 1e-13; po.solve = b.opt;
        const PDeltaResult pr = runPDelta(model, po);
        std::printf("PDSTATUS %d %d %d\n", pr.converged ? 1 : 0, pr.diverged ? 1 : 0, pr.iterations);
        r = pr.finalState;
    } else {
        r = solve(model, b.opt);
    }
    printState(model, r);
    if (b.nModes > 0) {
        const PreparedSystem ps = assembleAndFactor(model, b.opt);
        ModalOptions mo; mo.numModes = b.nModes;
        const ModalResult mr = solveModal(ps, mo);
        std::printf("FREQ %d", (int)mr.modes.size());
        for (const auto& md : mr.modes) std::printf(" %.12g", md.omega);
        std::printf("\n");
    }
}
}  // namespace

int main() {
    // Provenance to stderr; stdout is parsed by clients, so it stays clean.
    std::fprintf(stderr, "# frame_cli | build %s | compiled %s %s\n",
                 FRAMECORE_BUILD_SHA, __DATE__, __TIME__);

    Block b;
    bool inBlock = false;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream ss(line);
        std::string tag; if (!(ss >> tag)) continue;
        if (tag == "END") {
            runBlock(b);
            std::printf("EOR\n");
            std::fflush(stdout);     // daemon: a streaming client reads up to EOR per request
            b = Block{};             // reset for the next block
            inBlock = false;
            continue;
        }
        inBlock = true;
        parseLine(b, tag, ss);
    }
    // Lenient: a final block that reached EOF without a trailing END is still solved once.
    if (inBlock) { runBlock(b); std::printf("EOR\n"); std::fflush(stdout); }
    return 0;
}
