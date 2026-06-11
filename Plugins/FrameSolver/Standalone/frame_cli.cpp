// frame_cli — a tiny stdin/stdout driver around frame::solve for the OpenSees offline
// comparison harness (#14). It reads a line-based model description and prints the solved
// node displacements + member end forces at full precision, so Tools/opensees_compare.py
// can diff the SAME model solved by OpenSeesPy and by an analytic oracle. Engine-only;
// never linked into the game. Units: N, mm, MPa.
//
// INPUT (whitespace/line tokens; MAT/SMAT/SEC must precede the elements that index them):
//   MAT    E G rho                                    (beam material; nu unused -> 0)
//   SMAT   E nu G                                     (shell material; carries nu)
//   SEC    A Iy Iz J cy cz Asy Asz
//   NODE   id x y z  fUx fUy fUz fRx fRy fRz  [pUx pUy pUz pRx pRy pRz]
//          (6 fixed flags 0/1; optional 6 prescribed displacement values, default 0)
//   MEMBER id i j matIdx secIdx  refx refy refz  [active]   (optional 0/1, default 1)
//   SHELL  id n0 n1 n2 n3 matIdx t  [active]          (MITC4 flat-shell facet; optional 0/1, default 1)
//   NLOAD  node Fx Fy Fz Mx My Mz
//   UDL    member wx wy wz                            (local)
//   SPRESS shellId p                                  (transverse pressure)
//   HINGE  member dof Mp                              (plastic hinge: dof 4/5/10/11, signed Mp;
//                                                      the node-side moment is the caller's NLOAD)
//   OPT    enableReleases useTimoshenko pivotTol
//   PDELTA path                                       (second-order analysis: 0=frozen reuse,
//                                                      1=K_T reference; absent/<0 = linear solve.
//                                                      DISP/MF then report the second-order state;
//                                                      an extra "PDSTATUS conv div iters" line leads)
//   END
// MAT and SMAT append to ONE material pool in input order; matIdx indexes that pool.
// OUTPUT:
//   SINGULAR <0|1>
//   DISP nodeId ux uy uz rx ry rz                     (one per node, node order)
//   RF   nodeId Fx Fy Fz Mx My Mz                     (one per node, node order)
//   MF   id Ni Vyi Vzi Ti Myi Mzi Nj Vyj Vzj Tj Myj Mzj   (one per member)
//   SF   id Mxx Myy Mxy Qx Qy Nxx Nyy Nxy                  (one per shell)
#include "FrameCore/FrameSolver.h"
#include "FrameCore/ModalAnalysis.h"
#include "FrameCore/PDeltaAnalysis.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdio>

using namespace frame;

#ifndef FRAMECORE_BUILD_SHA
#define FRAMECORE_BUILD_SHA "unknown"   // overridden by the build script via /D (git short SHA)
#endif

namespace {
struct RawMat { real E, G, rho, nu; };
struct RawSec { real A, Iy, Iz, J, cy, cz, Asy, Asz; };
struct RawNode { int id; real x, y, z; int f[6]; real p[6]; };
struct RawMem { int id, i, j, mat, sec; real rx, ry, rz; int active; };
struct RawShell { int id, n[4], mat; real t; int active; };
struct RawNL { int node; real c[6]; };
struct RawUDL { int member; real wx, wy, wz; };
struct RawSP { int shell; real p; };
struct RawHinge { int member, dof; real Mp; };
}

int main() {
    // Provenance to stderr; stdout is parsed by the Python audit harness, so it stays clean.
    std::fprintf(stderr, "# frame_cli | build %s | compiled %s %s\n",
                 FRAMECORE_BUILD_SHA, __DATE__, __TIME__);
    std::vector<RawMat> mats; std::vector<RawSec> secs;
    std::vector<RawNode> nodes; std::vector<RawMem> mems;
    std::vector<RawShell> shes;
    std::vector<RawNL> nls; std::vector<RawUDL> udls; std::vector<RawSP> sps;
    std::vector<RawHinge> hins;
    SolveOptions opt;
    int nModes = 0;
    int pdelta = -1;   // -1 = linear solve; 0 = P-Delta frozen reuse path; 1 = P-Delta K_T reference

    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream ss(line);
        std::string tag; if (!(ss >> tag)) continue;
        if (tag == "MAT") { RawMat m{}; ss >> m.E >> m.G >> m.rho; m.nu = 0; mats.push_back(m); }
        else if (tag == "SMAT") { RawMat m{}; ss >> m.E >> m.nu >> m.G; m.rho = 0; mats.push_back(m); }
        else if (tag == "SEC") { RawSec s{}; ss >> s.A >> s.Iy >> s.Iz >> s.J >> s.cy >> s.cz >> s.Asy >> s.Asz; secs.push_back(s); }
        else if (tag == "NODE") { RawNode n{}; ss >> n.id >> n.x >> n.y >> n.z; for (int k=0;k<6;++k) ss >> n.f[k]; for (int k=0;k<6;++k) ss >> n.p[k]; nodes.push_back(n); }
        else if (tag == "MEMBER") { RawMem mm{}; mm.active = 1; ss >> mm.id >> mm.i >> mm.j >> mm.mat >> mm.sec >> mm.rx >> mm.ry >> mm.rz;
                                    int act; if (ss >> act) mm.active = act;   // optional token (a FAILED >> writes 0, so guard it)
                                    mems.push_back(mm); }
        else if (tag == "SHELL") { RawShell s{}; s.active = 1; ss >> s.id >> s.n[0] >> s.n[1] >> s.n[2] >> s.n[3] >> s.mat >> s.t;
                                   int act; if (ss >> act) s.active = act;     // optional token (a FAILED >> writes 0, so guard it)
                                   shes.push_back(s); }
        else if (tag == "NLOAD") { RawNL l{}; ss >> l.node; for (int k=0;k<6;++k) ss >> l.c[k]; nls.push_back(l); }
        else if (tag == "UDL") { RawUDL u{}; ss >> u.member >> u.wx >> u.wy >> u.wz; udls.push_back(u); }
        else if (tag == "SPRESS") { RawSP s{}; ss >> s.shell >> s.p; sps.push_back(s); }
        else if (tag == "HINGE") { RawHinge h{}; ss >> h.member >> h.dof >> h.Mp; hins.push_back(h); }
        else if (tag == "OPT") { int er=0, ut=0; real pt=1e-12; ss >> er >> ut >> pt; opt.enableReleases=er!=0; opt.useTimoshenko=ut!=0; opt.pivotTol=pt; }
        else if (tag == "EIGEN") { ss >> nModes; }
        else if (tag == "PDELTA") { ss >> pdelta; }
        else if (tag == "END") break;
    }

    FrameModel model;
    model.materials.reserve(mats.size());
    model.sections.reserve(secs.size());
    for (const auto& m : mats) { Material fm(m.E, m.G, m.rho); fm.nu = m.nu; fm.cap = Capacity::make(300,300,180); model.materials.push_back(fm); }
    for (const auto& s : secs) {
        Section fs; fs.A=s.A; fs.Iy=s.Iy; fs.Iz=s.Iz; fs.J=s.J; fs.cy=s.cy; fs.cz=s.cz; fs.Asy=s.Asy; fs.Asz=s.Asz;
        model.sections.push_back(fs);
    }
    model.nodes.reserve(nodes.size());
    for (const auto& n : nodes) {
        Node fn(n.id, n.x, n.y, n.z);
        for (int k=0;k<6;++k) { fn.fixed[k] = (n.f[k]!=0); fn.prescribed[k] = n.p[k]; }
        model.nodes.push_back(fn);
    }
    model.members.reserve(mems.size());
    for (const auto& mm : mems) {
        Member fmem(mm.id, mm.i, mm.j, mm.mat, mm.sec);   // mat/sec are already pool indices
        fmem.refVec = Vec3(mm.rx, mm.ry, mm.rz);
        fmem.active = (mm.active != 0);
        model.members.push_back(fmem);
    }
    model.shells.reserve(shes.size());
    for (const auto& s : shes) {
        ShellQuad sq(s.id, s.n[0], s.n[1], s.n[2], s.n[3], s.mat, s.t);
        sq.active = (s.active != 0);
        model.shells.push_back(sq);
    }
    for (const auto& l : nls) { NodalLoad nl; nl.node=l.node; for (int k=0;k<6;++k) nl.comp[k]=l.c[k]; model.nodalLoads.push_back(nl); }
    for (const auto& u : udls) { MemberUDL mu; mu.member=u.member; mu.w_local=Vec3(u.wx,u.wy,u.wz); model.memberUDLs.push_back(mu); }
    for (const auto& s : sps) { ShellPressure sp; sp.shell=s.shell; sp.p=s.p; model.shellPressures.push_back(sp); }
    for (const auto& h : hins) { model.hinges.push_back(PlasticHinge{ h.member, h.dof, h.Mp }); }

    SolveResult r;
    if (pdelta >= 0) {
        PDeltaOptions po; po.refactorPath = (pdelta != 0); po.maxIter = 5000; po.tolU = 1e-13; po.solve = opt;
        const PDeltaResult pr = runPDelta(model, po);
        std::printf("PDSTATUS %d %d %d\n", pr.converged ? 1 : 0, pr.diverged ? 1 : 0, pr.iterations);
        r = pr.finalState;
    } else {
        r = solve(model, opt);
    }

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
    if (nModes > 0) {
        const PreparedSystem ps = assembleAndFactor(model, opt);
        ModalOptions mo; mo.numModes = nModes;
        const ModalResult mr = solveModal(ps, mo);
        std::printf("FREQ %d", (int)mr.modes.size());
        for (const auto& md : mr.modes) std::printf(" %.12g", md.omega);   // rad/s
        std::printf("\n");
    }
    return 0;
}
