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
//   NODE   id x y z  fUx fUy fUz fRx fRy fRz          (fixed flags 0/1)
//   MEMBER id i j matIdx secIdx  refx refy refz
//   SHELL  id n0 n1 n2 n3 matIdx t                    (MITC4 flat-shell facet)
//   NLOAD  node Fx Fy Fz Mx My Mz
//   UDL    member wx wy wz                            (local)
//   SPRESS shellId p                                  (transverse pressure)
//   OPT    enableReleases useTimoshenko pivotTol
//   END
// MAT and SMAT append to ONE material pool in input order; matIdx indexes that pool.
// OUTPUT:
//   SINGULAR <0|1>
//   DISP id ux uy uz rx ry rz                         (one per node, node order)
//   MF   id Ni Vyi Vzi Ti Myi Mzi Nj Vyj Vzj Tj Myj Mzj   (one per member)
//   SF   id Mxx Myy Mxy Qx Qy Nxx Nyy Nxy                  (one per shell)
#include "FrameCore/FrameSolver.h"

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdio>

using namespace frame;

namespace {
struct RawMat { real E, G, rho, nu; };
struct RawSec { real A, Iy, Iz, J, cy, cz, Asy, Asz; };
struct RawNode { int id; real x, y, z; int f[6]; };
struct RawMem { int id, i, j, mat, sec; real rx, ry, rz; };
struct RawShell { int id, n[4], mat; real t; };
struct RawNL { int node; real c[6]; };
struct RawUDL { int member; real wx, wy, wz; };
struct RawSP { int shell; real p; };
}

int main() {
    std::vector<RawMat> mats; std::vector<RawSec> secs;
    std::vector<RawNode> nodes; std::vector<RawMem> mems;
    std::vector<RawShell> shes;
    std::vector<RawNL> nls; std::vector<RawUDL> udls; std::vector<RawSP> sps;
    SolveOptions opt;

    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream ss(line);
        std::string tag; if (!(ss >> tag)) continue;
        if (tag == "MAT") { RawMat m{}; ss >> m.E >> m.G >> m.rho; m.nu = 0; mats.push_back(m); }
        else if (tag == "SMAT") { RawMat m{}; ss >> m.E >> m.nu >> m.G; m.rho = 0; mats.push_back(m); }
        else if (tag == "SEC") { RawSec s{}; ss >> s.A >> s.Iy >> s.Iz >> s.J >> s.cy >> s.cz >> s.Asy >> s.Asz; secs.push_back(s); }
        else if (tag == "NODE") { RawNode n{}; ss >> n.id >> n.x >> n.y >> n.z; for (int k=0;k<6;++k) ss >> n.f[k]; nodes.push_back(n); }
        else if (tag == "MEMBER") { RawMem mm{}; ss >> mm.id >> mm.i >> mm.j >> mm.mat >> mm.sec >> mm.rx >> mm.ry >> mm.rz; mems.push_back(mm); }
        else if (tag == "SHELL") { RawShell s{}; ss >> s.id >> s.n[0] >> s.n[1] >> s.n[2] >> s.n[3] >> s.mat >> s.t; shes.push_back(s); }
        else if (tag == "NLOAD") { RawNL l{}; ss >> l.node; for (int k=0;k<6;++k) ss >> l.c[k]; nls.push_back(l); }
        else if (tag == "UDL") { RawUDL u{}; ss >> u.member >> u.wx >> u.wy >> u.wz; udls.push_back(u); }
        else if (tag == "SPRESS") { RawSP s{}; ss >> s.shell >> s.p; sps.push_back(s); }
        else if (tag == "OPT") { int er=0, ut=0; real pt=1e-12; ss >> er >> ut >> pt; opt.enableReleases=er!=0; opt.useTimoshenko=ut!=0; opt.pivotTol=pt; }
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
        for (int k=0;k<6;++k) fn.fixed[k] = (n.f[k]!=0);
        model.nodes.push_back(fn);
    }
    model.members.reserve(mems.size());
    for (const auto& mm : mems) {
        const Material* pm = &model.materials[(size_t)mm.mat];
        const Section*  ps = &model.sections[(size_t)mm.sec];
        Member fmem(mm.id, mm.i, mm.j, pm, ps);
        fmem.refVec = Vec3(mm.rx, mm.ry, mm.rz);
        model.members.push_back(fmem);
    }
    model.shells.reserve(shes.size());
    for (const auto& s : shes) {
        const Material* pm = &model.materials[(size_t)s.mat];
        model.shells.push_back(ShellQuad(s.id, s.n[0], s.n[1], s.n[2], s.n[3], pm, s.t));
    }
    for (const auto& l : nls) { NodalLoad nl; nl.node=l.node; for (int k=0;k<6;++k) nl.comp[k]=l.c[k]; model.nodalLoads.push_back(nl); }
    for (const auto& u : udls) { MemberUDL mu; mu.member=u.member; mu.w_local=Vec3(u.wx,u.wy,u.wz); model.memberUDLs.push_back(mu); }
    for (const auto& s : sps) { ShellPressure sp; sp.shell=s.shell; sp.p=s.p; model.shellPressures.push_back(sp); }

    const SolveResult r = solve(model, opt);

    std::printf("SINGULAR %d\n", r.singular ? 1 : 0);
    for (size_t k = 0; k < model.nodes.size(); ++k) {
        std::printf("DISP %d %.12g %.12g %.12g %.12g %.12g %.12g\n", (int)k,
                    r.disp((int)k,Ux), r.disp((int)k,Uy), r.disp((int)k,Uz),
                    r.disp((int)k,Rx), r.disp((int)k,Ry), r.disp((int)k,Rz));
    }
    for (size_t e = 0; e < r.memberForces.size(); ++e) {
        const MemberForcePair& mf = r.memberForces[e];
        std::printf("MF %d %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n", (int)e,
                    mf.endI.N, mf.endI.Vy, mf.endI.Vz, mf.endI.T, mf.endI.My, mf.endI.Mz,
                    mf.endJ.N, mf.endJ.Vy, mf.endJ.Vz, mf.endJ.T, mf.endJ.My, mf.endJ.Mz);
    }
    for (size_t e = 0; e < r.shellForces.size(); ++e) {
        const ShellElementForces& sf = r.shellForces[e];
        std::printf("SF %d %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g\n", (int)e,
                    sf.Mxx, sf.Myy, sf.Mxy, sf.Qx, sf.Qy, sf.Nxx, sf.Nyy, sf.Nxy);
    }
    return 0;
}
