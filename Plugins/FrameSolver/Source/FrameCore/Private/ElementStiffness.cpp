#include "ElementStiffness.h"
#include "FrameCore/MemberGeometry.h"
#include <algorithm>
#include <cmath>

namespace frame {

// Single shared source of the member local frame (Eigen-free public API).
void memberLocalAxes(const Vec3& pi, const Vec3& pj, const Vec3& refVec,
                     Vec3& outX, Vec3& outY, Vec3& outZ) {
    Vec3e x = toE(pj) - toE(pi);
    x /= x.norm();                             // local x = unit(j - i)

    Vec3e ref = toE(refVec);
    Vec3e z = x.cross(ref);                     // z perpendicular to x and ref
    if (z.norm() < 1e-8) {                      // refVec parallel to axis -> fallback +Y
        ref = Vec3e(0, 1, 0);
        z = x.cross(ref);
        if (z.norm() < 1e-8) { ref = Vec3e(1, 0, 0); z = x.cross(ref); }  // member along +Y
    }
    z.normalize();
    Vec3e y = z.cross(x);                        // right-handed: cross(x,y) = z
    y.normalize();

    outX = Vec3(x.x(), x.y(), x.z());
    outY = Vec3(y.x(), y.y(), y.z());
    outZ = Vec3(z.x(), z.y(), z.z());
}

Mat3 localAxes(const Vec3& pi, const Vec3& pj, const Vec3& refVec) {
    Vec3 ax, ay, az;
    memberLocalAxes(pi, pj, refVec, ax, ay, az);
    Mat3 R;                                     // rows are local axes in global coords
    R.row(0) = toE(ax).transpose();
    R.row(1) = toE(ay).transpose();
    R.row(2) = toE(az).transpose();
    return R;
}

Mat12 localStiffness12(real E, real G, real A, real Iy, real Iz, real J, real L) {
    Mat12 k = Mat12::Zero();
    const real EA = E * A / L;
    const real GJ = G * J / L;
    const real L2 = L * L, L3 = L * L * L;
    const real az = 12 * E * Iz / L3, bz = 6 * E * Iz / L2, cz = 4 * E * Iz / L, dz = 2 * E * Iz / L;
    const real ay = 12 * E * Iy / L3, by = 6 * E * Iy / L2, cy = 4 * E * Iy / L, dy = 2 * E * Iy / L;

    // axial  (u1=0, u2=6)
    k(0,0) = EA;  k(0,6) = -EA;  k(6,0) = -EA;  k(6,6) = EA;
    // torsion (rx1=3, rx2=9)
    k(3,3) = GJ;  k(3,9) = -GJ;  k(9,3) = -GJ;  k(9,9) = GJ;
    // bending about local z, x-y plane (v=1,7 ; rz=5,11)
    k(1,1)=az;  k(1,5)=bz;   k(1,7)=-az;  k(1,11)=bz;
    k(5,1)=bz;  k(5,5)=cz;   k(5,7)=-bz;  k(5,11)=dz;
    k(7,1)=-az; k(7,5)=-bz;  k(7,7)=az;   k(7,11)=-bz;
    k(11,1)=bz; k(11,5)=dz;  k(11,7)=-bz; k(11,11)=cz;
    // bending about local y, x-z plane (w=2,8 ; ry=4,10) — sign flips vs the z-block
    k(2,2)=ay;   k(2,4)=-by;  k(2,8)=-ay;  k(2,10)=-by;
    k(4,2)=-by;  k(4,4)=cy;   k(4,8)=by;   k(4,10)=dy;
    k(8,2)=-ay;  k(8,4)=by;   k(8,8)=ay;   k(8,10)=by;
    k(10,2)=-by; k(10,4)=dy;  k(10,8)=by;  k(10,10)=cy;

    return k;
}

Mat12 localStiffness12T(real E, real G, real A, real Iy, real Iz, real J, real L,
                        real Asy, real Asz) {
    Mat12 k = Mat12::Zero();
    const real EA = E * A / L;
    const real GJ = G * J / L;
    const real L2 = L * L;

    // bending about local z (x-y plane); shear in local y resisted by Asy
    const real Phiz = (Asy > 0) ? 12.0 * E * Iz / (G * Asy * L2) : 0.0;
    const real bz0  = E * Iz / ((1.0 + Phiz) * L * L2);   // = E Iz / ((1+Phi) L^3)
    const real az = 12.0 * bz0;
    const real bz = 6.0  * bz0 * L;
    const real cz = (4.0 + Phiz) * bz0 * L2;
    const real dz = (2.0 - Phiz) * bz0 * L2;

    // bending about local y (x-z plane); shear in local z resisted by Asz
    const real Phiy = (Asz > 0) ? 12.0 * E * Iy / (G * Asz * L2) : 0.0;
    const real by0  = E * Iy / ((1.0 + Phiy) * L * L2);
    const real ay = 12.0 * by0;
    const real by = 6.0  * by0 * L;
    const real cy = (4.0 + Phiy) * by0 * L2;
    const real dy = (2.0 - Phiy) * by0 * L2;

    // axial  (u1=0, u2=6)
    k(0,0) = EA;  k(0,6) = -EA;  k(6,0) = -EA;  k(6,6) = EA;
    // torsion (rx1=3, rx2=9)
    k(3,3) = GJ;  k(3,9) = -GJ;  k(9,3) = -GJ;  k(9,9) = GJ;
    // bending about local z, x-y plane (v=1,7 ; rz=5,11)
    k(1,1)=az;  k(1,5)=bz;   k(1,7)=-az;  k(1,11)=bz;
    k(5,1)=bz;  k(5,5)=cz;   k(5,7)=-bz;  k(5,11)=dz;
    k(7,1)=-az; k(7,5)=-bz;  k(7,7)=az;   k(7,11)=-bz;
    k(11,1)=bz; k(11,5)=dz;  k(11,7)=-bz; k(11,11)=cz;
    // bending about local y, x-z plane (w=2,8 ; ry=4,10) — sign flips vs the z-block
    k(2,2)=ay;   k(2,4)=-by;  k(2,8)=-ay;  k(2,10)=-by;
    k(4,2)=-by;  k(4,4)=cy;   k(4,8)=by;   k(4,10)=dy;
    k(8,2)=-ay;  k(8,4)=by;   k(8,8)=ay;   k(8,10)=by;
    k(10,2)=-by; k(10,4)=dy;  k(10,8)=by;  k(10,10)=cy;

    return k;
}

Mat12 transform12(const Mat3& R) {
    Mat12 T = Mat12::Zero();
    for (int b = 0; b < 4; ++b) T.block<3, 3>(3 * b, 3 * b) = R;
    return T;
}

Mat12 localGeometric12(real P, real L) {
    Mat12 g = Mat12::Zero();
    const real L2 = L * L;
    const real c = P / (30.0 * L);
    // bending about local z, x-y plane (v=1,7 ; rz=5,11) — canonical signs
    g(1,1)=36*c;    g(1,5)=3*L*c;   g(1,7)=-36*c;   g(1,11)=3*L*c;
    g(5,1)=3*L*c;   g(5,5)=4*L2*c;  g(5,7)=-3*L*c;  g(5,11)=-L2*c;
    g(7,1)=-36*c;   g(7,5)=-3*L*c;  g(7,7)=36*c;    g(7,11)=-3*L*c;
    g(11,1)=3*L*c;  g(11,5)=-L2*c;  g(11,7)=-3*L*c; g(11,11)=4*L2*c;
    // bending about local y, x-z plane (w=2,8 ; ry=4,10) — L-terms sign-flipped (mirror stiffness)
    g(2,2)=36*c;    g(2,4)=-3*L*c;  g(2,8)=-36*c;   g(2,10)=-3*L*c;
    g(4,2)=-3*L*c;  g(4,4)=4*L2*c;  g(4,8)=3*L*c;   g(4,10)=-L2*c;
    g(8,2)=-36*c;   g(8,4)=3*L*c;   g(8,8)=36*c;    g(8,10)=3*L*c;
    g(10,2)=-3*L*c; g(10,4)=-L2*c;  g(10,8)=3*L*c;  g(10,10)=4*L2*c;
    return g;
}

Mat12 localMass12(real rho, real A, real Iy, real Iz, real L) {
    Mat12 m = Mat12::Zero();
    const real L2 = L * L;
    const real ma = rho * A * L / 6.0;          // axial / torsion lumping factor
    const real mt = rho * (Iy + Iz) * L / 6.0;  // torsion uses the polar inertia Ip = Iy + Iz
    const real mb = rho * A * L / 420.0;        // bending consistent-mass factor

    // axial (u1=0, u2=6)
    m(0,0) = 2*ma; m(0,6) = ma;  m(6,0) = ma;  m(6,6) = 2*ma;
    // torsion (rx1=3, rx2=9)
    m(3,3) = 2*mt; m(3,9) = mt;  m(9,3) = mt;  m(9,9) = 2*mt;
    // bending about local z, x-y plane (v=1,7 ; rz=5,11) — canonical signs
    m(1,1)=156*mb;    m(1,5)=22*L*mb;   m(1,7)=54*mb;     m(1,11)=-13*L*mb;
    m(5,1)=22*L*mb;   m(5,5)=4*L2*mb;   m(5,7)=13*L*mb;   m(5,11)=-3*L2*mb;
    m(7,1)=54*mb;     m(7,5)=13*L*mb;   m(7,7)=156*mb;    m(7,11)=-22*L*mb;
    m(11,1)=-13*L*mb; m(11,5)=-3*L2*mb; m(11,7)=-22*L*mb; m(11,11)=4*L2*mb;
    // bending about local y, x-z plane (w=2,8 ; ry=4,10) — L-terms sign-flipped (mirror stiffness)
    m(2,2)=156*mb;    m(2,4)=-22*L*mb;  m(2,8)=54*mb;     m(2,10)=13*L*mb;
    m(4,2)=-22*L*mb;  m(4,4)=4*L2*mb;   m(4,8)=-13*L*mb;  m(4,10)=-3*L2*mb;
    m(8,2)=54*mb;     m(8,4)=-13*L*mb;  m(8,8)=156*mb;    m(8,10)=22*L*mb;
    m(10,2)=13*L*mb;  m(10,4)=-3*L2*mb; m(10,8)=22*L*mb;  m(10,10)=4*L2*mb;

    return m;
}

bool condenseReleases(Mat12& kl, Vec12& Qf, const std::array<bool, 12>& release) {
    int cidx[12]; int nc = 0;
    for (int i = 0; i < 12; ++i) if (release[i]) cidx[nc++] = i;
    if (nc == 0) return true;     // nothing released
    if (nc >= 12) return false;   // everything released -> no element left

    // released sub-block kcc (nc x nc) + singularity guard via rank-revealing LU
    Eigen::MatrixXd kcc(nc, nc);
    for (int a = 0; a < nc; ++a)
        for (int b = 0; b < nc; ++b)
            kcc(a, b) = kl(cidx[a], cidx[b]);

    Eigen::FullPivLU<Eigen::MatrixXd> lu(kcc);
    lu.setThreshold(real(1e-9));                  // RELATIVE to max pivot (Eigen scales it)
    if (lu.rank() < nc) return false;             // singular released block -> free mechanism
    const Eigen::MatrixXd kccInv = lu.inverse();

    // kRC = kl(:, c) (12 x nc); kl symmetric so kCR = kRC^T. Qc = Qf(c).
    Eigen::MatrixXd kRC(12, nc);
    for (int i = 0; i < 12; ++i)
        for (int a = 0; a < nc; ++a)
            kRC(i, a) = kl(i, cidx[a]);
    Eigen::VectorXd Qc(nc);
    for (int a = 0; a < nc; ++a) Qc(a) = Qf(cidx[a]);

    // Compute the full-size correction, subtract, then zero released rows/cols/entries
    // (zeroing overwrites the garbage the correction left in released positions).
    const Eigen::MatrixXd dK = kRC * kccInv * kRC.transpose();   // 12x12
    const Eigen::VectorXd dQ = kRC * (kccInv * Qc);             // 12
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) kl(i, j) -= dK(i, j);
        Qf(i) -= dQ(i);
    }
    for (int a = 0; a < nc; ++a) {
        const int idx = cidx[a];
        for (int j = 0; j < 12; ++j) { kl(idx, j) = 0; kl(j, idx) = 0; }
        Qf(idx) = 0;
    }
    return true;
}

} // namespace frame
