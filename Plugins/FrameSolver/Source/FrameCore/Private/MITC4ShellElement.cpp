#include "MITC4ShellElement.h"
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"
#include "FrameCore/SolveOptions.h"   // prepare() reads opts.useIncompatibleMembrane / useDKQPlate

#include <cmath>

namespace frame {

// ============================================================================
// MITC4 Reissner-Mindlin flat-shell facet (24 DOF = 4 nodes x 6). The local stiffness
// is three blocks mapped into the nodal DOFs [Ux,Uy,Uz,Rx,Ry,Rz]:
//   * plate BENDING (Uz,Rx,Ry) with the MITC4 assumed transverse shear (below),
//   * plane-stress MEMBRANE (Ux,Uy)        -- see membraneK / membraneToShellMap,
//   * Hughes-Brezzi DRILLING (Rz)          -- see Bdrill (gamma = G*t).
// The facet's local frame is built from the corner geometry and the whole 24x24 block
// is rotated into 3-D by T = blockdiag(R) (the solver then treats it like any element).
//
// The BENDING block is formulated in FIBER ROTATIONS (w, bx, by) with
//   u = z*bx,  v = z*by,  w = w(x,y)
//   curvatures  kappa = [ bx,x ; by,y ; bx,y + by,x ]
//   shear       gamma = [ w,x + bx ; w,y + by ]
// The transverse shear uses the MITC4 assumed natural strain (Bathe-Dvorkin 1985):
// the covariant shears g_xi and g_eta are sampled at the 4 edge midpoints (tying
// points A,B,C,D) and interpolated, which removes shear locking in the thin limit.
//
// The fiber rotations map to the engine's nodal rotations (Rx,Ry) by
//   bx = Ry,  by = -Rx,  w = Uz
// (a rigid rotation Ry about local y gives u = z*Ry  -> bx = Ry; Rx about x gives
//  v = -z*Rx -> by = -Rx). This sign map is applied through Pmap below and is
// pinned down by the constant-curvature patch test.
// ============================================================================

namespace {

// Bilinear Q4 shape functions and natural derivatives. Node order CCW:
// 1=(-1,-1) 2=(+1,-1) 3=(+1,+1) 4=(-1,+1).
inline void shapeN(real xi, real eta, real N[4]) {
    N[0] = 0.25 * (1 - xi) * (1 - eta);
    N[1] = 0.25 * (1 + xi) * (1 - eta);
    N[2] = 0.25 * (1 + xi) * (1 + eta);
    N[3] = 0.25 * (1 - xi) * (1 + eta);
}
inline void shapeDN(real xi, real eta, real dNxi[4], real dNeta[4]) {
    dNxi[0] = -0.25 * (1 - eta); dNxi[1] =  0.25 * (1 - eta);
    dNxi[2] =  0.25 * (1 + eta); dNxi[3] = -0.25 * (1 + eta);
    dNeta[0] = -0.25 * (1 - xi); dNeta[1] = -0.25 * (1 + xi);
    dNeta[2] =  0.25 * (1 + xi); dNeta[3] =  0.25 * (1 - xi);
}

// Jacobian J = [[x_xi, y_xi],[x_eta, y_eta]] at (xi,eta); returns detJ and Jinv.
inline real jacobian(const real xl[4], const real yl[4], real xi, real eta,
                     real Jinv[2][2]) {
    real dNxi[4], dNeta[4];
    shapeDN(xi, eta, dNxi, dNeta);
    real x_xi = 0, y_xi = 0, x_eta = 0, y_eta = 0;
    for (int i = 0; i < 4; ++i) {
        x_xi  += dNxi[i]  * xl[i];  y_xi  += dNxi[i]  * yl[i];
        x_eta += dNeta[i] * xl[i];  y_eta += dNeta[i] * yl[i];
    }
    const real det = x_xi * y_eta - y_xi * x_eta;
    const real inv = (det != 0.0) ? 1.0 / det : 0.0;
    Jinv[0][0] =  y_eta * inv; Jinv[0][1] = -y_xi * inv;
    Jinv[1][0] = -x_eta * inv; Jinv[1][1] =  x_xi * inv;
    return det;
}

using Mat3x12 = Eigen::Matrix<real, 3, 12>;
using Mat2x12 = Eigen::Matrix<real, 2, 12>;
using Row12   = Eigen::Matrix<real, 1, 12>;
using Mat3x4  = Eigen::Matrix<real, 3, 4>;    // QM6 incompatible-mode strain (3 strains x 4 bubble DOF)
using Mat12x4 = Eigen::Matrix<real, 12, 4>;   // QM6 compatible-incompatible coupling
using Mat4    = Eigen::Matrix<real, 4, 4>;    // QM6 internal bubble block
using Row4    = Eigen::Matrix<real, 1, 4>;    // QM6 bubble drilling contribution

// Bending strain-displacement (3x12) in plate DOF order [w,bx,by] x4.
// kappa = [bx,x ; by,y ; bx,y + by,x].
inline Mat3x12 Bbending(const real xl[4], const real yl[4], real xi, real eta) {
    real Jinv[2][2];
    jacobian(xl, yl, xi, eta, Jinv);
    real dNxi[4], dNeta[4];
    shapeDN(xi, eta, dNxi, dNeta);
    Mat3x12 B = Mat3x12::Zero();
    for (int i = 0; i < 4; ++i) {
        const real dNdx = Jinv[0][0] * dNxi[i] + Jinv[0][1] * dNeta[i];
        const real dNdy = Jinv[1][0] * dNxi[i] + Jinv[1][1] * dNeta[i];
        const int bx = 3 * i + 1, by = 3 * i + 2;
        B(0, bx) = dNdx;                 // kappa_x  = bx,x
        B(1, by) = dNdy;                 // kappa_y  = by,y
        B(2, bx) = dNdy;                 // kappa_xy = bx,y + by,x
        B(2, by) = dNdx;
    }
    return B;
}

// Covariant transverse-shear rows at a sample point, plate DOF order [w,bx,by] x4.
//   g_xi  = w,xi  + x_xi*bx + y_xi*by      (J row 1 . cartesian shear)
//   g_eta = w,eta + x_eta*bx + y_eta*by    (J row 2 . cartesian shear)
inline void covariantRows(const real xl[4], const real yl[4], real xi, real eta,
                          Row12& gxi, Row12& geta) {
    real N[4], dNxi[4], dNeta[4];
    shapeN(xi, eta, N);
    shapeDN(xi, eta, dNxi, dNeta);
    real x_xi = 0, y_xi = 0, x_eta = 0, y_eta = 0;
    for (int i = 0; i < 4; ++i) {
        x_xi  += dNxi[i]  * xl[i];  y_xi  += dNxi[i]  * yl[i];
        x_eta += dNeta[i] * xl[i];  y_eta += dNeta[i] * yl[i];
    }
    gxi.setZero(); geta.setZero();
    for (int i = 0; i < 4; ++i) {
        const int w = 3 * i, bx = 3 * i + 1, by = 3 * i + 2;
        gxi(w)  = dNxi[i];   gxi(bx)  = N[i] * x_xi;   gxi(by)  = N[i] * y_xi;
        geta(w) = dNeta[i];  geta(bx) = N[i] * x_eta;  geta(by) = N[i] * y_eta;
    }
}

// MITC4 assumed transverse-shear B (2x12, cartesian [gxz;gyz]) at (xi,eta).
// Tying: g_xi sampled at A(0,-1) & C(0,+1); g_eta at D(-1,0) & B(+1,0).
inline Mat2x12 Bshear(const real xl[4], const real yl[4], real xi, real eta) {
    Row12 rA, rC, rB, rD, tmp;
    covariantRows(xl, yl, 0.0, -1.0, rA, tmp);   // A: g_xi
    covariantRows(xl, yl, 0.0,  1.0, rC, tmp);   // C: g_xi
    covariantRows(xl, yl,  1.0, 0.0, tmp, rB);   // B: g_eta
    covariantRows(xl, yl, -1.0, 0.0, tmp, rD);   // D: g_eta

    const Row12 gxi  = 0.5 * (1 - eta) * rA + 0.5 * (1 + eta) * rC;
    const Row12 geta = 0.5 * (1 - xi)  * rD + 0.5 * (1 + xi)  * rB;

    real Jinv[2][2];
    jacobian(xl, yl, xi, eta, Jinv);             // [gxz;gyz] = Jinv * [gxi;geta]
    Mat2x12 B;
    B.row(0) = Jinv[0][0] * gxi + Jinv[0][1] * geta;
    B.row(1) = Jinv[1][0] * gxi + Jinv[1][1] * geta;
    return B;
}

// 2x2 Gauss rule.
const real kGP = 0.5773502691896258;   // 1/sqrt(3)
const real kG[2] = { -kGP, kGP };
const real kW[2] = { 1.0, 1.0 };
const real kShearCorr = 5.0 / 6.0;     // Reissner-Mindlin shear correction factor

// Plate-bending + MITC4-shear 12x12 stiffness in plate DOF order [w,bx,by] x4.
Mat12 plateK(const real xl[4], const real yl[4], real E, real nu, real G, real t) {
    Eigen::Matrix<real, 3, 3> Db = Eigen::Matrix<real, 3, 3>::Zero();
    const real Dfac = E * t * t * t / (12.0 * (1.0 - nu * nu));
    Db(0, 0) = Dfac;       Db(0, 1) = nu * Dfac;
    Db(1, 0) = nu * Dfac;  Db(1, 1) = Dfac;
    Db(2, 2) = Dfac * (1.0 - nu) * 0.5;
    const real Ds = kShearCorr * G * t;          // isotropic transverse shear modulus*t

    Mat12 K = Mat12::Zero();
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            const real xi = kG[a], eta = kG[b];
            real Jinv[2][2];
            const real detJ = jacobian(xl, yl, xi, eta, Jinv);
            const real wgt = kW[a] * kW[b] * detJ;
            const Mat3x12 Bb = Bbending(xl, yl, xi, eta);
            const Mat2x12 Bs = Bshear(xl, yl, xi, eta);
            K += wgt * (Bb.transpose() * Db * Bb);
            K += wgt * (Bs.transpose() * (Ds * Bs));
        }
    return K;
}

// ---- Membrane (plane stress) + drilling ----------------------------------------
// Membrane plane-stress B (3x12) in DOF order [u,v,thz] x4 (the thz columns are 0:
// the in-plane normal strain comes only from u,v).  eps = [u,x ; v,y ; u,y + v,x].
inline Mat3x12 Bmembrane(const real xl[4], const real yl[4], real xi, real eta) {
    real Jinv[2][2];
    jacobian(xl, yl, xi, eta, Jinv);
    real dNxi[4], dNeta[4];
    shapeDN(xi, eta, dNxi, dNeta);
    Mat3x12 B = Mat3x12::Zero();
    for (int i = 0; i < 4; ++i) {
        const real dNdx = Jinv[0][0] * dNxi[i] + Jinv[0][1] * dNeta[i];
        const real dNdy = Jinv[1][0] * dNxi[i] + Jinv[1][1] * dNeta[i];
        const int u = 3 * i, v = 3 * i + 1;
        B(0, u) = dNdx;                  // eps_x  = u,x
        B(1, v) = dNdy;                  // eps_y  = v,y
        B(2, u) = dNdy; B(2, v) = dNdx;  // gamma_xy = u,y + v,x
    }
    return B;
}

// Drilling B (1x12): the Hughes-Brezzi residual (thz - omega), omega = 0.5(v,x - u,y).
// Penalizing its energy gives the drilling DOF real stiffness tied to the membrane
// rotation, so it vanishes in constant-strain states (clean patch test) yet removes
// the in-plane rotational zero-energy mode that would otherwise make a coplanar shell
// patch singular.
inline Row12 Bdrill(const real xl[4], const real yl[4], real xi, real eta) {
    real Jinv[2][2];
    jacobian(xl, yl, xi, eta, Jinv);
    real N[4], dNxi[4], dNeta[4];
    shapeN(xi, eta, N);
    shapeDN(xi, eta, dNxi, dNeta);
    Row12 B; B.setZero();
    for (int i = 0; i < 4; ++i) {
        const real dNdx = Jinv[0][0] * dNxi[i] + Jinv[0][1] * dNeta[i];
        const real dNdy = Jinv[1][0] * dNxi[i] + Jinv[1][1] * dNeta[i];
        const int u = 3 * i, v = 3 * i + 1, tz = 3 * i + 2;
        B(tz) =  N[i];                   // + thz
        B(u)  =  0.5 * dNdy;             // + 0.5 u,y
        B(v)  = -0.5 * dNdx;             // - 0.5 v,x
    }
    return B;
}

// Membrane + drilling 12x12 stiffness in DOF order [u,v,thz] x4.
Mat12 membraneK(const real xl[4], const real yl[4], real E, real nu, real G, real t) {
    Eigen::Matrix<real, 3, 3> Dm = Eigen::Matrix<real, 3, 3>::Zero();
    const real f = E / (1.0 - nu * nu);
    Dm(0, 0) = f;       Dm(0, 1) = nu * f;
    Dm(1, 0) = nu * f;  Dm(1, 1) = f;
    Dm(2, 2) = f * (1.0 - nu) * 0.5;
    const real gamma = G * t;            // Hughes-Brezzi drilling penalty (per area)

    Mat12 K = Mat12::Zero();
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            const real xi = kG[a], eta = kG[b];
            real Jinv[2][2];
            const real detJ = jacobian(xl, yl, xi, eta, Jinv);
            const real wgt = kW[a] * kW[b] * detJ;
            const Mat3x12 Bm = Bmembrane(xl, yl, xi, eta);
            const Row12   Bd = Bdrill(xl, yl, xi, eta);
            K += wgt * t * (Bm.transpose() * Dm * Bm);
            K += wgt * gamma * (Bd.transpose() * Bd);
        }
    return K;
}

// ---- QM6 incompatible-mode membrane (S8-8a, opt-in) -----------------------------
// Wilson Q6 (1973) + Taylor (1976) correction: add 2 bubble modes phi1=1-xi^2, phi2=1-eta^2
// on each of u,v (4 element-internal DOF a=[a1,a2,a3,a4]), then static-condense them out.
//   u_inc = (1-xi^2) a1 + (1-eta^2) a2,   v_inc = (1-xi^2) a3 + (1-eta^2) a4
// The bubble strain B_inc uses the CENTRE Jacobian J0 (Taylor's correction): integral(B_inc)
// over the element = 0, so a constant-stress state stores no bubble energy and the distorted-
// mesh constant-strain patch test passes. This defeats the in-plane (membrane) bending lock
// that cripples the bilinear Q4 membrane on coarse meshes. Drilling is left untouched.
inline Mat3x4 BincQM6(const real J0inv[2][2], real xi, real eta) {
    const real dP1dxi  = -2.0 * xi;    // d(1-xi^2)/dxi   (P1 has no eta dependence)
    const real dP2deta = -2.0 * eta;   // d(1-eta^2)/deta (P2 has no xi dependence)
    // dP/dx = Jinv00 dP/dxi + Jinv01 dP/deta ; dP/dy = Jinv10 dP/dxi + Jinv11 dP/deta, with J0.
    const real dP1dx = J0inv[0][0] * dP1dxi;
    const real dP1dy = J0inv[1][0] * dP1dxi;
    const real dP2dx = J0inv[0][1] * dP2deta;
    const real dP2dy = J0inv[1][1] * dP2deta;
    Mat3x4 B; B.setZero();
    B(0, 0) = dP1dx; B(0, 1) = dP2dx;                                      // eps_x  = du_inc/dx
    B(1, 2) = dP1dy; B(1, 3) = dP2dy;                                      // eps_y  = dv_inc/dy
    B(2, 0) = dP1dy; B(2, 1) = dP2dy; B(2, 2) = dP1dx; B(2, 3) = dP2dx;    // gamma_xy
    return B;
}

// Drilling contribution of the QM6 bubbles (1x4): the bubble part of the Hughes-Brezzi
// drilling strain thz + 0.5 u,y - 0.5 v,x, i.e. 0.5 du_inc/dy - 0.5 dv_inc/dx (centre J0).
// CRITICAL: without this, the drilling penalty keeps using the COMPATIBLE rotation (no bubble),
// so a pure-bending state -- where the membrane needs v ~ x^2 (a bubble mode) to kill parasitic
// shear -- stores spurious drilling energy that RE-LOCKS what the membrane bubbles released
// (observed: QM6 unlocked only to -38% instead of ~-1% of the Euler-Bernoulli tip). Letting the
// drilling rotation track the bubble-enhanced membrane rotation restores the full unlock.
inline Row4 BdrillIncQM6(const real J0inv[2][2], real xi, real eta) {
    const real dP1dxi = -2.0 * xi, dP2deta = -2.0 * eta;
    const real dP1dx = J0inv[0][0] * dP1dxi, dP1dy = J0inv[1][0] * dP1dxi;
    const real dP2dx = J0inv[0][1] * dP2deta, dP2dy = J0inv[1][1] * dP2deta;
    Row4 B;
    B(0, 0) =  0.5 * dP1dy;  B(0, 1) =  0.5 * dP2dy;    // +0.5 du_inc/dy  (a1,a2 = P1,P2 on u)
    B(0, 2) = -0.5 * dP1dx;  B(0, 3) = -0.5 * dP2dx;    // -0.5 dv_inc/dx  (a3,a4 = P1,P2 on v)
    return B;
}

// QM6 membrane + drilling 12x12 in DOF order [u,v,thz] x4. The compatible part is exactly
// membraneK (reused -> no formula drift); the incompatible bubbles -- coupled to BOTH the
// membrane strain AND the drilling rotation -- are condensed via
//   K* = Kc - Kca Kaa^-1 Kca^T   (Kaa is SPD, so the inverse is safe).
Mat12 membraneK_QM6(const real xl[4], const real yl[4], real E, real nu, real G, real t) {
    const Mat12 Kc = membraneK(xl, yl, E, nu, G, t);    // compatible Q4 membrane + drilling

    Eigen::Matrix<real, 3, 3> Dm = Eigen::Matrix<real, 3, 3>::Zero();
    const real f = E / (1.0 - nu * nu);
    Dm(0, 0) = f;       Dm(0, 1) = nu * f;
    Dm(1, 0) = nu * f;  Dm(1, 1) = f;
    Dm(2, 2) = f * (1.0 - nu) * 0.5;
    const real gamma = G * t;                           // Hughes-Brezzi drilling penalty

    real J0inv[2][2];
    jacobian(xl, yl, 0.0, 0.0, J0inv);                  // centre Jacobian (constant) for bubbles

    Mat12x4 Kca = Mat12x4::Zero();
    Mat4    Kaa = Mat4::Zero();
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            const real xi = kG[a], eta = kG[b];
            real Jinv[2][2];
            const real detJ = jacobian(xl, yl, xi, eta, Jinv);
            const real w = kW[a] * kW[b] * detJ;
            const Mat3x12 Bm  = Bmembrane(xl, yl, xi, eta);
            const Row12   Bd  = Bdrill(xl, yl, xi, eta);
            const Mat3x4  Bi  = BincQM6(J0inv, xi, eta);
            const Row4    Bdi = BdrillIncQM6(J0inv, xi, eta);
            Kca += w * t     * (Bm.transpose() * Dm * Bi);   // membrane coupling
            Kca += w * gamma * (Bd.transpose() * Bdi);        // drilling coupling
            Kaa += w * t     * (Bi.transpose() * Dm * Bi);
            Kaa += w * gamma * (Bdi.transpose() * Bdi);
        }
    return Kc - Kca * Kaa.inverse() * Kca.transpose();
}

// ---- DKQ discrete-Kirchhoff THIN-plate bending (S8-8b, opt-in) -------------------------
// Batoz & Tahar (1982). A 12-DOF [w,bx,by] x4 plate-bending block (same DOF as plateK -> the
// same plateToShellMap reuses it), but with NO transverse-shear DOF: the rotation field beta
// is an 8-node serendipity interpolation whose 4 mid-side rotations are eliminated by enforcing
// ZERO tangential shear along each edge, so beta is expressed from the corner (w,bx,by). Thin
// plates only (t/L < ~1/20); mid/thick plates must use the MITC4 plateK. Curvature
// chi = [bx,x; by,y; bx,y+by,x] matches plateK; Batoz's native DOF (w,theta_x,theta_y) map to
// mine by theta_x = by, theta_y = -bx (the beta-vs-slope sign squares out of the stiffness and
// is pinned for moment recovery by the constant-curvature patch test F49a).

// 8-node serendipity values + natural derivatives. Corners 1..4 at (+-1,+-1); mid nodes
// 5=(0,-1) edge(1,2), 6=(+1,0) edge(2,3), 7=(0,+1) edge(3,4), 8=(-1,0) edge(4,1).
inline void serendipity8(real xi, real eta, real N[8], real dNx[8], real dNe[8]) {
    const real xc[4] = { -1, 1, 1, -1 }, ec[4] = { -1, -1, 1, 1 };
    for (int i = 0; i < 4; ++i) {
        const real xp = xi * xc[i], ep = eta * ec[i];
        N[i]   = 0.25 * (1 + xp) * (1 + ep) * (xp + ep - 1);
        dNx[i] = 0.25 * xc[i] * (1 + ep) * (2 * xp + ep);
        dNe[i] = 0.25 * ec[i] * (1 + xp) * (xp + 2 * ep);
    }
    N[4] = 0.5 * (1 - xi * xi) * (1 - eta);    dNx[4] = -xi * (1 - eta);        dNe[4] = -0.5 * (1 - xi * xi);
    N[5] = 0.5 * (1 + xi) * (1 - eta * eta);   dNx[5] =  0.5 * (1 - eta * eta); dNe[5] = -eta * (1 + xi);
    N[6] = 0.5 * (1 - xi * xi) * (1 + eta);    dNx[6] = -xi * (1 + eta);        dNe[6] =  0.5 * (1 - xi * xi);
    N[7] = 0.5 * (1 - xi) * (1 - eta * eta);   dNx[7] = -0.5 * (1 - eta * eta); dNe[7] = -eta * (1 - xi);
}

// DKQ edge coefficients per edge k=0..3 connecting corners k and (k+1)%4 (Batoz-Tahar).
inline void dkqEdgeCoeffs(const real xl[4], const real yl[4],
                          real a[4], real b[4], real c[4], real d[4], real e[4]) {
    for (int k = 0; k < 4; ++k) {
        const int i = k, j = (k + 1) % 4;
        const real xij = xl[i] - xl[j], yij = yl[i] - yl[j];
        const real L2 = xij * xij + yij * yij;
        a[k] = -xij / L2;
        b[k] = 0.75 * xij * yij / L2;
        c[k] = (0.25 * xij * xij - 0.5 * yij * yij) / L2;
        d[k] = -yij / L2;
        e[k] = (0.25 * yij * yij - 0.5 * xij * xij) / L2;
    }
}

// Build Hx,Hy (each 12, Batoz DOF order (w,theta_x,theta_y) per node) from a serendipity
// value array s (s = N, or N,xi, or N,eta -> same formula gives Hx, Hx,xi, Hx,eta...).
inline void dkqH(const real s[8], const real a[4], const real b[4], const real c[4],
                 const real d[4], const real e[4], real Hx[12], real Hy[12]) {
    for (int n = 0; n < 4; ++n) {
        const int en = n;            // "next" edge (node n -> n+1), mid node s[4+en]
        const int ep = (n + 3) % 4;  // "prev" edge (node n-1 -> n), mid node s[4+ep]
        const real sn = s[4 + en], sp = s[4 + ep];
        Hx[3 * n + 0] = 1.5 * (a[en] * sn - a[ep] * sp);
        Hx[3 * n + 1] = b[en] * sn + b[ep] * sp;
        Hx[3 * n + 2] = s[n] - c[en] * sn - c[ep] * sp;
        Hy[3 * n + 0] = 1.5 * (d[en] * sn - d[ep] * sp);
        Hy[3 * n + 1] = -s[n] + e[en] * sn + e[ep] * sp;
        Hy[3 * n + 2] = -b[en] * sn - b[ep] * sp;
    }
}

// DKQ curvature B (3x12) in MY plate DOF [w,bx,by] x4 at (xi,eta). Built in Batoz (w,tx,ty)
// then column-remapped: B_mine(:,w)=B_B(:,w), B_mine(:,bx)=-B_B(:,ty), B_mine(:,by)=B_B(:,tx).
inline Mat3x12 BdkqMine(const real xl[4], const real yl[4], real xi, real eta) {
    real N[8], dNx[8], dNe[8];
    serendipity8(xi, eta, N, dNx, dNe);
    real a[4], b[4], c[4], d[4], e[4];
    dkqEdgeCoeffs(xl, yl, a, b, c, d, e);
    real Hxx[12], Hyx[12], Hxe[12], Hye[12];          // d/dxi (..x) and d/deta (..e)
    dkqH(dNx, a, b, c, d, e, Hxx, Hyx);
    dkqH(dNe, a, b, c, d, e, Hxe, Hye);
    real Jinv[2][2];
    jacobian(xl, yl, xi, eta, Jinv);
    Mat3x12 Bb;                                       // Batoz-order curvature
    for (int k = 0; k < 12; ++k) {
        const real Hx_x = Jinv[0][0] * Hxx[k] + Jinv[0][1] * Hxe[k];
        const real Hx_y = Jinv[1][0] * Hxx[k] + Jinv[1][1] * Hxe[k];
        const real Hy_x = Jinv[0][0] * Hyx[k] + Jinv[0][1] * Hye[k];
        const real Hy_y = Jinv[1][0] * Hyx[k] + Jinv[1][1] * Hye[k];
        Bb(0, k) = Hx_x;            // chi_x  = bx,x
        Bb(1, k) = Hy_y;           // chi_y  = by,y
        Bb(2, k) = Hx_y + Hy_x;    // chi_xy = bx,y + by,x
    }
    Mat3x12 B;                                        // remap Batoz (w,tx,ty) -> mine (w,bx,by)
    for (int n = 0; n < 4; ++n) {                     // bx=theta_y, by=-theta_x (pinned by F49a)
        B.col(3 * n + 0) =  Bb.col(3 * n + 0);        // w
        B.col(3 * n + 1) =  Bb.col(3 * n + 2);        // bx <-  theta_y
        B.col(3 * n + 2) = -Bb.col(3 * n + 1);        // by <- -theta_x
    }
    return B;
}

// DKQ plate-bending 12x12 in plate DOF [w,bx,by] x4 (no shear; 2x2 Gauss). G unused (Kirchhoff).
Mat12 plateK_DKQ(const real xl[4], const real yl[4], real E, real nu, real /*G*/, real t) {
    Eigen::Matrix<real, 3, 3> Db = Eigen::Matrix<real, 3, 3>::Zero();
    const real Dfac = E * t * t * t / (12.0 * (1.0 - nu * nu));
    Db(0, 0) = Dfac;       Db(0, 1) = nu * Dfac;
    Db(1, 0) = nu * Dfac;  Db(1, 1) = Dfac;
    Db(2, 2) = Dfac * (1.0 - nu) * 0.5;

    Mat12 K = Mat12::Zero();
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            const real xi = kG[a], eta = kG[b];
            real Jinv[2][2];
            const real detJ = jacobian(xl, yl, xi, eta, Jinv);
            const real wgt = kW[a] * kW[b] * detJ;
            const Mat3x12 Bp = BdkqMine(xl, yl, xi, eta);
            K += wgt * (Bp.transpose() * Db * Bp);
        }
    return K;
}

// 24x24 consistent mass in the facet's local DOF order [u,v,w,thx,thy,thz] x4. Translational
// inertia rho*t for u,v,w; rotary inertia rho*t^3/12 for thx,thy,thz (the drilling thz gets the
// same small rotary term so the mass matrix is positive-definite -> the modal eigenproblem is
// well-posed). `rho` is already in consistent units (tonne/mm^3).
Eigen::Matrix<real, 24, 24> shellMass24(const real xl[4], const real yl[4], real rho, real t) {
    Eigen::Matrix<real, 24, 24> M = Eigen::Matrix<real, 24, 24>::Zero();
    const real mTrans = rho * t;                 // per unit area, translational
    const real mRot   = rho * t * t * t / 12.0;  // per unit area, rotary
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            const real xi = kG[a], eta = kG[b];
            real Jinv[2][2];
            const real detJ = jacobian(xl, yl, xi, eta, Jinv);
            real N[4]; shapeN(xi, eta, N);
            const real wdet = kW[a] * kW[b] * detJ;
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j) {
                    const real NN = wdet * N[i] * N[j];
                    for (int d = 0; d < 3; ++d) M(6 * i + d, 6 * j + d) += mTrans * NN;  // u,v,w
                    for (int d = 3; d < 6; ++d) M(6 * i + d, 6 * j + d) += mRot * NN;    // thx,thy,thz
                }
        }
    return M;
}

} // anonymous namespace

// Map plate DOF [w,bx,by] x4 -> shell-local DOF [u,v,w,Rx,Ry,Rz] x4 (12x24):
//   w  = Uz   (local index 6k+2)
//   bx = Ry   (local index 6k+4)
//   by = -Rx  (local index 6k+3)
static Eigen::Matrix<real, 12, 24> plateToShellMap() {
    Eigen::Matrix<real, 12, 24> P = Eigen::Matrix<real, 12, 24>::Zero();
    for (int k = 0; k < 4; ++k) {
        P(3 * k + 0, 6 * k + 2) =  1.0;   // w   <- Uz
        P(3 * k + 1, 6 * k + 4) =  1.0;   // bx  <- Ry
        P(3 * k + 2, 6 * k + 3) = -1.0;   // by  <- -Rx
    }
    return P;
}

// Map membrane DOF [u,v,thz] x4 -> shell-local DOF [u,v,w,Rx,Ry,Rz] x4 (12x24):
//   u -> Ux (6k+0),  v -> Uy (6k+1),  thz -> Rz (6k+5).
static Eigen::Matrix<real, 12, 24> membraneToShellMap() {
    Eigen::Matrix<real, 12, 24> P = Eigen::Matrix<real, 12, 24>::Zero();
    for (int k = 0; k < 4; ++k) {
        P(3 * k + 0, 6 * k + 0) = 1.0;    // u   <- Ux
        P(3 * k + 1, 6 * k + 1) = 1.0;    // v   <- Uy
        P(3 * k + 2, 6 * k + 5) = 1.0;    // thz <- Rz
    }
    return P;
}

bool MITC4ShellElement::prepare(const FrameModel& model, const SolveOptions& opts, std::string& why) {
    const ShellQuad& sh = model.shells[static_cast<size_t>(s_)];
    id_ = sh.id;
    t_  = sh.t;
    const Material& smat = model.materials[static_cast<size_t>(sh.matIdx)];
    E_  = smat.E;
    nu_ = smat.nu;
    G_  = smat.G;

    int idx[4];
    for (int k = 0; k < 4; ++k) {
        idx[k] = model.nodeIndex(sh.n[k]);
        if (idx[k] < 0) { why = "shell references missing node"; return false; }
    }
    for (int k = 0; k < 4; ++k)
        for (int d = 0; d < 6; ++d)
            dofs_[6 * k + d] = gdof(idx[k], d);

    // ---- facet local frame: e1,e2 in-plane, n = normal. Rows of R_ are the local
    // axes in global coords, so R_ * v_global = v_local. ----
    const Vec3 P0 = model.nodes[idx[0]].pos, P1 = model.nodes[idx[1]].pos;
    const Vec3 P2 = model.nodes[idx[2]].pos, P3 = model.nodes[idx[3]].pos;
    Vec3 n = cross(P2 - P0, P3 - P1);
    const real nlen = norm(n);
    if (nlen <= 0) { why = "degenerate shell quad (zero normal)"; return false; }
    n = n * (1.0 / nlen);
    Vec3 e1 = P1 - P0;
    e1 = e1 - n * dot(e1, n);                 // project edge onto facet plane
    const real e1len = norm(e1);
    if (e1len <= 0) { why = "degenerate shell quad (edge parallel to normal)"; return false; }
    e1 = e1 * (1.0 / e1len);
    const Vec3 e2 = cross(n, e1);

    R_(0, 0) = e1.x; R_(0, 1) = e1.y; R_(0, 2) = e1.z;
    R_(1, 0) = e2.x; R_(1, 1) = e2.y; R_(1, 2) = e2.z;
    R_(2, 0) = n.x;  R_(2, 1) = n.y;  R_(2, 2) = n.z;

    // Corner coordinates in the facet 2D frame (origin at P0; planar projection).
    const Vec3 Pg[4] = { P0, P1, P2, P3 };
    for (int k = 0; k < 4; ++k) {
        const Vec3 r = Pg[k] - P0;
        xl_[k] = dot(r, e1);
        yl_[k] = dot(r, e2);
    }

    // ---- transform T_ = blockdiag(R_) x8 (disp + rot block per node). ----
    T_.setZero();
    for (int blk = 0; blk < 8; ++blk)
        T_.block(3 * blk, 3 * blk, 3, 3) = R_;

    // ---- 24x24 local stiffness: plate bending mapped into the bending DOFs +
    // membrane (plane stress) + drilling mapped into the in-plane DOFs. ----
    useQM6_ = opts.useIncompatibleMembrane;   // 8a: QM6 incompatible-mode membrane (opt-in)
    useDKQ_ = opts.useDKQPlate;               // 8b: DKQ discrete-Kirchhoff thin plate (opt-in)
    const Mat12 Kp = useDKQ_ ? plateK_DKQ(xl_, yl_, E_, nu_, G_, t_)
                             : plateK(xl_, yl_, E_, nu_, G_, t_);
    const Mat12 Km = useQM6_ ? membraneK_QM6(xl_, yl_, E_, nu_, G_, t_)
                             : membraneK(xl_, yl_, E_, nu_, G_, t_);
    const Eigen::Matrix<real, 12, 24> Pb = plateToShellMap();
    const Eigen::Matrix<real, 12, 24> Pm = membraneToShellMap();
    kl_.setZero();
    kl_ += Pb.transpose() * Kp * Pb;
    kl_ += Pm.transpose() * Km * Pm;

    // consistent mass (rho kg/m^3 -> tonne/mm^3 via 1e-12) for modal analysis
    ml_ = shellMass24(xl_, yl_, smat.rho * 1.0e-12, t_);

    // ---- transverse pressure -> consistent nodal loads on the local-w (Uz) DOFs.
    // Qf_ holds the LOCAL equivalent nodal load; addEquivalentNodalLoads rotates it
    // to global via F += T^T Qf_. ----
    Qf_.setZero();
    real pTot = 0;
    for (const auto& sp : model.shellPressures)
        if (sp.shell == id_) pTot += sp.p;
    if (pTot != 0.0) {
        for (int a = 0; a < 2; ++a)
            for (int b = 0; b < 2; ++b) {
                const real xi = kG[a], eta = kG[b];
                real Jinv[2][2];
                const real detJ = jacobian(xl_, yl_, xi, eta, Jinv);
                real N[4]; shapeN(xi, eta, N);
                const real wgt = kW[a] * kW[b] * detJ;
                for (int i = 0; i < 4; ++i)
                    Qf_(6 * i + 2) += N[i] * pTot * wgt;   // load along local +z
            }
    }
    return true;
}

void MITC4ShellElement::assemble(std::vector<Triplet>& trips) const {
    const Mat24 kg = T_.transpose() * kl_ * T_;
    for (int a = 0; a < 24; ++a)
        for (int b = 0; b < 24; ++b)
            if (kg(a, b) != 0.0) trips.emplace_back(dofs_[a], dofs_[b], kg(a, b));
}

void MITC4ShellElement::assembleMass(std::vector<Triplet>& trips) const {
    const Mat24 mg = T_.transpose() * ml_ * T_;
    for (int a = 0; a < 24; ++a)
        for (int b = 0; b < 24; ++b)
            if (mg(a, b) != 0.0) trips.emplace_back(dofs_[a], dofs_[b], mg(a, b));
}

void MITC4ShellElement::addEquivalentNodalLoads(VecX& F) const {
    if (Qf_.isZero(0)) return;
    const Vec24 fg = T_.transpose() * Qf_;     // local nodal load -> global
    for (int a = 0; a < 24; ++a) F(dofs_[a]) += fg(a);
}

void MITC4ShellElement::recover(const VecX& u, SolveResult& R) const {
    // Gather global element displacements, rotate to local, project to plate DOFs.
    Vec24 ug;
    for (int a = 0; a < 24; ++a) ug(a) = u(dofs_[a]);
    const Vec24 ul = T_ * ug;
    const Eigen::Matrix<real, 12, 24> Pb = plateToShellMap();
    const Eigen::Matrix<real, 12, 24> Pmem = membraneToShellMap();
    const Eigen::Matrix<real, 12, 1> dp = Pb * ul;     // plate (w,bx,by)
    const Eigen::Matrix<real, 12, 1> dm = Pmem * ul;   // membrane (u,v,thz)

    // Constitutive at the element centre.
    Eigen::Matrix<real, 3, 3> Db = Eigen::Matrix<real, 3, 3>::Zero();
    const real Dfac = E_ * t_ * t_ * t_ / (12.0 * (1.0 - nu_ * nu_));
    Db(0, 0) = Dfac;       Db(0, 1) = nu_ * Dfac;
    Db(1, 0) = nu_ * Dfac; Db(1, 1) = Dfac;
    Db(2, 2) = Dfac * (1.0 - nu_) * 0.5;
    const real Ds = kShearCorr * G_ * t_;

    Eigen::Matrix<real, 3, 3> Dm = Eigen::Matrix<real, 3, 3>::Zero();
    const real f = E_ / (1.0 - nu_ * nu_);
    Dm(0, 0) = f;       Dm(0, 1) = nu_ * f;
    Dm(1, 0) = nu_ * f; Dm(1, 1) = f;
    Dm(2, 2) = f * (1.0 - nu_) * 0.5;

    const Mat3x12 Bb = useDKQ_ ? BdkqMine(xl_, yl_, 0.0, 0.0) : Bbending(xl_, yl_, 0.0, 0.0);
    const Mat3x12 Bm = Bmembrane(xl_, yl_, 0.0, 0.0);
    const Eigen::Matrix<real, 3, 1> M = Db * (Bb * dp);
    Eigen::Matrix<real, 2, 1> Q;
    if (useDKQ_) {
        Q.setZero();                                   // DKQ is Kirchhoff: no transverse shear
    } else {
        const Mat2x12 Bs = Bshear(xl_, yl_, 0.0, 0.0);
        Q = Ds * (Bs * dp);
    }
    const Eigen::Matrix<real, 3, 1> N = t_ * (Dm * (Bm * dm));

    ShellElementForces sf;
    sf.shell = id_;
    sf.Mxx = M(0); sf.Myy = M(1); sf.Mxy = M(2);
    sf.Qx  = Q(0); sf.Qy  = Q(1);
    sf.Nxx = N(0); sf.Nyy = N(1); sf.Nxy = N(2);
    // Per-corner bending moments (natural coords of the 4 CCW corners) for design peak
    // recovery. Linear field -> these combine/envelope correctly; a constant-moment field
    // gives every corner == the centre value above.
    // DKQ CAVEAT: the discrete-Kirchhoff B has large serendipity derivatives at the corners, so
    // BdkqMine at a corner is NOT a pointwise curvature estimator (a constant-curvature field does
    // NOT give corner == centre for DKQ). The per-corner values stay LINEAR in dp (combine/envelope
    // remain valid), but for a DKQ design PEAK use the centre value (Gauss-point extrapolation is
    // NOT IMPLEMENTED) -- do not read DKQ MxxC[k] as the moment AT that corner.
    const real cxi[4]  = { -1.0, 1.0, 1.0, -1.0 };
    const real ceta[4] = { -1.0, -1.0, 1.0, 1.0 };
    for (int k = 0; k < 4; ++k) {
        const Mat3x12 Bbk = useDKQ_ ? BdkqMine(xl_, yl_, cxi[k], ceta[k]) : Bbending(xl_, yl_, cxi[k], ceta[k]);
        const Eigen::Matrix<real, 3, 1> Mk = Db * (Bbk * dp);
        sf.MxxC[k] = Mk(0); sf.MyyC[k] = Mk(1); sf.MxyC[k] = Mk(2);
    }
    R.shellForces[static_cast<size_t>(s_)] = sf;
}

} // namespace frame
