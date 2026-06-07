#include "MITC4ShellElement.h"
#include "FrameCore/FrameModel.h"
#include "FrameCore/SolveResult.h"

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

bool MITC4ShellElement::prepare(const FrameModel& model, const SolveOptions& /*opts*/, std::string& why) {
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
    const Mat12 Kp = plateK(xl_, yl_, E_, nu_, G_, t_);
    const Mat12 Km = membraneK(xl_, yl_, E_, nu_, G_, t_);
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

    const Mat3x12 Bb = Bbending(xl_, yl_, 0.0, 0.0);
    const Mat2x12 Bs = Bshear(xl_, yl_, 0.0, 0.0);
    const Mat3x12 Bm = Bmembrane(xl_, yl_, 0.0, 0.0);
    const Eigen::Matrix<real, 3, 1> M = Db * (Bb * dp);
    const Eigen::Matrix<real, 2, 1> Q = Ds * (Bs * dp);
    const Eigen::Matrix<real, 3, 1> N = t_ * (Dm * (Bm * dm));

    ShellElementForces sf;
    sf.shell = id_;
    sf.Mxx = M(0); sf.Myy = M(1); sf.Mxy = M(2);
    sf.Qx  = Q(0); sf.Qy  = Q(1);
    sf.Nxx = N(0); sf.Nyy = N(1); sf.Nxy = N(2);
    // Per-corner bending moments (natural coords of the 4 CCW corners) for design peak
    // recovery. Linear field -> these combine/envelope correctly; a constant-moment field
    // gives every corner == the centre value above.
    const real cxi[4]  = { -1.0, 1.0, 1.0, -1.0 };
    const real ceta[4] = { -1.0, -1.0, 1.0, 1.0 };
    for (int k = 0; k < 4; ++k) {
        const Mat3x12 Bbk = Bbending(xl_, yl_, cxi[k], ceta[k]);
        const Eigen::Matrix<real, 3, 1> Mk = Db * (Bbk * dp);
        sf.MxxC[k] = Mk(0); sf.MyyC[k] = Mk(1); sf.MxyC[k] = Mk(2);
    }
    R.shellForces[static_cast<size_t>(s_)] = sf;
}

} // namespace frame
