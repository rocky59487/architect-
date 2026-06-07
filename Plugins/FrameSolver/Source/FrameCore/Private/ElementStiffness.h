#pragma once
#include "FrameEigen.h"
#include "FrameCore/Member.h"
#include "FrameCore/Node.h"

namespace frame {

// Rows of the returned matrix are the local axis unit vectors (x, y, z) in global
// coordinates, so R * v_global = v_local. Guards the vertical-member degeneracy:
// when refVec is (near) parallel to the member axis it falls back to global +Y.
Mat3 localAxes(const Vec3& pi, const Vec3& pj, const Vec3& refVec);

// 12x12 local Euler-Bernoulli beam-column stiffness in the canonical DOF order
// [u1 v1 w1 rx1 ry1 rz1 u2 v2 w2 rx2 ry2 rz2].
Mat12 localStiffness12(real E, real G, real A, real Iy, real Iz, real J, real L);

// Timoshenko variant: adds shear flexibility through the effective shear areas
// Asy (resists local-y shear, pairs with the Iz bending block) and Asz (local-z
// shear, pairs with Iy). The closed-form "interdependent interpolation" stiffness
// is exact (no shear locking) and reduces to localStiffness12 as Asy,Asz -> inf
// (shear parameter Phi -> 0). Axial and torsion are identical to the EB element.
Mat12 localStiffness12T(real E, real G, real A, real Iy, real Iz, real J, real L,
                        real Asy, real Asz);

// 12x12 transform T = blockdiag(R, R, R, R); u_local = T * u_global, k_g = T^T k_l T.
Mat12 transform12(const Mat3& R);

// 12x12 CONSISTENT mass matrix (Euler-Bernoulli, translational; no cross-section rotary
// inertia, matching the analytic frequency omega_n = (n.pi/L)^2 sqrt(EI/rho.A)). Same DOF
// order + sign convention as localStiffness12. `rho` MUST already be in the engine's
// consistent units (tonne/mm^3 = rho[kg/m^3] * 1e-12). Torsion uses the polar Ip = Iy + Iz.
Mat12 localMass12(real rho, real A, real Iy, real Iz, real L);

// Static condensation of released local DOFs out of the element stiffness AND the
// fixed-end forces. release[k]=true => local DOF k transmits no force (hinge/pin):
//   kl*_rr = kl_rr - kl_rc kcc^-1 kl_cr ;  Qf*_r = Qf_r - kl_rc kcc^-1 Qf_c ;  rest = 0.
// Modifies kl and Qf IN PLACE. Returns false (leaving kl/Qf UNCHANGED) when the
// released sub-block kcc is singular -- a free mechanism (e.g. both torsional ends
// released) -- so the caller can flag the model singular with a clear diagnostic.
bool condenseReleases(Mat12& kl, Vec12& Qf, const std::array<bool, 12>& release);

} // namespace frame
