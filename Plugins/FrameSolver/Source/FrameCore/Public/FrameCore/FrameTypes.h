#pragma once
//
// FrameCore — authoritative conventions (engine-agnostic, no Eigen / no Unreal here).
//
// Units (this milestone): N, mm, MPa  (MPa = N/mm^2). All fixture checks are
// relative, so the unit choice does not affect tolerances.
//
// DOF order PER NODE (authoritative):   [ Ux, Uy, Uz, Rx, Ry, Rz ]  = 0..5
// Global DOF index:                      6 * nodeIndex + localDof
// Element DOF vector (12):               [ node_i: 6 ][ node_j: 6 ]
//
// Sign conventions:
//   * Axial N is COMPRESSION-POSITIVE (a vertical column under gravity reads N > 0).
//   * Member end forces are reported in LOCAL coordinates.
//   * Reactions are reported in GLOBAL coordinates.
//   * Local axes: local x = unit(j - i); the member refVec defines the local x-y
//     plane (default global +Z, fallback +Y when parallel to the member).
//
#include <cmath>

// Cross-module export (dual build). In a UE build FRAMECORE_UE is defined and
// UnrealBuildTool provides FRAMECORE_API (dllexport when building FrameCore,
// dllimport for consumers). Standalone has no DLL boundary -> define it empty.
#ifndef FRAMECORE_UE
  #define FRAMECORE_API
#endif

namespace frame {

using real = double;                 // MUST be double (see plan Gotcha #3)

enum Dof : int { Ux = 0, Uy = 1, Uz = 2, Rx = 3, Ry = 4, Rz = 5 };
inline constexpr int DOF_PER_NODE = 6;
inline constexpr int ELEM_DOF     = 12;   // [node_i 6][node_j 6]

using NodeId   = int;
using MemberId = int;

inline int gdof(int nodeIndex, int localDof) { return DOF_PER_NODE * nodeIndex + localDof; }

// Minimal POD 3-vector — keeps the public API free of Eigen and Unreal types.
struct Vec3 {
    real x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(real X, real Y, real Z) : x(X), y(Y), z(Z) {}
};

inline Vec3 operator-(const Vec3& a, const Vec3& b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
inline Vec3 operator+(const Vec3& a, const Vec3& b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
inline Vec3 operator*(const Vec3& a, real s)        { return { a.x * s, a.y * s, a.z * s }; }
inline Vec3 operator*(real s, const Vec3& a)        { return { s * a.x, s * a.y, s * a.z }; }
inline real dot(const Vec3& a, const Vec3& b)       { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline real norm(const Vec3& a)                     { return std::sqrt(dot(a, a)); }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
}

inline constexpr real kPi   = real(3.14159265358979323846);
inline constexpr real twoPi = real(2) * kPi;

} // namespace frame
