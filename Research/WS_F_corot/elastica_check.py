"""elastica_check.py — WS-F large-deflection cantilever oracle table (scratch, NOT engine).

Cantilever, point load P perpendicular to the undeformed axis at the tip.
Nondimensional load alpha = P L^2 / EI. With theta(s) = slope angle (positive toward the
load) and s in [0,1] (arc length / L):

    theta''(s) = -alpha * cos(theta(s)),   theta(0) = 0,  theta'(1) = 0
    (moment M(s) = P (x_tip - x(s)), differentiated; tip moment = 0)

Shooting on theta'(0) with high-accuracy RK integration. Outputs delta_v/L (load-direction
tip deflection), delta_h/L (axial shortening) and theta_tip for alpha in 0.5..10 — the
future F-fixture oracle table for the co-rotational stage (S9). Cross-check against
published elastica tables (e.g. Mattiasson 1981) happens in the research doc, not here.

Self-validation: the table is printed at two integrator tolerances; matching digits are
the trustworthy ones.
"""

import math

try:
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


def integrate(alpha, m0, rtol):
    """Integrate theta'' = -alpha*cos(theta) from s=0..1; return theta'(1), and the tip
    integrals x(1)=int cos, w(1)=int sin."""
    if HAVE_SCIPY:
        def rhs(s, y):
            return [y[1], -alpha * math.cos(y[0]), math.cos(y[0]), math.sin(y[0])]
        sol = solve_ivp(rhs, (0.0, 1.0), [0.0, m0, 0.0, 0.0],
                        rtol=rtol, atol=rtol * 1e-3, dense_output=False, method="DOP853")
        y = sol.y[:, -1]
        return y[1], y[2], y[3]
    # fallback: classic RK4, fixed step
    n = int(4.0 / rtol ** 0.25)  # crude: rtol 1e-10 -> ~ 1.3e3 steps... bump:
    n = max(n, 200000)
    h = 1.0 / n
    th, m, x, w = 0.0, m0, 0.0, 0.0

    def f(th, m):
        return m, -alpha * math.cos(th)

    for _ in range(n):
        k1t, k1m = f(th, m)
        k2t, k2m = f(th + 0.5 * h * k1t, m + 0.5 * h * k1m)
        k3t, k3m = f(th + 0.5 * h * k2t, m + 0.5 * h * k2m)
        k4t, k4m = f(th + h * k3t, m + h * k3m)
        x += h / 6.0 * (math.cos(th) + 2 * math.cos(th + 0.5 * h * k1t)
                        + 2 * math.cos(th + 0.5 * h * k2t) + math.cos(th + h * k3t))
        w += h / 6.0 * (math.sin(th) + 2 * math.sin(th + 0.5 * h * k1t)
                        + 2 * math.sin(th + 0.5 * h * k2t) + math.sin(th + h * k3t))
        th += h / 6.0 * (k1t + 2 * k2t + 2 * k3t + k4t)
        m += h / 6.0 * (k1m + 2 * k2m + 2 * k3m + k4m)
    return m, x, w


def solve_alpha(alpha, rtol):
    lo, hi = 1e-9, alpha  # theta'(0) = alpha*x_tip <= alpha
    f_lo = integrate(alpha, lo, rtol)[0]
    f_hi = integrate(alpha, hi, rtol)[0]
    if f_lo * f_hi > 0:
        hi = alpha * 1.5
        f_hi = integrate(alpha, hi, rtol)[0]
    if HAVE_SCIPY:
        m0 = brentq(lambda m: integrate(alpha, m, rtol)[0], lo, hi, xtol=1e-14, rtol=8.9e-16)
    else:
        for _ in range(200):  # bisection
            mid = 0.5 * (lo + hi)
            fm = integrate(alpha, mid, rtol)[0]
            if f_lo * fm <= 0:
                hi, f_hi = mid, fm
            else:
                lo, f_lo = mid, fm
        m0 = 0.5 * (lo + hi)
    mres, x1, w1 = integrate(alpha, m0, rtol)
    theta_tip = None
    # tip angle: integrate once more storing theta(1)
    if HAVE_SCIPY:
        def rhs(s, y):
            return [y[1], -alpha * math.cos(y[0])]
        sol = solve_ivp(rhs, (0.0, 1.0), [0.0, m0], rtol=rtol, atol=rtol * 1e-3, method="DOP853")
        theta_tip = sol.y[0, -1]
    return w1, 1.0 - x1, theta_tip, mres


def main():
    print(f"# elastica cantilever tip-load table (scipy={HAVE_SCIPY})")
    print("# alpha=PL^2/EI   deltaV/L      deltaH/L      thetaTip[rad]  residual")
    for rtol in (1e-10, 1e-12):
        print(f"## rtol={rtol:g}")
        for alpha in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0):
            dv, dh, tt, res = solve_alpha(alpha, rtol)
            tts = f"{tt:.10f}" if tt is not None else "n/a"
            print(f"{alpha:5.2f}  {dv:.10f}  {dh:.10f}  {tts}  {res:.2e}")
    # small-load sanity: alpha -> 0 gives the linear answer deltaV/L = alpha/3
    alpha = 1e-4
    dv, dh, tt, res = solve_alpha(alpha, 1e-12)
    print(f"# linear-limit check: alpha={alpha:g} deltaV/L={dv:.12e} vs alpha/3={alpha/3:.12e} "
          f"rel={(abs(dv - alpha/3)/(alpha/3)):.2e}")


if __name__ == "__main__":
    main()
