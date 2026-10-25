# teme.py - TEME -> J2000 (MEME/EME2000) for CircuitPython
#
# Input time is a UTC Julian date split as jd + fr (jd the day part, normally
# a half-integer such as 2460370.5; fr the fraction of a day), the same split
# SGP4 uses.
#
# TEME -> TOD (equation of equinoxes) -> MOD (nutation, IAU-1980,
# truncated) -> J2000 (precession, IAU-1976). Works
# the same for position and velocity vectors.
#
# note: CircuitPython floats are often 30-bit (~7 digits). Anything
# that would lose precision (the Julian date itself, multi-revolution
# angles) is done with Python integers, which are arbitrary precision.

import math

_ARCSEC = math.pi / (180.0 * 3600.0)
_NS_PER_CENTURY = 36525 * 86400 * 10**9

# Fundamental arguments (IAU-1980): constant deg, linear in integer
# micro-degrees/century, T^2 deg, T^3 deg.
_ARGS = (
    (134.96298139, 477198867398, 0.0086972, 1.78e-5),    # l   (Moon anomaly)
    (357.52772333, 35999050340, -0.0001603, -3.3e-6),    # l'  (Sun anomaly)
    (93.27191028, 483202017538, -0.0036825, 3.1e-6),     # F
    (297.85036306, 445267111480, -0.0019142, 5.3e-6),    # D
    (125.04452222, -1934136261, 0.0020708, 2.2e-6),      # Omega
)

# Largest IAU-1980 terms: l, l', F, D, Om, dpsi0, dpsi1, deps0, deps1
# (coefficients in 0.0001 arcsec, dpsi1/deps1 per century). Truncated series:
# error is roughly a meter at LEO, well under SGP4's own error.
_NUT = (
    (0, 0, 0, 0, 1, -171996, -174.2, 92025, 8.9),
    (0, 0, 2, -2, 2, -13187, -1.6, 5736, -3.1),
    (0, 0, 2, 0, 2, -2274, -0.2, 977, -0.5),
    (0, 0, 0, 0, 2, 2062, 0.2, -895, 0.5),
    (0, 1, 0, 0, 0, 1426, -3.4, 54, -0.1),
    (1, 0, 0, 0, 0, 712, 0.1, -7, 0.0),
    (0, 1, 2, -2, 2, -517, 1.2, 224, -0.6),
    (0, 0, 2, 0, 1, -386, -0.4, 200, 0.0),
    (1, 0, 2, 0, 2, -301, 0.0, 129, -0.1),
    (0, -1, 2, -2, 2, 217, -0.5, -95, 0.3),
)


def tt_ns_from_jd(jd, fr, leap_seconds=37):
    """Integer ns of TT since J2000.0 from a UTC Julian date split as jd + fr.

    jd is snapped to the nearest half day and handled as an integer, so the
    large number never goes through float arithmetic; fr only needs ~1e-7 day
    (about 10 ms) of precision. leap_seconds is TAI-UTC (37 s since
    2017-01-01 - check IERS if that may have changed). If your jd/fr is
    already TT, pass leap_seconds=-32.184.
    """
    half_days = int(round(jd * 2.0)) - 2 * 2451545   # exact for x.0 / x.5
    ns = half_days * 43200 * 10**9
    ns += int(round(fr * 86400e9))
    ns += int(round((leap_seconds + 32.184) * 1e9))
    return ns


def _fund_args(ns, T):
    out = []
    for c0, c1, c2, c3 in _ARGS:
        # linear term with exact integer math, reduced mod 360 deg
        lin = ((c1 * ns) // _NS_PER_CENTURY) % 360000000
        deg = c0 + lin / 1e6 + T * T * (c2 + T * c3)
        out.append(math.radians(deg))
    return out


def _angles(ns_tt):
    """Return (eps_mean, eps_true, dpsi, eqeq, zeta, theta, z) in radians."""
    T = ns_tt / _NS_PER_CENTURY  # float is fine for slowly varying terms

    # Mean obliquity (IAU-1976)
    eps_mean = (84381.448 - 46.8150 * T - 0.00059 * T * T + 0.001813 * T * T * T) * _ARCSEC

    # Nutation
    args = _fund_args(ns_tt, T)
    dpsi = 0.0
    deps = 0.0
    for nl, nlp, nf, nd, nom, a0, a1, b0, b1 in _NUT:
        arg = nl * args[0] + nlp * args[1] + nf * args[2] + nd * args[3] + nom * args[4]
        dpsi += (a0 + a1 * T) * math.sin(arg)
        deps += (b0 + b1 * T) * math.cos(arg)
    dpsi *= 1e-4 * _ARCSEC
    deps *= 1e-4 * _ARCSEC
    eps_true = eps_mean + deps

    # Equation of the equinoxes (with the post-1997 complementary terms)
    om = args[4]
    eqeq = dpsi * math.cos(eps_mean) + (0.00264 * math.sin(om) + 0.000063 * math.sin(2 * om)) * _ARCSEC

    # Precession angles (IAU-1976)
    zeta = (2306.2181 * T + 0.30188 * T * T + 0.017998 * T**3) * _ARCSEC
    theta = (2004.3109 * T - 0.42665 * T * T - 0.041833 * T**3) * _ARCSEC
    z = (2306.2181 * T + 1.09468 * T * T + 0.018203 * T**3) * _ARCSEC

    return eps_mean, eps_true, dpsi, eqeq, zeta, theta, z


# Elementary frame rotations (passive, Vallado's ROT1/ROT2/ROT3)
def _r1(v, a):
    c, s = math.cos(a), math.sin(a)
    return (v[0], c * v[1] + s * v[2], -s * v[1] + c * v[2])


def _r2(v, a):
    c, s = math.cos(a), math.sin(a)
    return (c * v[0] - s * v[2], v[1], s * v[0] + c * v[2])


def _r3(v, a):
    c, s = math.cos(a), math.sin(a)
    return (c * v[0] + s * v[1], -s * v[0] + c * v[1], v[2])


def teme_to_j2000(vec, jd, fr, leap_seconds=37):
    """Rotate a TEME vector (km or km/s) to J2000.

    vec : (x, y, z)
    jd, fr : UTC Julian date split (e.g. 2460370.5, 0.5 for 2024-03-01 12:00)
    """
    eps_mean, eps_true, dpsi, eqeq, zeta, theta, z = _angles(
        tt_ns_from_jd(jd, fr, leap_seconds)
    )
    r = tuple(vec)
    # TEME -> TOD: true equinox sits eqeq behind the mean equinox
    r = _r3(r, -eqeq)
    # TOD -> MOD: undo nutation
    r = _r1(r, eps_true)
    r = _r3(r, dpsi)
    r = _r1(r, -eps_mean)
    # MOD -> J2000: undo precession
    r = _r3(r, z)
    r = _r2(r, -theta)
    r = _r3(r, zeta)
    return r


if __name__ == "__main__":
    print(teme_to_j2000((7000.0, 0.0, 0.0), 2460370.5, 0.5))
