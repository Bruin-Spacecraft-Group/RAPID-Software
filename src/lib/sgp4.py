"""
Simplified General Perturbations Model 4 Implementation for Orbit Propagation
Original code by Brandon-Rhodes python-sgp4

As outlined in https://celestrak.org/publications/AIAA/2008-6770/AIAA-2008-6770.pdf
"""

# Disable some linter issues across the *otherwise verified* file
# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals,too-many-statements
# pylint: disable=multiple-statements,unnecessary-semicolon

import math

try:
    import ulab.numpy as np  # For CircuitPython
except ImportError:
    import numpy as np  # For GitHub Actions / PC testing

# Error codes
ECCENTRICITY = 1 # eccentricity is not within 0-1
MOTION = 2 # error in propagating mean motion
SEMIRECT = 4 # apoapsis, periapsis characteristics error
DECAY = 6 # orbit has decayed

TAU = 2.0 * math.pi
X2O3 = 2.0 / 3.0

class SGP4Error(Exception):
    """
    Exception class for SGP4 errors as described with # Error codes at the top of the file
    """
    def __init__(self, code):
        super().__init__(code)
        self.code = code

def _initl(xke, j2,
            ecc, epoch, inclo, no):
    #  ----------------------- earth constants ----------------------
    #  sgp4fix identify constants and allow alternate values
    #  only xke and j2 are used here so pass them in directly
    #  tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = whichconst
    x2o3   = 2.0 / 3.0
    tau = np.pi*2

    #  ------------- calculate auxillary epoch quantities ----------
    eccsq  = ecc * ecc
    omeosq = 1.0 - eccsq
    rteosq = np.sqrt(omeosq)
    cosio  = np.cos(inclo)
    cosio2 = cosio * cosio

    #  ------------------ un-kozai the mean motion -----------------
    ak    = pow(xke / no, x2o3)
    d1    = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq)
    del_  = d1 / (ak * ak)
    adel  = ak * (1.0 - del_ * del_ - del_ *
            (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0))
    del_  = d1/(adel * adel)
    no    = no / (1.0 + del_)

    ao    = pow(xke / no, x2o3)
    sinio = np.sin(inclo)
    po    = ao * omeosq
    con42 = 1.0 - 5.0 * cosio2
    con41 = -con42-cosio2-cosio2
    posq  = po * po
    rp    = ao * (1.0 - ecc)

    ts70  = epoch - 7305.0
    ds70 = (ts70 + 1.0e-8) // 1.0
    tfrac = ts70 - ds70
    #  find greenwich location at epoch
    c1    = 1.72027916940703639e-2
    thgr70= 1.7321343856509374
    fk5r  = 5.07551419432269442e-15
    c1p2p = c1 + tau
    gsto  = (thgr70 + c1*ds70 + c1p2p*tfrac + ts70*ts70*fk5r) % tau
    if gsto < 0.0:
        gsto = gsto + tau
    # gsto = _gstime(epoch + 2433281.5)

    return (no,
       ao,    con41,  con42, cosio,
       cosio2, omeosq, posq,
       rp,    rteosq,sinio , gsto
    )


# -- SGP4_Update helper functions --
def _secular_update(satrec, t):
    """Secular gravity and drag. Returns (mm, argpm, nodem, tempa, tempe, templ)."""
    t2 = t * t
    xmdf = satrec.mo + satrec.mdot * t
    argpdf = satrec.argp + satrec.argpdot * t
    nodem = satrec.raan + satrec.nodedot * t + satrec.nodecf * t2
    tempa = 1.0 - satrec.cc1 * t
    tempe = satrec.bstar * satrec.cc4 * t
    templ = satrec.t2cof * t2
    mm = xmdf
    argpm = argpdf

    if satrec.isimp != 1:
        d_mm, d_argp, d_a, d_e, d_l = _higher_order_drag(satrec, xmdf, t, t2)
        mm += d_mm
        argpm += d_argp
        tempa += d_a
        tempe += d_e
        templ += d_l

    return mm, argpm, nodem, tempa, tempe, templ


def _higher_order_drag(satrec, xmdf, t, t2):
    """Extra terms for non-simplified (isimp != 1) orbits.
    Returns increments (d_mm, d_argp, d_a, d_e, d_l)."""
    delomg = satrec.omgcof * t
    delmtemp = 1.0 + satrec.eta * math.cos(xmdf)
    delm = satrec.xmcof * (delmtemp * delmtemp * delmtemp - satrec.delmo)
    temp = delomg + delm

    t3 = t2 * t
    t4 = t3 * t

    d_a = -(satrec.d2 * t2 + satrec.d3 * t3 + satrec.d4 * t4)
    d_e = satrec.bstar * satrec.cc5 * (math.sin(xmdf + temp) - satrec.sinmao)
    d_l = satrec.t3cof * t3 + t4 * (satrec.t4cof + t * satrec.t5cof)
    return temp, -temp, d_a, d_e, d_l


def _mean_elements(satrec, mm, argpm, nodem, tempa, tempe, templ):
    """Validate and wrap mean elements. Returns (am, em, nm, mm, argpm, nodem)."""
    if satrec.n <= 0.0:
        raise SGP4Error(MOTION)

    am = math.pow(satrec.xke / satrec.n, X2O3) * tempa * tempa
    nm = satrec.xke / math.pow(am, 1.5)
    em = satrec.ecc - tempe

    if em >= 1.0 or em < 0.0:
        raise SGP4Error(ECCENTRICITY)
    em = max(em, 1e-06)

    mm = mm + satrec.n * templ
    xlm = (mm + argpm + nodem) % TAU
    nodem = nodem % TAU
    argpm = argpm % TAU
    mm = (xlm - argpm - nodem) % TAU
    return am, em, nm, mm, argpm, nodem


def _store_mean(satrec, am, em, nm, mm, argpm, nodem):
    satrec.am = am
    satrec.em = em
    satrec.im = satrec.inclo
    satrec.Om = nodem
    satrec.om = argpm
    satrec.mm = mm
    satrec.nm = nm


def _long_period_periodics(satrec, am, em, mm, argpm, nodem):
    axnl = em * math.cos(argpm)
    temp = 1.0 / (am * (1.0 - em * em))
    aynl = em * math.sin(argpm) + temp * satrec.aycof
    xl = mm + argpm + nodem + temp * satrec.xlcof * axnl
    return axnl, aynl, xl


def _solve_kepler(xl, nodep, axnl, aynl, tol=1.0e-12, max_iter=10):
    """Newton iteration. Returns (sineo1, coseo1) from the last iteration's
    input, matching the reference implementation exactly."""
    u = (xl - nodep) % TAU
    eo1 = u
    sineo1 = coseo1 = 0.0

    for _ in range(max_iter):
        sineo1 = math.sin(eo1)
        coseo1 = math.cos(eo1)
        step = 1.0 - coseo1 * axnl - sineo1 * aynl
        step = (u - aynl * coseo1 + axnl * sineo1 - eo1) / step
        if step >= 0.95:
            step = 0.95
        elif step <= -0.95:
            step = -0.95
        eo1 += step
        if abs(step) < tol:
            break

    return sineo1, coseo1


def _short_period_prelim(am, axnl, aynl, sineo1, coseo1):
    ecose = axnl * coseo1 + aynl * sineo1
    esine = axnl * sineo1 - aynl * coseo1
    el2 = axnl * axnl + aynl * aynl
    pl = am * (1.0 - el2)
    if pl < 0.0:
        raise SGP4Error(SEMIRECT)

    rl = am * (1.0 - ecose)
    rdotl = math.sqrt(am) * esine / rl
    rvdotl = math.sqrt(pl) / rl
    betal = math.sqrt(1.0 - el2)
    temp = esine / (1.0 + betal)
    sinu = am / rl * (sineo1 - aynl - axnl * temp)
    cosu = am / rl * (coseo1 - axnl + aynl * temp)
    su = math.atan2(sinu, cosu)
    sin2u = (cosu + cosu) * sinu
    cos2u = 1.0 - 2.0 * sinu * sinu
    return rl, rdotl, rvdotl, betal, su, sin2u, cos2u, pl


def _short_period_corrections(satrec, nm, inclm, nodep,
                              rl, rdotl, rvdotl, betal, su, sin2u, cos2u, pl):
    sinip = math.sin(inclm)
    cosip = math.cos(inclm)

    temp = 1.0 / pl
    temp1 = 0.5 * satrec.j2 * temp
    temp2 = temp1 * temp

    mrt = (rl * (1.0 - 1.5 * temp2 * betal * satrec.con41)
           + 0.5 * temp1 * satrec.x1mth2 * cos2u)
    if mrt < 1.0:
        raise SGP4Error(DECAY)

    su = su - 0.25 * temp2 * satrec.x7thm1 * sin2u
    xnode = nodep + 1.5 * temp2 * cosip * sin2u
    xinc = inclm + 1.5 * temp2 * cosip * sinip * cos2u
    mvt = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / satrec.xke
    rvdot = rvdotl + nm * temp1 * (
        satrec.x1mth2 * cos2u + 1.5 * satrec.con41) / satrec.xke
    return mrt, su, xnode, xinc, mvt, rvdot


def _orientation_vectors(su, xnode, xinc):
    """Returns (ux, uy, uz, vx, vy, vz)."""
    sinsu = math.sin(su)
    cossu = math.cos(su)
    snod = math.sin(xnode)
    cnod = math.cos(xnode)
    sini = math.sin(xinc)
    cosi = math.cos(xinc)
    xmx = -snod * cosi
    xmy = cnod * cosi
    return (xmx * sinsu + cnod * cossu,
            xmy * sinsu + snod * cossu,
            sini * sinsu,
            xmx * cossu - cnod * sinsu,
            xmy * cossu - snod * sinsu,
            sini * cossu)


def _to_cartesian(satrec, mrt, mvt, rvdot, basis):
    ux, uy, uz, vx, vy, vz = basis
    vkmpersec = satrec.radiusearthkm * satrec.xke / 60.0
    mr = mrt * satrec.radiusearthkm
    r = (mr * ux, mr * uy, mr * uz)
    v = ((mvt * ux + rvdot * vx) * vkmpersec,
         (mvt * uy + rvdot * vy) * vkmpersec,
         (mvt * uz + rvdot * vz) * vkmpersec)
    return r, v


def sgp4_update(satrec, tsince):
    """Propagate satrec to tsince (minutes). Returns (r, v) in km, km/s,
    or (False, False) with satrec.error set on failure."""
    satrec.t = tsince
    try:
        mm, argpm, nodem, tempa, tempe, templ = _secular_update(satrec, tsince)
        am, em, nm, mm, argpm, nodem = _mean_elements(
            satrec, mm, argpm, nodem, tempa, tempe, templ)
        _store_mean(satrec, am, em, nm, mm, argpm, nodem)

        axnl, aynl, xl = _long_period_periodics(satrec, am, em, mm, argpm, nodem)
        sineo1, coseo1 = _solve_kepler(xl, nodem, axnl, aynl)
        rl, rdotl, rvdotl, betal, su, sin2u, cos2u, pl = _short_period_prelim(
            am, axnl, aynl, sineo1, coseo1)
        mrt, su, xnode, xinc, mvt, rvdot = _short_period_corrections(
            satrec, nm, satrec.inclo, nodem,
            rl, rdotl, rvdotl, betal, su, sin2u, cos2u, pl)
        basis = _orientation_vectors(su, xnode, xinc)
        satrec.error = 0
        return _to_cartesian(satrec, mrt, mvt, rvdot, basis)
    except SGP4Error as e:
        satrec.error = e.code
        return False, False


def sgp4_init(satrec, satn,   epoch,
                bstar, dn, ddn, ecc, argp,
                inclo, mo,
                raan
                ):
    """
    (Math Heavy) initialisation of sgp4 which instantiates all important
    variables in the satrec object
    """

    temp4    =   1.5e-12

    # Near Earth Variables
    satrec.isimp   = 0;   satrec.aycof  = 0.0
    satrec.con41   = 0.0; satrec.cc1    = 0.0; satrec.cc4      = 0.0
    satrec.cc5     = 0.0; satrec.d2     = 0.0; satrec.d3       = 0.0
    satrec.d4      = 0.0; satrec.delmo  = 0.0; satrec.eta      = 0.0
    satrec.argpdot = 0.0; satrec.omgcof = 0.0; satrec.sinmao   = 0.0
    satrec.t       = 0.0; satrec.t2cof  = 0.0; satrec.t3cof    = 0.0
    satrec.t4cof   = 0.0; satrec.t5cof  = 0.0; satrec.x1mth2   = 0.0
    satrec.x7thm1  = 0.0; satrec.mdot   = 0.0; satrec.nodedot  = 0.0
    satrec.xlcof   = 0.0; satrec.xmcof  = 0.0; satrec.nodecf   = 0.0

    # Earth Constants
    satrec.mu     = 398600.5           #  in km3 / s2
    satrec.radiusearthkm = 6378.137      #  km
    satrec.xke    = 0.07436685317
    satrec.tumin  = 13.44685108
    satrec.j2     =  0.00108262998905
    satrec.j3     = -0.00000253215306
    satrec.j4     = -0.00000161098761
    satrec.j3oj2  = -0.002338890559

    ss = 1.012229276
    qzms2ttemp = 0.00658499496
    qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp
    x2o3   =  2.0 / 3.0

    # -- initialisation markers
    satrec.t	 = 0.0

    # --
    satrec.norad = satn
    satrec.classification = 'U'

    # --
    satrec.bstar   = bstar
    satrec.dn    = dn
    satrec.ddn   = ddn
    satrec.ecc    = ecc
    satrec.argp   = argp
    satrec.inclo   = inclo
    satrec.mo	     = mo
    satrec.raan   = raan

    # single averaged mean elements
    satrec.am = 0.0
    satrec.em = 0.0
    satrec.im = 0.0
    satrec.Om = 0.0
    satrec.mm = 0.0
    satrec.nm = 0.0

    satrec.error = 0

    # --
    (
    satrec.n,
    ao,    satrec.con41,  con42, cosio,
    cosio2, omeosq, posq,
    rp,    rteosq,sinio , satrec.gsto,
    ) = _initl(
        satrec.xke, satrec.j2, satrec.ecc, epoch, satrec.inclo, satrec.n
        )
    satrec.a    = pow( satrec.n*satrec.tumin , (-2.0/3.0) )
    satrec.alta = satrec.a*(1.0 + satrec.ecc) - 1.0
    satrec.altp = satrec.a*(1.0 - satrec.ecc) - 1.0

    sfour  = ss
    qzms24 = qzms2t

    if omeosq >= 0.0 or satrec.n >= 0.0:
        satrec.isimp = 0
        if rp < 220.0 / satrec.radiusearthkm + 1.0:
            satrec.isimp = 1
        perige = (rp - 1.0) * satrec.radiusearthkm;

        #  - for perigees below 156 km, s and qoms2t are altered -
        if perige < 156.0:

            sfour = perige - 78.0;
            if perige < 98.0:
                sfour = 20.0;
            #  sgp4fix use multiply for speed instead of pow
            qzms24temp =  (120.0 - sfour) / satrec.radiusearthkm;
            qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
            sfour  = sfour / satrec.radiusearthkm + 1.0;
            print(perige)

    pinvsq = 1.0 / posq

    tsi  = 1.0 / (ao - sfour)
    satrec.eta  = ao * satrec.ecc * tsi
    etasq = satrec.eta * satrec.eta
    eeta  = satrec.ecc * satrec.eta
    psisq = np.fabs(1.0 - etasq)
    coef  = qzms24 * pow(tsi, 4.0);
    coef1 = coef / pow(psisq, 3.5)
    cc2   = coef1 * satrec.n * (ao * (1.0 + 1.5 * etasq + eeta *
                (4.0 + etasq)) + 0.375 * satrec.j2 * tsi / psisq * satrec.con41 *
                (8.0 + 3.0 * etasq * (8.0 + etasq)))
    satrec.cc1   = satrec.bstar * cc2
    cc3   = 0.0
    if satrec.ecc > 1.0e-4:
        cc3 = -2.0 * coef * tsi * satrec.j3oj2 * satrec.n * sinio / satrec.ecc
    satrec.x1mth2 = 1.0 - cosio2
    satrec.cc4    = 2.0* satrec.n * coef1 * ao * omeosq * \
                    (satrec.eta * (2.0 + 0.5 * etasq) + satrec.ecc *
                    (0.5 + 2.0 * etasq) - satrec.j2 * tsi / (ao * psisq) *
                    (-3.0 * satrec.con41 * (1.0 - 2.0 * eeta + etasq *
                    (1.5 - 0.5 * eeta)) + 0.75 * satrec.x1mth2 *
                    (2.0 * etasq - eeta * (1.0 + etasq)) * np.cos(2.0 * satrec.argp)))
    satrec.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 *
                (etasq + eeta) + eeta * etasq)
    cosio4 = cosio2 * cosio2
    temp1  = 1.5 * satrec.j2 * pinvsq * satrec.n
    temp2  = 0.5 * temp1 * satrec.j2 * pinvsq
    temp3  = -0.46875 * satrec.j4 * pinvsq * pinvsq * satrec.n
    satrec.mdot     = satrec.n + 0.5 * temp1 * rteosq * satrec.con41 + 0.0625 * \
                    temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4)
    satrec.argpdot  = (-0.5 * temp1 * con42 + 0.0625 * temp2 *
                        (7.0 - 114.0 * cosio2 + 395.0 * cosio4) +
                        temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4))
    xhdot1            = -temp1 * cosio
    satrec.nodedot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * cosio2) +
                        2.0 * temp3 * (3.0 - 7.0 * cosio2)) * cosio
    satrec.omgcof   = satrec.bstar * cc3 * np.cos(satrec.argp)
    satrec.xmcof    = 0.0
    if satrec.ecc > 1.0e-4:
        satrec.xmcof = -x2o3 * coef * satrec.bstar / eeta
    satrec.nodecf = 3.5 * omeosq * xhdot1 * satrec.cc1
    satrec.t2cof   = 1.5 * satrec.cc1

    if np.fabs(cosio+1.0) > 1.5e-12:
        satrec.xlcof = -0.25 * satrec.j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio)
    else:
        satrec.xlcof = -0.25 * satrec.j3oj2 * sinio * (3.0 + 5.0 * cosio) / temp4
    satrec.aycof   = -0.5 * satrec.j3oj2 * sinio

    delmotemp = 1.0 + satrec.eta * np.cos(satrec.mo)
    satrec.delmo   = delmotemp * delmotemp * delmotemp
    satrec.sinmao  = np.sin(satrec.mo)
    satrec.x7thm1  = 7.0 * cosio2 - 1.0

    if satrec.isimp != 1:
        cc1sq          = satrec.cc1 * satrec.cc1
        satrec.d2    = 4.0 * ao * tsi * cc1sq
        temp           = satrec.d2 * tsi * satrec.cc1 / 3.0
        satrec.d3    = (17.0 * ao + sfour) * temp
        satrec.d4    = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * \
                        satrec.cc1
        satrec.t3cof = satrec.d2 + 2.0 * cc1sq
        satrec.t4cof = 0.25 * (3.0 * satrec.d3 + satrec.cc1 *
                        (12.0 * satrec.d2 + 10.0 * cc1sq))
        satrec.t5cof = 0.2 * (3.0 * satrec.d4 +
                        12.0 * satrec.cc1 * satrec.d3 +
                        6.0 * satrec.d2 * satrec.d2 +
                        15.0 * cc1sq * (2.0 * satrec.d2 + cc1sq))

    # propagate to 0
    sgp4_update(satrec, 0)

    return True
