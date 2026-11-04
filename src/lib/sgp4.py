"""
Simplified General Perturbations Model 4 Implementation for Orbit Propagation

As outlined in https://celestrak.org/publications/AIAA/2008-6770/AIAA-2008-6770.pdf
"""

# Error codes
ECCENTRICITY = 1 # eccentricity is not within 0-1
MOTION = 2 # error in propagating mean motion
SEMIRECT = 4 # apoapsis, periapsis characteristics error
DECAY = 6 # orbit has decayed

try:
    import ulab.numpy as np  # For CircuitPython
except ImportError:
    import numpy as np  # For GitHub Actions / PC testing

def _gstime(jdut1):
    deg2rad = np.pi / 180.0
    tau = np.pi*2
    tut1 = (jdut1 - 2451545.0) / 36525.0
    temp = -6.2e-6* tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + \
            (876600.0*3600 + 8640184.812866) * tut1 + 67310.54841  #  sec
    temp = (temp * deg2rad / 240.0) % tau # 360/86400 = 1/240, to deg, to rad

    #  ------------------------ check quadrants ---------------------
    if temp < 0.0:
        temp += tau

    return temp

def _initl(xke, j2,
            ecc, epoch, inclo, no):
    #  ----------------------- earth constants ----------------------
    #  sgp4fix identify constants and allow alternate values
    #  only xke and j2 are used here so pass them in directly
    #  tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2 = whichconst
    x2o3   = 2.0 / 3.0;
    tau = np.pi*2

    #  ------------- calculate auxillary epoch quantities ----------
    eccsq  = ecc * ecc;
    omeosq = 1.0 - eccsq;
    rteosq = np.sqrt(omeosq);
    cosio  = np.cos(inclo);
    cosio2 = cosio * cosio;

    #  ------------------ un-kozai the mean motion -----------------
    ak    = pow(xke / no, x2o3);
    d1    = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
    del_  = d1 / (ak * ak);
    adel  = ak * (1.0 - del_ * del_ - del_ *
            (1.0 / 3.0 + 134.0 * del_ * del_ / 81.0));
    del_  = d1/(adel * adel)
    no    = no / (1.0 + del_)

    ao    = pow(xke / no, x2o3)
    sinio = np.sin(inclo)
    po    = ao * omeosq
    con42 = 1.0 - 5.0 * cosio2
    con41 = -con42-cosio2-cosio2
    posq  = po * po
    rp    = ao * (1.0 - ecc)

    ts70  = epoch - 7305.0;
    ds70 = (ts70 + 1.0e-8) // 1.0;
    tfrac = ts70 - ds70;
    #  find greenwich location at epoch
    c1    = 1.72027916940703639e-2;
    thgr70= 1.7321343856509374;
    fk5r  = 5.07551419432269442e-15;
    c1p2p = c1 + tau;
    gsto  = (thgr70 + c1*ds70 + c1p2p*tfrac + ts70*ts70*fk5r) % tau
    if gsto < 0.0:
        gsto = gsto + tau
    # gsto = _gstime(epoch + 2433281.5)

    return (no,
       ao,    con41,  con42, cosio,
       cosio2, omeosq, posq,
       rp,    rteosq,sinio , gsto
    )

def sgp4_update(satrec, tsince):
    """
    Transforms to cartesian (r, v) from a Satrec (TLE 7+) input, then performs matrix operations
    to propagate the system 

    Current implementation uses an analytical approach to compute the Jacobian A matrix
    Possible alternative is numerical stepping (for lower accuracy) of associated partials
    """
    # -- Set mathematical constants
    x2o3  = 2.0 / 3.0
    tau = 2.0 * np.pi
    vkmpersec = satrec.radiusearthkm * satrec.xke/60.0;

    satrec.t = tsince

    # -- Update for secular gravity and atmospheric drag
    xmdf    = satrec.mo + satrec.mdot * satrec.t
    argpdf  = satrec.argp + satrec.argpdot * satrec.t
    nodedf  = satrec.raan + satrec.nodedot * satrec.t
    argpm   = argpdf
    mm      = xmdf
    t2      = satrec.t * satrec.t
    nodem   = nodedf + satrec.nodecf * t2
    tempa   = 1.0 - satrec.cc1 * satrec.t
    tempe   = satrec.bstar * satrec.cc4 * satrec.t
    templ   = satrec.t2cof * t2

    if satrec.isimp != 1:

        delomg = satrec.omgcof * satrec.t
        #  sgp4fix use mutliply for speed instead of pow
        delmtemp =  1.0 + satrec.eta * np.cos(xmdf)
        delm   = satrec.xmcof * \
                (delmtemp * delmtemp * delmtemp -
                satrec.delmo)
        temp   = delomg + delm
        mm     = xmdf + temp
        argpm  = argpdf - temp
        t3     = t2 * satrec.t
        t4     = t3 * satrec.t
        tempa  = tempa - satrec.d2 * t2 - satrec.d3 * t3 - \
                        satrec.d4 * t4
        tempe  = tempe + satrec.bstar * satrec.cc5 * (np.sin(mm) -
                        satrec.sinmao)
        templ  = templ + satrec.t3cof * t3 + t4 * (satrec.t4cof +
                        satrec.t * satrec.t5cof)

    nm    = satrec.n
    em    = satrec.ecc
    inclm = satrec.inclo

    if nm <= 0.0:
        satrec.error = satrec.MOTION
        return False, False

    am = pow((satrec.xke / nm),x2o3) * tempa * tempa
    nm = satrec.xke / pow(am, 1.5)
    em = em - tempe

    if em >= 1.0 or em < 0.0:
        satrec.error = satrec.ECCENTRICITY

        return False, False

    if em < 1.0e-6:
        em  = 1.0e-6
    mm     = mm + satrec.n * templ
    xlm    = mm + argpm + nodem
    emsq   = em * em
    temp   = 1.0 - emsq

    nodem  = nodem % tau
    argpm  = argpm % tau
    xlm    = xlm % tau
    mm     = (xlm - argpm - nodem) % tau

    satrec.am = am
    satrec.em = em
    satrec.im = inclm
    satrec.Om = nodem
    satrec.om = argpm
    satrec.mm = mm
    satrec.nm = nm

    #  ----------------- compute extra mean quantities -------------
    sinim = np.sin(inclm)
    cosim = np.cos(inclm)

    #  -------------------- add lunar-solar periodics --------------
    ep     = em
    xincp  = inclm
    argpp  = argpm
    nodep  = nodem
    mp     = mm
    sinip  = sinim
    cosip  = cosim

    #  -------------------- long period periodics ------------------
    axnl = ep * np.cos(argpp)
    temp = 1.0 / (am * (1.0 - ep * ep))
    aynl = ep* np.sin(argpp) + temp * satrec.aycof
    xl   = mp + argpp + nodep + temp * satrec.xlcof * axnl

    #  --------------------- solve kepler's equation ---------------
    u    = (xl - nodep) % tau
    eo1  = u
    tem5 = 9999.9
    ktr = 1

    while np.fabs(tem5) >= 1.0e-12 and ktr <= 10:

        sineo1 = np.sin(eo1)
        coseo1 = np.cos(eo1)
        tem5   = 1.0 - coseo1 * axnl - sineo1 * aynl
        tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5
        if np.fabs(tem5) >= 0.95:
            tem5 = 0.95 if tem5 > 0.0 else -0.95
        eo1    = eo1 + tem5
        ktr = ktr + 1

    #  ------------- short period preliminary quantities -----------
    ecose = axnl*coseo1 + aynl*sineo1
    esine = axnl*sineo1 - aynl*coseo1
    el2   = axnl*axnl + aynl*aynl
    pl    = am*(1.0-el2)
    if pl < 0.0:
        satrec.error = satrec.SEMIRECT

        return False, False

    else:

        rl     = am * (1.0 - ecose)
        rdotl  = np.sqrt(am) * esine/rl
        rvdotl = np.sqrt(pl) / rl
        betal  = np.sqrt(1.0 - el2)
        temp   = esine / (1.0 + betal)
        sinu   = am / rl * (sineo1 - aynl - axnl * temp)
        cosu   = am / rl * (coseo1 - axnl + aynl * temp)
        su     = np.atan2(sinu, cosu)
        sin2u  = (cosu + cosu) * sinu
        cos2u  = 1.0 - 2.0 * sinu * sinu
        temp   = 1.0 / pl
        temp1  = 0.5 * satrec.j2 * temp
        temp2  = temp1 * temp

        #  -------------- update for short period periodics ------------
        mrt   = rl * (1.0 - 1.5 * temp2 * betal * satrec.con41) + \
                0.5 * temp1 * satrec.x1mth2 * cos2u
        su    = su - 0.25 * temp2 * satrec.x7thm1 * sin2u
        xnode = nodep + 1.5 * temp2 * cosip * sin2u
        xinc  = xincp + 1.5 * temp2 * cosip * sinip * cos2u
        mvt   = rdotl - nm * temp1 * satrec.x1mth2 * sin2u / satrec.xke
        rvdot = rvdotl + nm * temp1 * (satrec.x1mth2 * cos2u +
                1.5 * satrec.con41) / satrec.xke

        #  --------------------- orientation vectors -------------------
        sinsu =  np.sin(su)
        cossu =  np.cos(su)
        snod  =  np.sin(xnode)
        cnod  =  np.cos(xnode)
        sini  =  np.sin(xinc)
        cosi  =  np.cos(xinc)
        xmx   = -snod * cosi
        xmy   =  cnod * cosi
        ux    =  xmx * sinsu + cnod * cossu
        uy    =  xmy * sinsu + snod * cossu
        uz    =  sini * sinsu
        vx    =  xmx * cossu - cnod * sinsu
        vy    =  xmy * cossu - snod * sinsu
        vz    =  sini * cossu

        #  --------- position and velocity (in km and km/sec) ----------
        _mr = mrt * satrec.radiusearthkm
        r = (_mr * ux, _mr * uy, _mr * uz)
        v = ((mvt * ux + rvdot * vx) * vkmpersec,
            (mvt * uy + rvdot * vy) * vkmpersec,
            (mvt * uz + rvdot * vz) * vkmpersec)

    if mrt < 1.0:
            satrec.error = satrec.DECAY
            return False, False

    return r, v

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

    if omeosq >= 0.0 or satrec.n >= 0.0:
        satrec.isimp = 0
        if rp < 220.0 / satrec.radiusearthkm + 1.0:
            satrec.isimp = 1
        sfour  = ss
        qzms24 = qzms2t
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
