"""
Sun position algorithm adapted from https://aa.usno.navy.mil/faq/sun_approx

References:
[1] U.S. Naval Observatory, "Computing Approximate Solar Coordinates",
    https://aa.usno.navy.mil/faq/sun_approx
[2] The Astronomical Almanac, "Low precision formulas for the Sun", p. C5
    (the source of [1]); gives apparent coordinates of the Sun referred to the
    equinox and ecliptic of date, as summarised in
    https://en.wikipedia.org/wiki/Position_of_the_Sun
[3] D. A. Vallado, "Fundamentals of Astrodynamics and Applications", 4th ed.,
    errata for pp. 277-281: the Almanac Sun algorithm is true-of-date (TOD),
    https://celestrak.org/software/vallado/ErrataVer4.pdf
[4] O. Montenbruck and E. Gill, "Satellite Orbits: Models, Methods and
    Applications", Springer, 2000. Sec. 3.3.1, Eq. (3.35): satellite to body
    vector s - r from geocentric positions. Sec. 3.3.2, Eq. (3.46): Sun position
    vector; precession moves the equinox 1.3972 deg per century relative to J2000
[5] IAU 2012 Resolution B2: 1 au = 149 597 870 700 m exactly,
    https://iauarchive.eso.org/static/resolutions/IAU2012_English.pdf
[6] D. A. Vallado, "Fundamentals of Astrodynamics and Applications", 4th ed.,
    Microcosm Press, 2013. Sec. 4.4.1, Eq. (4-1): position vector from right
    ascension and declination. Sec. 5.1.1, Algorithm 29: Sun position vector
[7] D. A. Vallado, P. Crawford, R. Hujsak, and T. S. Kelso, "Revisiting Spacetrack
    Report #3", AIAA 2006-6753: SGP4 output is in the True Equator Mean Equinox
    (TEME) frame, https://celestrak.org/publications/AIAA/2006-6753/
[8] python-sgp4 documentation: Satrec.sgp4 returns position in km from the
    center of the Earth in TEME; jday returns (jd, fr), https://pypi.org/project/sgp4/
"""

import numpy as np

J2000 = 2451545.0 #For Julian Date January 1st 2000 12:00 UTC 
AU_KM = 149597870.7 #Kilometers in 1 au, exact by definition [5]

def vec_sun_from_earth(jd, fr):
    """
    Calculating the Earth to Sun distance and direction from Julian Date using [1]

    Args:
        jd: Julian Date (UT), e.g. jd as returned by sgp4.api.jday
        fr: Fraction of a day added to jd, e.g. fr as returned by sgp4.api.jday

    Returns:
        r_sun: Earth to Sun distance (AU)
        alpha: Right ascension of the Sun (rad), in (-pi, pi]
        delta: Declination of the Sun (rad)

    alpha and delta are geocentric (Earth to Sun) apparent equatorial coordinates
    referred to the equator and equinox of date [1][2], not J2000. Mean vs true of
    date [2][3] differ by nutation (< 20 arcsec), below the ~1 arcmin accuracy of [1]
    within two centuries of 2000. The equinox of date is ~0.37 deg from J2000 in
    2026 [4, Sec. 3.3.2].
    """
    
    #Days since epoch J2000.0 [1]: D = JD - 2451545.0
    d = (jd + fr) - J2000
    print("days since j2000: ", d) 

    #Mean anomaly of the Sun in degrees [1]: g = 357.529 + 0.98560028 D
    g_deg = 357.529 + 0.98560028 * d 
    g = np.deg2rad(g_deg)
    print("g new: ", g)

    #Mean longitude of the Sun in degrees [1]: q = 280.459 + 0.98564736 D
    q_deg = 280.459 + 0.98564736 * d
    q = np.deg2rad(q_deg)
    print("q new: ", q)
    
    #Geocentric apparent ecliptic longitude of the Sun (adjusted for aberration) in degrees [1]:
    #L = q + 1.915 sin(g) + 0.020 sin(2g)
    #Reducing g, q, and L to 0-360 deg is optional [1] since sin/cos are periodic
    L_deg = q_deg + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g)
    L = np.deg2rad(L_deg)
    print("L new: ", L) 

    #Ecliptic latitude of the Sun, approximated as 0 [1]
    b = 0

    #Earth to Sun distance in AU [1]: R = 1.00014 - 0.01671 cos(g) - 0.00014 cos(2g)
    r_sun = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)
    print("r sun new: ", r_sun) 
    
    #Mean obliquity of the ecliptic in degrees [1]: e = 23.439 - 0.00000036 D
    epsilon_deg = 23.439 - 0.00000036 * d
    epsilon = np.deg2rad(epsilon_deg) 
    print("ep new: ", epsilon) 

    #Right ascension of the Sun [1]: tan(RA) = cos(e) sin(L) / cos(L)
    #atan2 keeps RA in the same quadrant as L [1]
    aarg = np.cos(epsilon) * np.sin(L)
    alpha = np.arctan2(aarg, np.cos(L))
    print("alpha new: ", alpha, aarg) 
          
    #Declination of the Sun [1]: sin(d) = sin(e) sin(L)
    darg = np.sin(epsilon) * np.sin(L)
    delta = np.arcsin(darg)
    print("delta new: ", delta, darg)
                
    return r_sun, alpha, delta


def vec_satellite_to_sun(jd, fr, r_sat):
    """
    Calculating vector from satellite to Sun (in km) from Julian Date and the
    satellite position vector returned by SGP4 (Earth-centered, in km)

    Args:
        jd: Julian Date (UT), e.g. jd as returned by sgp4.api.jday [8]
        fr: Fraction of a day added to jd, e.g. fr as returned by sgp4.api.jday [8]
        r_sat: Earth to satellite position (km) in TEME, as returned by Satrec.sgp4 [7][8]

    Returns:
        v_sat_sun: Satellite to Sun vector (km) in the equator and equinox of date
                   frame of vec_sun_from_earth. TEME differs from that frame only by
                   nutation terms (< 20 arcsec) [3][7], below the ~1 arcmin accuracy
                   of [1], so r_sat is used without a frame rotation.
    """

    r_sun, alpha, delta = vec_sun_from_earth(jd, fr)

    #Earth to Sun vector (km) from distance, RA, and Dec [6, Eq. (4-1)]:
    #r * [cos(delta)cos(alpha), cos(delta)sin(alpha), sin(delta)]
    #x points to the equinox of date, z to the north celestial pole [2]
    #Equivalent to r * [cos(L), cos(e)sin(L), sin(e)sin(L)] [4, Eq. (3.46)][6, Alg. 29]
    r_earth_sun = r_sun * AU_KM * np.array([
        np.cos(delta) * np.cos(alpha),
        np.cos(delta) * np.sin(alpha),
        np.sin(delta)
    ])

    #Satellite to Sun = Earth to Sun - Earth to satellite [4, Eq. (3.35)]
    v_sat_sun = r_earth_sun - np.array(r_sat)

    return v_sat_sun

if __name__ == '__main__':
        jd = 2461221
        fr = 10446
        print("Returned: " + str(vec_sun_from_earth(jd, fr)))
