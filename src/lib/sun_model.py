"""
Sun position algorithm adapted from https://aa.usno.navy.mil/faq/sun_approx
"""

import numpy as np

J2000 = 2451545.0 #For Julian Date January 1st 2000 12:00 UTC 
AU_KM = 149597870.7 #Kilometers in 1 AU

def dist_sun_from_earth(jd, fr):
    """
    Calculating Earth to Sun vector from Julian Date NOTE: Is it sun to earth or earth to sun...
    """
    
    #Days Since Epoch J2000: D = JD - 2451545.0
    d = (jd + fr) - J2000
    print("days since j2000: ", d) 

    #Mean Anomaly of Sun Equation in degrees: g = 357.529 + 0.98560028D 
    g_deg = 357.529 + 0.98560028 * d 
    g = np.deg2rad(g_deg)
    print("g new: ", g)

    #Mean Longitude of Sun in degrees: q = 280.459 + 0.98564736D
    q_deg = 280.459 + 0.98564736 * d
    q = np.deg2rad(q_deg)
    print("q new: ", q)
    
    #Mean Eliptical Longitude of Sun in degrees: L = q + 1.915 + sin(g) + 0.020 * sin(2g)
    #?NOTE: Add correction for angles not in between 0 to 360 / 0 to 2pi 
    L_deg = q_deg + 1.915 * np.sin(g) + 0.020 * np.sin(2 * g)
    L = np.deg2rad(L_deg)
    print("L new: ", L) 

    #Approx Eliptical Latitude 
    b = 0

    #Distance of Sun from Earth using degrees (in AU): r_sun = 1.0014 - 0.01671cos(g) - 0.00014cos(2g)
    r_sun = 1.0014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)
    print("r sun new: ", r_sun) 
    
    #Mean obliquity of ecliptic in degrees: epsilon = 23.439 - 0.00000036D
    epsilon_deg = 23.439 - 0.00000036 * d
    epsilon = np.deg2rad(epsilon_deg) 
    print("ep new: ", epsilon) 

    #Right ascension of Sun: tan(alpha) = cos(epsilon) * sin(L) / cos(L) 
    aarg = np.cos(epsilon) * np.sin(L)
    alpha = np.arctan2(aarg, np.cos(L))
    print("alpha new: ", alpha, aarg) 
          
    #Declination of sun: sin(delta) = sin(epsilon) * sin(L) 
    darg = np.sin(epsilon) * np.sin(L)
    delta = np.arcsin(darg)
    print("delta new: ", delta, darg)
                
    return r_sun, alpha, delta


def vec_satellite_to_sun(jd, fr, r_sat):
    """
    Calculating vector from satellite to Sun (in km) from Julian Date and the
    satellite position vector returned by SGP4 (Earth-centered, in km)
    """

    r_sun, alpha, delta = dist_sun_from_earth(jd, fr)

    #Earth to Sun vector from RA/Dec: r * [cos(delta)cos(alpha), cos(delta)sin(alpha), sin(delta)]
    r_earth_sun = r_sun * AU_KM * np.array([
        np.cos(delta) * np.cos(alpha),
        np.cos(delta) * np.sin(alpha),
        np.sin(delta)
    ])

    #Satellite to Sun = Earth to Sun - satellite to earth
    v_sat_sun = r_earth_sun + np.array(r_sat)

    return v_sat_sun

if __name__ == '__main__':
        jd = 2461221
        fr = 10446
        print("Returned: " + str(dist_sun_from_earth(jd, fr)))
