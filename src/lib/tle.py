"""
Functions and Variables to update and use TLE (two-line element data)

Adapted from the TLE-tools library by @FedericoStra on GitHub
"""

from sgp4 import sgp4_init, sgp4_update
from datastores.adcs import Time

try:
    import ulab.numpy as np  # For CircuitPython
except ImportError:
    import numpy as np  # For GitHub Actions / PC testing

MIN_PER_DAY = 1440
EPOCH0 = 2433281.5 # jan 0 1950

def _day_of_year_to_month_day(day_of_year, is_leap):
    """Core logic for turning days into months, for easy testing."""
    february_bump = (2 - is_leap) * (day_of_year >= 60 + is_leap)
    august = day_of_year >= 215
    month, day = divmod(2 * (day_of_year - 1 + 30 * august + february_bump), 61)
    month += 1 - august
    day //= 2
    day += 1
    return month, day

def _days2mdhms(year, days, round_to_microsecond=6):
    """Convert a float point number of days into the year into date and time.

    >>> days2mdhms(2000, 32.0)  # February 1
    (2, 1, 0, 0, 0.0)
    >>> days2mdhms(2000, 366.0)  # December 31, since 2000 was a leap year
    (12, 31, 0, 0, 0.0)

    The floating point seconds are rounded to an even number of
    microseconds if ``round_to_microsecond`` is true.
    """
    day_of_year, day_fraction = divmod(days, 1.0)

    second = day_fraction * 86400.0
    if round_to_microsecond:
        second = round(second, round_to_microsecond)

    minute, second = divmod(second, 60.0)
    if round_to_microsecond:
        second = round(second, round_to_microsecond)

    minute = int(minute)
    hour, minute = divmod(minute, 60)
    hour = int(hour)

    is_leap = year % 400 == 0 or (year % 4 == 0 and year % 100 != 0)
    month, day = _day_of_year_to_month_day(int(day_of_year), is_leap)
    if month == 13:  # behave like the original in case of overflow
        month = 12
        day += 31

    return month, day, hour, minute, second

def _conv_year(s):
    """Interpret a two-digit year string."""
    if isinstance(s, int):
        return s
    y = int(s)
    return y + (1900 if y >= 57 else 2000)

def _parse_decimal(s):
    """Parse a floating point with implicit leading dot.

    >>> _parse_decimal('378')
    0.378
    """
    return float('.' + s)

def _parse_float(s):
    """Parse a floating point with implicit dot and exponential notation.

    >>> _parse_float(' 12345-3')
    0.00012345
    >>> _parse_float('+12345-3')
    0.00012345
    >>> _parse_float('-12345-3')
    -0.00012345
    """
    return float(s[0] + '.' + s[1:6] + 'e' + s[6:8])

def _sgp4_jday(time):
    """
    Converts jdsatepoch into a compatible number for sgp4
    """

    year, mon, day, hr, minute, sec = time

    return (367.0 * year -
            7.0 * (year + ((mon + 9.0) // 12.0)) * 0.25 // 1.0 +
            275.0 * mon // 9.0 +
            day + 1721013.5 +
            ((sec / 60.0 + minute) / 60.0 + hr) / 24.0
            )

def jday(time):
    """
    From a date, return a Julian date in its date + fractional form

    Used to build the jd, fr parameters in the function call of sgp4_update
    """

    year, mon, day, hr, minute, sec = time

    jd = (367.0 * year
         - 7 * (year + ((mon + 9) // 12.0)) * 0.25 // 1.0
	   + 275 * mon / 9.0 // 1.0
	   + day
         + 1721013.5)
    fr = (sec + minute * 60.0 + hr * 3600.0) / 86400.0
    return jd, fr

class Satrec:
    """
    Satellite record object

    In this implementation, built from TLE data

    Two line-elements (TLEs) are unpacked from both given and propagated data.
    This implementation uses Keplerian orbital parameters

    All the attributes parsed from the TLE are expressed in the same units that
    are used in the TLE format.

    :str name:
        Name of the satellite.
    :str norad:
        NORAD catalog number (https://en.wikipedia.org/wiki/Satellite_Catalog_Number).
    :str classification:
        'U', 'C', 'S' for unclassified, classified, secret.
    :str int_desig:
        International designator (https://en.wikipedia.org/wiki/International_Designator),
    :int epoch_year:
        Year of the epoch.
    :float epoch_day:
        Day of the year plus fraction of the day.
    :float dn:
        First time derivative of the mean motion (divided by 2 in TLE)
    :float ddn:
        Second time derivative of the mean motion (divided by 6 in TLE).
    :float bstar:
        BSTAR coefficient (https://en.wikipedia.org/wiki/BSTAR).
    :int set_num:
        Element set number.
    :float inc:
        Inclination.
    :float raan:
        Right ascension of the ascending node.
    :float ecc:
        Eccentricity.
    :float argp:
        Argument of perigee.
    :float mo:
        Mean anomaly.
    :float n:
        Mean motion.
    :int rev_num:
        Revolution number.
    """

    def __init__(self):
        self.name = ""

        self.norad = ""
        self.classification = ""
        self.int_desig = ""

        self.epoch_year = 0
        self.epoch_day = 0.0
        self.dn = 0.0
        self.ddn = 0.0
        self.bstar = 0.0
        self.ephtype = ""
        self.set_num = 0

        self.inc = 0.0
        self.raan = 0.0
        self.ecc = 0.0
        self.argp = 0.0
        self.mo = 0.0
        self.n = 0.0
        self.rev_num = 0

        self.tle_str = ""
        self.error = ""

        self.jdsatepoch = 0.0
        self.jdsatepoch_f = 0.0

    @classmethod
    def from_tle_lines(cls, name, line1, line2):
        """Parse a TLE from its constituent lines.

        All the attributes parsed from the TLE are expressed in the same units that
        are used in the TLE format.
        """

        sat = cls()
        sat.name = name
        sat.norad = line1[2:7]
        sat.classification = line1[7] or 'U'
        sat.int_desig = line1[9:17].strip()
        sat.epoch_year = _conv_year(line1[18:20])
        sat.epoch_day = float(line1[20:32])
        sat.dn = float(line1[33:43])
        sat.ddn = _parse_float(line1[44:52])
        sat.bstar = _parse_float(line1[53:61])
        sat.ephtype = int(line1[62])
        sat.set_num = int(line1[64:68])
        sat.inc = float(line2[8:16])
        sat.raan = float(line2[17:25])
        sat.ecc = _parse_decimal(line2[26:33])
        sat.argp = float(line2[34:42])
        sat.mo = float(line2[43:51])
        sat.n = float(line2[52:63])
        sat.rev_num = int(line2[63:68])
        sat.tle_str = name + line1 + line2
        return sat

    @classmethod
    def from_tle_file(cls, filename):
        """Load TLE from a file."""
        with open(filename, encoding="utf-8") as fp:
            return cls.from_tle_str(fp.read())

    @classmethod
    def from_tle_str(cls, string):
        """Load TLE from a string."""
        lines = [line.strip() for line in string.splitlines() if line.strip()]
        if len(lines) == 3:
            name, line1, line2 = lines
        elif len(lines) == 2:
            name = ""
            line1, line2 = lines
        else:
            raise ValueError("TLE data must contain two element lines and an optional name")

        return cls.from_tle_lines(name, line1, line2)

    @classmethod
    def sgp4_init(cls, tle):
        """
        Creates a satrec object specifically modified from TLE to be used in sgp4.

        Changes units, activates certain new parameters, etc. 
        """

        self = tle

        # constants for unit change
        deg2rad  =  np.pi / 180.0         #    0.0174532925199433
        xpdotp   =  1440.0 / (2.0 *np.pi)  #  229.1831180523293

        #  ---- convert to sgp4 units ----
        self.n = self.n / xpdotp
        self.dn = self.dn  / (xpdotp*1440.0)
        self.ddn= self.ddn / (xpdotp*1440.0*1440)

        #  ---- find standard orbital elements ----
        self.inc = self.inc  * deg2rad
        self.raan = self.raan  * deg2rad
        self.argp = self.argp  * deg2rad
        self.mo    = self.mo     * deg2rad

        year = self.epoch_year

        mon, day, hr, minute, sec = _days2mdhms(year, self.epoch_day)
        time = Time.from_ymdhms(year, mon, day, hr, minute, sec)
        jd_full = _sgp4_jday(time)

        # Split into two-part representation
        self.jdsatepoch   = np.floor(jd_full-0.5) + 0.5   # noon-to-noon JD boundary
        self.jdsatepoch_f = jd_full - self.jdsatepoch

        sgp4_init(self, self.set_num, jd_full-EPOCH0, self.bstar,
                  self.dn, self.ddn, self.ecc, self.argp, self.inc, self.mo,
                  self.raan)

        return self

    def to_array(self):
        """
        Return 2D array of TLE values

        Indexed as 
        [line, col]

        n is mean motion, d suggests time derivative

        name: [0,0] 

        norad: [1,0] classification: [1,1] int_desig: [1,2] epoch_year: [1,3] day: [1,4] 
        dn: [1,5] ddn: [1,6] bstar: [1,7] set_num: [1,8]

        inclination: [2,0] RAAN: [2,1] eccentricity: [2,2] arg_perigee: [2,3] Mean Anomaly: [2,4] 
        n: [2,5] rev_num: [2,6]
        """

        return [
            [self.name], # Line 0
            [self.norad, self.classification, self.int_desig, # line 1 ID
              # line 1 time-derivative
             self.epoch_year, self.epoch_day, self.dn, self.ddn, self.bstar, self.set_num],
              # line 2 orbital params
            [self.inc, self.raan, self.ecc, self.argp, self.mo, self.n, self.rev_num]
        ]

    def sgp4_update(self, jd, fr):
        """
        For a julian date (jd) and its fractional representation (fr), 
        propagate Satrec using sgp4
        """

        tsince = ((jd - self.jdsatepoch) * MIN_PER_DAY +
                  (fr - self.jdsatepoch_f) * MIN_PER_DAY)

        r, v = sgp4_update(self, tsince)

        return self.error, r, v

    def error_message(self):
        """
        Return error message from current self.error (when polled)
        """
        if self.error == 2:
            return (f'mean motion {0:f} is less than zero').format(self.n)

        return ""

def test():
    """
    Test script to test the sgp4 script using the ISS TLE
    """
    iss_tle = """ISS (ZARYA)
    \n1 25544U 98067A   26121.81277072  .00006771  00000-0  13067-3 0  9997 
    \n2 25544  51.6311 169.6452 0007227  12.7206 347.3964 15.49051775564564"""

    sat: Satrec = Satrec.from_tle_str(
        iss_tle
    )

    sgp4_obj : Satrec= Satrec.sgp4_init(sat)
    time = Time.from_ymdhms(2026, 5, 3, 0, 0, 0)
    error, r, v = sgp4_obj.sgp4_update(*jday(time))
    print(error, r, v)

if __name__ == "__main__":
    test()
