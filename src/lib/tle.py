"""
Functions and Variables to update and use TLE (two-line element data)

Adapted from the TLE-tools library by @FedericoStra on GitHub
"""

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
        self.n = 0.0 # mean motion
        self.rev_num = 0

        self.tle_str = ""

    @classmethod
    def from_tle_lines(cls, name, line1, line2):
        """Parse a TLE from its constituent lines.

        All the attributes parsed from the TLE are expressed in the same units that
        are used in the TLE format.
        """

        cls()
        cls.name=name
        cls.norad=line1[2:7]
        cls.classification=line1[7] or 'U'
        cls.int_desig=line1[9:17]
        cls.epoch_year=_conv_year(line1[18:20])
        cls.epoch_day=float(line1[20:32])
        cls.dn=float(line1[33:43])
        cls.bstar=_parse_float(line1[53:61])
        cls.ddn=_parse_float(line1[44:52])
        cls.ephtype = line1[62]
        cls.set_num=line1[64:68]
        cls.inc=float(line2[8:16])
        cls.raan=float(line2[17:25])
        cls.ecc=_parse_decimal(line2[26:33])
        cls.argp=float(line2[34:42])
        cls.mo=float(line2[43:51])
        cls.n=float(line2[52:63])
        cls.rev_num=line2[63:68]
        cls.tle_str=name+line1+line2

    @classmethod
    def from_tle_file(cls, filename):
        """Load TLE from a file."""
        if isinstance(filename, str):
            with open(filename, encoding="utf-8") as fp:
                return cls.from_tle_lines(*fp.readlines[:2])

        return None

    @classmethod
    def from_tle_str(cls, string):
        """Load TLE from a string."""
        return cls.from_tle_lines(*string.split('\n')[:3])

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
