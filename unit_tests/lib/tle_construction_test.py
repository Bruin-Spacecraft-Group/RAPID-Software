# one line of tle may not get through
# check the exponents (self.ddn and self.bstar)
# the check sum isnt where its suppoed to be ( migth be shifted so check column 69)
# decimal parsing might be off (put decimal in the right place)
# create a tle string for each case (take current ISS on)





import unittest
import tle

TLE_tests = [
    """
    """, # exponent 1
    """
    """,
]


class TLEConstruction_Test(unittest.TestCase):

    def size_check(tle_text): 
        lines = [l for l in tle_text.strip("\n").splitlines()]
        
        if len(lines) == 3:
            lines = lines[1:]
        if len(lines) != 2:
            return False
        l1, l2 = lines
        return (
            len(l1) == 69 and len(l2) == 69
            and l1.startswith("1 ") and l2.startswith("2 ")
        )
        
        
    def test_parsing(self):
        #Test TLE data parsing, like inclination, eccentricity, mean motion etc.
        
        line1 = "1 25544U 98067A   20029.54791435  .00001264  00000-0  29621-4 0  9993"
        line2 = "2 25544  51.6435  21.0930 0007417  99.9999 260.0000 15.49112339211667"

        satellite = tle.Satrec.from_tle_lines("ISS (ZARYA)", line1, line2)

        """Test TLE data parsing, like inclination, eccentricity, mean motion etc."""
        
        line1 = "1 25544U 98067A   20029.54791435  .00001264  00000-0  29621-4 0  9993"
        line2 = "2 25544  51.6435  21.0930 0007417  99.9999 260.0000 15.49112339211667"

        satellite = tle.Satrec()
        satellite.tle.Satrec.from_tle_lines("ISS (ZARYA)", line1, line2)

        self.assertEqual(satellite.norad, '98067A') # NORAD catalog number
        self.assertEqual(satellite.classification, 'U') # Classification
        self.assertEqual(satellite.int_desig, '98067A') # International Designator

        self.assertEqual(satellite.epoch_day, 029.54791435) # Days after year as flaot
        self.assertEqual(satellite.epoch_year, 2020) # Year
        self.assertAlmostEqual(satellite.dn, 15.49112339) # first time deriv. mean motion
        self.assertAlmostEqual(satellite.ddn,  00000-0) # Second time deriv. mean motion
        self.assertAlmostEqual(satellite.bstar, 0.000029621) #BSTAR coeff
        self.assertAlmostEqual(satellite.ephtype, 0) # Ephemeris type
        self.assertAlmostEqual(satellite.set_num, 999) # Set num(increment on update, usually fixed 999)

        
        self.assertAlmostEqual(satellite.inc, 51.6435) # Inclination
        self.assertAlmostEqual(satellite.raan, 21.0930) # Right ascension of ascending node
        self.assertAlmostEqual(satellite.ecc, 0.0007417) # eccentricity 
        self.assertAlmostEqual(satellite.argp, 99.9999) # Argument of perigee
        self.assertAlmostEqual(satellite.mo, 260.0000) # Mean anomaly
        self.assertEqual(satellite.n, 25544) # Mean Motion
        self.assertEqual(satellite.rev_num, 21166) # Revolution




