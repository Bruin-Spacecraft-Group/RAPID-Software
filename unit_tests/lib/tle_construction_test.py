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

    def eccentricity_decimal_check(self):
        

    
        
    def test_parsing(self):
        #Test TLE data parsing, like inclination, eccentricity, mean motion etc.
        
        line1 = "1 25544U 98067A   20029.54791435  .00001264  00000-0  29621-4 0  9993"
        line2 = "2 25544  51.6435  21.0930 0007417  99.9999 260.0000 15.49112339211667"

        satellite = tle.TLE(line1, line2)

        self.assertAlmostEqual(satellite.inclination, 51.6435)
        self.assertAlmostEqual(satellite.eccentricity, 0.0007417)
        self.assertAlmostEqual(satellite.mean_motion, 15.49112339)
        self.assertAlmostEqual(satellite.bstar, 0.000029621)
        self.assertAlmostEqual(satellite.mean_anomaly, 99.9999)
        self.assertAlmostEqual(satellite.argument_of_perigee, 260.0000)
        self.assertAlmostEqual(satellite.right_ascension, 21.0930)
        self.assertEqual(satellite.satellite_number, 25544)
        self.assertEqual(satellite.classification, 'U')
        self.assertEqual(satellite.international_designator, '98067A')
        self.assertEqual(satellite.epoch, 20029.54791435)
        self.assertEqual(satellite.revolution_number, 21166)




