#TODO: Epoch is not beingg parsed/set correctly, need to implement in tle.py then test in here
#
import unittest
import os
import tle

# TLE data for the International Space Station (ISS) used in the tests
TLE_DATA = """ISS (ZARYA)
1 25544U 98067A   20029.54791435  .00001264  00000-0  29621-4 0  9998
2 25544  51.6435  21.0930 0007417  99.9999 260.0000 15.49112339211666"""
TLE_NAME, TLE_LINE1, TLE_LINE2 = TLE_DATA.splitlines()


class TLEConstruction_Test(unittest.TestCase):
    """Checks if the TLE string is valid."""

    def size_check(self, tle_text): 
        """This method checks if both lines of the TLE are 69 characters long."""
        lines = [line.strip() for line in tle_text.strip().splitlines() if line.strip()]
        data_lines = [l for l in lines if l.startswith(('1 ', '2 '))]

        for line in data_lines:

            if len(line) < 69:
                return False

        return True

    def is_check_sum_valid(self, tle_text):
        """This method goes through both lines of the TLE, adds their numbers up and checks if the remainder
        of the sum of all the numbers is equal to the check sum. If it is, it returns True. If it isn't, it 
        return"""
        lines = [line.strip() for line in tle_text.strip().splitlines() if line.strip()]
        data_lines = [l for l in lines if l.startswith(('1 ', '2 '))]

        if not data_lines:
            return False

        for line in data_lines:

            if len(line) < 69:
                return False

            data_part = line[:68]
            expected_checksum = int(line[68])

            total_sum = 0
            for char in data_part:
                if char.isdigit():
                    total_sum += int(char)
                elif char == '-':
                    total_sum += 1
            if total_sum % 10 != expected_checksum:
                return False
        return True
   
    def test_tle_format(self):
        self.assertTrue(self.size_check(TLE_DATA))
        self.assertTrue(self.is_check_sum_valid(TLE_DATA))
        
    def test_parse_float(self):

        # Negative values with negative exponent
        self.assertAlmostEqual(
            tle._parse_float("-12345-3"), -0.00012345
        )

        # Positive exponent
        self.assertAlmostEqual(
            tle._parse_float(" 12345+2"), 12.345
        )
            
    def assert_parsing(self, sat):
        """Helper Function to test different TLE parsing methods based on the ISS TLE"""

        self.assertEqual(sat.name, 'ISS (ZARYA)')
        self.assertEqual(sat.norad, '25544') # NORAD catalog number
        self.assertEqual(sat.classification, 'U') # Classification
        self.assertEqual(sat.int_desig, '98067A') # International Designator

        self.assertEqual(sat.epoch, (2020, 1, 29, 13, 8, 59, 799840)) # Epoch (year, month, day, hour, minute, second, microsecond)
        self.assertAlmostEqual(sat.dn, 0.00001264) # first time deriv. mean motion
        self.assertAlmostEqual(sat.ddn, 0.0) # Second time deriv. mean motion
        self.assertAlmostEqual(sat.bstar, 0.000029621) #BSTAR coeff
        self.assertAlmostEqual(sat.ephtype, 0) # Ephemeris type
        self.assertAlmostEqual(sat.set_num, 999) # Set num(increment on update, usually fixed 999)

        
        self.assertAlmostEqual(sat.inc, 51.6435) # Inclination
        self.assertAlmostEqual(sat.raan, 21.0930) # Right ascension of ascending node
        self.assertAlmostEqual(sat.ecc, 0.0007417) # eccentricity 
        self.assertAlmostEqual(sat.argp, 99.9999) # Argument of perigee
        self.assertAlmostEqual(sat.mo, 260.0000) # Mean anomaly
        self.assertAlmostEqual(sat.n, 15.49112339) # Mean Motion
        self.assertEqual(sat.rev_num, 21166) # Revolution

    def test_from_tle_lines(self):
        """Test parsing TLE from individual lines."""
        sat = tle.Satrec.from_tle_lines(TLE_NAME, TLE_LINE1, TLE_LINE2)
        self.assert_parsing(sat)

    def test_from_tle_str(self):
        """Test parsing TLE from a string."""
        sat = tle.Satrec.from_tle_str(TLE_DATA)
        self.assert_parsing(sat)
    
    def test_from_tle_file(self):
        """Test parsing TLE from a file."""
        tle_path = "iss_test.tle"
        file_created = False
        try:
            with open(tle_path, "w", encoding="utf-8") as tle_file:
                file_created = True
                tle_file.write(TLE_DATA)

            sat = tle.Satrec.from_tle_file(tle_path)
            self.assert_parsing(sat)
        finally:
            if file_created:
                os.remove(tle_path)

    def test_to_array(self):
        """Test TLE fields are returned in the documented 2D array layout."""
        sat = tle.Satrec.from_tle_lines(TLE_NAME, TLE_LINE1, TLE_LINE2)

        self.assertEqual(sat.to_array(), [
            [TLE_NAME],
            ["25544", "U", "98067A", 2020, 29.54791435, 0.00001264,
             0.0, 0.000029621, 999],
            [51.6435, 21.0930, 0.0007417, 99.9999, 260.0, 15.49112339, 21166],
        ])

    #I am ending the liveshare because I am leaving
    def test_epoch_parsing(self):
        """Test epoch dates, leap years, and invalid days."""
        return
        
    def make_epoch_line(epoch):
        line1 = TLE_LINE1[:18] + epoch + TLE_LINE1[32:]
        return line1[:68] + str(tle._compute_checksum(line1[:68]))
        
if __name__ == "__main__":
    unittest.main()