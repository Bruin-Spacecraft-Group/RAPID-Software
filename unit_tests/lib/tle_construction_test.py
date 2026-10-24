# one line of tle may not get through
# check the exponents (self.ddn and self.bstar)
# the check sum isnt where its supposed to be ( might be shifted so check column 69)
# decimal parsing might be off (put decimal in the right place)
# create a tle string for each case (take current ISS one)

import unittest
import os
import tempfile
import tle

class TLEConstruction_Test(unittest.TestCase):
    """Checks if the TLE string is valid."""

    def size_check(tle_text): 
        """This method checks if both lines of the TLE are 69 characters long."""
        lines = [line.strip() for line in tle_text.strip().splitlines() if line.strip()]
        data_lines = [l for l in lines if l.startswith(('1 ', '2 '))]

        for line in data_lines:

            if len(line) < 69:
                return False

        return True

    def is_check_sum_valid(tle_text):
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
        

    def test_check_exponents_from_tle(self):
        """Check DDN and BSTAR in a complete TLE."""
        
        line1 = "1 25544U 98067A   20029.54791435  .00001264  00000-0  29621-4 0  9993"
        line2 = "2 25544  51.6435  21.0930 0007417  99.9999 260.0000 15.49112339211667"

        # Case 1: Original ISS TLE
        sat = tle.Satrec.from_tle_lines("ISS (ZARYA)", line1, line2)


        self.assertAlmostEqual(sat.ddn, 0.0)
        self.assertAlmostEqual(sat.bstar, 0.000029621)

        # Case 2: Negative exponent
        line1_case2 = (
            line1[:44] + " 12345-3"
            + line1[52:53] + " 12345-3"
            + line1[61:]
        )

        sat = tle.Satrec.from_tle_lines("ISS (ZARYA)", line1_case2, line2)

        self.assertAlmostEqual(sat.ddn, 0.00012345)
        self.assertAlmostEqual(sat.bstar, 0.00012345)

        # Case 3: Negative value
        self.assertAlmostEqual(
            tle._parse_float("-12345-3"), -0.00012345
        )

        # Case 4: Positive exponent
        self.assertAlmostEqual(
            tle._parse_float(" 12345+2"), 12.345
        )
            
    def assert_parsing(self, sat):
        """Helper Function to test different TLE parsing methods based on the ISS TLE"""

        self.assertEqual(sat.name, 'ISS (ZARYA)')
        self.assertEqual(sat.norad, '25544') # NORAD catalog number
        self.assertEqual(sat.classification, 'U') # Classification
        self.assertEqual(sat.int_desig, '98067A') # International Designator

        self.assertEqual(sat.epoch_day, 29.54791435) # Days after year as float
        self.assertEqual(sat.epoch_year, 2020) # Year
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
        line1 = "1 25544U 98067A   20029.54791435  .00001264  00000-0  29621-4 0  9993"
        line2 = "2 25544  51.6435  21.0930 0007417  99.9999 260.0000 15.49112339211667"
        sat = tle.Satrec.from_tle_lines("ISS (ZARYA)", line1, line2)
        self.assert_parsing(sat)

    def test_from_tle_str(self):
        """Test parsing TLE from a string."""
        tle_string = """ISS (ZARYA)
1 25544U 98067A   20029.54791435  .00001264  00000-0  29621-4 0  9993
2 25544  51.6435  21.0930 0007417  99.9999 260.0000 15.49112339211667"""
        sat = tle.Satrec.from_tle_str(tle_string)
        self.assert_parsing(sat)
    
    def test_from_tle_file(self):
        """Test parsing TLE from a file."""
        tle_data = """ISS (ZARYA)
1 25544U 98067A   20029.54791435  .00001264  00000-0  29621-4 0  9993
2 25544  51.6435  21.0930 0007417  99.9999 260.0000 15.49112339211667"""
        with tempfile.TemporaryDirectory() as temp_dir:
            tle_path = os.path.join(temp_dir, "iss.tle")
            with open(tle_path, "w", encoding="utf-8") as tle_file:
                tle_file.write(tle_data)

            sat = tle.Satrec.from_tle_file(tle_path)
            self.assert_parsing(sat)

if __name__ == "__main__":
    unittest.main()