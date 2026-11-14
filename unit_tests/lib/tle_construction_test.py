#TODO: Verify that the tests are correct
import unittest
import os
import math
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
        """Test the TLE format validation methods."""
        self.assertTrue(self.size_check(TLE_DATA))
        short_line = TLE_DATA.replace(TLE_LINE1, TLE_LINE1[:68])
        self.assertFalse(self.size_check(short_line))

        self.assertTrue(self.is_check_sum_valid(TLE_DATA))
        bad_checksum = TLE_DATA[:-1] + (
            "0" if TLE_DATA[-1] != "0" else "1"
        )
        self.assertFalse(self.is_check_sum_valid(bad_checksum))
            
    def assert_parsing(self, sat):
        """Helper Function to test different TLE parsing methods based on the ISS TLE"""

        self.assertEqual(sat.name, 'ISS (ZARYA)') #Name of the satellite
        self.assertNotEqual(sat.name, 'OTHER')
        self.assertEqual(sat.norad, '25544') # NORAD catalog number
        self.assertNotEqual(sat.norad, '25545')
        self.assertEqual(sat.classification, 'U') # Classification
        self.assertNotEqual(sat.classification, 'C')
        self.assertEqual(sat.int_desig, '98067A') # International Designator
        self.assertNotEqual(sat.int_desig, '98067B')

        self.assertEqual(sat.epoch, (2020, 1, 29, 13, 8, 59, 799840)) # Epoch (year, month, day, hour, minute, second, microsecond)
        self.assertNotEqual(sat.epoch, (2020, 1, 29, 13, 8, 59, 799841))
        self.assertAlmostEqual(sat.dn, 0.00001264) # first time deriv. mean motion
        self.assertNotAlmostEqual(sat.dn, 1.00001264)
        self.assertAlmostEqual(sat.ddn, 0.0) # Second time deriv. mean motion
        self.assertNotAlmostEqual(sat.ddn, 1.0)
        self.assertAlmostEqual(sat.bstar, 0.000029621) #BSTAR coeff
        self.assertNotAlmostEqual(sat.bstar, 1.000029621)
        self.assertAlmostEqual(sat.ephtype, 0) # Ephemeris type
        self.assertNotAlmostEqual(sat.ephtype, 1)
        self.assertAlmostEqual(sat.set_num, 999) # Set num(increment on update, usually fixed 999)
        self.assertNotAlmostEqual(sat.set_num, 1000)
        
        self.assertAlmostEqual(sat.inc, 51.6435) # Inclination
        self.assertNotAlmostEqual(sat.inc, 52.6435)
        self.assertAlmostEqual(sat.raan, 21.0930) # Right ascension of ascending node
        self.assertNotAlmostEqual(sat.raan, 22.0930)
        self.assertAlmostEqual(sat.ecc, 0.0007417) # eccentricity 
        self.assertNotAlmostEqual(sat.ecc, 1.0007417)
        self.assertAlmostEqual(sat.argp, 99.9999) # Argument of perigee
        self.assertNotAlmostEqual(sat.argp, 100.9999)
        self.assertAlmostEqual(sat.mo, 260.0000) # Mean anomaly
        self.assertNotAlmostEqual(sat.mo, 261.0)
        self.assertAlmostEqual(sat.n, 15.49112339) # Mean Motion
        self.assertNotAlmostEqual(sat.n, 16.49112339)
        self.assertEqual(sat.rev_num, 21166) # Revolution
        self.assertNotEqual(sat.rev_num, 21167)

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

        expected = [
            [TLE_NAME],
            ["25544", "U", "98067A", 2020, 29.54791435, 0.00001264,
             0.0, 0.000029621, 999],
            [51.6435, 21.0930, 0.0007417, 99.9999, 260.0, 15.49112339, 21166],
        ]
        self.assertEqual(sat.to_array(), expected)
        wrong_expected = [row[:] for row in expected]
        wrong_expected[1][0] = "25545"
        self.assertNotEqual(sat.to_array(), wrong_expected)

    def test_epoch_parsing(self):
        """Test the TLE epoch fields are parsed correctly."""

        sat = tle.Satrec.from_tle_lines(TLE_NAME, TLE_LINE1, TLE_LINE2)

        self.assertEqual(sat.epoch_year, 2020)
        self.assertNotEqual(sat.epoch_year, 2021)
        self.assertEqual(sat.epoch_day_int, 29)
        self.assertNotEqual(sat.epoch_day_int, 30)
        self.assertEqual(sat.epoch_frac, 54791435)
        self.assertNotEqual(sat.epoch_frac, 54791436)
        self.assertAlmostEqual(sat.epoch_day, 29.54791435)
        self.assertNotAlmostEqual(sat.epoch_day, 30.54791435)

    def test_parse_helpers(self):
        """Test the parsing helpers like _parse_float(), _num(), _parse_int(), etc."""

        self.assertEqual(tle._parse_int("0042"), 42) # Leading 0 test
        self.assertNotEqual(tle._parse_int("0042"), 43)
        self.assertEqual(tle._parse_int("  "), 0) # Empty string test
        self.assertNotEqual(tle._parse_int("  "), 1)
        self.assertEqual(tle._parse_int(" 12 "), 12) # Leading/trailing whitespace test
        with self.assertRaises(ValueError): # Invalid integer string test
            tle._parse_int("12x")
        for invalid_integer in ("-1", "+1", "1.0"):
            with self.subTest(invalid_integer=invalid_integer):
                with self.assertRaises(ValueError):
                    tle._parse_int(invalid_integer)

        self.assertAlmostEqual(tle._parse_decimal("378"), 0.378) # Decimal parsing test
        self.assertNotAlmostEqual(tle._parse_decimal("378"), 0.379)
        self.assertEqual(tle._parse_decimal("000"), 0.0) # Leading zeros test
        self.assertAlmostEqual(tle._parse_decimal("0007417"), 0.0007417)
        for invalid_decimal in ("", "abc", "-378"):
            with self.subTest(invalid_decimal=invalid_decimal):
                with self.assertRaises(ValueError):
                    tle._parse_decimal(invalid_decimal)

        self.assertAlmostEqual(tle._num(" -12.5 "), -12.5) # Leading/trailing whitespace and negative number test
        self.assertNotAlmostEqual(tle._num(" -12.5 "), 12.5)
        self.assertAlmostEqual(tle._num("+1.25"), 1.25) # Leading sign
        self.assertAlmostEqual(tle._num(".5"), 0.5)
        for invalid_number in (
            "", " ", ".", "+-", "nan", "inf", "1e5", "1_0", "9" * 400,
        ):
            with self.subTest(invalid_number=invalid_number[:10]):
                with self.assertRaises(ValueError):
                    tle._num(invalid_number)

        # Negative values with negative exponent
        self.assertAlmostEqual(
            tle._parse_float("-12345-3"), -0.00012345
        )
        self.assertNotAlmostEqual(
            tle._parse_float("-12345-3"), 1.0
        )

        # Positive exponent
        self.assertAlmostEqual(
            tle._parse_float(" 12345+2"), 12.345
        )
        self.assertNotAlmostEqual(
            tle._parse_float(" 12345+2"), 13.345
        )
        self.assertAlmostEqual(tle._parse_float("+12345-3"), 0.00012345)
        for invalid_tle_float in ("123", " 12345x3"):
            with self.subTest(invalid_tle_float=invalid_tle_float):
                with self.assertRaises(ValueError):
                    tle._parse_float(invalid_tle_float)

    def test_calendar_helper_functions(self):
        """Test calendar and Julian-date helper functions."""
        # Test _conv_year()
        for year_text, expected in (("00", 2000), ("56", 2056), ("57", 1957),
                                    ("99", 1999)):
            with self.subTest(year_text=year_text):
                self.assertEqual(tle._conv_year(year_text), expected)
                self.assertNotEqual(tle._conv_year(year_text), expected + 1)
        self.assertEqual(tle._conv_year(1999), 1999)
        with self.assertRaises(ValueError):
            tle._conv_year("xx")
            
        sat = tle.Satrec.from_tle_lines(TLE_NAME, TLE_LINE1, TLE_LINE2)
        self.assertEqual(sat.epoch_unix, (1580303339, 799840))
        self.assertNotEqual(sat.epoch_unix, (1580303339, 799841))

        # Test _days_from_civil()
        for year, month, day, expected in (
            (1970, 1, 1, 0),
            (1969, 12, 31, -1),
            (1970, 1, 2, 1),
            (2000, 1, 1, 10957),
            (2000, 2, 29, 11016),
            (2020, 2, 29, 18321),
        ):
            with self.subTest(civil_date=(year, month, day)):
                actual = tle._days_from_civil(year, month, day)
                self.assertEqual(actual, expected)
                self.assertNotEqual(actual, expected + 1)

        # Test _days2mdhms()
        for year, day_of_year, expected in (
            (2000, 1.0, (1, 1, 0, 0, 0.0)),
            (2000, 32.0, (2, 1, 0, 0, 0.0)),
            (2000, 60.0, (2, 29, 0, 0, 0.0)),
            (1900, 60.0, (3, 1, 0, 0, 0.0)),
            (2023, 365.0, (12, 31, 0, 0, 0.0)),
            (2020, 366.0, (12, 31, 0, 0, 0.0)),
            (2020, 29.54791435, (1, 29, 13, 8, 59.79984)),
            (2000, 1.0 + 59.9999996 / 86400, (1, 1, 0, 1, 0.0)),
        ):
            with self.subTest(year=year, day_of_year=day_of_year):
                actual = tle._days2mdhms(year, day_of_year)
                self.assertEqual(actual, expected)
                self.assertNotEqual(actual, expected[:-1] + (expected[-1] + 1,))

        unrounded = tle._days2mdhms(2000, 1.0 + 1.23456789 / 86400, 0)
        self.assertAlmostEqual(unrounded[-1], 1.23456789, places=8)
        self.assertNotAlmostEqual(unrounded[-1], 1.33456789, places=8)

        # Test _is_leap()
        for year, expected in (
            (2000, True), (1900, False), (2024, True),
            (2023, False), (2400, True), (2100, False),
        ):
            with self.subTest(leap_year=year):
                self.assertEqual(tle._is_leap(year), expected)
                self.assertNotEqual(tle._is_leap(year), not expected)

        # Test epoch_unix property of Satrec
        epoch = tle.Satrec()
        epoch.epoch_year = 1970
        epoch.epoch_day_int = 1
        epoch.epoch_frac = 0
        self.assertEqual(epoch.epoch_unix, (0, 0))
        self.assertNotEqual(epoch.epoch_unix, (1, 0))
        epoch.epoch_year = 1969
        epoch.epoch_day_int = 365
        self.assertEqual(epoch.epoch_unix, (-86400, 0))
        self.assertNotEqual(epoch.epoch_unix, (-86400, 1))
        for year, day, fraction in ((1900, 366, 0), (2020, 1, 100000000)):
            epoch.epoch_year = year
            epoch.epoch_day_int = day
            epoch.epoch_frac = fraction
            with self.subTest(epoch=(year, day, fraction)):
                with self.assertRaises(ValueError):
                    epoch.epoch_unix

        # Test jday() and _sgp4_jday() functions
        for ymdhms, expected_jd, expected_fraction in (
            ((2000, 1, 1, 0, 0, 0), 2451544.5, 0.0),
            ((2000, 1, 1, 12, 0, 0), 2451544.5, 0.5),
            ((2000, 1, 1, 12, 0, 0.5), 2451544.5, 43200.5 / 86400),
            ((2000, 1, 2, 0, 0, 0), 2451545.5, 0.0),
            ((2000, 2, 29, 23, 59, 59), 2451603.5, 86399 / 86400),
        ):
            with self.subTest(julian_input=ymdhms):
                time = tle.Time.from_ymdhms(*ymdhms)
                jd, fraction = tle.jday(time)
                self.assertAlmostEqual(jd, expected_jd)
                self.assertNotAlmostEqual(jd, expected_jd + 1)
                self.assertAlmostEqual(fraction, expected_fraction)
                self.assertNotAlmostEqual(fraction, expected_fraction + 0.1)
                self.assertAlmostEqual(
                    tle._sgp4_jday(time), expected_jd + expected_fraction
                )
                self.assertNotAlmostEqual(
                    tle._sgp4_jday(time), expected_jd + expected_fraction + 1
                )

    def test_validation_helpers(self):
        """Tests tle_checksum(), and _fix_name() functions."""

        # Test _fix_name() with various inputs
        self.assertEqual(tle._fix_name(0, [TLE_LINE1, TLE_LINE2]), "")
        self.assertNotEqual(tle._fix_name(0, [TLE_LINE1, TLE_LINE2]), TLE_NAME)
        self.assertEqual(
            tle._fix_name(1, [TLE_NAME, TLE_LINE1, TLE_LINE2]), TLE_NAME
        )
        self.assertEqual(
            tle._fix_name(1, ["0 " + TLE_NAME, TLE_LINE1, TLE_LINE2]), TLE_NAME
        )
        for start, lines in (
            (-1, [TLE_NAME, TLE_LINE1, TLE_LINE2]),
            (2, [TLE_NAME, TLE_LINE1, TLE_LINE2]),
            (0, [TLE_LINE1]),
            (1, [TLE_NAME, TLE_LINE1]),
        ):
            with self.subTest(name_start=start, line_count=len(lines)):
                with self.assertRaises(ValueError):
                    tle._fix_name(start, lines)

        # Test tle_checksum() with various (mostly nonsense) inputs to check its calculations
        self.assertEqual(tle.tle_checksum(TLE_LINE1[:68]), 8)
        self.assertEqual(tle.tle_checksum(TLE_LINE2[:68]), 6)
        self.assertEqual(tle.tle_checksum("123-+."), 7)
        self.assertNotEqual(tle.tle_checksum("123-+."), 8)
        self.assertEqual(tle.tle_checksum("9" * 68), 2)

    def test_validate(self):
        """Test valid TLE boundaries and reject malformed fields."""

        # Helper function to create modified TLE lines with specific field changes
        def changed_pair(fields):
            lines = [TLE_LINE1, TLE_LINE2]
            for line_number, start, end, value in fields:
                self.assertEqual(len(value), end - start)
                index = line_number - 1
                body = lines[index][:start] + value + lines[index][end:68]
                lines[index] = body + str(tle.tle_checksum(body))
            return tuple(lines)

        valid_boundary_pairs = (
            changed_pair((
                (1, 2, 7, "A1234"),
                (1, 18, 20, "20"),
                (1, 20, 32, "366.99999999"),
                (1, 64, 68, "    "),
                (2, 2, 7, "A1234"),
                (2, 8, 16, "180.0000"),
                (2, 17, 25, "359.9999"),
                (2, 26, 33, "9999999"),
                (2, 34, 42, "359.9999"),
                (2, 43, 51, "359.9999"),
                (2, 52, 63, "99.99999999"),
                (2, 63, 68, "     "),
            )),
            changed_pair((
                (1, 2, 7, "A1234"),
                (1, 7, 8, " "),
                (1, 9, 17, "        "),
                (1, 18, 20, "20"),
                (1, 20, 32, "001.00000000"),
                (1, 64, 68, "    "),
                (2, 2, 7, "A1234"),
                (2, 8, 16, "  0.0000"),
                (2, 17, 25, "  0.0000"),
                (2, 26, 33, "0000000"),
                (2, 34, 42, "  0.0000"),
                (2, 43, 51, "  0.0000"),
                (2, 52, 63, " 0.00000001"),
                (2, 63, 68, "     "),
            )),
        )
        for pair in valid_boundary_pairs:
            with self.subTest(valid_boundary_epoch=pair[0][18:32]):
                self.assertIsNone(tle._validate(*pair))

        invalid_pairs = [
            ("short line", (TLE_LINE1[:-1], TLE_LINE2)),
            ("long line", (TLE_LINE1 + " ", TLE_LINE2)),
            ("non-string line", (None, TLE_LINE2)),
            ("bad line 1 number", changed_pair(((1, 0, 1, "2"),))),
            ("bad line 2 number", changed_pair(((2, 0, 1, "1"),))),
            ("non-digit checksum", (TLE_LINE1[:68] + "x", TLE_LINE2)),
            (
                "checksum mismatch",
                (
                    TLE_LINE1[:68] + str((int(TLE_LINE1[68]) + 1) % 10),
                    TLE_LINE2,
                ),
            ),
            ("catalog mismatch", changed_pair(((2, 2, 7, "25545"),))),
            (
                "invalid Alpha-5 catalog",
                changed_pair(((1, 2, 7, "I1234"), (2, 2, 7, "I1234"))),
            ),
            ("classification", changed_pair(((1, 7, 8, "X"),))),
            ("international designator", changed_pair(((1, 9, 17, "98?67A  "),))),
            ("epoch year", changed_pair(((1, 18, 20, "+1"),))),
            ("epoch day syntax", changed_pair(((1, 20, 32, "029x54791435"),))),
            ("epoch day zero", changed_pair(((1, 20, 32, "000.00000000"),))),
            ("epoch day 367", changed_pair(((1, 20, 32, "367.00000000"),))),
            (
                "day 366 in non-leap year",
                changed_pair(((1, 18, 20, "19"), (1, 20, 32, "366.00000000"))),
            ),
            ("epoch fraction syntax", changed_pair(((1, 20, 32, "029.1234567x"),))),
            ("mean motion derivative", changed_pair(((1, 33, 43, " .00001x64"),))),
            ("second derivative", changed_pair(((1, 44, 52, " 0000x-0"),))),
            ("BSTAR", changed_pair(((1, 53, 61, " 29621x4"),))),
            ("ephemeris type", changed_pair(((1, 62, 63, "X"),))),
            ("element set number", changed_pair(((1, 64, 68, "A999"),))),
            ("inclination syntax", changed_pair(((2, 8, 16, " 51x6435"),))),
            ("inclination range", changed_pair(((2, 8, 16, "181.0000"),))),
            ("RAAN range", changed_pair(((2, 17, 25, "360.0000"),))),
            ("eccentricity syntax", changed_pair(((2, 26, 33, "00074x7"),))),
            ("argument of perigee range", changed_pair(((2, 34, 42, "360.0000"),))),
            ("mean anomaly range", changed_pair(((2, 43, 51, "360.0000"),))),
            ("mean motion syntax", changed_pair(((2, 52, 63, "15.4911233x"),))),
            ("zero mean motion", changed_pair(((2, 52, 63, "00.00000000"),))),
            ("revolution number", changed_pair(((2, 63, 68, "21x66"),))),
        ]

        for label, pair in invalid_pairs:
            with self.subTest(invalid_tle=label):
                with self.assertRaises(ValueError):
                    tle._validate(*pair)

    def test_sgp4_methods(self):
        """Check SGP4 against published Vallado verification vectors."""
        # Vanguard 1, case 00005 from the Vallado verification data:
        # https://github.com/brandon-rhodes/python-sgp4/blob/master/sgp4/SGP4-VER.TLE
        # https://github.com/brandon-rhodes/python-sgp4/blob/master/sgp4/tcppver.out
        # https://celestrak.org/publications/AIAA/2008-6770/AIAA-2008-6770.pdf
        line1 = (
            "1 00005U 58002B   00179.78495062  .00000023  00000-0  "
            "28098-4 0  4753"
        )
        line2 = (
            "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 "
            "10.82419157413667"
        )
        sat = tle.Satrec.from_tle_lines("VANGUARD 1", line1, line2)
        original = (sat.n, sat.dn, sat.ddn, sat.inc, sat.raan, sat.argp, sat.mo)
        epoch_seconds, epoch_microseconds = sat.epoch_unix
        expected_epoch_jd = 2440587.5 + (
            epoch_seconds + epoch_microseconds / 1_000_000
        ) / 86400

        # Init tests
        self.assertIs(tle.Satrec.sgp4_init(sat), sat)
        self.assertEqual(sat.error, 0)
        radians_per_revolution = 2 * math.pi
        converted_n = original[0] * radians_per_revolution / 1440
        self.assertGreater(sat.n, 0)
        self.assertAlmostEqual(sat.n, converted_n, delta=converted_n * 0.001)
        converted_dn = original[1] * radians_per_revolution / 1440**2
        self.assertAlmostEqual(
            sat.dn, converted_dn, delta=max(abs(converted_dn) * 1e-10, 1e-15)
        )
        self.assertEqual(sat.ddn, original[2] * radians_per_revolution / 1440**3)
        initialized_angles = (sat.inc, sat.raan, sat.argp, sat.mo)
        for actual, degrees in zip(initialized_angles, original[3:]):
            self.assertAlmostEqual(actual, degrees * math.pi / 180, places=12)
        self.assertAlmostEqual(
            sat.jdsatepoch + sat.jdsatepoch_f, expected_epoch_jd, places=8
        )

        # Update tests
        error, position, velocity = sat.sgp4_update(
            sat.jdsatepoch, sat.jdsatepoch_f
        )
        self.assertEqual(error, 0)
        self.assertEqual(sat.error, 0)
        self.assertAlmostEqual(sat.t, 0.0)
        # tcppver.out uses different Earth constants; allow their small offset.
        reference_position = (7022.46529266, -1400.08296755, 0.03995155)
        reference_velocity = (1.893841015, 6.405893759, 4.534807250)
        for actual, expected in zip(position, reference_position):
            self.assertAlmostEqual(actual, expected, delta=0.02)
        for actual, expected in zip(velocity, reference_velocity):
            self.assertAlmostEqual(actual, expected, delta=0.00002)

        # WGS-84 reference values from python-sgp4's assert_wgs84() test:
        # https://github.com/brandon-rhodes/python-sgp4/blob/master/sgp4/tests.py
        minutes = 309.67110720001529
        day_offset, epoch_fraction = divmod(
            sat.jdsatepoch_f + minutes / 1440, 1
        )
        error, position, velocity = sat.sgp4_update(
            sat.jdsatepoch + day_offset, epoch_fraction
        )
        self.assertEqual(error, 0)
        self.assertEqual(sat.error, 0)
        self.assertAlmostEqual(sat.t, minutes, places=6)
        reference_position = (
            -3754.2437675772426,
            7876.3549956188945,
            4719.227897029576,
        )
        for actual, expected in zip(position, reference_position):
            self.assertAlmostEqual(actual, expected, delta=0.0002)

        # Error checks for invalid inputs to SGP4 update
        for field, invalid_value, expected_error in (
            ("ecc", 1.0, 1),
            ("n", -1.0, 2),
        ):
            with self.subTest(invalid_field=field):
                invalid_sat = tle.Satrec.from_tle_lines("VANGUARD 1", line1, line2)
                tle.Satrec.sgp4_init(invalid_sat)
                setattr(invalid_sat, field, invalid_value)
                error, position, velocity = invalid_sat.sgp4_update(
                    invalid_sat.jdsatepoch, invalid_sat.jdsatepoch_f
                )
                self.assertEqual(error, expected_error)
                self.assertEqual(invalid_sat.error, expected_error)
                self.assertFalse(position)
                self.assertFalse(velocity)

if __name__ == "__main__":
    unittest.main()
