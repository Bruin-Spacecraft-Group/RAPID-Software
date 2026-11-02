import os
from decimal import Decimal
from unittest import TestCase
from unittest.mock import DEFAULT, mock_open, patch
import datetime

from geomag import (
    BlackoutZoneException,
    CautionZoneException,
    GeoMag,
    GeoMagResult,
    GeoMagUncertaintyResult,
)

from wmm.wmm_2025 import WMM_2025

TEST_STYLE_0 = 0
TEST_STYLE_1 = 1
TEST_STYLE_2 = 2


def get_test_filename(filename):
    return os.path.join(os.path.dirname(__file__), filename)


def get_os_based_test_path(path):
    return path.replace("/", os.sep)


class TestGeoMagResult(TestCase):
    def test_calculate_uncertainty(self):
        geo_mag = GeoMag(coefficients_file="wmm/WMM_2025.COF")
        result = geo_mag.calculate(80, 0, 0, 2025)
        uncertainty = result.calculate_uncertainty()
        self.assertIsInstance(uncertainty, GeoMagUncertaintyResult)
        self.assertAlmostEqual(uncertainty.d, 0.87, 2)

    def test_d_property(self):
        result = GeoMagResult(0, 0, 0, 0)
        result.d = 5
        self.assertEqual(result.d, result.dec)

    def test_f_properties(self):
        result = GeoMagResult(0, 0, 0, 0)
        result.f = 5
        self.assertEqual(result.f, result.ti)
        self.assertEqual(result.f, result.total_intensity)

    def test_i_property(self):
        result = GeoMagResult(0, 0, 0, 0)
        result.i = 5
        self.assertEqual(result.i, result.dip)
        self.assertEqual(result.i, result.inclination)


class TestGeoMagUncertaintyResult(TestCase):
    def test_static_values_wmm2025(self):
        geo_mag = GeoMag(coefficients_file="wmm/WMM_2025.COF")
        result = geo_mag.calculate(80, 0, 0, 2025)
        self.assertFalse(result.is_high_resolution)
        uncertainty = GeoMagUncertaintyResult(result)
        self.assertEqual(uncertainty.x, 137.0)
        self.assertEqual(uncertainty.y, 89.0)
        self.assertEqual(uncertainty.z, 141.0)
        self.assertEqual(uncertainty.h, 133.0)
        self.assertEqual(uncertainty.f, 138.0)
        self.assertEqual(uncertainty.i, 0.20)

    # def test_static_values_wmmhr2025(self):
    #     geo_mag = GeoMag(coefficients_file="wmm/WMMHR_2025.COF", high_resolution=True)
    #     result = geo_mag.calculate(80, 0, 0, 2025)
    #     self.assertTrue(result.is_high_resolution)
    #     uncertainty = GeoMagUncertaintyResult(result)
    #     self.assertEqual(uncertainty.x, 135.0)
    #     self.assertEqual(uncertainty.y, 85.0)
    #     self.assertEqual(uncertainty.z, 134.0)
    #     self.assertEqual(uncertainty.h, 130.0)
    #     self.assertEqual(uncertainty.f, 134.0)
    #     self.assertEqual(uncertainty.i, 0.19)

    def test_time_out_of_supported_range(self):
        geo_mag = GeoMag(coefficients_file="wmm/WMM_2025.COF")
        result = geo_mag.calculate(80, 0, 0, 2031.0, allow_date_outside_lifespan=True)
        with self.assertRaisesRegex(
            ValueError, "GeoMagResult outside of known uncertainty estimates."
        ):
            GeoMagUncertaintyResult(result)


class TestGeoMagCoefficients(TestCase):
    def get_test_values(self, test_parameter, style):
        if style == TEST_STYLE_0:
            print(len(test_parameter.split()))
            time, alt, glat, glon, d, i, h, x, y, z, f, gv, _, _, _, _, _, _, _ = (
                t(s)
                for t, s in zip(
                    (
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                    ),
                    test_parameter.split(),
                )
            )

        elif style == TEST_STYLE_1:
            time, alt, glat, glon, x, y, z, h, f, i, d, gv, _, _, _, _, _, _, _ = (
                t(s)
                for t, s in zip(
                    (
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                    ),
                    test_parameter.split(),
                )
            )

        else:
            time, alt, glat, glon, d, i, h, x, y, z, f, _, _, _, _, _, _, _ = (
                t(s)
                for t, s in zip(
                    (
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                        float,
                    ),
                    test_parameter.split(),
                )
            )
            gv = None

        return time, alt, glat, glon, x, y, z, h, f, i, d, gv

    def run_tests(self, geo_mag, test_filename, style):
        with open(get_test_filename(test_filename)) as test_values_file:
            for row, test_parameter in enumerate(test_values_file):
                if test_parameter[0] == "#":
                    continue

                time, alt, glat, glon, x, y, z, h, f, i, d, gv = self.get_test_values(
                    test_parameter, style
                )

                result = geo_mag.calculate(glat, glon, alt, time)
                gv_test = -999 if result.gv is None else result.gv

                self.assertAlmostEqual(
                    x, result.x, 1, f"Row {row}: X (nT) expected {x}, result {result.x}"
                )
                self.assertAlmostEqual(
                    y, result.y, 1, f"Row {row}: Y (nT) expected {y}, result {result.y}"
                )
                self.assertAlmostEqual(
                    z, result.z, 1, f"Row {row}: Z (nT) expected {z}, result {result.z}"
                )
                self.assertAlmostEqual(
                    h, result.h, 1, f"Row {row}: H (nT) expected {h}, result {result.h}"
                )
                self.assertAlmostEqual(
                    f, result.f, 1, f"Row {row}: F (nT) expected {f}, result {result.f}"
                )
                self.assertAlmostEqual(
                    i,
                    result.i,
                    2,
                    f"Row {row}: I (Deg) expected {i}, result {result.i}",
                )
                self.assertAlmostEqual(
                    d,
                    result.d,
                    2,
                    f"Row {row}: D (Deg) expected {d}, result {result.d}",
                )
                if style != TEST_STYLE_2:
                    self.assertAlmostEqual(
                        gv,
                        gv_test,
                        2,
                        f"Row {row}: GV (Deg) expected {gv}, result {result.gv}",
                    )

    def test_calculate_declination_from_2025_wmm_style_2_file(self):
        self.run_tests(
            GeoMag(coefficients_file="wmm/WMM_2025.COF"),
            "wmm/WMM2025_TEST_VALUES.txt",
            TEST_STYLE_2,
        )

    def test_calculate_declination_from_2025_wmm_style_2_data(self):
        self.run_tests(
            GeoMag(coefficients_data=WMM_2025),
            "wmm/WMM2025_TEST_VALUES.txt",
            TEST_STYLE_2,
        )

    def test_calculate_declination_from_2025_wmm_style_2_year(self):
        self.run_tests(
            GeoMag(base_year=2025),
            "wmm/WMM2025_TEST_VALUES.txt",
            TEST_STYLE_2,
        )

    # def test_calculate_declination_from_2025_wmm_hr_style_2_file(self):
    #     self.run_tests(
    #         GeoMag(coefficients_file="wmm/WMMHR_2025.COF", high_resolution=True),
    #         "test_values/WMMHR2025_TEST_VALUES.txt",
    #         TEST_STYLE_2,
    #     )

    # def test_calculate_declination_from_2025_wmm_hr_style_2_data(self):
    #     self.run_tests(
    #         GeoMag(coefficients_data=WMMHR_2025, high_resolution=True),
    #         "test_values/WMMHR2025_TEST_VALUES.txt",
    #         TEST_STYLE_2,
    #     )

    # def test_calculate_declination_from_2025_wmm_hr_style_2_year(self):
    #     self.run_tests(
    #         GeoMag(base_year=2025, high_resolution=True),
    #         "test_values/WMMHR2025_TEST_VALUES.txt",
    #         TEST_STYLE_2,
    #     )


class TestGeoMag(TestCase):
    def test_calculate(self):
        geo_mag = GeoMag(coefficients_file="wmm/WMM_2025.COF")
        result = geo_mag.calculate(89, -121, 28, 2025)
        self.assertAlmostEqual(result.d, -99.77, 2)
        self.assertAlmostEqual(result.i, 88.47, 2)
        self.assertAlmostEqual(result.h, 1504.298146, 4)
        self.assertAlmostEqual(result.x, -255.388723, 4)
        self.assertAlmostEqual(result.y, -1482.460628, 4)
        self.assertAlmostEqual(result.z, 56194.288771, 4)
        self.assertAlmostEqual(result.f, 56214.419888, 4)
        self.assertFalse(result.is_high_resolution)

    # def test_calculate_high_resolution(self):
    #     geo_mag = GeoMag(coefficients_file="wmm/WMMHR_2025.COF", high_resolution=True)
    #     result = geo_mag.calculate(89, -121, 28, 2025)
    #     self.assertAlmostEqual(result.d, -101.96, 2)
    #     self.assertAlmostEqual(result.i, 88.44, 2)
    #     self.assertAlmostEqual(result.h, 1526.369169, 4)
    #     self.assertAlmostEqual(result.x, -316.211505, 4)
    #     self.assertAlmostEqual(result.y, -1493.255880, 4)
    #     self.assertAlmostEqual(result.z, 56107.758937, 4)
    #     self.assertAlmostEqual(result.f, 56128.516957, 4)
    #     self.assertTrue(result.is_high_resolution)

    def test_create_list(self):
        self.assertEqual(GeoMag._create_list(2), [None, None])
        self.assertEqual(GeoMag._create_list(3, 0), [0, 0, 0])
        self.assertNotEqual(GeoMag._create_list(4, 0), [1, 1, 1, 1])

    def test_create_matrix(self):
        self.assertEqual(GeoMag._create_matrix(2, 2), [[None, None], [None, None]])
        self.assertEqual(GeoMag._create_matrix(2, 3, 0), [[0, 0, 0], [0, 0, 0]])
        self.assertNotEqual(GeoMag._create_matrix(2, 4), [[1, 1, 1, 1], [1, 1, 1, 1]])

    def test_exception_get_coefficients_year_pre_2025(self):
        with self.assertRaisesRegex(
            ValueError,
            "There are no coefficients for the year 2024",
        ):
            GeoMag._get_coefficients_year(2024)

    def test_exception_get_coefficients_year_post_2029(self):
        with self.assertRaisesRegex(
            ValueError,
            "There are no coefficients for the year 2030",
        ):
            GeoMag._get_coefficients_year(2030)

    def test_get_coefficients_year(self):
        self.assertEqual(GeoMag._get_coefficients_year(2025), "2025")
        self.assertEqual(GeoMag._get_coefficients_year(2029), "2025")
        self.assertEqual(GeoMag._get_coefficients_year(2025.0), "2025")
        self.assertEqual(GeoMag._get_coefficients_year(Decimal("2025.0")), "2025")
        self.assertEqual(
            GeoMag._get_coefficients_year(datetime.datetime(2025, 1, 1)), "2025"
        )
        self.assertEqual(
            GeoMag._get_coefficients_year(datetime.date(2025, 1, 1)), "2025"
        )

    def test_get_model_filename_default(self):
        geo_mag = GeoMag()
        model_filename = geo_mag._get_model_filename()
        expected = get_os_based_test_path("wmm/WMM.COF")
        self.assertTrue(
            model_filename.endswith(expected),
            f"{model_filename!r} does not end with {expected!r}",
        )

    # def test_get_model_filename_default_with_high_resolution(self):
    #     geo_mag = GeoMag(high_resolution=True)
    #     model_filename = geo_mag._get_model_filename()
    #     self.assertEqual(
    #         model_filename[-22:], get_os_based_test_path("pygeomag/wmm/WMMHR.COF")
    #     )

    def test_get_model_filename_default_not_in_wmm_path(self):
        m = mock_open(read_data="")
        m.side_effect = [OSError, DEFAULT]
        with patch("geomag.open", m):
            geo_mag = GeoMag()
            model_filename = geo_mag._get_model_filename()
            expected = get_os_based_test_path("wmm/WMM.COF")
            self.assertTrue(
                model_filename.endswith(expected),
                f"{model_filename!r} does not end with {expected!r}",
            )
        self.assertEqual(m.call_count, 2)

    def test_get_model_filename_default_when_neither_file_exists(self):
        m = mock_open()
        m.side_effect = OSError
        with patch("geomag.open", m):
            geo_mag = GeoMag()
            model_filename = geo_mag._get_model_filename()
            expected = get_os_based_test_path("wmm/WMM.COF")
            self.assertTrue(
                model_filename.endswith(expected),
                f"{model_filename!r} does not end with {expected!r}",
            )
        self.assertEqual(m.call_count, 2)

    def test_get_model_filename_different(self):
        geo_mag = GeoMag(coefficients_file="wmm/WMM_NEW.COF")
        model_filename = geo_mag._get_model_filename()
        self.assertEqual(model_filename[-15:], "wmm/WMM_NEW.COF")

    def test_get_model_filename_fullpath(self):
        geo_mag = GeoMag(coefficients_file="/wmm/WMM_NEW.COF")
        model_filename = geo_mag._get_model_filename()
        self.assertEqual(model_filename, "/wmm/WMM_NEW.COF")

    # def test_load_coefficients_invalid_header(self):
    #     geo_mag = GeoMag(coefficients_file="../test/test_files/INVALID_HEADER.COF")
    #     with self.assertRaisesRegex(ValueError, "Invalid header in model file"):
    #         geo_mag.calculate(0, 80, 0, 2030)

    # def test_load_coefficients_invalid_row(self):
    #     geo_mag = GeoMag(coefficients_file="../test/test_files/INVALID_ROW.COF")
    #     with self.assertRaisesRegex(ValueError, "Corrupt record in model file"):
    #         geo_mag.calculate(0, 80, 0, 2030)

    # def test_load_coefficients_invalid_row_data(self):
    #     geo_mag = GeoMag(coefficients_file="../test/test_files/INVALID_ROW_DATA.COF")
    #     with self.assertRaisesRegex(ValueError, "Corrupt record in model file"):
    #         geo_mag.calculate(0, 80, 0, 2030)

    def test_load_coefficients_maxord(self):
        maxord_11_value = -3.1868
        maxord_12_value = -3.1796
        self.assertNotEqual(maxord_11_value, maxord_12_value)
        geo_mag = GeoMag(coefficients_file="wmm/WMM_2025.COF")
        geo_mag._maxord = 11
        result = geo_mag.calculate(0, 80, 0, 2025)
        self.assertAlmostEqual(result.d, maxord_11_value, 4)
        geo_mag = GeoMag(coefficients_file="wmm/WMM_2025.COF")
        geo_mag._maxord = 12
        result = geo_mag.calculate(0, 80, 0, 2025)
        self.assertAlmostEqual(result.d, maxord_12_value, 4)

    def test_load_coefficients_missing_file(self):
        geo_mag = GeoMag(coefficients_file="missing.cof")
        with self.assertRaisesRegex(FileNotFoundError, "No such file or directory"):
            geo_mag.calculate(0, 80, 0, 2030)

    def test_multiple_parameters_supplied_raises(self):
        with self.assertRaisesRegex(
            ValueError,
            "Only one of coefficients_file, coefficients_data, base_year can be set.",
        ):
            GeoMag(coefficients_file="wmm/WMM_2025.COF", coefficients_data=WMM_2025)
        with self.assertRaisesRegex(
            ValueError,
            "Only one of coefficients_file, coefficients_data, base_year can be set.",
        ):
            GeoMag(coefficients_file="wmm/WMM_2025.COF", base_year=2025)
        with self.assertRaisesRegex(
            ValueError,
            "Only one of coefficients_file, coefficients_data, base_year can be set.",
        ):
            GeoMag(coefficients_data=WMM_2025, base_year=2025)
        with self.assertRaisesRegex(
            ValueError,
            "Only one of coefficients_file, coefficients_data, base_year can be set.",
        ):
            GeoMag(
                coefficients_file="wmm/WMM_2025.COF",
                coefficients_data=WMM_2025,
                base_year=2025,
            )

    def test_property_life_span(self):
        geo_mag = GeoMag(coefficients_file="wmm/WMM_2025.COF")
        self.assertTupleEqual(geo_mag.life_span, (2025.0, 2030.0))

    def test_property_model(self):
        geo_mag = GeoMag(coefficients_file="wmm/WMM_2025.COF")
        self.assertEqual(geo_mag.model, "WMM-2025")
