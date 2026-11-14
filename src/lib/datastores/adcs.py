"""
Module of objects and classes which store data pertaining to the current state of the
ADCS system for all tasks to update and use. Readings that have not yet been initialized
are set to `None` throughout this module.
"""

class Datastore:
    """
    Datastore class for adcs processes. Holds time, sensor, and attitude data to be used system-wide
    """

    # TLE String
    TLE = """ISS (ZARYA)\n
1 25544U 98067A   26166.51237796  .00007685  00000-0  14626-3 0  9999\n
2 25544  51.6337 308.3821 0004850 189.0196 171.0706 15.49243792571497"""

    # Action types
    DETUMBLE = 0
    POINT_TO_SUN = 1
    POINT_TO_EARTH = 2
    NOMINAL_PROCESSES = 3

    # Constant matrices for MEKF
    CV_MATRIX = (
        None
    )
    GYRO_NOISE = (
        None
    )
    MEAS_NOISE = (
        None
    )

    def __init__(self):
        self.time: Time = Time()
        self.sensor: SensorData = SensorData()
        self.quaternion = (
            None  # Quaternion representing attitude from body frame to inertial frame
        )
        self.mode = self.DETUMBLE
        # import tle
        # self.satrec: tle.Satrec = tle.Satrec.from_tle_str(self.TLE)

class Time:
    """
    Time helper class
    """

    def __init__(self):
        self.ymdhms = (0, 0, 0, 0, 0, 0) # year month day hour minute second
        self.current_time = None
        self.last_cdh_update = None
        self.update_interval = 1.0  # secondsgit
        self.time_since_last_mekf = 0.0  # dt from tasks/mekf

    @classmethod
    # pylint: disable = too-many-positional-arguments
    def from_ymdhms(cls, year, mon, day, hr, minute, sec):
        """
        Creating a time object from just year, month, day, hour, minute, second
        Minimal class to be used in time_sgp4 that returns a tuple with said values
        """
        obj = cls()
        obj.ymdhms = (year, mon, day, hr, minute, sec)
        return obj

    def time_sgp4(self):
        """
        Time function to return the time in a format that is minimal for SGP4
        Returns a tuple of year, month, day, hour, minute, second
        """
        return self.ymdhms

class SensorData:
    """
    Sensor helper class
    """
    def __init__(self):
        self.sun = None
        self.magnetometer = None
        self.gyroscope = None
