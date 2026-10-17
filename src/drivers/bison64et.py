"""
Driver module for the BiSon64-ET Sun Sensor module.
"""

import digitalio


class SunSensor:

    def __init__(self, q1, q2, q3, q4):
        self.q1 = digitalio.DigitalInOut(q1)
        self.q1 = digitalio.Direction.INPUT
        self.q2 = digitalio.DigitalInOut(q2)
        self.q2 = digitalio.Direction.INPUT
        self.q3 = digitalio.DigitalInOut(q3)
        self.q3 = digitalio.Direction.INPUT
        self.q4 = digitalio.DigitalInOut(q4)
        self.q4 = digitalio.Direction.INPUT

    def get_angles(self):
        pass

    def get_calibrated_angles(self):
        pass
