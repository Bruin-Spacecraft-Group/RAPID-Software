"""
Driver module for the DRV8876 H-Bridge Motor Driver

H-bridge is set to operate purely in EN/PH mode
"""

import digitalio

class HBridge:
    """
    Class for the H-Bridge Motor Driver
    """

    def __init__(self, en, ph, nsleep):
        # in en/ph operation
        self.en = digitalio.DigitalInOut(en)
        self.en.direction = digitalio.Direction.OUTPUT

        self.ph = digitalio.DigitalInOut(ph)
        self.ph.direction = digitalio.Direction.OUTPUT

        self.nsleep = digitalio.DigitalInOut(nsleep)
        self.nsleep.direction = digitalio.Direction.OUTPUT

        # env variables
        self.manual_wake : bool = False

    def wake(self):
        """
        Both internal and manual wake function.
        """

        self.nsleep.value = True
        self.manual_wake = False

    def sleep(self, until_manual_wake = False):
        """
        Sleep function to set the h-bridge into low power mode, in light sleep, wake up using any movement function

        :param until_manual_wake: deep sleep (until the wake function is directly called) if True, light sleep if false
        """
        self.nsleep.value = False
        self.manual_wake = until_manual_wake

    # these functions might look more redundant but I think this will ultimately just depend on who codes with it.
    # may delete later

    def forward(self):
        """
        Full speed fwd
        """
        if not self.manual_wake:
            self.wake()

            self.en.value = True
            self.ph.value = True

    def reverse(self):
        """
        Full speed rev
        """
        if not self.manual_wake:
            self.wake()

            self.en.value = True
            self.ph.value = False

    def brake(self):
        """
        Brake function, en/ph mode does NOT have coast.
        """
        if not self.manual_wake:
            self.wake()

            self.en.value = False

    def spin(self, d):
        """
        Spin Function that combines forward, reverse, and brake into one.

        :param d: direction that is either forward (1) brake (0) or reverse (-1)
        """
        if not self.manual_wake:
            self.wake()
        else:
            return

        if d == 1:
            self.en.value = True
            self.ph.value = True
        elif d == 0:
            self.en.value = False
        elif d == -1:
            self.en.value = True
            self.ph.value = False
