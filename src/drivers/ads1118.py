"""
Driver module for the ADS1118 SPI analog-to-digital converter from Texas Instruments.
"""

import time
import asyncio

import digitalio

import pin_manager


class MuxSelection:
    """
    Data points  that can be read by an ADS1118. Includes all single-ended and
    differential voltage inputs, as well as the internal temperature sensor.

    Values correspond to the MUX bitfield in the Config Register, except for the
    `TEMPERATURE` selection, which is handled through a separate TS_MODE field
    and is given the special flag value 0xFF.
    """

    CH0_SINGLE_END = 4
    CH1_SINGLE_END = 5
    CH2_SINGLE_END = 6
    CH3_SINGLE_END = 7
    CH0_CH1_DIFF = 0
    CH0_CH3_DIFF = 1
    CH1_CH3_DIFF = 2
    CH2_CH3_DIFF = 3
    TEMPERATURE = 255


class InputRange:
    """
    Internal amplifier programmable gains that can be used on an ADS1118.

    Values correspond to the PGA bitfield in the Config Register. For the 0.256V FSR PGA
    setting (FSR_0_256V), the constant value 0b101 is arbitrarily selected from the valid
    set of settings [0b101, 0b110, 0b111]. Note that changing the internal amplifier gain
    for the ADS1118 affects the signal-to-noise ratio of the measurements themselves.
    """

    FSR_6_144V = 0
    FSR_4_096V = 1
    FSR_2_048V = 2
    FSR_1_024V = 3
    FSR_0_512V = 4
    FSR_0_256V = 5  # 6 and 7 are also valid here


class SamplingRate:
    """
    Sampling rates, in samples per second, that can be used on an ADS1118.

    Values correspond to the DR bitfield in the Config Register. Note that changing the sample
    rate for the ADS1118 affects the signal-to-noise ratio of the measurements themselves.
    """

    RATE_8 = 0
    RATE_16 = 1
    RATE_32 = 2
    RATE_64 = 3
    RATE_128 = 4
    RATE_250 = 5
    RATE_475 = 6
    RATE_860 = 7


ADS1118_LSB_SIZES = dict(
    [
        (InputRange.FSR_6_144V, 187.5e-6),
        (InputRange.FSR_4_096V, 125e-6),
        (InputRange.FSR_2_048V, 62.5e-6),
        (InputRange.FSR_1_024V, 31.25e-6),
        (InputRange.FSR_0_512V, 15.625e-6),
        (InputRange.FSR_0_256V, 7.8125e-6),
    ]
)
ADS1118_SPS_DELAYS = dict(
    [
        (SamplingRate.RATE_8, 0.125),
        (SamplingRate.RATE_16, 0.063),
        (SamplingRate.RATE_32, 0.032),
        (SamplingRate.RATE_64, 0.016),
        (SamplingRate.RATE_128, 0.008),
        (SamplingRate.RATE_250, 0.004),
        (SamplingRate.RATE_475, 0.003),
        (SamplingRate.RATE_860, 0.002),
    ]
)
ADS1118_SPI_RESET_TIME = 0.030  # ideally 28ms, but give it some wiggle room


class Ads1118:
    """
    Driver class for the ADS1118 SPI analog-to-digital converter from Texas Instruments.
    """

    def __init__(self, sck, mosi, miso, ss):
        pm = pin_manager.PinManager.get_instance()
        self.spi_bus = pm.create_spi(sck, mosi, miso)
        self.drdy_gpio = pm.create_digital_in_out(miso)
        self.ss_gpio = pm.create_digital_in_out(ss)
        with self.ss_gpio as ss_gpio:
            ss_gpio.direction = digitalio.Direction.OUTPUT
            ss_gpio.value = True

    # Returns either the voltage in volts, or the temperature in degrees Celsius
    async def take_sample(
        self,
        channel,
        input_range=InputRange.FSR_4_096V,
        sample_rate=SamplingRate.RATE_128,
    ):
        """
        Asynchronous coroutine to sample the ADC on with a given set of settings.

        First, a single-shot conversion command is sent to the ADC with the settings specified.
        Then, the coroutine yields for other tasks for a fixed amount of time based on the
        selected sampling rate to allow the ADC to process the data. Finally, the ADC is polled
        for data readiness and, once readiness is confirmed, the data is read from the ADC.

        In the case that an issue with the ADS1118 prevents data from being ready on schedule,
        the driver will begin the process of resetting the ADC. This process takes 28ms, so in
        the unlikely case that this process needs to occur, overall performance may be degraded.

        The measurement defaults to a full-scale range of 4.096V and a sampling rate of 128
        samples per second if not otherwise specified. The channel selected must be specified
        explicitly and has no default.
        """
        Ads1118._check_sampling_params(channel, input_range, sample_rate)
        transmit_buffer = Ads1118._build_config_register_bytearray(
            channel, input_range, sample_rate
        )
        receive_buffer = bytearray([0, 0])

        data_ready = False

        while not data_ready:

            # send data-getting command
            with self.spi_bus as spi, self.ss_gpio as ss:
                spi.try_lock()
                spi.configure(baudrate=1000000, polarity=0, phase=1)
                ss.direction = digitalio.Direction.OUTPUT
                ss.value = False
                spi.write_readinto(transmit_buffer, receive_buffer)
                ss.value = True
                spi.unlock()

            # wait for data to be ready
            await asyncio.sleep(ADS1118_SPS_DELAYS[sample_rate])

            # check if data is ready
            with self.drdy_gpio as drdy, self.ss_gpio as ss:
                ss.direction = digitalio.Direction.OUTPUT
                ss.value = False
                # busy-wait for CS to DRDY propogation time, unless it takes so long that the
                # ADC resets its SPI peripheral (at which point, quit looking for DRDY and retry
                # the whole transaction)
                # don't do this async because we're currently holding onto hardware
                # we've already waited the delay time above so this should be near-instant
                t0 = time.monotonic_ns()
                while (not data_ready) and (
                    (time.monotonic_ns() - t0) < ADS1118_SPI_RESET_TIME
                ):
                    data_ready = not drdy.value
                ss.value = True

        transmit_buffer[0] = transmit_buffer[0] & 0x7F
        with self.spi_bus as spi, self.ss_gpio as ss:
            spi.try_lock()
            spi.configure(baudrate=1000000, polarity=0, phase=1)
            ss.direction = digitalio.Direction.OUTPUT
            ss.value = False
            spi.write_readinto(transmit_buffer, receive_buffer)
            ss.value = True
            spi.unlock()

        return (
            Ads1118._temperature_from_bytes(receive_buffer)
            if (channel == MuxSelection.TEMPERATURE)
            else Ads1118._voltage_from_bytes(receive_buffer, input_range)
        )

    @staticmethod
    def _check_channel_param(channel):
        if channel is not MuxSelection.TEMPERATURE:
            assert isinstance(channel, int)
            assert channel >= 0
            assert channel < 8

    @staticmethod
    def _check_fsr_param(input_range):
        assert isinstance(input_range, int)
        assert input_range >= 0
        assert input_range < 8

    @staticmethod
    def _check_sps_param(sample_rate):
        assert isinstance(sample_rate, int)
        assert sample_rate >= 0
        assert sample_rate < 8

    @staticmethod
    def _check_sampling_params(channel, input_range, sample_rate):
        Ads1118._check_channel_param(channel)
        Ads1118._check_fsr_param(input_range)
        Ads1118._check_sps_param(sample_rate)

    @staticmethod
    def _build_config_register_bytearray(channel, input_range, sample_rate):
        return bytearray(
            [
                (0b1 << 7)
                | ((channel & 0b111) << 4)
                | ((input_range & 0b111) << 1)
                | 1,
                ((sample_rate & 0b111) << 5)
                | ((channel == MuxSelection.TEMPERATURE) << 4)
                | 0b1010,
            ]
        )

    @staticmethod
    def _int_from_two_bytes_signed_be(buffer: bytearray):
        output = 0
        if buffer[0] & 0x80:
            output -= 65536
        output += 256 * buffer[0]
        output += buffer[1]
        return output

    @staticmethod
    def _temperature_from_bytes(receive_buffer):
        reading = Ads1118._int_from_two_bytes_signed_be(receive_buffer) >> 2
        return reading * 0.03125

    @staticmethod
    def _voltage_from_bytes(receive_buffer, fsr):
        lsb_size = ADS1118_LSB_SIZES[fsr]
        return Ads1118._int_from_two_bytes_signed_be(receive_buffer) * lsb_size
