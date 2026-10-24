"""
Testing module for the hbridge magnetorquer controller
"""


import microcontroller as mc

import h_bridge as hb

# en, nsleep, ph
mt = hb.HBridge(mc.pin.PA00, mc.pin.PA01, mc.pin.PA04)

if __name__ == "__main__":
    mode = int(input("input mode: "))
    print(mode)

    while True:
        if mode == 0:
            mt.forward()
        elif mode == 1:
            mt.reverse()
        elif mode == 2:
            mt.brake()
        elif mode == 3:
            mt.sleep()
