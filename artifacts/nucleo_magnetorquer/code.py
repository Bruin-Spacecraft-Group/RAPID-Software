"""
Testing module for the hbridge magnetorquer controller
"""


import microcontroller as mc

import drivers.h_bridge as hb # remove drivers. after development

# en, nsleep, ph on the z-axis magnetorquer
mt_x = hb.HBridge(mc.pin.PA00, mc.pin.PA01, mc.pin.PA04)

if __name__ == "__main__":
    
    pass
