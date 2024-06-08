#!/usr/bin/env python3
# -*- coding: utf-8 -*

import pickle
from pathlib import Path

if __name__ == '__main__':
    print()
    print(" ************************************************ ")
    print(" *        Starting BroadbandSlip Program        * ")
    print(" ************************************************ ")
    print()

    # Number of broadband scenarios
    n_scen = 4
    # Characteristic asperity size (Km)
    asp_size = 5.0

    # Corner wave number
    kc = 5.0

    # Velocity Model H(Km),Vp(Km/s),Vs(Km/s) and Rho(Kg/m3) Bianchi, Ameri*
    H = [1.00, 1.00, 3.00, 22.0, 15.0, 99.0]
    Vp = [3.16, 4.83, 5.76, 6.51, 7.00, 7.80]
    Vs = [1.70, 2.60, 3.10, 3.50, 3.80, 4.20]
    Rho = [2500, 2840, 2940, 3150, 3260, 3500]
