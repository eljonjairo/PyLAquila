#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
    Fault Generator
    Params:
        dh_f: float = Size of output square subfaults
        .mat = Structure matlab file
        name: string = Name of structure matlab file
    Return:
        Fault object with slip, rupture time, and rise time distributions
        over the interpolated fault
"""
from scipy.io import loadmat

from Classes.Fault import Fault
from pathlib import Path

if __name__ == '__main__':

    print()
    print(" ************************************************ ")
    print(" *        Starting Generate Fault Program       * ")
    print(" ************************************************ ")
    print()

    # Fault Name
    name = 'LAquilaCirella03_eff'

    # Output subfaults size in km
    dh_f = 1.0  # Output subfaults size in km

    # Input mat file name
    name_mat = 's2009LAQUIL03CIRE'

    # Input dir mat file
    input_mat = 'Inputs/'

    # Output dir for object file
    out_dir = 'Outputs/Fault/'

    # Creates a new instance of the Fault class
    LAquilaFault = Fault(name, dh_f)
    print(LAquilaFault)

    mat_dir = Path('')
    file_mat = name_mat + '.mat'
    input_mat = mat_dir / input_mat / file_mat

    print(input_mat)

    inFault = loadmat(input_mat)
    LAquilaFault.load_mat_file(inFault[name_mat])
    LAquilaFault.plot_input_slip()
    LAquilaFault.interpolate_coords()
    LAquilaFault.interpolate_slip()
    LAquilaFault.interpolate_rise_time()
    LAquilaFault.interpolate_rupture_time()
    LAquilaFault.plot_slip()
    LAquilaFault.set_effective_size_fault()
    LAquilaFault.plot_slip()

    LAquilaFault.triangulate_fault()
    LAquilaFault.plot_triangulation()

    LAquilaFault.write_uni_vector(out_dir)
    LAquilaFault.save_fault(out_dir)

