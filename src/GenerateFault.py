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

from src.Classes import Fault

if __name__ == '__main__':

    print()
    print(" ************************************************ ")
    print(" *        Starting Generate Fault Program       * ")
    print(" ************************************************ ")
    print()

    # Fault Name
    name = 'LAquilaCirella03_eff'

    # Output dir for topo file
    out_dir = 'Fault/Outputs/'

    dh_f = 1.0  # Output subfaults size in km

    # Creates a new instance of the Fault class
    LAquilaFault = Fault.Fault(name, dh_f)
    print(LAquilaFault)

    LAquilaFault.load_mat_file('s2009LAQUIL03CIRE')
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















