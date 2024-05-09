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

from src.Fault import Fault

if __name__ == '__main__':

    print()
    print(" ************************************************ ")
    print(" *        Starting Generate Fault Program       * ")
    print(" ************************************************ ")
    print()

    # Fault Name
    name = 'LAquilaCirella03_eff'

    # Output dir for topo file
    out_dir = 'Outputs/'

    dh_f = 0.5  # Output subfaults size in km

    # Creates a new instance of the Fault class
    LAquilaFault = Fault.Fault(name, dh_f)
    print(LAquilaFault)

    LAquilaFault.load_mat_file('s2009LAQUIL03CIRE')
    LAquilaFault.plot_fault_input_slip_2d()
    LAquilaFault.plot_fault_input_slip_3d()
    LAquilaFault.interpolate_coords()
    LAquilaFault.interpolate_slip()
    LAquilaFault.interpolate_rise_time()
    LAquilaFault.interpolate_rupture_time()
    LAquilaFault.plot_fault_xyz_slip_3d()
    LAquilaFault.set_effective_size_fault()

    LAquilaFault.triangulate_fault()
    # LAquilaFault.add_nodes_above_below()
    # LAquilaFault.plot_triangulation()
    #
    # LAquilaFault.write_univector(out_dir)
    #
    # LAquilaFault.save_fault(out_dir)















