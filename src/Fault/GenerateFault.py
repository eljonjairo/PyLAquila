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

    dh_f = 500  # Output subfaults size in m

    # Creates a new instance of the Fault class
    LAquilaFault = Fault.Fault(name, dh_f)
    print(LAquilaFault)

    LAquilaFault.load_mat_file('s2009LAQUIL03CIRE')
    LAquilaFault.plot_fault_inputs()
    #LAquilaFault.set_full_fault()
    #LAquilaFault.set_effective_fault()
    #LAquilaFault.plot_xyz_slipin()
    #LAquilaFault.interpolate_xyz_coords()
    #LAquilaFault.interpolate_slip()
    #LAquilaFault.interpolate_rise_time()
    #LAquilaFault.interpolate_rupt_time()
 #
 #   LAquilaFault.plot_xyz_slip()
 #   LAquilaFault.plot_xyz_model_slip()
 #   LAquilaFault.compare_xyz_slip()
    
    # df = px.data.tips()
    # fig = px.scatter(df, x="total_bill", y="tip", color="size",
    #              title="Numeric 'size' values mean continuous color")

    # fig.show()
    
    # LAquilaFault.triangulate_fault()
    # LAquilaFault.add_nodes_above_below()
    # LAquilaFault.plot_triangulation()
    #
    # LAquilaFault.write_univector(out_dir)
    #
    # LAquilaFault.save_fault(out_dir)















