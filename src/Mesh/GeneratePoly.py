#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Generate model's output files
#
#
# John Diaz June 2023
#

import Model


if __name__ == '__main__':

    print()
    print(" ************************************************ ")
    print(" *        Starting Generate Model Program       * ")
    print(" ************************************************ ")
    print()

    # Fault Name
    name = 'LAquilaCirella03_eff_dhF1000'

    # Output dir for topo file
    out_dir = 'Outputs/'

    # X and Y UTM limits (Km) for Central Italy Model 13 Stations
    x_min = 300.0
    x_max = 420.0
    y_min = 4610.0
    y_max = 4750.0
    z_min = -60.0


    # Creates a new instance of the Model class
    LAquilaModel = Model.Model(name, x_min, x_max, y_min, y_max)
    print(LAquilaModel)
    LAquilaModel.generate_model_topo(zone_number=33, zone_letter='N')
#    LAquilaModel.save(out_dir)





