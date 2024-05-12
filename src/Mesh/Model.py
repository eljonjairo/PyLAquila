#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Model class
# 

# TO DO LIST:
# Generate Topo
# Generate Velo full
# Generate velo fault
#
# John Diaz June 2023
#

import utm
import pygmt
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.cm
import pickle
import numpy as np
from scipy.interpolate import griddata

class Model:
    def __init__(self, name, x_min, x_max, y_min, y_max):
        self.name = name
        self.x_min = x_min
        self.x_max = x_max
        self.y_min = y_min
        self.y_max = y_max
        self.X_topo_in = None
        self.Y_topo_in = None
        self.Z_topo_in = None
        self.X_topo = None
        self.Y_topo = None
        self.Z_topo = None

    def __str__(self):
        model = " Model's name: " + str(self.name) + "\n"
        model += (" X utm limits (km) [" + str(self.x_min)
                  + ',' + str(self.x_max) + "]\n")
        model += (" Y utm limits (km) [" + str(self.y_min)
                  + ',' + str(self.y_max) + "]\n")

        return model

    def generate_region_topo(self, zone_number, zone_letter):
        print(" Generate model topography of the region around the model ")

        # Transform from UTM to lat lon coordinates
        lat_min, lon_min = utm.to_latlon(self.x_min * 1000, self.y_min * 1000,
                                         zone_number, zone_letter)
        lat_max, lon_max = utm.to_latlon(self.x_max * 1000, self.y_max * 1000,
                                         zone_number, zone_letter)

        lat_min = round(lat_min * 10 - 2) / 10
        lat_max = round(lat_max * 10 + 2) / 10
        lon_min = round(lon_min * 10 - 2) / 10
        lon_max = round(lon_max * 10 + 2) / 10

        # Load Earth relief data for the entire globe and a subset region
        region = [lon_min, lon_max, lat_min, lat_max]
        grid = pygmt.datasets.load_earth_relief(resolution="03s", region=region)

        self.Z_topo_in = grid.data.flatten(order='F').transpose()/1000
        LON_topo, LAT_topo = np.meshgrid(grid.lon.data, grid.lat.data)

        fig = pygmt.Figure()
        pygmt.makecpt(cmap="geo", series=[-6000, 6000])
        fig.grdimage(grid=grid, region=region,
                     frame=['WSrt+t" Topografia Italia Central"',
                            "xa0.4", "ya0.4"])
        fig.colorbar(frame=["xa2000f500+lElevación ", "y+lm"])
        fig.show()

        # (x, y) utm coords (Km)
        X_topo, Y_topo, tmp1, tmp2 = utm.from_latlon(LAT_topo, LON_topo,
                                                     zone_number, zone_letter)
        self.X_topo_in = X_topo.flatten(order='F').transpose()/1000
        self.Y_topo_in = Y_topo.flatten(order='F').transpose()/1000

        model = " Topo downloaded with limits: \n"
        model += (f" X utm (km) [ {np.min(self.X_topo_in):.3f}, "
                 f"{np.max(self.X_topo_in):.3f}] \n")
        model += (f" Y utm (km) [{np.min(self.Y_topo_in):.3f}, "
                  f"{np.max(self.Y_topo_in):.3f}] \n")
        print(model)

    def generate_model_topo(self, dh_topo):

        x_model_topo = np.arange(self.x_min, self.x_max+dh_topo, dh_topo)
        y_model_topo = np.arange(self.y_min, self.y_max+dh_topo, dh_topo)

        X_topo, Y_topo = np.meshgrid(x_model_topo, y_model_topo,
                                     indexing='xy')
        xy_topo = np.array((X_topo, Y_topo)).transpose()
        xy_topo_in = np.array((self.X_topo_in, self.Y_topo_in)).transpose()

        Z_topo = griddata(xy_topo_in, self.Z_topo_in, xy_topo, method='cubic')
        X_topo = X_topo.flatten(order='F').transpose()
        Y_topo = Y_topo.flatten(order='F').transpose()
        Z_topo = Z_topo.flatten(order='F').transpose()

        fig = plt.figure()
        ax = fig.subplots(1,1)
        mp = ax.pcolormesh(X_topo, Y_topo, Z_topo, cmap='terrain', vmin=-1, vmax=4)
        plt.colorbar(mp, location='bottom', label="Elevación (m)", shrink=.6)
        ax.set_aspect('equal', adjustable='box')
        ax.set_title('Topografía Italia Central')
        fig.show()

    def save(self, out_dir):
        outTopofile = out_dir + self.out_name + ".pickle"
        TopoObj = open(outTopofile, 'wb')
        pickle.dump(self, TopoObj)
        TopoObj.close()

        print(f" Topo info saved in:  {outTopofile} ")
        print()

        #cvals  = [ -100, 0, 100, 1000, 2000, 3000, 4000, ]
        #colors = ["blue","darkgreen","green","forestgreen","yellowgreen",
        #          "orange","maroon","sienna","brown","white"]
#
#        norm = plt.Normalize(min(cvals),max(cvals))
#        tuples = list(zip(map(norm,cvals), colors))
#        cTopo = matplotlib.colors.LinearSegmentedColormap.from_list("", tuples)
