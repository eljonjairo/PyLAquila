import utm
import pygmt
import numpy as np
import plotly.io as io
import pickle
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
from scipy.spatial import distance
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import plotly.express as px
from tqdm import tqdm

class Model:
    def __init__(self, name, x_min, x_max, y_min, y_max):
        self.name = name
        self.x_min = x_min  # x minimum border of the model
        self.x_max = x_max  # x max border of the model
        self.y_min = y_min  # y minimum border of the model
        self.y_max = y_max  # y maximum border of the model
        self.x_topo = None  # x coordinate of topo points vector
        self.y_topo = None  # y coordinate of topo points vector
        self.z_topo = None  # z coordinate of topo points vector

    def __str__(self):
        model = " Model's name: " + str(self.name) + "\n"
        return model

    def generate_topo_model(self, zone_number, zone_letter, dh_topo):
        print(" Generate model topography x, y z points ")

        # Transform from UTM to lat lon coordinates
        lat_min, lon_min = utm.to_latlon(self.x_min * 1000, self.y_min * 1000,
                                         zone_number, zone_letter)
        lat_max, lon_max = utm.to_latlon(self.x_max * 1000, self.y_max * 1000,
                                         zone_number, zone_letter)

        lat_min_region = round(lat_min * 10 - 2) / 10
        lat_max_region = round(lat_max * 10 + 2) / 10
        lon_min_region = round(lon_min * 10 - 2) / 10
        lon_max_region = round(lon_max * 10 + 2) / 10

        # Load Earth relief data for the entire globe and a subset region
        region = [lon_min_region, lon_max_region, lat_min_region, lat_max_region]
        grid = pygmt.datasets.load_earth_relief(resolution="03s", region=region)

        z_region = grid.data.flatten(order='F').transpose() / 1000
        LON_topo, LAT_topo = np.meshgrid(grid.lon.data, grid.lat.data)

        fig = pygmt.Figure()
        pygmt.makecpt(cmap="geo", series=[-5000, 4000])
        fig.grdimage(grid=grid, region=region,
                     frame=['WSrt+t" Topografia Italia Central"', "xa0.4",
                            "ya0.4"])
        fig.colorbar(frame=["xa2000f500+lElevación ", "y+lm"])
        fig.plot(data=np.array([[lon_min, lat_min, lon_max, lat_max]]),
                 style='r+s', pen="2p,black")
        fig.show()

        # (x, y) utm coords (Km) of downloaded region topo
        X_region, Y_region, t1, t2 = utm.from_latlon(LAT_topo, LON_topo,
                                                     zone_number, zone_letter)
        # Convert matrices to vectors
        x_region = X_region.flatten(order='F').transpose() / 1000
        y_region = Y_region.flatten(order='F').transpose() / 1000

        model = " \n Topo downloaded with limits: \n"
        model += (f" X utm (km) [ {np.min(x_region):.3f}, "
                  f"{np.max(x_region):.3f}] \n")
        model += (f" Y utm (km) [{np.min(y_region):.3f}, "
                  f"{np.max(y_region):.3f}] \n")
        print(model)

        # (x, y) utm coords (Km) of model topo
        x_model_topo = np.arange(self.x_min, self.x_max + dh_topo, dh_topo)
        y_model_topo = np.arange(self.y_min, self.y_max + dh_topo, dh_topo)

        X_topo, Y_topo = np.meshgrid(x_model_topo, y_model_topo,
                                     indexing='xy')
        x_topo = X_topo.flatten(order='F').transpose()
        y_topo = Y_topo.flatten(order='F').transpose()
        xy_topo = np.array((x_topo, y_topo)).transpose()
        xy_topo_in = np.array((x_region, y_region)).transpose()

        z_topo = griddata(xy_topo_in, z_region, xy_topo, method='cubic')
        self.x_topo = x_topo.flatten(order='F').transpose()
        self.y_topo = y_topo.flatten(order='F').transpose()
        self.z_topo = z_topo.flatten(order='F').transpose()

        model = " Model topo created with limits: \n"
        model += (f" X utm (km) [ {np.min(self.x_topo):.3f}, "
                  f"{np.max(self.x_topo):.3f}] \n")
        model += (f" Y utm (km) [{np.min(self.y_topo):.3f}, "
                  f"{np.max(self.y_topo):.3f}] \n")
        print(model)

    def generate_rectangular_model(self, dh, z_layer, in_fault, dist_to_fault):
        pass

    # *************************************************************************
    # *                                                                       *
    # *                 Plots methods section                                 *
    # *                                                                       *
    # *************************************************************************
    def plot_model_topo(self):
        """ Plot the input model topo in x, y, z coordinates (Km)"""
        colorbar = dict(lenmode='fraction', len=0.5, thickness=20,
                        bordercolor="black", title="<b> z (Km) </b>", x=-0.1)

        tickfont = dict(color="black", size=16, family="Arial Black")
        xaxis = dict(title="<b> x (Km) </b>", tickfont=tickfont)
        yaxis = dict(title="<b> y (Km) </b>", tickfont=tickfont)

        data = go.Heatmap(z=self.z_topo, x=self.x_topo, y=self.y_topo,
                          hoverongaps=False, colorscale='viridis',
                          colorbar=colorbar)
        fig = go.Figure(data=data)
        fig.update_layout(xaxis=xaxis, yaxis=yaxis)
        fig.show(renderer="browser")
    def save(self, out_dir):
        outTopofile = out_dir + self.out_name + ".pickle"
