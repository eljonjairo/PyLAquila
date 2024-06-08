import copy

import utm
import pygmt
import numpy as np
from scipy.interpolate import griddata
from scipy.spatial import Delaunay
from scipy.spatial import distance
import plotly.graph_objs as go
from tqdm import tqdm
import plotly.figure_factory as ff


class Model:
    def __init__(self, name, x_min, x_max, y_min, y_max):
        self.name = name
        self.x_min = x_min  # x minimum border of the model
        self.x_max = x_max  # x max border of the model
        self.y_min = y_min  # y minimum border of the model
        self.y_max = y_max  # y maximum border of the model
        self.z_min = None  # z bottom border of the model
        self.z_max = None  # z top border of the model
        self.x_topo = None  # x coordinate of topo points vector
        self.y_topo = None  # y coordinate of topo points vector
        self.z_topo = None  # z coordinate of topo points vector
        self.x_rect = None  # x coordinate of rectangular model points vector
        self.y_rect = None  # x coordinate of rectangular model points vector
        self.z_rect = None  # x coordinate of rectangular model points vector
        self.x_model = None  # x coordinate of the whole 3d model
        self.y_model = None  # y coordinate of the whole 3d model
        self.z_model = None  # z coordinate of the whole 3d model
        self.tri_global = None  # Triangulation of the model facets
        self.trib_marker = None  # Markers of the model facets

    def __str__(self):
        model = " Model's name: " + str(self.name) + "\n"
        return model

    def generate_topo_model(self, zone_number, zone_letter, dh_topo):
        print(" Generate model topography x, y, and z points ")

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

        x_layer = np.arange(self.x_min, self.x_max + dh, dh)
        y_layer = np.arange(self.y_min, self.y_max + dh, dh)

        X, Y = np.meshgrid(x_layer, y_layer, indexing='xy')
        x_layer = X.flatten(order='F').transpose()
        y_layer = Y.flatten(order='F').transpose()

        # First layer at z = 0
        x_rect = copy.deepcopy(x_layer)
        y_rect = copy.deepcopy(y_layer)
        z_rect = np.zeros_like(x_rect)

        for izl in z_layer:
            z_tmp = -np.ones_like(x_layer) * izl
            x_rect = np.concatenate((x_rect, x_layer))
            y_rect = np.concatenate((y_rect, y_layer))
            z_rect = np.concatenate((z_rect, z_tmp))

        # Fault coordinates
        fcoor = np.array([in_fault.xf, in_fault.yf, in_fault.zf]).transpose()

        print("Adding nodes with distance to the fault > dist_to_fault")

        for icoor in tqdm(fcoor):
            xyz = np.array([x_rect, y_rect, z_rect]).transpose()
            nx = int(xyz.size / 3)
            dist = np.zeros(nx, )
            for ix in range(nx):
                dist[ix] = distance.euclidean(xyz[ix], icoor)
            i_del = np.where(dist < dist_to_fault)
            x_rect = np.delete(x_rect, i_del)
            y_rect = np.delete(y_rect, i_del)
            z_rect = np.delete(z_rect, i_del)

        self.z_min = np.min(z_rect)
        self.x_rect = x_rect
        self.y_rect = y_rect
        self.z_rect = z_rect

        print(f" Rectangular model created with {x_rect.size} points")

    def build_triangulate_3d_model(self, in_fault):
        # Start building the model adding the topo and its triangulation
        print(" Adding topo nodes")
        x_model = copy.deepcopy(self.x_topo)
        y_model = copy.deepcopy(self.y_topo)
        z_model = copy.deepcopy(self.z_topo)

        print(" Adding topo triangulation")
        XY_topo = np.array([self.x_topo, self.y_topo]).transpose()
        tri_global = Delaunay(XY_topo).simplices
        n_tri = int(tri_global.size / 3)
        trib_marker = np.zeros(n_tri, )

        print(" Adding fault nodes")
        nx = len(x_model)
        x_model = np.concatenate((x_model, in_fault.xf))
        x_model = np.concatenate((x_model, in_fault.xf_add))
        y_model = np.concatenate((y_model, in_fault.yf))
        y_model = np.concatenate((y_model, in_fault.yf_add))
        z_model = np.concatenate((z_model, in_fault.zf))
        z_model = np.concatenate((z_model, in_fault.zf_add))

        print(" Adding fault triangulation")
        tri_global = np.concatenate((tri_global, in_fault.tri + nx))
        trib_marker = np.concatenate((trib_marker, in_fault.trib_marker))

        print(" Adding rectangular model nodes")
        x_model = np.concatenate((x_model, self.x_rect))
        y_model = np.concatenate((y_model, self.y_rect))
        z_model = np.concatenate((z_model, self.z_rect))

        print(" Adding rectangular model triangulation")
        x_min = np.min(x_model)
        x_max = np.max(x_model)
        y_min = np.min(y_model)
        y_max = np.max(y_model)
        z_min = np.min(z_model)

        tri_w = facet_x(x_model, y_model, z_model, x_min)
        tri_global = np.concatenate((tri_global, tri_w))
        tri_bw = np.zeros(int(tri_w.size / 3), )
        trib_marker = np.concatenate((trib_marker, tri_bw))
        print(f" Facets at x: {x_min} Km")

        tri_e = facet_x(x_model, y_model, z_model, x_max)
        tri_global = np.concatenate((tri_global, tri_e))
        tri_be = np.zeros(int(tri_e.size / 3), )
        trib_marker = np.concatenate((trib_marker, tri_be))
        print(f" Facets at x: {x_max} Km")

        tri_s = facet_y(x_model, y_model, z_model, y_min)
        tri_global = np.concatenate((tri_global, tri_s))
        tri_bs = np.zeros(int(tri_s.size / 3), )
        trib_marker = np.concatenate((trib_marker, tri_bs))
        print(f" Facets at y: {y_min} Km")

        tri_n = facet_y(x_model, y_model, z_model, y_max)
        tri_global = np.concatenate((tri_global, tri_n))
        tri_bn = np.zeros(int(tri_n.size / 3), )
        trib_marker = np.concatenate((trib_marker, tri_bn))
        print(f" Facets at y: {y_max} Km")

        tri_z = facet_z(x_model, y_model, z_model, z_min)
        tri_global = np.concatenate((tri_global, tri_z))
        tri_bz = np.zeros(int(tri_z.size / 3), )
        trib_marker = np.concatenate((trib_marker, tri_bz))
        print(f" Facets at z: {z_min} Km")

        self.x_model = x_model
        self.y_model = y_model
        self.z_model = z_model
        self.tri_global = tri_global
        self.trib_marker = trib_marker

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

    def plot_model_rect_3d(self):

        data = [go.Scatter3d(x=self.x_rect, y=self.y_rect, z=self.z_rect,
                             mode='markers', marker=dict(size=12,
                                                         color=self.z_rect,
                                                         colorscale='Viridis',
                                                         opacity=0.8))]


        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                      eye=dict(x=-1.7, y=-2.0, z=1.2))

        margin = dict(r=30, l=10, b=30, t=20)
        tickfont = dict(color="black", size=16, family="Arial Black")

        title = dict(text="<b> Rectangular model </b>",
                     font_family="Arial Black", font_color="black", x=0.5,
                     y=0.85)

        xaxis = dict(title="<b>X (km)</b>", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=4, range=[self.x_min, self.x_max],
                     tickfont=tickfont, showbackground=True)

        yaxis = dict(title="<b>Y (km)</b> ", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=4, range=[self.y_min, self.y_max],
                     tickfont=tickfont, showbackground=True)

        zaxis = dict(title="<b>Z (Km )</b>", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=4, range=[self.z_min, self.z_max],
                     tickfont=tickfont, showbackground=True)

        scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
                     aspectmode='cube')

        layout = go.Layout(scene=scene, margin=margin, width=1200, height=800,
                           title=title)

        fig = go.Figure(data=data)
        fig.show(renderer="browser")

    def plot_model_rect(self):
        tickfont = dict(color="black", size=16, family="Arial Black")
        title = dict(text="<b> Rectangular model </b>",
                     font_family="Arial Black", font_color="black", x=0.5,
                     y=0.85)

        xaxis = dict(title="<b>X (km)</b>", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=4, range=[self.x_min, self.x_max],
                     tickfont=tickfont, showbackground=True)

        yaxis = dict(title="<b>Y (km)</b> ", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=4, range=[self.y_min, self.y_max],
                     tickfont=tickfont, showbackground=True)

        zaxis = dict(title="<b>Z (Km )</b>", showgrid=True, showticklabels=True,
                     gridcolor="black", nticks=4, range=[self.z_min, self.z_max],
                     tickfont=tickfont, showbackground=True)

        scene = dict(xaxis=xaxis, yaxis=yaxis, zaxis=zaxis, aspectmode='cube')
        margin = dict(r=30, l=10, b=30, t=20)

        layout = go.Layout(scene=scene, margin=margin, width=1200, height=800,
                           title=title)
        fig = go.Figure(data=[go.Scatter3d(x=self.x_rect, y=self.y_rect,
                                           z=self.z_rect, mode='markers',
                                           marker=dict(size=2,
                                                       color=self.z_rect,
                                                       colorscale='Viridis',
                                                       opacity=0.9))], layout=layout)
        fig.show(renderer="browser")

    def plot_triangulation(self):
        tickfont = dict(color="black", size=19, family="Arial Black")

        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                      eye=dict(x=-1.0, y=-1.0, z=0.7))

        label_z = dict(text="<b>z (Km)</b>", font_family="Arial Black",
                       font_color="black", font_size=26)

        xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True,
                     gridcolor="white", showticklabels=True, nticks=2,
                     range=[self.x_min, self.x_max], tickfont=tickfont,
                     showbackground=True)

        yaxis = dict(title="<b> yUTM (Km) </b>", showgrid=True,
                     showticklabels=True,
                     gridcolor="white", nticks=2, range=[self.y_min, self.y_max],
                     tickfont=tickfont, showbackground=True)

        zaxis = dict(title=label_z, showgrid=True, showticklabels=True,
                     gridcolor="white", nticks=4, range=[self.z_min, self.z_max],
                     tickfont=tickfont, showbackground=True)

        margin = dict(r=10, l=20, b=20, t=20)

        scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
                     aspectmode='cube')

        title = dict(text="<b>Fault Triangulation</b>", font_family="Arial Black",
                     font_color="black", font_size=24, x=0.5, y=0.97)

        layout = go.Layout(scene=scene, margin=margin, width=1600,
                           height=1200, title=title)

        fig = ff.create_trisurf(x=self.x_model, y=self.y_model, z=self.z_model,
                                colormap=[(1.0, 1.0, 0.6), (0.95, 0.95, 0.6)],
                                simplices=self.tri_global, plot_edges=True,
                                show_colorbar=False)

        fig.update_layout(layout)
        fig.update(layout_coloraxis_showscale=False)
        fig.show(renderer="browser")

    def save(self, out_dir):
        outTopofile = out_dir + self.out_name + ".pickle"


def facet_x(x, y, z, xb):
    ix = np.where(x == xb)
    ye = y[ix]
    ze = z[ix]
    yz = np.array([ye, ze]).transpose()
    tri = Delaunay(yz).simplices
    n_tri = int(tri.size / 3)
    tri_x = np.zeros((n_tri, 3))
    jx = np.asarray(ix).transpose()

    for itri in range(0, n_tri):
        xtri = tri[itri, 0]
        ytri = tri[itri, 1]
        ztri = tri[itri, 2]
        tri_x[itri, 0] = jx[xtri]
        tri_x[itri, 1] = jx[ytri]
        tri_x[itri, 2] = jx[ztri]

    return tri_x


def facet_y(x, y, z, yb):
    iy = np.where(y == yb)
    xe = x[iy]
    ze = z[iy]
    xz = np.array([xe, ze]).transpose()
    tri = Delaunay(xz).simplices
    n_tri = int(tri.size / 3)
    tri_y = np.zeros((n_tri, 3))
    jy = np.asarray(iy).transpose()

    for itri in range(0, n_tri):
        xtri = tri[itri, 0]
        ytri = tri[itri, 1]
        ztri = tri[itri, 2]
        tri_y[itri, 0] = jy[xtri]
        tri_y[itri, 1] = jy[ytri]
        tri_y[itri, 2] = jy[ztri]

    return tri_y


def facet_z(x, y, z, zb):
    iz = np.where(z == zb)
    xe = x[iz]
    ye = y[iz]
    xy = np.array([xe, ye]).transpose()
    tri = Delaunay(xy).simplices
    n_tri = int(tri.size / 3)
    tri_z = np.zeros((n_tri, 3))
    jz = np.asarray(iz).transpose()

    for itri in range(0, n_tri):
        xtri = tri[itri, 0]
        ytri = tri[itri, 1]
        ztri = tri[itri, 2]
        tri_z[itri, 0] = jz[xtri]
        tri_z[itri, 1] = jz[ytri]
        tri_z[itri, 2] = jz[ztri]

    return tri_z
