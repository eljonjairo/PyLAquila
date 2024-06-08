"""
 Fault Class
"""
import pickle

import numpy as np
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objs as go
import scipy
import utm
from scipy.spatial import distance
from plotly.subplots import make_subplots


class Fault:
    def __init__(self, name: str, dh_f: float):

        self.dh_f = dh_f
        self.name = name + "_dhF" + str(int(self.dh_f * 1000)) + "m"

        self.in_X = None  # Input matriz with the x coordinate of the fault
        self.in_Y = None  # Input matriz with the y coordinate of the fault
        self.in_Z = None  # Input matriz with the z coordinate of the fault
        self.in_LON = None  # Input longitude coord (°)
        self.in_LAT = None  # Input latitude coord (°)
        self.in_SLIP = None  # Input slip distribution (cm)
        self.in_RISE = None  # Input rise time distribution
        self.in_RUPT = None  # Input rupture time distribution
        self.in_stk_vec = None  # Input strike vector
        self.in_dip_vec = None  # Input dip vector
        self.in_dstk_vec = None  # Input delta strike vector
        self.in_ddip_vec = None  # Input delta dip vector
        self.in_ddip = None  # Input size of subfaults in dip direction
        self.in_dstk = None  # Input size of subfaults in strike direction
        self.XF = None  # Interpolated matriz x coordinate matriz
        self.YF = None  # Interpolated matriz y coordinate matriz
        self.ZF = None  # Interpolated matriz z coordinate matriz
        self.xf = None  # Interpolated vector with x coordinate of the fault
        self.yf = None  # Interpolated vector with y coordinate of the fault
        self.zf = None  # Interpolated vector with z coordinate of the fault
        self.xf_add = None  # Extra nodes over and under the subfaults
        self.yf_add = None  # Extra nodes over and under the subfaults
        self.zf_add = None  # Extra nodes over and under the subfaults
        self.SLIP = None  # Interpolated slip distribution
        self.RUPT = None  # Interpolated rupture time distribution
        self.RISE = None  # Interpolated rise time distribution
        self.stk_len = None  # Length in the strike direction
        self.dip_len = None  # Length in the dip direction
        self.stk_vec = None  # Strike vector
        self.dip_vec = None  # Dip vector
        self.stk_uni = None  # Strike vector
        self.dip_uni = None  # Dip vector
        self.nstk = None  # Number of subfaults in strike direction
        self.ndip = None  # Number of subfaults in dip direction
        self.hypo_lon = None  # Hypocenter longitude coordinate
        self.hypo_lat = None  # Hypocenter latitude coordinate
        self.hypo_x = None  # Hypocenter x coordinate
        self.hypo_y = None  # Hypocenter y coordinate
        self.hypo_z = None  # Hypocenter z coordinate
        self.hypo_idip = None  # Hypocenter index in dip direction
        self.hypo_istk = None  # Hypocenter index in strike direction
        self.ddip = None  # Size of subfaults in dip direction
        self.dstk = None  # Size of subfaults in strike direction
        self.n_tri = None  # number of triangles of Delaunay
        self.tri = None  # delaunay triangulation array
        self.facet_norm = None  # normal, strike and dip vector at each facet
        self.trib_marker = None  # Triangulation facets marker

    def __repr__(self):
        return "Fault name: " + str(self.name) + "dhF: " + str(self.dh_f) + " Km"

    def load_mat_file(self, in_fault):
        """
        Load matlab structure file
        :param in_fault: matlab structure with input fault.
        : return:
        Fault object with input file attributes
        """

        # Load original fault coordinates, slip, rise time and rupture time
        self.in_Z = -in_fault['geoZ'][-1][0]
        self.in_LAT = in_fault['geoLAT'][-1][0]
        self.in_LON = in_fault['geoLON'][-1][0]
        self.in_SLIP = in_fault['slipSPL'][-1][0]
        self.in_RISE = in_fault['riseSPL'][-1][0]
        self.in_RUPT = in_fault['timeSPL'][-1][0]

        # Hypocenter coordinates (Km) 
        self.hypo_lon = in_fault['evLON'][-1][0][0][0]
        self.hypo_lat = in_fault['evLAT'][-1][0][0][0]
        self.hypo_z = -in_fault['evDPT'][-1][0][0][0]

        # Hypocenter (x, y) utm coords (Km)
        self.hypo_x, self.hypo_y, tmp1, tmp2 = utm.from_latlon(self.hypo_lat,
                                                               self.hypo_lon,
                                                               33, 'T')
        self.hypo_x = self.hypo_x / 1000
        self.hypo_y = self.hypo_y / 1000

        # (x, y) utm coords (Km)
        self.in_X, self.in_Y, tmp1, tmp2 = utm.from_latlon(self.in_LAT,
                                                           self.in_LON,
                                                           33, 'T')
        self.in_X = self.in_X / 1000
        self.in_Y = self.in_Y / 1000

        print()
        print(f" Fault input loaded successfully ")
        print()

        # Calculate magnitude and unitary delta vector in strike direction
        in_dstk_vec = np.array([self.in_X[0, 1] - self.in_X[0, 0],
                                self.in_Y[0, 1] - self.in_Y[0, 0],
                                self.in_Z[0, 1] - self.in_Z[0, 0]])
        self.in_dstk = np.linalg.norm(in_dstk_vec)
        self.stk_uni = in_dstk_vec / self.in_dstk

        # Calculate magnitude and unitary vector of delta in dip direction
        in_ddip_vec = np.array([self.in_X[1, 0] - self.in_X[0, 0],
                                self.in_Y[1, 0] - self.in_Y[0, 0],
                                self.in_Z[1, 0] - self.in_Z[0, 0]])
        self.in_ddip = np.linalg.norm(in_ddip_vec)
        self.dip_uni = in_ddip_vec / self.in_ddip

        # Calculate length and number of points in strike and dip direction
        in_ndip, in_nstk = self.in_Z.shape
        self.stk_len = round((in_nstk - 1) * self.in_dstk)
        self.dip_len = round((in_ndip - 1) * self.in_ddip)

        # Create input vectors in strike and dip directions
        self.in_stk_vec = np.linspace(0, self.stk_len, in_nstk)
        self.in_dip_vec = np.linspace(0, self.dip_len, in_ndip)

        print()
        print(f" Fault's input dimensions: ")
        print(f" Strike (Km): {self.stk_len:.3f} nstk: {in_nstk} "
              f" dstk (Km): {self.in_dstk:.3f} ")
        print(f" Dip (Km): {self.dip_len:.3f} ndip: {in_ndip} "
              f" ddip (Km): {self.in_ddip:.3f} ")
        print()
        print(f" Hypocenter coordinates: ")
        print(f" hypo_x (km): {self.hypo_x:.3f} hypo_y (km): {self.hypo_y:.3f}"
              f" hypo_z (km): {self.hypo_z:.3f} ")
        print()

    def set_effective_size_fault(self):
        """ Calculate the effective length f the fault following: Source Scaling
        Properties from Finite-Fault-Rupture Models (Mai & Beroza 2000)
        :param: Slip distribution over the entire fault
        :return: Slip distribution over the effective length of the fault
        """

        slip = self.SLIP

        # Cumulative slip in strike and dip direction
        sum_slip_stk = np.sum(slip, axis=0)
        sum_slip_dip = np.sum(slip, axis=1)

        # Calculate the integral, and the effective length
        acorr_stk = np.correlate(sum_slip_stk, sum_slip_stk, mode='full')
        x_stk = np.arange(0, len(acorr_stk), 1, dtype=int) * self.dstk
        Wacf_stk = scipy.integrate.simpson(acorr_stk, x_stk) / max(acorr_stk)

        acorr_dip = np.correlate(sum_slip_dip, sum_slip_dip, mode='full')
        x_dip = np.arange(0, len(acorr_dip), 1, dtype=int) * self.ddip
        Wacf_dip = scipy.integrate.simpson(acorr_dip, x_dip) / max(acorr_dip)

        print(f" Effective length in strike direction: {Wacf_stk:.2f} Km.")
        print(f" Effective length in dip direction: {Wacf_dip:.2f} Km.")

        # Calculate effective number of nodes
        nstk_eff = round(Wacf_stk / self.dstk) + 1
        ndip_eff = round(Wacf_dip / self.ddip) + 1

        # Select the interval in the strike direction to maximize the slip sum
        istk = 0
        jstk = nstk_eff
        istk_eff = istk
        sum_slip_eff = np.sum(sum_slip_stk[istk:jstk])
        while jstk < self.nstk:
            istk += 1
            jstk += 1
            sum_slip_tmp = np.sum(sum_slip_stk[istk:jstk])
            if sum_slip_tmp > sum_slip_eff:
                istk_eff = istk
                sum_slip_eff = sum_slip_tmp

        # Select the interval in the dip direction to maximize the slip sum
        idip = 0
        jdip = ndip_eff

        idip_eff = idip
        sum_slip_eff = np.sum(sum_slip_dip[idip:jdip])
        while jdip < self.ndip:
            idip += 1
            jdip += 1
            sum_slip_tmp = np.sum(sum_slip_dip[idip:jdip])
            if sum_slip_tmp > sum_slip_eff:
                idip_eff = idip
                sum_slip_eff = sum_slip_tmp

        slip_eff = slip[idip_eff:idip_eff + ndip_eff, istk_eff:istk_eff + nstk_eff]

        # Replace fault attributes using the effective size
        self.nstk = nstk_eff
        self.ndip = ndip_eff
        self.SLIP = self.SLIP[idip_eff:idip_eff + ndip_eff,
                    istk_eff:istk_eff + nstk_eff]
        self.RUPT = self.RUPT[idip_eff:idip_eff + ndip_eff,
                    istk_eff:istk_eff + nstk_eff]
        self.RISE = self.RISE[idip_eff:idip_eff + ndip_eff,
                    istk_eff:istk_eff + nstk_eff]
        self.XF = self.XF[idip_eff:idip_eff + ndip_eff,
                  istk_eff:istk_eff + nstk_eff]
        self.YF = self.YF[idip_eff:idip_eff + ndip_eff,
                  istk_eff:istk_eff + nstk_eff]
        self.ZF = self.ZF[idip_eff:idip_eff + ndip_eff,
                  istk_eff:istk_eff + nstk_eff]
        self.stk_vec = self.stk_vec[istk_eff:istk_eff + nstk_eff]
        self.dip_vec = self.dip_vec[idip_eff:idip_eff + ndip_eff]

        self.xf = self.XF.flatten(order='F').transpose()
        self.yf = self.YF.flatten(order='F').transpose()
        self.zf = self.ZF.flatten(order='F').transpose()

        print()
        print(f" Original slip matrix dimensions: {slip.shape} ")
        print(f" Effective slip matrix dimensions: {slip_eff.shape} ")
        print()

    # *************************************************************************
    # *                                                                       *
    # *              Interpolation methods section                            *
    # *                                                                       *
    # *************************************************************************
    def interpolate_coords(self):
        """
        Interpolate (X, Y, Z) coordinates from subfaults size (in_dstk, in_dstk)
        in strike and dip directions to square subfaults of size dh_f
        :param: self object input (in_X, in_Y, in_Z) coordinates
        :return: self object output (X, Y, Z) coordinates
        """

        # Calculate output fault's size
        self.nstk = int(self.stk_len / self.dh_f) + 1
        self.ndip = int(self.dip_len / self.dh_f) + 1

        # Create arrays for interpolated coordinates
        self.XF = np.zeros((self.ndip, self.nstk))
        self.YF = np.zeros((self.ndip, self.nstk))
        self.ZF = np.zeros((self.ndip, self.nstk))

        # Initiate arrays to calculate the index of hypocenter in interpolated
        # XYZ coordinates
        HYPO_DIST = np.zeros((self.ndip, self.nstk))  # Hypocenter distance

        # Vectors of size dh_f in strike and dip directions
        dstk_vec = self.stk_uni * self.dh_f
        ddip_vec = self.dip_uni * self.dh_f

        # Effective size of subfaults in strike and dip directions
        self.dstk = np.linalg.norm(dstk_vec)
        self.ddip = np.linalg.norm(ddip_vec)

        # Create interpolated vectors in strike and dip directions
        self.stk_vec = np.linspace(0, self.stk_len, self.nstk)
        self.dip_vec = np.linspace(0, self.dip_len, self.ndip)

        for istk in range(0, self.nstk):
            delta_stk = istk * dstk_vec
            for idip in range(0, self.ndip):
                delta_dip = idip * ddip_vec
                self.XF[idip, istk] = (self.in_X[0, 0] + delta_stk[0]
                                       + delta_dip[0])
                self.YF[idip, istk] = (self.in_Y[0, 0] + delta_stk[1]
                                       + delta_dip[1])
                self.ZF[idip, istk] = (self.in_Z[0, 0] + delta_stk[2]
                                       + delta_dip[2])

                # Calculate hypocenter distance in XY coords.
                xy = [self.XF[idip, istk], self.YF[idip, istk]]
                HYPO_DIST[idip, istk] = distance.euclidean(xy,
                                                           [self.hypo_x, self.hypo_y])
        # Find the indexes of the hypocenter
        self.hypo_idip, self.hypo_istk = np.where(HYPO_DIST == np.min(HYPO_DIST))

        # Make vertical correction ensuring that hypocenter relays on the fault
        self.ZF += self.hypo_z - self.ZF[self.hypo_idip, self.hypo_istk]

        # From matrix to column vector following fortran
        self.xf = self.XF.flatten(order='F').transpose()
        self.yf = self.YF.flatten(order='F').transpose()
        self.zf = self.ZF.flatten(order='F').transpose()

        print(f" Fault's output dimensions: ")
        print(f" Strike (Km): {self.stk_len:.3f} nstk: {self.nstk} "
              f" dstk (Km): {self.dstk:.3f} ")
        print(f" Dip (Km): {self.dip_len:.3f} ndip: {self.ndip} "
              f" ddip (Km): {self.ddip:.3f} ")
        print()
        print(" (x, y, z) coordinates has been successfully interpolated ")

    def interpolate_slip(self):
        """
        Interpolate slip distribution from subfaults size (in_dstk, in_dstk)
        in strike and dip directions to square subfaults of size dh_f
        :param: self object input in_SLIP, in_stk_vec, in_dip_vec, stk_vec,
                dip_vec
        :return: self object output SLIP
        """
        SlipF = scipy.interpolate.interp2d(self.in_stk_vec, self.in_dip_vec,
                                           self.in_SLIP, kind="cubic")
        self.SLIP = SlipF(self.stk_vec, self.dip_vec)

        print()
        print(" Slip distribution has been successfully interpolated ")

    def interpolate_rise_time(self):
        """
        Interpolate rise time distribution from subfaults size (in_dstk, in_dstk)
        in strike and dip directions to square subfaults of size dh_f
        :param: self object input in_RISE, in_stk_vec, in_dip_vec, stk_vec,
                dip_vec
        :return: self object output RISE
        """
        rise_fun = scipy.interpolate.interp2d(self.in_stk_vec, self.in_dip_vec,
                                              self.in_RISE, kind="cubic")
        self.RISE = rise_fun(self.stk_vec, self.dip_vec)
        print()
        print(" Rise time distribution has been successfully interpolated ")

    def interpolate_rupture_time(self):
        """
        Interpolate slip distribution from subfaults size (in_dstk, in_dstk)
        in strike and dip directions to square subfaults of size dh_f
        :param: self object input in_RUPT, in_stk_vec, in_dip_vec, stk_vec,
                dip_vec
        :return: self object output RUPT
        """
        rupt_fun = scipy.interpolate.interp2d(self.in_stk_vec, self.in_dip_vec,
                                              self.in_RUPT, kind="cubic")
        self.RUPT = rupt_fun(self.stk_vec, self.dip_vec)
        print()
        print(" Rupture time distribution has been interpolated successfully")

    # *************************************************************************
    # *                                                                       *
    # *                 Plots methods section                                 *
    # *                                                                       *
    # *************************************************************************

    def plot_fault_xyz_slip_3d(self):
        """ Plot the interpolated slip distribution in (x, y, z) coordinates """

        colorbar_a = dict(lenmode='fraction', len=0.5, thickness=20,
                          bordercolor="black", title="<b> slip (cm) </b>", x=0.1)

        data_a = [go.Surface(x=self.XF, y=self.YF, z=self.ZF,
                             surfacecolor=self.SLIP,
                             colorscale=px.colors.sequential.Viridis,
                             colorbar=colorbar_a, showscale=True,
                             lighting=dict(ambient=0.7))]

        tickfont_a = dict(color="black", size=16, family="Arial Black")

        camera_a = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                        eye=dict(x=-1.7, y=-2.0, z=1.2))

        margin_a = dict(r=30, l=10, b=30, t=20)

        title_a = dict(text="<b> Interpolated Slip </b>",
                       font_family="Arial Black",
                       font_color="black", x=0.5, y=0.85)

        xaxis_a = dict(title="<b>X (km)</b>", showgrid=True,
                       showticklabels=True,
                       gridcolor="black", nticks=4, range=[360, 400],
                       tickfont=tickfont_a, showbackground=True)

        yaxis_a = dict(title="<b> Y (km)</b> ", showgrid=True,
                       showticklabels=True,
                       gridcolor="black", nticks=4, range=[4670, 4710],
                       tickfont=tickfont_a, showbackground=True)

        zaxis_a = dict(title="<b>z (Km )</b>", showgrid=True,
                       showticklabels=True, gridcolor="black", nticks=4,
                       range=[-15, 0], tickfont=tickfont_a, showbackground=True)

        scene_a = dict(camera=camera_a, xaxis=xaxis_a, yaxis=yaxis_a,
                       zaxis=zaxis_a, aspectmode='cube')

        layout_a = go.Layout(scene=scene_a, margin=margin_a, width=1200,
                             height=800, title=title_a)

        fig = go.Figure(data=data_a, layout=layout_a)
        fig.show(renderer="browser")

    def plot_input_slip(self):
        """ Plot the input slip distribution in lat-lon-z coordinates """

        colorbar = dict(lenmode='fraction', len=0.5, thickness=20,
                        bordercolor="black", title="<b> slip (cm) </b>", x=-0.1)

        data_a = go.Surface(x=self.in_LON, y=self.in_LAT, z=self.in_Z,
                            surfacecolor=self.in_SLIP,
                            colorscale=px.colors.sequential.Viridis,
                            colorbar=colorbar, showscale=False,
                            lighting=dict(ambient=0.7))

        tickfont = dict(color="black", size=18, family="Arial Black")

        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                      eye=dict(x=-1.7, y=-2.0, z=1.2))

        margin = dict(r=10, l=10, b=10, t=70)

        xaxis_a = dict(title="<b> Longitude (°) </b>", showgrid=True,
                       showticklabels=True,
                       gridcolor="black", nticks=4, range=[13.2, 13.7],
                       tickfont=tickfont, showbackground=True)

        yaxis_a = dict(title="<b> Latitude (°) </b> ", showgrid=True,
                       showticklabels=True,
                       gridcolor="black", nticks=4, range=[42.1, 42.5],
                       tickfont=tickfont, showbackground=True)

        zaxis_a = dict(title="<b> z (Km ) </b>", showgrid=True,
                       showticklabels=True, gridcolor="black", nticks=4,
                       range=[-15, 0], tickfont=tickfont, showbackground=True)

        scene_a = dict(camera=camera, xaxis=xaxis_a, yaxis=yaxis_a,
                       zaxis=zaxis_a, aspectmode='cube')

        title = dict(text="<b> Input Slip </b>", font_family="Arial Black",
                     font_color="black", x=0.5, y=0.85)

        fig = make_subplots(rows=1, cols=2, specs=[[{"type": "surface"},
                                                    {"type": "heatmap"}]])

        fig.add_trace(data_a, row=1, col=1)
        fig.update_layout(scene=scene_a, margin=margin)

        xaxis = dict(title="<b> strike (Km) </b>", tickfont=tickfont)
        yaxis = dict(title="<b> dip (Km) </b>", tickfont=tickfont,
                     autorange="reversed")

        colorscale = [[0, 'rgb(250, 250, 250)'], [1.0, 'rgb(255, 255, 255)']]
        contours = dict(coloring='lines',
                        labelfont=dict(size=14, color='white', ),
                        showlabels=True)

        data = go.Heatmap(z=self.in_SLIP, x=self.in_stk_vec, y=self.in_dip_vec,
                          hoverongaps=False, colorscale='viridis',
                          colorbar=colorbar)
        fig.add_trace(data, row=1, col=2)
        fig.update_layout(xaxis=xaxis, yaxis=yaxis, margin=margin, title=title)

        data_c = go.Contour(z=self.in_RUPT, x=self.in_stk_vec, y=self.in_dip_vec,
                            contours=contours, showscale=False,
                            colorscale=colorscale)
        fig.add_trace(data_c, row=1, col=2)
        fig.show(renderer="browser")

    def plot_slip(self):
        """ Plot the interpolated slip distribution in (x, y, z) coordinates """

        colorbar = dict(lenmode='fraction', len=0.5, thickness=20,
                        bordercolor="black", title="<b> slip (cm) </b>", x=-0.1)

        data_a = go.Surface(x=self.XF, y=self.YF, z=self.ZF,
                            surfacecolor=self.SLIP,
                            colorscale=px.colors.sequential.Viridis,
                            colorbar=colorbar, showscale=False,
                            lighting=dict(ambient=0.7))

        tickfont = dict(color="black", size=18, family="Arial Black")

        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                      eye=dict(x=-1.7, y=-2.0, z=1.2))

        margin = dict(r=10, l=10, b=10, t=70)

        xaxis_a = dict(title="<b> Longitude (°) </b>", showgrid=True,
                       showticklabels=True,
                       gridcolor="black", nticks=4, range=[360, 390],
                       tickfont=tickfont, showbackground=True)

        yaxis_a = dict(title="<b> Latitude (°) </b> ", showgrid=True,
                       showticklabels=True,
                       gridcolor="black", nticks=4, range=[4670, 4710],
                       tickfont=tickfont, showbackground=True)

        zaxis_a = dict(title="<b> z (Km ) </b>", showgrid=True,
                       showticklabels=True, gridcolor="black", nticks=4,
                       range=[-15, 0], tickfont=tickfont, showbackground=True)

        scene_a = dict(camera=camera, xaxis=xaxis_a, yaxis=yaxis_a,
                       zaxis=zaxis_a, aspectmode='cube')

        title = dict(text="<b> Input Slip </b>", font_family="Arial Black",
                     font_color="black", x=0.5, y=0.85)

        fig = make_subplots(rows=1, cols=2, specs=[[{"type": "surface"},
                                                    {"type": "heatmap"}]])

        fig.add_trace(data_a, row=1, col=1)
        fig.update_layout(scene=scene_a, margin=margin)

        xaxis = dict(title="<b> strike (Km) </b>", tickfont=tickfont)
        yaxis = dict(title="<b> dip (Km) </b>", tickfont=tickfont,
                     autorange="reversed")

        colorscale = [[0, 'rgb(250, 250, 250)'], [1.0, 'rgb(255, 255, 255)']]
        contours = dict(coloring='lines',
                        labelfont=dict(size=14, color='white', ),
                        showlabels=True)

        data = go.Heatmap(z=self.SLIP, x=self.stk_vec, y=self.dip_vec,
                          hoverongaps=False, colorscale='viridis',
                          colorbar=colorbar)
        fig.add_trace(data, row=1, col=2)
        fig.update_layout(xaxis=xaxis, yaxis=yaxis, margin=margin, title=title)

        data_c = go.Contour(z=self.RUPT, x=self.stk_vec, y=self.dip_vec,
                            contours=contours, showscale=False,
                            colorscale=colorscale)
        fig.add_trace(data_c, row=1, col=2)
        fig.show(renderer="browser")

    def plot_fault_input_slip_2d(self):
        tickfont = dict(color="black", size=12, family="Arial Black")
        xaxis = dict(title="<b> strike (Km) </b>", tickfont=tickfont)
        yaxis = dict(title="<b> dip (Km) </b>", tickfont=tickfont,
                     autorange="reversed")

        fig = go.Figure(data=go.Heatmap(z=self.in_SLIP, x=self.in_stk_vec,
                                        y=self.in_dip_vec, hoverongaps=False))

        fig.update_layout(xaxis=xaxis, yaxis=yaxis)
        colorscale = [[0, 'rgb(250, 250, 250)'], [1.0, 'rgb(255, 255, 255)']]
        contours = dict(coloring='lines', showlabels=True)
        figc = go.Figure(data=go.Contour(z=self.in_RUPT, x=self.in_stk_vec,
                                         y=self.in_dip_vec, contours=contours,
                                         line_width=2, showscale=False,
                                         colorscale=colorscale))
        fig.add_trace(figc.data[0])
        fig.show(renderer="browser")

    def plot_triangulation(self):
        tickfont = dict(color="black", size=19, family="Arial Black")

        camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                      eye=dict(x=-1.0, y=-1.0, z=0.7))

        label_z = dict(text="<b>z (Km)</b>", font_family="Arial Black",
                       font_color="black", font_size=26)

        xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True,
                     gridcolor="white", showticklabels=True, nticks=2,
                     range=[360, 390], tickfont=tickfont, showbackground=True)

        yaxis = dict(title="<b> yUTM (Km) </b>", showgrid=True,
                     showticklabels=True,
                     gridcolor="white", nticks=2, range=[4670, 4710],
                     tickfont=tickfont, showbackground=True)

        zaxis = dict(title=label_z, showgrid=True, showticklabels=True,
                     gridcolor="white", nticks=4, range=[-15, -1],
                     tickfont=tickfont, showbackground=True)

        margin = dict(r=10, l=20, b=20, t=20)

        scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
                     aspectmode='cube')

        title = dict(text="<b>Fault Triangulation</b>", font_family="Arial Black",
                     font_color="black", font_size=24, x=0.5, y=0.97)

        layout = go.Layout(scene=scene, margin=margin, width=1600,
                           height=1200, title=title)

        fig = ff.create_trisurf(x=self.xf, y=self.yf, z=self.zf,
                                colormap=[(1.0, 1.0, 0.6), (0.95, 0.95, 0.6)],
                                simplices=self.tri, plot_edges=True,
                                show_colorbar=False)

        fig.update_layout(layout)
        fig.update(layout_coloraxis_showscale=False)
        fig.show(renderer="browser")

    # *************************************************************************
    # *                                                                       *
    # *         Triangulation methods section                                *
    # *                                                                       *
    # *************************************************************************

    def triangulate_fault(self):
        """
        Triangulate square subfaults
        :param:
        :return:
        """
        print()
        print(" Fault Triangulation")

        index_fault = np.arange(0, self.xf.size).reshape((self.ndip, self.nstk),
                                                         order='F')

        # Calculate number of triangles
        self.n_tri = (self.nstk - 1) * (self.ndip - 1) * 2

        # Initiate triangulation array
        self.tri = np.zeros([self.n_tri, 3], dtype=int)

        # Build a Delaunay triangulation
        jtri = -1
        for istk in range(0, self.nstk - 1):
            for idip in range(0, self.ndip - 1):
                jtri += 1
                self.tri[jtri, 0] = index_fault[idip, istk]
                self.tri[jtri, 1] = index_fault[idip, istk + 1]
                self.tri[jtri, 2] = index_fault[idip + 1, istk + 1]
                jtri += 1
                self.tri[jtri, 0] = index_fault[idip, istk]
                self.tri[jtri, 1] = index_fault[idip + 1, istk + 1]
                self.tri[jtri, 2] = index_fault[idip + 1, istk]

        print(f" Number of facets: {self.n_tri}")
        self.trib_marker = np.ones(self.n_tri, )

        # Calculate unitary normal, strike and dip vector at each facet
        self.facet_norm = np.zeros((self.n_tri, 9))

        # Vector normal to earth surface
        n_surf = np.array([0, 0, -1])

        for itri in range(0, self.n_tri):
            iv0 = self.tri[itri, 0]
            iv1 = self.tri[itri, 1]
            iv2 = self.tri[itri, 2]

            v0 = np.array([self.xf[iv0], self.yf[iv0], self.zf[iv0]])
            v1 = np.array([self.xf[iv1], self.yf[iv1], self.zf[iv1]])
            v2 = np.array([self.xf[iv2], self.yf[iv2], self.zf[iv2]])

            vec_normal = np.cross(v1 - v0, v2 - v0)
            vec_normal = vec_normal / np.linalg.norm(vec_normal)
            vec_strike = np.cross(vec_normal, n_surf)
            vec_strike = vec_strike / np.linalg.norm(vec_strike)
            vec_dip = np.cross(vec_strike, vec_normal)
            vec_dip = vec_dip / np.linalg.norm(vec_dip)

            self.facet_norm[itri, 0:3] = vec_normal
            self.facet_norm[itri, 3:6] = vec_strike
            self.facet_norm[itri, 6:9] = vec_dip

        # Create nodes above and below to the fault
        x_above = self.xf + vec_normal[0] * self.dh_f
        y_above = self.yf + vec_normal[1] * self.dh_f
        z_above = self.zf + vec_normal[2] * self.dh_f
        x_below = self.xf - vec_normal[0] * self.dh_f
        y_below = self.yf - vec_normal[1] * self.dh_f
        z_below = self.zf - vec_normal[2] * self.dh_f

        print()
        print(f" {len(x_above)} nodes added above the fault ")
        print(f" {len(x_below)} nodes added below the fault ")
        print()

        self.xf_add = np.concatenate((x_above, x_below), axis=None)
        self.yf_add = np.concatenate((y_above, y_below), axis=None)
        self.zf_add = np.concatenate((z_above, z_below), axis=None)

    # *************************************************************************
    # *                                                                       *
    # *                     save files methods section                        *
    # *                                                                       *
    # *************************************************************************
    def save_fault(self, dir_):
        print()
        out_file = dir_ + self.name + ".pickle"
        object_file = open(out_file, 'wb')
        pickle.dump(self, object_file)
        object_file.close()
        print(f" Fault object saved in: {out_file} ")

    def write_uni_vector(self, dir_):
        print()
        # Write .vector file
        f_vector_header = f"{self.n_tri:d}"
        f_vector = dir_ + self.name + ".vector"
        with open(f_vector, 'wb') as f:
            np.savetxt(f, self.facet_norm, header=f_vector_header,
                       comments=' ', fmt='%10.6f')
        f.close()

        print(f" vector file saved in: {f_vector} ")

    def write_fcoor(self, dir_, fault_id, nt, dt):
        print()
        # Write fcoor file
        fcoorHeader = "%d  %d %4.2f " % (self.nstk * self.ndip, nt, dt)

        fcoorName = dir_ + self.name + "_ID_" + fault_id + ".in"
        fcoor = np.array([self.xf, self.yf, self.zf])

        with open(fcoorName, 'wb') as f:
            np.savetxt(f, fcoor, header=fcoorHeader, comments=' ', fmt='%9.4f')

        print(f" Coordinates file saved in: {fcoorName} ")
        print(f" Number of subfaults: {self.nstk * self.ndip} ")
        print(f" Number of time steps: {nt} ")
        print(f" time step: {dt} ")

        f.close()
