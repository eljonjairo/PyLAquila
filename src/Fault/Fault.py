#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import scipy
#
# Fault class
# 

import utm
import numpy as np
import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots
from scipy.io import loadmat
from scipy.spatial import distance

class Fault:
    def __init__(self, name: str, dh_f: float):

        self.dh_f = dh_f
        self.name = name + "_dhF" + str(int(self.dh_f)) + "m"

        self.in_X = None  # Input X coordinate of the fault
        self.in_Y = None  # Input Y coordinate matriz
        self.in_Z = None  # Input X coordinate of the fault
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
        self.X = None  # Interpolated X coordinate matriz
        self.Y = None  # Interpolated Y coordinate matriz
        self.Z = None  # Interpolated Z coordinate matriz
        self.SLIP = None  # Interpolated slip distribution
        self.RUPT = None  # Interpolated rupture time distribution
        self.RISE = None  # Interpolated rise time distribution
        self.stk_len = None  # Length in the strike direction
        self.dip_len = None  # Length in the dip direction
        self.stk_vec = None  # Strike vector
        self.dip_vec = None  # Dip vector
        self.stk_univec = None  # Strike vector
        self.dip_univec = None  # Dip vector
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

    def __repr__(self):
        return "Fault name: " + str(self.name) + "dhF: " + str(self.dh_f) + " Km"

    def load_mat_file(self, input_mat):
        """
        Load matlab structure file
        :param input_mat: matlab structure file path:
        : return:
        Fault object with input file attributes
        """

        print()
        print(f" Loading matlab file from {input_mat} ")

        # Load Matlab Finite Fault input file
        inFault = loadmat("".join(["Inputs/", input_mat, ".mat"]))

        # Load original fault coordinates, slip, rise time and rupture time
        self.in_Z = -inFault[input_mat]['geoZ'][-1][0]
        self.in_LAT = inFault[input_mat]['geoLAT'][-1][0]
        self.in_LON = inFault[input_mat]['geoLON'][-1][0]
        self.in_SLIP = inFault[input_mat]['slipSPL'][-1][0]
        self.in_RISE = inFault[input_mat]['riseSPL'][-1][0]
        self.in_RUPT = inFault[input_mat]['timeSPL'][-1][0]

        # Hypocenter coordinates (Km) 
        self.hypo_lon = inFault[input_mat]['evLON'][-1][0][0][0]
        self.hypo_lat = inFault[input_mat]['evLAT'][-1][0][0][0]
        self.hypo_z = -inFault[input_mat]['evDPT'][-1][0][0][0]

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
        self.stk_univec = in_dstk_vec / self.in_dstk

        # Calculate magnitude and unitary vector of delta in dip direction
        in_ddip_vec = np.array([self.in_X[1, 0] - self.in_X[0, 0],
                                self.in_Y[1, 0] - self.in_Y[0, 0],
                                self.in_Z[1, 0] - self.in_Z[0, 0]])
        self.in_ddip = np.linalg.norm(in_ddip_vec)
        self.dip_univec = in_ddip_vec / self.in_ddip

        # Calculate length and number of points in strike an dip direction
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

        # Coordinates of the first fault point
        ini_vec = np.array([self.in_X[0, 0], self.in_Y[0, 0], self.in_Z[0, 0]])

        # Calculate output fault's size
        self.nstk = int(self.stk_len / self.dh_f) + 1
        self.ndip = int(self.dip_len / self.dh_f) + 1

        # Create arrays for interpolated coordinates
        self.X = np.zeros((self.ndip, self.nstk))
        self.Y = np.zeros((self.ndip, self.nstk))
        self.Z = np.zeros((self.ndip, self.nstk))

        # Initiate arrays to calculate the index of hypocenter in interpolated
        # XYZ coordinates
        HYPO_DIST = np.zeros((self.ndip, self.nstk))  # Hypocenter distance
        hypo_xy = [self.hypo_x, self.hypo_y]

        # Vectors of size dh_f in strike and dip directions
        dstk_vec = self.stk_univec * self.dh_f
        ddip_vec = self.dip_univec * self.dh_f

        # Effective size of subfaults in strike and dip directions
        self.dstk = np.linalg.norm(dstk_vec)
        self.ddip = np.linalg.norm(ddip_vec)

        for istk in range(0, self.nstk):
            delta_stk = istk * dstk_vec
            for idip in range(0, self.ndip):
                delta_dip = idip * ddip_vec
                self.X[idip, istk] = (self.in_X[0, 0] + delta_stk[0]
                                      + delta_dip[0])
                self.Y[idip, istk] = (self.in_Y[0, 0] + delta_stk[1]
                                      + delta_dip[1])
                self.Z[idip, istk] = (self.in_Z[0, 0] + delta_stk[2]
                                      + delta_dip[2])

                # Calculate hypocenter distance in XY coords.
                XYvec = [self.X[idip, istk], self.Y[idip, istk]]
                HYPO_DIST[idip, istk] = distance.euclidean(XYvec,
                                                           [self.hypo_x, self.hypo_y])
        # Find the indexes of the hypocenter
        self.hypo_idip, self.hypo_istk = np.where(HYPO_DIST == np.min(HYPO_DIST))

        # Make vertical correction ensuring that hypocenter relays on the fault
        self.Z += self.hypo_z - self.Z[self.hypo_idip, self.hypo_istk]

        print(f" Fault's output dimensions: ")
        print(f" Strike (Km): {self.stk_len:.3f} nstk: {self.nstk} "
              f" dstk (Km): {self.dstk:.3f} ")
        print(f" Dip (Km): {self.dip_len:.3f} ndip: {self.ndip} "
              f" ddip (Km): {self.ddip:.3f} ")

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









    # *************************************************************************
    # *                                                                       *
    # *                 Plots methods section                                 *
    # *                                                                       *
    # *************************************************************************
    def plot_fault_inputs(self):
        """ Plot the input slip distribution in lat-lon-z coordinates """

        colorbar_a = dict(lenmode='fraction', len=0.5, thickness=20,
                          bordercolor="black", title="<b> slip (cm) </b>", x=0.1)

        data_a = [go.Surface(x=self.in_LON, y=self.in_LAT, z=self.in_Z,
                             surfacecolor=self.in_SLIP,
                             colorscale=px.colors.sequential.Viridis,
                             colorbar=colorbar_a, showscale=True,
                             lighting=dict(ambient=0.7))]

        tickfont_a = dict(color="black", size=16, family="Arial Black")

        camera_a = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
                        eye=dict(x=-1.7, y=-2.0, z=1.2))

        margin_a = dict(r=30, l=10, b=30, t=20)

        title_a = dict(text="<b>Input Slip </b>", font_family="Arial Black",
                       font_color="black", x=0.5, y=0.85)

        xaxis_a = dict(title="<b>Longitude (°)</b>", showgrid=True,
                       showticklabels=True,
                       gridcolor="black", nticks=4, range=[13.2, 13.7],
                       tickfont=tickfont_a, showbackground=True)

        yaxis_a = dict(title="<b>Latitude (°)</b> ", showgrid=True,
                       showticklabels=True,
                       gridcolor="black", nticks=4, range=[42.1, 42.5],
                       tickfont=tickfont_a, showbackground=True)

        zaxis_a = dict(title="<b>z (Km )</b>", showgrid=True,
                       showticklabels=True, gridcolor="black", nticks=4,
                       range=[-15, 0], tickfont=tickfont_a, showbackground=True)

        scene_a = dict(camera=camera_a, xaxis=xaxis_a, yaxis=yaxis_a,
                       zaxis=zaxis_a, aspectmode='cube')

        layout_a = go.Layout(scene=scene_a, margin=margin_a, width=1200,
                             height=800, title=title_a)
        fig = make_subplots(rows=1, cols=3)
        fig.add_trace(go.Figure(data=data_a, layout=layout_a), row=1, col=1)
        fig.show(renderer="browser")

    # def set_full_fault(self):
    #
    #     # Xutm, Yutm coord Km
    #     self.XFinMat, self.YFinMat, tmp1, tmp2 = \
    #         utm.from_latlon(self.LatFinMat, self.LonFinMat, 33, 'T')
    #     self.XFin3D = self.XFinMat.flatten(order='F').transpose()
    #     self.YFin3D = self.YFinMat.flatten(order='F').transpose()
    #     self.ZFin3D = self.ZFinMat.flatten(order='F').transpose()
    #     self.Slipin3D = self.SlipinMat.flatten(order='F').transpose()
    #     self.RiseTin3D = self.RiseTinMat.flatten(order='F').transpose()
    #     self.RupTin3D = self.RiseTinMat.flatten(order='F').transpose()
    #
    #     # Calculate magnitude and unitary vector of delta in strike direction
    #     vec_dstk = np.array([self.XFinMat[0, 1] - self.XFinMat[0, 0],
    #                          self.YFinMat[0, 1] - self.YFinMat[0, 0],
    #                          self.ZFinMat[0, 1] - self.ZFinMat[0, 0]])
    #     dstk_in = np.linalg.norm(vec_dstk)
    #     self.univec_dstk = vec_dstk / dstk_in
    #
    #     # Calculate magnitude and unitary vector of delta in dip direction
    #     vec_ddip = np.array([self.XFinMat[1, 0] - self.XFinMat[0, 0],
    #                          self.YFinMat[1, 0] - self.YFinMat[0, 0],
    #                          self.ZFinMat[1, 0] - self.ZFinMat[0, 0]])
    #     ddip_in = np.linalg.norm(vec_ddip)
    #     self.univec_ddip = vec_ddip / ddip_in
    #
    #     # Calculate output delta vector in strike and dip directions
    #     self.dstkVec = self.univec_dstk * self.dhF
    #     self.ddipVec = self.univec_ddip * self.dhF
    #     self.dstk = np.linalg.norm(self.dstkVec)
    #     self.ddip = np.linalg.norm(self.ddipVec)
    #
    #     # Calculate output fault's size
    #     self.nstk = int(self.stk_len_in / self.dhF) + 1
    #     self.ndip = int(self.dip_len_in / self.dhF) + 1
    #
    #     # Define input and output strike and dip vectors
    #     self.dipinVec = np.linspace(0, self.dip_len_in, self.ndip_in)
    #     self.stkinVec = np.linspace(0, self.stk_len_in, self.nstk_in)
    #     self.stkVec = np.linspace(0, self.stk_len_in, self.nstk)
    #     self.dipVec = np.linspace(0, self.dip_len_in, self.ndip)
    #
    #     stk_len = round((self.nstk - 1) * self.dstk)
    #     dip_len = round((self.ndip - 1) * self.ddip)
    #
    #     print(" Output Fault Dimensions:")
    #     print(" Strike (Km): %6.2f nstk: %d dstk (Km): %6.2f "
    #           % (stk_len / 1000, self.nstk, self.dstk / 1000))
    #     print(" Dip (Km): %6.2f ndip: %d ddip (Km): %6.2f"
    #           % (dip_len / 1000, self.ndip, self.ddip / 1000))
    #
    # def effective_size(self):
    #
    #     slip = self.SlipinMat
    #
    #     # Cumulative slip in strike and dip direction
    #     sum_slip_stk = np.sum(slip, axis=0)
    #     sum_slip_dip = np.sum(slip, axis=1)
    #
    #     # Calculate the integral, and the effective length following
    #     # Source Scaling Properties from Finite-Fault-Rupture Models (Mai & Beroza 2000)
    #     acorr_stk = np.correlate(sum_slip_stk, sum_slip_stk, mode='full')
    #     x_stk = np.arange(0, len(acorr_stk), 1, dtype=int) * self.dstk_in
    #     Wacf_stk = scipy.integrate.simpson(acorr_stk, x_stk) / max(acorr_stk)
    #
    #     acorr_dip = np.correlate(sum_slip_dip, sum_slip_dip, mode='full')
    #     x_dip = np.arange(0, len(acorr_dip), 1, dtype=int) * self.ddip_in
    #     Wacf_dip = scipy.integrate.simpson(acorr_dip, x_dip) / max(acorr_dip)
    #
    #     print(" Effective length in strike direction: %5.2f Km." % (Wacf_stk))
    #     print(" Effective length in dip direction: %5.2f Km." % (Wacf_dip))
    #
    #     # Calculate effective number of nodes
    #     nstk_eff = round(Wacf_stk / self.dstk_in) + 1
    #     ndip_eff = round(Wacf_dip / self.ddip_in) + 1
    #
    #     # Select the interval with length Leff in the strike dir which
    #     # maximize the slip sum
    #     istk = 0
    #     jstk = nstk_eff
    #     istk_eff = istk
    #     sum_slip_eff = np.sum(sum_slip_stk[istk:jstk])
    #     while (jstk < self.nstk_in):
    #         istk += 1
    #         jstk += 1
    #         sum_slip_tmp = np.sum(sum_slip_stk[istk:jstk])
    #         if (sum_slip_tmp > sum_slip_eff):
    #             istk_eff = istk
    #             sum_slip_eff = sum_slip_tmp
    #
    #     # Select the interval with length Leff in the dip dir which
    #     # maximize the slip sum
    #     idip = 0
    #     jdip = ndip_eff
    #
    #     idip_eff = idip
    #     sum_slip_eff = np.sum(sum_slip_dip[idip:jdip])
    #     while (jdip < self.ndip_in):
    #         idip += 1
    #         jdip += 1
    #         print(jdip)
    #         sum_slip_tmp = np.sum(sum_slip_dip[idip:jdip])
    #         if (sum_slip_tmp > sum_slip_eff):
    #             idip_eff = idip
    #             sum_slip_eff = sum_slip_tmp
    #
    #     slip_eff = slip[idip_eff:idip_eff + ndip_eff, istk_eff:istk_eff + nstk_eff]
    #     print()
    #     print(f" Original slip matrix dimensions: {slip.shape} ")
    #     print(f" Effective slip matrix dimensions: {slip_eff.shape} ")
    #     print()
    #     return np.array([idip_eff, idip_eff + ndip_eff, istk_eff, istk_eff + nstk_eff])
    #
    # def set_effective_fault(self):
    #
    #     # Xutm, Yutm coord Km
    #     self.XFinMat, self.YFinMat, tmp1, tmp2 = \
    #         utm.from_latlon(self.LatFinMat, self.LonFinMat, 33, 'T')
    #
    #     # Calculate effective dimensions
    #
    #     eff_dim = self.effective_size()
    #
    #     # Extracting effective part of the fault.
    #     self.XFinMat = self.XFinMat[eff_dim[0]:eff_dim[1],
    #                    eff_dim[2]:eff_dim[3]]
    #     self.YFinMat = self.YFinMat[eff_dim[0]:eff_dim[1],
    #                    eff_dim[2]:eff_dim[3]]
    #     self.ZFinMat = self.ZFinMat[eff_dim[0]:eff_dim[1],
    #                    eff_dim[2]:eff_dim[3]]
    #     self.SlipinMat = self.SlipinMat[eff_dim[0]:eff_dim[1],
    #                      eff_dim[2]:eff_dim[3]]
    #     self.RiseTinMat = self.RiseTinMat[eff_dim[0]:eff_dim[1],
    #                       eff_dim[2]:eff_dim[3]]
    #     self.RupTinMat = self.RupTinMat[eff_dim[0]:eff_dim[1],
    #                      eff_dim[2]:eff_dim[3]]
    #
    #     self.ndip_in, self.nstk_in = self.XFinMat.shape
    #     # Calculate length and output number of points in strike an dip direction
    #     self.stk_len_in = round((self.nstk_in - 1) * self.dstk_in)
    #     self.dip_len_in = round((self.ndip_in - 1) * self.ddip_in)
    #
    #     self.XFin3D = self.XFinMat.flatten(order='F').transpose()
    #     self.YFin3D = self.YFinMat.flatten(order='F').transpose()
    #     self.ZFin3D = self.ZFinMat.flatten(order='F').transpose()
    #     self.Slipin3D = self.SlipinMat.flatten(order='F').transpose()
    #     self.RiseTin3D = self.RiseTinMat.flatten(order='F').transpose()
    #     self.RupTin3D = self.RiseTinMat.flatten(order='F').transpose()
    #
    #     # Calculate magnitude and unitary vector of delta in strike direction
    #     vec_dstk = np.array([self.XFinMat[0, 1] - self.XFinMat[0, 0],
    #                          self.YFinMat[0, 1] - self.YFinMat[0, 0],
    #                          self.ZFinMat[0, 1] - self.ZFinMat[0, 0]])
    #     dstk_in = np.linalg.norm(vec_dstk)
    #     self.univec_dstk = vec_dstk / dstk_in
    #
    #     # Calculate magnitude and unitary vector of delta in dip direction
    #     vec_ddip = np.array([self.XFinMat[1, 0] - self.XFinMat[0, 0],
    #                          self.YFinMat[1, 0] - self.YFinMat[0, 0],
    #                          self.ZFinMat[1, 0] - self.ZFinMat[0, 0]])
    #     ddip_in = np.linalg.norm(vec_ddip)
    #     self.univec_ddip = vec_ddip / ddip_in
    #
    #     # Calculate output delta vector in strike and dip directions
    #     self.dstkVec = self.univec_dstk * self.dhF
    #     self.ddipVec = self.univec_ddip * self.dhF
    #     self.dstk = np.linalg.norm(self.dstkVec)
    #     self.ddip = np.linalg.norm(self.ddipVec)
    #
    #     # Calculate output fault's size
    #     self.nstk = int(self.stk_len_in / self.dhF) + 1
    #     self.ndip = int(self.dip_len_in / self.dhF) + 1
    #
    #     # Define input and output strike and dip vectors
    #     self.dipinVec = np.linspace(0, self.dip_len_in, self.ndip_in)
    #     self.stkinVec = np.linspace(0, self.stk_len_in, self.nstk_in)
    #
    #     self.stkVec = np.linspace(0, self.stk_len_in, self.nstk)
    #     self.dipVec = np.linspace(0, self.dip_len_in, self.ndip)
    #
    #     stk_len = round((self.nstk - 1) * self.dstk)
    #     dip_len = round((self.ndip - 1) * self.ddip)
    #
    #     print(" Output Fault Dimensions:")
    #     print(" Strike (Km): %6.2f nstk: %d dstk (Km): %6.2f "
    #           % (stk_len / 1000, self.nstk, self.dstk / 1000))
    #     print(" Dip (Km): %6.2f ndip: %d ddip (Km): %6.2f"
    #           % (dip_len / 1000, self.ndip, self.ddip / 1000))
    #
    #
    #
    #
    # # *************************************************************************
    # # *                                                                       *
    # # *         Triangulatiopn methods section                                *
    # # *                                                                       *
    # # *************************************************************************
    # def triangulate_fault(self):
    #     print()
    #     print(" Fault Triangulation")
    #     index_fault = np.arange(0, self.XF3D.size).reshape((self.ndip, self.nstk)
    #                                                        , order='F')
    #     ntri = (self.nstk - 1) * (self.ndip - 1) * 2
    #     self.tri = np.zeros([ntri, 3], dtype=int)
    #     #xy_3D = np.array((self.XF3D,self.YF3D)).transpose()
    #
    #     # Delaunay triangulation
    #     # tri = Delaunay(xy_3D).simplices
    #     self.ntri = int(self.tri.size / 3)
    #     jtri = -1
    #     for istk in range(0, self.nstk - 1):
    #         for idip in range(0, self.ndip - 1):
    #             jtri += 1
    #             self.tri[jtri, 0] = index_fault[idip, istk]
    #             self.tri[jtri, 1] = index_fault[idip, istk + 1]
    #             self.tri[jtri, 2] = index_fault[idip + 1, istk + 1]
    #             jtri += 1
    #             self.tri[jtri, 0] = index_fault[idip, istk]
    #             self.tri[jtri, 1] = index_fault[idip + 1, istk + 1]
    #             self.tri[jtri, 2] = index_fault[idip + 1, istk]
    #
    #     print(f" Number of facets: {ntri}")
    #     self.trib_marker = np.ones(self.ntri, )
    #     # Calculate unitary normal, strike and dip vector at each facet
    #     self.univector = np.zeros((self.ntri, 9))
    #     # Vector normal to earth surface
    #     nsurf = np.array([0, 0, -1])
    #
    #     for itri in range(0, self.ntri):
    #         iv0 = self.tri[itri, 0]
    #         iv1 = self.tri[itri, 1]
    #         iv2 = self.tri[itri, 2]
    #         v0 = np.array([self.XF3D[iv0], self.YF3D[iv0], self.ZF3D[iv0]])
    #         v1 = np.array([self.XF3D[iv1], self.YF3D[iv1], self.ZF3D[iv1]])
    #         v2 = np.array([self.XF3D[iv2], self.YF3D[iv2], self.ZF3D[iv2]])
    #         self.vec_normal = np.cross(v1 - v0, v2 - v0)
    #         self.vec_normal = self.vec_normal / np.linalg.norm(self.vec_normal)
    #         self.vec_strike = np.cross(self.vec_normal, nsurf)
    #         self.vec_strike = self.vec_strike / np.linalg.norm(self.vec_strike)
    #         self.vec_dip = np.cross(self.vec_strike, self.vec_normal)
    #         self.vec_dip = self.vec_dip / np.linalg.norm(self.vec_dip)
    #
    #         self.univector[itri, 0:3] = self.vec_normal
    #         self.univector[itri, 3:6] = self.vec_strike
    #         self.univector[itri, 6:9] = self.vec_dip
    #
    # def add_nodes_above_below(self):
    #     # Add nodes above an below to the fault
    #     x_above = self.XF3D + self.vec_normal[0] * self.dhF
    #     y_above = self.YF3D + self.vec_normal[1] * self.dhF
    #     z_above = self.ZF3D + self.vec_normal[2] * self.dhF
    #     x_below = self.XF3D - self.vec_normal[0] * self.dhF
    #     y_below = self.YF3D - self.vec_normal[1] * self.dhF
    #     z_below = self.ZF3D - self.vec_normal[2] * self.dhF
    #     print()
    #     print(f" {len(x_above)} nodes added above the fault ")
    #     print(f" {len(x_below)} nodes added below the fault ")
    #     print()
    #     self.XF3D_add = np.concatenate((x_above, x_below), axis=None)
    #     self.YF3D_add = np.concatenate((y_above, y_below), axis=None)
    #     self.ZF3D_add = np.concatenate((z_above, z_below), axis=None)
    #
    # def plot_xyz_slipin(self):
    #
    #     colorbar = dict(lenmode='fraction', len=0.75, thickness=20,
    #                     bordercolor="black", title="<b> slip(m) </b>", x=0.2)
    #
    #     data = [go.Surface(x=self.XFinMat / 1000, y=self.YFinMat / 1000,
    #                        z=self.ZFinMat / 1000, surfacecolor=self.SlipinMat,
    #                        colorscale=px.colors.sequential.Viridis,
    #                        colorbar=colorbar, showscale=True,
    #                        lighting=dict(ambient=0.7))]
    #
    #     tickfont = dict(color="black", size=16, family="Arial Black")
    #
    #     camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
    #                   eye=dict(x=-1.5, y=-2.0, z=1.2))
    #
    #     xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True, showticklabels=True,
    #                  gridcolor="black", nticks=6, range=[360, 385],
    #                  tickfont=tickfont, showbackground=True)
    #     yaxis = dict(title="<b> yUTM (Km) </b> ", showgrid=True, showticklabels=True,
    #                  gridcolor="black", nticks=6, range=[4670, 4720],
    #                  tickfont=tickfont, showbackground=True)
    #     zaxis = dict(title="<b> z (Km ) </b>", showgrid=True, showticklabels=True,
    #                  gridcolor="black", nticks=6, range=[-18, 2], tickfont=tickfont,
    #                  showbackground=True)
    #
    #     margin = dict(r=30, l=10, b=30, t=20)
    #
    #     title = dict(text="<b>Input Slip </b>", font_family="Arial Blak",
    #                  font_color="black", x=0.5, y=0.85)
    #     scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
    #                  aspectmode='cube')
    #
    #     layout = go.Layout(scene=scene, margin=margin, width=800, height=350,
    #                        title=title)
    #
    #     fig = go.Figure(data=data, layout=layout)
    #     fig.show()
    #
    # def plot_xyz_slip(self):
    #
    #     colorbar = dict(lenmode='fraction', len=0.75, thickness=20, bordercolor="black",
    #                     title="<b> slip(m) </b>", x=0.2)
    #
    #     data = [go.Surface(x=self.XFMat / 1000, y=self.YFMat / 1000, z=self.ZFMat / 1000,
    #                        surfacecolor=self.SlipMat,
    #                        colorscale=px.colors.sequential.Viridis, colorbar=colorbar, showscale=True,
    #                        lighting=dict(ambient=0.7))]
    #
    #     tickfont = dict(color="black", size=16, family="Arial Black")
    #
    #     camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
    #                   eye=dict(x=-1.5, y=-2.0, z=1.2))
    #
    #     xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True, showticklabels=True,
    #                  gridcolor="black", nticks=6, range=[360, 385],
    #                  tickfont=tickfont, showbackground=True)
    #     yaxis = dict(title="<b> yUTM (Km) </b> ", showgrid=True, showticklabels=True,
    #                  gridcolor="black", nticks=6, range=[4670, 4720],
    #                  tickfont=tickfont, showbackground=True)
    #     zaxis = dict(title="<b> z (Km ) </b>", showgrid=True, showticklabels=True,
    #                  gridcolor="black", nticks=6, range=[-18, 2], tickfont=tickfont,
    #                  showbackground=True)
    #
    #     margin = dict(r=30, l=10, b=30, t=20)
    #
    #     title = dict(text="<b>Interpolated Slip </b>", font_family="Arial Blak",
    #                  font_color="black", x=0.5, y=0.85)
    #
    #     scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
    #                  aspectmode='cube')
    #
    #     layout = go.Layout(scene=scene, margin=margin, width=800, height=350,
    #                        title=title)
    #
    #     fig = go.Figure(data=data, layout=layout)
    #     fig.show()
    #
    # def plot_xyz_model_slip(self):
    #
    #     colorbar = dict(lenmode='fraction', len=0.75, thickness=20, bordercolor="black",
    #                     title="<b> slip(m) </b>", x=0.2)
    #
    #     data = [go.Surface(x=self.XFMat / 1000, y=self.YFMat / 1000, z=self.ZFMat / 1000,
    #                        surfacecolor=self.SlipMat,
    #                        colorscale="Viridis", colorbar=colorbar, showscale=True,
    #                        lighting=dict(ambient=0.7))]
    #
    #     tickfont = dict(color="black", size=16, family="Arial Black")
    #
    #     camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
    #                   eye=dict(x=-0.5, y=-2.0, z=1.5))
    #
    #     xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True, showticklabels=True,
    #                  gridcolor="black", nticks=10, range=[300, 420],
    #                  tickfont=tickfont, showbackground=True)
    #     yaxis = dict(title="<b> yUTM (Km) </b> ", showgrid=True, showticklabels=True,
    #                  gridcolor="black", nticks=10, range=[4610, 4750],
    #                  tickfont=tickfont, showbackground=True)
    #     zaxis = dict(title="<b> z (Km ) </b>", showgrid=True, showticklabels=True,
    #                  gridcolor="black", nticks=6, range=[-60, 2], tickfont=tickfont,
    #                  showbackground=True)
    #
    #     margin = dict(r=30, l=50, b=20, t=20)
    #
    #     title = dict(text="<b>Interpolated Slip </b>", font_family="Arial Blak",
    #                  font_color="black", x=0.5, y=0.85)
    #
    #     scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
    #                  aspectmode='cube')
    #
    #     layout = go.Layout(scene=scene, margin=margin, width=800, height=350,
    #                        title=title)
    #
    #     fig = go.Figure(data=data, layout=layout)
    #     fig.show()
    #
    # def compare_xyz_slip(self):
    #
    #     colorbar = dict(lenmode='fraction', len=0.9, thickness=10,
    #                     bordercolor="black", title="<b> slip(m) </b>")
    #
    #     data_a = go.Surface(x=self.XFinMat / 1000, y=self.YFinMat / 1000, z=self.ZFinMat / 1000,
    #                         surfacecolor=self.SlipinMat,
    #                         colorbar=colorbar,
    #                         showscale=False, lighting=dict(ambient=0.9))
    #
    #     data_b = go.Surface(x=self.XFMat / 1000, y=self.YFMat / 1000, z=self.ZFMat / 1000,
    #                         surfacecolor=self.SlipMat, colorbar=colorbar,
    #                         showscale=True, lighting=dict(ambient=0.9))
    #
    #     tickfont = dict(color="black", size=12, family="Arial Black")
    #
    #     camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
    #                   eye=dict(x=-1.5, y=-2.0, z=1.2))
    #
    #     xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True, gridcolor="white",
    #                  showticklabels=True, nticks=6, range=[360, 385],
    #                  tickfont=tickfont, showbackground=True)
    #     yaxis = dict(title="<b> yUTM (Km) </b>", showgrid=True, showticklabels=True,
    #                  gridcolor="white", nticks=6, range=[4670, 4720],
    #                  tickfont=tickfont, showbackground=True)
    #     zaxis = dict(title="<b> z (Km) </b>", showgrid=True, showticklabels=True,
    #                  gridcolor="white", nticks=6, range=[-18, 2], tickfont=tickfont,
    #                  showbackground=True)
    #     margin = dict(r=5, l=5, b=10, t=20)
    #     scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
    #                  aspectmode='cube')
    #
    #     fig = tls.make_subplots(rows=1, cols=2, specs=[[{'is_3d': True},
    #                                                     {'is_3d': True}]],
    #                             horizontal_spacing=0,
    #                             subplot_titles=["<b>Input Slip </b>",
    #                                             "<b>Interpolated Slip </b>"])
    #     fig.append_trace(data_a, 1, 1)
    #     fig.append_trace(data_b, 1, 2)
    #     fig.update_scenes(scene, row=1, col=2)
    #     fig.update_scenes(scene, row=1, col=1)
    #     fig.update_layout(height=310, width=800, margin=margin)
    #     fig.show()
    #
    # def plot_triangulation(self):
    #     tickfont = dict(color="black", size=19, family="Arial Black")
    #
    #     camera = dict(up=dict(x=0, y=0, z=1), center=dict(x=0, y=0, z=0),
    #                   eye=dict(x=-1.0, y=-1.0, z=0.7))
    #     label_z = dict(text="<b>z (Km)</b>", font_family="Arial Blak",
    #                    font_color="black", font_size=26)
    #     xaxis = dict(title="<b> xUTM (Km) </b>", showgrid=True, gridcolor="white",
    #                  showticklabels=True, nticks=2, range=[360, 390],
    #                  tickfont=tickfont, showbackground=True)
    #     yaxis = dict(title="<b> yUTM (Km) </b>", showgrid=True, showticklabels=True,
    #                  gridcolor="white", nticks=2, range=[4670, 4710],
    #                  tickfont=tickfont, showbackground=True)
    #     zaxis = dict(title=label_z, showgrid=True, showticklabels=True,
    #                  gridcolor="white", nticks=4, range=[-15, -1], tickfont=tickfont,
    #                  showbackground=True)
    #     margin = dict(r=10, l=20, b=20, t=20)
    #
    #     scene = dict(camera=camera, xaxis=xaxis, yaxis=yaxis, zaxis=zaxis,
    #                  aspectmode='cube')
    #
    #     title = dict(text="<b>Fault Triangulation</b>", font_family="Arial Blak",
    #                  font_color="black", font_size=24, x=0.5, y=0.97)
    #
    #     layout = go.Layout(scene=scene, margin=margin, width=1600,
    #                        height=1200, title=title)
    #
    #     fig = ff.create_trisurf(x=self.XF3D / 1000, y=self.YF3D / 1000, z=self.ZF3D / 1000,
    #                             colormap=[(1.0, 1.0, 0.6), (0.95, 0.95, 0.6)],
    #                             simplices=self.tri, plot_edges=True)
    #     fig.update_layout(layout)
    #     fig.update(layout_coloraxis_showscale=False)
    #
    #     fig.show()
    #
    #     # *************************************************************************
    #
    # # *                                                                       *
    # # *                     save files methods section                        *
    # # *                                                                       *
    # # *************************************************************************
    # def save_fault(self, dir):
    #     print()
    #     out_file = dir + self.name + ".pickle"
    #     object_file = open(out_file, 'wb')
    #     pickle.dump(self, object_file)
    #     object_file.close()
    #     print(f" Fault object saved in: {out_file} ")
    #
    # def write_univector(self, dir):
    #     print()
    #     # Write .vector file
    #     fvector_header = "%d" % (self.ntri)
    #     fvector = dir + self.name + ".vector"
    #     with open(fvector, 'wb') as f:
    #         np.savetxt(f, self.univector, header=fvector_header,
    #                    comments=' ', fmt='%10.6f')
    #     f.close()
    #
    #     print(f" vector file saved in: {fvector} ")
    #
    # def write_fcoor(self, dir, id, nt, dt):
    #     print()
    #     # Write fcoor file
    #     fcoorHeader = "%d  %d %4.2f " % (self.n_sub_faults, nt, dt)
    #
    #     fcoorName = dir + self.name + "_ID_" + id + ".in"
    #
    #     with open(fcoorName, 'wb') as f:
    #         np.savetxt(f, self.fcoor, header=fcoorHeader, comments=' ', fmt='%9.4f')
    #
    #     print(f" Coordinates file saved in: {fcoorName} ")
    #     print(f" Number of subfaults: {self.n_sub_faults} ")
    #     print(f" Number of time steps: {nt} ")
    #     print(f" time step: {dt} ")
    #
    #     f.close()
