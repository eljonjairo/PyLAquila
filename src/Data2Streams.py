#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Generate Obspy Stream per Station with Seismic Data
#
# John Diaz December 2022

import numpy as np
from os import scandir
import obspy
from obspy.core import UTCDateTime
from Classes.LAquilaData import Data

def file2trace(in_file: str):
#    def file2trace(in_file: str, tfinal: float):

    with open(in_file, 'r') as f:
        lines = f.readlines()

    tr = obspy.Trace()
    stime = lines[26].split(':')[1].split('_')
    time = stime[0]+'T'+stime[1]
    tr.stats.starttime = UTCDateTime(time)
    tr.stats.delta = float(lines[28].split(':')[1])
    tr.stats.npts = int(lines[29].split(':')[1])
    tr.data = np.array(lines[64:-1], dtype=float)/100  # Velocity [m/s]

    return tr


if __name__ == '__main__':
    # Data folder
    data_folder = "../Data/ASCII/"

    # Object Folder
    output_folder = "../Data/"

    stations = ['AQK', 'AQU', 'AQV', 'AQA', 'AQG', 'GSA', 'MTR', 'FMG',
                'ANT', 'AVZ', 'CSO1', 'LSS', 'SUL']

    # 1.0 Hz simulations
    nt = 816
    dt = 0.0490204

    # Cutoff freqs
    freqmin = 0.02  # Cirella
    freqmax = 1.0

    # Event time
    event_time = "2009-04-06T01:32:40"

    # Add 40 seconds of simulation
    final_time = "2009-04-06T01:33:20"

    # Field DISP: displacement, VEL: velocity, ACC: acceleration
    field = "VEL"

    nstats = len(stations)
    time = np.linspace(0, dt*nt, nt)

    real_data = Data('real_data', nstats, nt, dt, freqmin, freqmax)
    real_data.time = time
    real_data.data = {}

    velo_x = np.zeros((nstats, nt))
    velo_y = np.zeros((nstats, nt))
    velo_z = np.zeros((nstats, nt))

    for istat, stat in enumerate(stations):
        # Real Data Processing
        print(f" Station {stat} ")

        # East component
        data_file = [f for f in scandir(data_folder + stat + '/') if f.is_file()
                     and f.name.__contains__(field)
                     and (f.name.__contains__("HNE")
                     or f.name.__contains__("HLE"))]

        print(f" East component file: {data_file[0].name}")
        trace = file2trace(data_file[0].path)
        ti = UTCDateTime(event_time)
        tf = UTCDateTime(final_time)
        trace_vx = trace.slice(ti, tf)
        trace_vx = trace_vx.resample(sampling_rate=1/dt)
        trace_vx = trace_vx.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
        # trace.plot()
        if trace_vx.stats.npts < nt:
            trace_vx = np.pad(trace_vx.data, (0, nt-trace_vx.stats.npts),
                              'constant', constant_values=(0, 0))

        velo_x[istat, :] = trace_vx.data

        # North component
        data_file = [f for f in scandir(data_folder + stat + '/') if f.is_file()
                     and f.name.__contains__(field)
                     and (f.name.__contains__("HNN")
                          or f.name.__contains__("HLN"))]

        print(f" North component file: {data_file[0].name}")
        trace = file2trace(data_file[0].path)
        trace_vy = trace.slice(ti, tf)
        trace_vy = trace_vy.resample(sampling_rate=1/dt)
        trace_vy = trace_vy.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
        # trace.plot()
        if trace_vy.stats.npts < nt:
            trace_vy.data = np.pad(trace_vy.data, (0, nt-trace_vy.stats.npts),
                                   'constant', constant_values=(0, 0))
        velo_y[istat, :] = trace_vy.data

        data_file = [f for f in scandir(data_folder + stat + '/') if f.is_file()
                     and f.name.__contains__(field)
                     and (f.name.__contains__("HNZ")
                          or f.name.__contains__("HLZ"))]

        print(f" Vertical component file: {data_file[0].name}")
        trace = file2trace(data_file[0].path)
        trace_vz = trace.slice(ti, tf)
        trace_vz = trace_vz.resample(sampling_rate=1/dt)
        trace_vz = trace_vz.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
        # trace.plot()
        if trace_vz.stats.npts < nt:
            trace_vz.data = np.pad(trace_vz.data, (0, nt-trace_vz.stats.npts),
                                   'constant', constant_values=(0, 0))
        velo_z[istat, :] = trace_vz.data

    # Set velo traces
    real_data.set_velo_acc(velo_x, velo_y, velo_z, stations)
    real_data.set_velo_spec()
    real_data.set_arias()
    real_data.plot_data()
    real_data.save_simulation_data(output_folder=output_folder)
