#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Generate Obspy Stream per Station with Seismic Data
#
# John Diaz December 2022


import numpy as np
from os import scandir
from pathlib import Path
import obspy
from obspy.core import UTCDateTime
from obspy.core.stream import Stream
from tqdm import tqdm

if __name__ == '__main__':
    # DG folder
    DGFolder = "../DGrun_1.0Hz/"

    # Output Folder
    SynFolder = "Outputs/SyntheticVelo_1.0Hz/"

    # Event time
    teven = "2009-04-06T01:32:40"

    # Stations name in the same order as GEODG3D.acquisition DG file.
    Stats = ["AQK", "AQU", "AQV", "AQA", "AQG", "GSA", "MTR", "FMG", "ANT", "AVZ",
             "CSO1", "LSS", "SUL"]

    nstats = int(len(Stats))
    # 0.5 Hz simulations
    # nt = 816
    # dt = 0.0489741

    # 1.0 Hz simulations
    nt = 816
    dt = 0.0490204

    print("  ")
    print(" START PROGRAM ")
    print("  ")

    ini = 0
    subfolders = [f.path for f in scandir(DGFolder) if f.is_dir()]
    for sub_f in tqdm(subfolders[ini:]):
        file = Path(sub_f)
        file_vx = file / "VX_1"
        if file_vx.exists():
            file_vy = file / "VY_1"
            file_vz = file / "VZ_1"
            print("Loading DG velocity synthetic files:")
            print(f"  {file_vx}")
            print(f"  {file_vy}")
            print(f"  {file_vz}")

            velo_syn_x = np.reshape(np.fromfile(file_vx, dtype=np.float32),
                                    (nstats, nt), order='F')
            velo_syn_y = np.reshape(np.fromfile(file_vy, dtype=np.float32),
                                    (nstats, nt), order='F')
            velo_syn_z = np.reshape(np.fromfile(file_vz, dtype=np.float32),
                                    (nstats, nt), order='F')

            for istat in range(0, nstats):
                # DG Data Processing
                velox_DG = velo_syn_x[istat, :]
                veloy_DG = velo_syn_y[istat, :]
                veloz_DG = velo_syn_z[istat, :]

                # Check if there is nan values
                if not np.isnan([velox_DG, veloy_DG, veloz_DG]).any():

                    trace_vx = obspy.Trace()
                    trace_vx.stats.station = Stats[istat]
                    trace_vx.stats.starttime = UTCDateTime(teven)
                    trace_vx.stats.delta = dt
                    trace_vx.stats.npts = nt
                    trace_vx.stats.channel = "DGx"
                    trace_vx.data = velox_DG
                    trace_vx.filter('lowpass', freq=1.0, corners=2, zerophase=True)
                    # trace_vx.plot()

                    trace_vy = obspy.Trace()
                    trace_vy.stats.station = Stats[istat]
                    trace_vy.stats.starttime = UTCDateTime(teven)
                    trace_vy.stats.delta = dt
                    trace_vy.stats.npts = nt
                    trace_vy.stats.channel = "DGy"

                    trace_vy.data = veloy_DG
                    trace_vy.filter('lowpass', freq=1.0, corners=2, zerophase=True)
                    # trace_vy.plot()

                    trace_vz = obspy.Trace()
                    trace_vz.stats.station = Stats[istat]
                    trace_vz.stats.starttime = UTCDateTime(teven)
                    trace_vz.stats.delta = dt
                    trace_vz.stats.npts = nt
                    trace_vz.stats.channel = "DGz"
                    trace_vz.data = veloz_DG
                    trace_vz.filter('lowpass', freq=1.0, corners=2, zerophase=True)
                    # trace_vz.plot()

                    st = Stream([trace_vx, trace_vy, trace_vz])
                    st_file = SynFolder + (sub_f.split('/')[-1] + '_DGVEL_' + Stats[istat]
                                           + '.pickle')

                    print(f" Writing Stream in file: {st_file} ")
                    st.write(st_file, format="PICKLE")

        else:
            print(f" file {file_vx} not found")
