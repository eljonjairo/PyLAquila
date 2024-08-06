#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Save original source fields (slip, rise time, rupture time) into pickle file
#
# John Diaz july 2024

from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.io import loadmat
from os import scandir
import pickle
import numpy as np

if __name__ == '__main__':
    # Input folder
    input_folder = "Inputs/"
    files = [f.path for f in scandir(input_folder) if f.is_file()
             and f.name.endswith('GALL.mat')]
    dstk_out = 1.0
    nstk_out = 25
    ddip_out = 1.0
    ndip_out = 15
    for file in files:
        f_name = file.split('.')[0].split('/')[-1]
        input_source = loadmat(file)
        # Load original src fields
        slip = input_source[f_name]['slipSPL'][-1][0]
        rise = input_source[f_name]['riseSPL'][-1][0]
        rupt_t = input_source[f_name]['timeSPL'][-1][0]
        dhF = input_source[f_name]['invDzDx'][-1][0]
        ddip = dhF[0][0]
        dstk = dhF[0][1]
        # Original strike and dip vectors
        n_dip, n_stk = np.shape(slip)
        ori_stk = np.arange(0, n_stk*dstk, dstk)
        ori_dip = np.arange(0, n_dip*ddip, ddip)
        # Output strike and dip vectors
        out_stk = np.arange(0, nstk_out*dstk_out, dstk_out)
        out_dip = np.arange(0, ndip_out*ddip_out, ddip_out)
        X, Y = np.meshgrid(out_dip, out_stk, indexing='ij')

        # Interpolate original src fields
        SlipF = RGI((ori_dip, ori_stk), slip, bounds_error=False,
                    fill_value=None)
        slip = SlipF((X, Y))
        RiseF = RGI((ori_dip, ori_stk), rise, bounds_error=False,
                    fill_value=None)
        rise = RiseF((X, Y))
        RuptF = RGI((ori_dip, ori_stk), rupt_t, bounds_error=False,
                    fill_value=None)
        rupt_t = RuptF((X, Y) )

        slip_fact = np.ceil(np.max(slip))
        rise_fact = np.ceil(np.max(rise))
        rupt_t_fact = np.ceil(np.max(rupt_t))
        slip_norm = slip / slip_fact
        rise_norm = rise / rise_fact
        rupt_t_norm = rupt_t / rupt_t_fact
        n_dip, n_stk = np.shape(slip)

        source = {'name': f_name, 'n_stk': n_stk, 'n_dip': n_dip,
                  'SLIP': slip, 'slip_fact': slip_fact, 'SLIP_NORM': slip_norm,
                  'RISE': rise, 'rise_fact': rise_fact, 'RISE_NORM': rise_norm,
                  'RUPT_T': rupt_t, 'rupt_t_fact': rupt_t_fact,
                  'RUPT_T_NORM': rupt_t_norm}

        out_file = input_folder + f_name + ".pickle"
        object_file = open(out_file, 'wb')
        pickle.dump(source, object_file)
        object_file.close()
        print(f" source fields saved in: {out_file} ")