#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# Save each test fields (slip, rise time, rupture time) into pickle file
#
# John Diaz july 2024

from scipy.io import loadmat
from os import scandir
import pickle
import numpy as np
from tqdm import tqdm

if __name__ == '__main__':
    # Input folder
    input_folder = "/home/jon/PhD/codes/LAquila_2023/LAquila/Source/SourceGen/srate_files/Sources_mat/"

    # Output folder
    output_folder = "Outputs/Fields/"

    model_id = "LAquilaCirella2012bdhF500m_ID_11"

    files = [f.path for f in scandir(input_folder) if f.is_file()
             and f.name.startswith(model_id)]

    for file in tqdm(files):
        f_name = file.split('/')[-1]
        f_model = "LAquila_ID" + f_name.split('_')[-1].split('.')[0]
        print(f" Loading file: {f_name}")
        input_source = loadmat(file)
        slip = input_source['SourceParam']['slip'][-1][0]
        rise = input_source['SourceParam']['riseT'][-1][0]
        rupt_t = input_source['SourceParam']['rupTime'][-1][0]
        rupt_v = input_source['SourceParam']['Vr'][-1][0]
        dstk = input_source['SourceParam']['dstk_STF'][-1][0]
        ddip = input_source['SourceParam']['ddip_STF'][-1][0]
        slip_fact = np.ceil(np.max(slip))
        rise_fact = np.ceil(np.max(rise))
        rupt_t_fact = np.ceil(np.max(rupt_t))
        rupt_v_fact = np.ceil(np.max(rupt_v))
        slip_norm = slip / slip_fact
        rise_norm = rise / rise_fact
        rupt_t_norm = rupt_t / rupt_t_fact
        rupt_v_norm = rupt_v / rupt_v_fact
        n_dip, n_stk = np.shape(slip)

        source = {'name': f_model,
                  'n_stk': n_stk, 'n_dip': n_dip, 'dstk': dstk, 'ddip': ddip,
                  'SLIP': slip, 'slip_fact': slip_fact, 'SLIP_NORM': slip_norm,
                  'RISE': rise, 'rise_fact': rise_fact, 'RISE_NORM': rise_norm,
                  'RUPT_T': rupt_t, 'rupt_t_fact': rupt_t_fact,
                  'RUPT_T_NORM': rupt_t_norm,
                  'RUPT_V': rupt_v, 'rupt_v_fact': rupt_v_fact,
                  'RUPT_V_NORM': rupt_v_norm}

        out_file = output_folder + f_model + ".pickle"
        object_file = open(out_file, 'wb')
        pickle.dump(source, object_file)
        object_file.close()
        print(f" source fields saved in: {out_file} ")
