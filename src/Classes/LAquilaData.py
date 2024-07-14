"""
 Fault Class
"""
import pickle

import numpy as np
from os import scandir
from os import remove
import obspy
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px


class Data:
    def __init__(self, name: str):
        self.name = name
        self.dg_data = {}
        self.src = None

    def set_dg_data(self, stations, data_folder):
        for stat in stations.keys():
            self.dg_data[stat] = {}
            self.dg_data[stat]['trim_t'] = stations[stat]

            dg_data = self.extract_data(data_folder, stat)
            self.dg_data[stat]['ori_vx'] = dg_data['ori_vx']
            self.dg_data[stat]['ori_vy'] = dg_data['ori_vy']
            self.dg_data[stat]['ori_vz'] = dg_data['ori_vz']
            self.dg_data[stat]['ori_time'] = dg_data['ori_time']
            self.dg_data[stat]['vx'] = dg_data['vx']
            self.dg_data[stat]['vy'] = dg_data['vy']
            self.dg_data[stat]['vz'] = dg_data['vz']
            self.dg_data[stat]['time'] = dg_data['time']

    def extract_data(self, data_folder, stat):
        # Scan for simulation obspy files
        velo = [f for f in scandir(data_folder) if f.is_file()
                and f.name.startswith(self.name)
                and f.name.__contains__(stat)]

        # print(velo[0].name)
        # for velo in velo_files:
        #     # Load synthetic data for each station
        #     stat = velo.name.split('_')[-1].split('.')[0]
        print(f" station: {stat} ")
        st = obspy.read(velo[0].path)
        ti = st.traces[0].stats.starttime
        tf = st.traces[0].stats.endtime
        # st.plot()

        # Save original data
        ori_vx = st.traces[0].data
        ori_vy = st.traces[1].data
        ori_vz = st.traces[2].data
        ori_time = st.traces[2].times()
        print(ori_vx[:5])

        ti = st.traces[0].stats.starttime
        tf = st.traces[0].stats.endtime
        st.trim(ti + self.dg_data[stat]['trim_t'], tf)
        st.normalize()

        # Save normalized and time trim data
        vx = st.traces[0].data
        vy = st.traces[1].data
        vz = st.traces[2].data
        time = st.traces[2].times()

        # Create and add a dict with each station synthetic data
        dg_data = {'ori_vx': ori_vx, 'ori_vy': ori_vy,
                   'ori_vz': ori_vz, 'ori_time': ori_time,
                   'vx': vx, 'vy': vy, 'vz': vz, 'time': time}
        return dg_data

    def extract_src_fields(self, sim_path):
        with open(sim_path, 'rb') as f:
            model_src = pickle.load(f)

        # Load and add a dict with src fields
        self.src = model_src

    def plot_simulation(self):
        stk = np.arange(self.src['n_stk'])
        dip = np.arange(self.src['n_dip'])

        tickfont = dict(color="black", size=12, family="Arial Black")
        xaxis = dict(title="<b> strike (Km) </b>", tickfont=tickfont)
        yaxis = dict(title="<b> dip (Km) </b>", tickfont=tickfont,
                     autorange="reversed")
        title_a = dict(text="<b> Interpolated Slip </b>",
                       font_family="Arial Black",
                       font_color="black", x=0.5, y=0.85)
        colorbar_a = dict(lenmode='fraction', len=0.55, thickness=20,
                          bordercolor="black", title="<b> slip (m) </b>",
                          x=1.0, y=0.8)
        data_a = go.Heatmap(z=self.src['SLIP'], x=stk, y=dip, hoverongaps=False,
                            colorbar=colorbar_a)
        colorscale = [[0, 'rgb(250, 250, 250)'], [1.0, 'rgb(255, 255, 255)']]
        contours = dict(coloring='lines', showlabels=True)
        data_b = go.Contour(z=self.src['RUPT_T'], x=stk, y=dip, contours=contours,
                            line_width=2, showscale=False, colorscale=colorscale)
        colorbar_c = dict(lenmode='fraction', len=0.55, thickness=20,
                          bordercolor="black", title="<b> rise time (s) </b>",
                          x=1.0, y=0.25)
        data_c = go.Heatmap(z=self.src['RISE'], x=stk, y=dip, hoverongaps=False,
                            colorbar=colorbar_c)

        fig_A = make_subplots(rows=2, cols=1,
                              subplot_titles=('slip with rupture time',
                                              'rise time'),
                              shared_xaxes=True, vertical_spacing=0.1)

        fig_A.add_trace(data_a, row=1, col=1)
        fig_A.add_trace(data_b, row=1, col=1)
        fig_A.add_trace(data_c, row=2, col=1)
        fig_A.update_yaxes(row=1, col=1, autorange="reversed", title='dip',
                           tickfont=tickfont, scaleanchor="x", scaleratio=1)
        fig_A.update_yaxes(row=2, col=1, autorange="reversed", title='dip',
                           tickfont=tickfont, scaleanchor="x", scaleratio=1)
        fig_A.update_xaxes(row=2, col=1, title='strike', tickfont=tickfont,
                           range=[np.min(stk), np.max(stk)])

        fig_A.show(renderer="browser")

        fig_B = go.Figure()
        for istat, stat in enumerate(self.dg_data.keys()):
            fig_B.add_trace(go.Scatter(x=self.dg_data[stat]['time'],
                                       y=self.dg_data[stat]['vx']+2.2*istat,
                                       mode='lines'))

        fig_B.show(renderer="browser")

    def save_simulation_data(self, output_folder):
        out_file = output_folder + self.name + ".pickle"
        object_file = open(out_file, 'wb')
        pickle.dump(self, object_file)
        object_file.close()
        print(f" Saving simulation data in: {out_file}")
