"""
 Fault Class
"""
import pickle

import numpy as np
from os import scandir
import obspy
from plotly.subplots import make_subplots
import plotly.graph_objects as go


class Data:
    def __init__(self, name: str):
        self.name = name
        self.src = None
        self.dg_data = {}

    def set_dg_data(self, stations, data_folder, low_freq=None, high_freq=None):

        for stat in stations.keys():

            #  Extract DG data for each station in a dict
            self.dg_data[stat] = self.extract_data(data_folder, stat,
                                                   low_freq=low_freq,
                                                   high_freq=high_freq)

    def extract_data(self, data_folder, stat, low_freq=None, high_freq=None):
        # Scan obspy files for each simulation (self.name)
        # and each station (stat) in folder (data_folder)

        velo = [f for f in scandir(data_folder) if f.is_file()
                and f.name.startswith(self.name)
                and f.name.__contains__(stat)]

        dg_data = {}
        if len(velo) > 0:  # Check if there are dg velo files for the simulation
            dg_data['low_freq'] = low_freq
            dg_data['high_freq'] = high_freq
            st = obspy.read(velo[0].path)
            # st.plot()

            # Save original data
            ori_st = st.copy()
            dg_data['ori_vx'] = ori_st.traces[0].data
            dg_data['ori_vy'] = ori_st.traces[1].data
            dg_data['ori_vz'] = ori_st.traces[2].data

            # Cut the initial zeros of the traces
            # ti = st.traces[0].stats.starttime
            # tf = st.traces[0].stats.endtime
            # st.trim(ti + self.dg_data[stat]['trim_t'], tf)

            # Save trimmed data
            cut_st = st.copy()
            dg_data['vx'] = cut_st.traces[0].data
            dg_data['vy'] = cut_st.traces[1].data
            dg_data['vz'] = cut_st.traces[2].data

            cut_norm_st = cut_st.copy()
            cut_norm_st.normalize()
            dg_data['norm_vx'] = cut_norm_st.traces[0].data
            dg_data['norm_vy'] = cut_norm_st.traces[1].data
            dg_data['norm_vz'] = cut_norm_st.traces[2].data

            # Only low pass filtering
            if low_freq and not high_freq:
                st.filter('lowpass', freq=low_freq, corners=2, zerophase=True)
                cut_low_filt_st = st.copy()
                dg_data['low_vx'] = cut_low_filt_st.traces[0].data
                dg_data['low_vy'] = cut_low_filt_st.traces[1].data
                dg_data['low_vz'] = cut_low_filt_st.traces[2].data
                st.normalize()
                cut_low_filt_norm_st = st.copy()
                dg_data['low_norm_vx'] = cut_low_filt_norm_st.traces[0].data
                dg_data['low_norm_vy'] = cut_low_filt_norm_st.traces[1].data
                dg_data['low_norm_vz'] = cut_low_filt_norm_st.traces[2].data

            # Only high pass filtering
            if not low_freq and high_freq:
                st.filter('highpass', freq=high_freq, corners=2, zerophase=True)
                cut_high_filt_st = st.copy()
                dg_data['high_vx'] = cut_high_filt_st.traces[0].data
                dg_data['high_vy'] = cut_high_filt_st.traces[1].data
                dg_data['high_vz'] = cut_high_filt_st.traces[2].data
                st.normalize()
                cut_high_filt_norm_st = st.copy()
                dg_data['high_norm_vx'] = cut_high_filt_norm_st.traces[0].data
                dg_data['high_norm_vy'] = cut_high_filt_norm_st.traces[1].data
                dg_data['high_norm_vz'] = cut_high_filt_norm_st.traces[2].data

            # Bandpass filtering
            if low_freq and high_freq:
                st.filter('bandpass', freqmin=low_freq, freqmax=high_freq,
                          corners=2, zerophase=True)
                cut_band_filt_st = st.copy()
                dg_data['band_vx'] = cut_band_filt_st.traces[0].data
                dg_data['band_vy'] = cut_band_filt_st.traces[1].data
                dg_data['band_vz'] = cut_band_filt_st.traces[2].data
                st.normalize()
                cut_band_filt_norm_st = st.copy()
                dg_data['band_norm_vx'] = cut_band_filt_norm_st.traces[0].data
                dg_data['band_norm_vy'] = cut_band_filt_norm_st.traces[1].data
                dg_data['band_norm_vz'] = cut_band_filt_norm_st.traces[2].data

        return dg_data

    def extract_src_fields(self, sim_path):
        with open(sim_path, 'rb') as f:
            model_src = pickle.load(f)

        # Load and add a dict with src fields
        self.src = model_src

    def plot_simulation(self):
        stk = np.arange(self.src['n_stk'])
        dip = np.arange(self.src['n_dip'])

        fig = make_subplots(rows=3, cols=4, subplot_titles=('slip', 'vx', 'vy',
                                                            'vz', 'rupt_time',
                                                            'cut_vx', 'cut_vy',
                                                            'cut_vz',
                                                            'rise time', 'filt_vx',
                                                            'filt_vy', 'filt_vz'))

        tickfont = dict(color="black", size=12, family="Arial Black")
        xaxis = dict(title="<b> strike (Km) </b>", tickfont=tickfont)
        yaxis = dict(title="<b> dip (Km) </b>", tickfont=tickfont,
                     autorange="reversed")
        title_a = dict(text="<b> Interpolated Slip </b>",
                       font_family="Arial Black",
                       font_color="black", x=0.5, y=0.85)

        # slip distribution with rupture time contours
        colorbar_slip = dict(lenmode='fraction', len=0.21, thickness=15,
                             bordercolor="black", orientation='h',
                             x=0.11, y=0.66)
        data_slip = go.Heatmap(z=self.src['SLIP'], x=stk, y=dip, hoverongaps=False,
                               colorscale='cividis', colorbar=colorbar_slip)
        colorscale = [[0, 'rgb(250, 250, 250)'], [1.0, 'rgb(255, 255, 255)']]
        contours = dict(coloring='lines', showlabels=True)
        data_slip_rupt = go.Contour(z=self.src['RUPT_T'], x=stk, y=dip,
                                    contours=contours, line_width=2,
                                    showscale=False, colorscale=colorscale)

        # rupture time distribution
        colorbar_rupt = dict(lenmode='fraction', len=0.21, thickness=15,
                             bordercolor="black", orientation='h',
                             x=0.11, y=0.28)
        data_rupt = go.Heatmap(z=self.src['RUPT_T'], x=stk, y=dip, hoverongaps=False,
                               colorbar=colorbar_rupt)

        # rupture time distribution
        colorbar_rise = dict(lenmode='fraction', len=0.21, thickness=15,
                             bordercolor="black", orientation='h',
                             x=0.11, y=-0.13)
        data_rise = go.Heatmap(z=self.src['RISE'], x=stk, y=dip, hoverongaps=False,
                               colorbar=colorbar_rise)

        fig.add_trace(data_slip, row=1, col=1)
        fig.add_trace(data_slip_rupt, row=1, col=1)
        fig.add_trace(data_rupt, row=2, col=1)
        fig.add_trace(data_rise, row=3, col=1)

        fig.update_yaxes(row=1, col=1, autorange="reversed", title='dip',
                         tickfont=tickfont, scaleanchor="x", scaleratio=1)
        fig.update_yaxes(row=2, col=1, autorange="reversed", title='dip',
                         tickfont=tickfont, scaleanchor="x", scaleratio=1)
        fig.update_yaxes(row=3, col=1, autorange="reversed", title='dip',
                         tickfont=tickfont, scaleanchor="x", scaleratio=1)
        fig.update_xaxes(row=1, col=1, range=[np.min(stk), np.max(stk)])
        fig.update_xaxes(row=2, col=1, range=[np.min(stk), np.max(stk)])
        fig.update_xaxes(row=3, col=1, range=[np.min(stk), np.max(stk)])

        for istat, stat in enumerate(self.dg_data.keys()):

            vx = self.dg_data[stat]['ori_vx']
            vy = self.dg_data[stat]['ori_vy']
            vz = self.dg_data[stat]['ori_vz']
            max_vx = np.max(np.abs(vx))
            max_vy = np.max(np.abs(vy))
            max_vz = np.max(np.abs(vz))
            max_stat_vx = f"{stat} vx: {max_vx:.3f}"
            max_stat_vy = f"{stat} vy: {max_vy:.3f}"
            max_stat_vz = f"{stat} vz: {max_vz:.3f}"
            fig.add_trace(go.Scatter(y=vx + istat * max_vx, mode='lines',
                                     name=max_stat_vx,
                                     line=dict(color='red', width=1)),
                          row=1, col=2)
            fig.add_trace(go.Scatter(y=vx + istat * max_vy, mode='lines',
                                     name=max_stat_vy,
                                     line=dict(color='black', width=1)),
                          row=1, col=3)
            fig.add_trace(go.Scatter(y=vx + istat * max_vz, mode='lines',
                                     name=max_stat_vz,
                                     line=dict(color='blue', width=1)),
                          row=1, col=4)

            vx = self.dg_data[stat]['vx']
            vy = self.dg_data[stat]['vy']
            vz = self.dg_data[stat]['vz']
            max_vx = np.max(np.abs(vx))
            max_vy = np.max(np.abs(vy))
            max_vz = np.max(np.abs(vz))
            max_stat_vx = f"{stat} vx: {max_vx:.3f}"
            max_stat_vy = f"{stat} vy: {max_vy:.3f}"
            max_stat_vz = f"{stat} vz: {max_vz:.3f}"

            fig.add_trace(go.Scatter(y=vx + istat * max_vx, mode='lines',
                                     name=max_stat_vx,
                                     line=dict(color='red', width=1)),
                          row=2, col=2)
            fig.add_trace(go.Scatter(y=vx + istat * max_vy, mode='lines',
                                     name=max_stat_vy,
                                     line=dict(color='black', width=1)),
                          row=2, col=3)
            fig.add_trace(go.Scatter(y=vx + istat * max_vz, mode='lines',
                                     name=max_stat_vz,
                                     line=dict(color='blue', width=1)),
                          row=2, col=4)

            low_freq = self.dg_data[stat]['low_freq']
            high_freq = self.dg_data[stat]['high_freq']

            if low_freq and not high_freq:
                vx = self.dg_data[stat]['low_norm_vx']
                vy = self.dg_data[stat]['low_norm_vy']
                vz = self.dg_data[stat]['low_norm_vz']
                max_vx = np.max(np.abs(vx))
                max_vy = np.max(np.abs(vy))
                max_vz = np.max(np.abs(vz))

                fig.add_trace(go.Scatter(y=vx + istat * max_vx, mode='lines',
                                         showlegend=False,
                                         line=dict(color='red', width=1)),
                              row=3, col=2)
                fig.add_trace(go.Scatter(y=vx + istat * max_vy, mode='lines',
                                         showlegend=False,
                                         line=dict(color='black', width=1)),
                              row=3, col=3)
                fig.add_trace(go.Scatter(y=vx + istat * max_vz, mode='lines',
                                         showlegend=False,
                                         line=dict(color='blue', width=1)),
                              row=3, col=4)

            if not low_freq and high_freq:
                vx = self.dg_data[stat]['high_norm_vx']
                vy = self.dg_data[stat]['high_norm_vy']
                vz = self.dg_data[stat]['high_norm_vz']
                max_vx = np.max(np.abs(vx))
                max_vy = np.max(np.abs(vy))
                max_vz = np.max(np.abs(vz))

                fig.add_trace(go.Scatter(y=vx + istat * max_vx, mode='lines',
                                         showlegend=False,
                                         line=dict(color='red', width=1)),
                              row=3, col=2)
                fig.add_trace(go.Scatter(y=vx + istat * max_vy, mode='lines',
                                         showlegend=False,
                                         line=dict(color='black', width=1)),
                              row=3, col=3)
                fig.add_trace(go.Scatter(y=vx + istat * max_vz, mode='lines',
                                         showlegend=False,
                                         line=dict(color='blue', width=1)),
                              row=3, col=4)

            if low_freq and high_freq:
                vx = self.dg_data[stat]['band_norm_vx']
                vy = self.dg_data[stat]['band_norm_vy']
                vz = self.dg_data[stat]['band_norm_vz']
                max_vx = np.max(np.abs(vx))
                max_vy = np.max(np.abs(vy))
                max_vz = np.max(np.abs(vz))

                fig.add_trace(go.Scatter(y=vx + istat * max_vx, mode='lines',
                                         showlegend=False,
                                         line=dict(color='red', width=1)),
                              row=3, col=2)
                fig.add_trace(go.Scatter(y=vx + istat * max_vy, mode='lines',
                                         showlegend=False,
                                         line=dict(color='black', width=1)),
                              row=3, col=3)
                fig.add_trace(go.Scatter(y=vx + istat * max_vz, mode='lines',
                                         showlegend=False,
                                         line=dict(color='blue', width=1)),
                              row=3, col=4)

            if not low_freq and not high_freq:
                vx = self.dg_data[stat]['norm_vx']
                vy = self.dg_data[stat]['norm_vy']
                vz = self.dg_data[stat]['norm_vz']
                max_vx = np.max(np.abs(vx))
                max_vy = np.max(np.abs(vy))
                max_vz = np.max(np.abs(vz))

                fig.add_trace(go.Scatter(y=vx + istat * max_vx, mode='lines',
                                         line=dict(color='red', width=1),
                              showlegend=False),
                              row=3, col=2)
                fig.add_trace(go.Scatter(y=vx + istat * max_vy, mode='lines',
                                         line=dict(color='black', width=1),
                              showlegend=False),
                              row=3, col=3)
                fig.add_trace(go.Scatter(y=vx + istat * max_vz, mode='lines',
                                         line=dict(color='blue', width=1),
                                         showlegend=False),
                              row=3, col=4)

            fig.update_layout(title=dict(text=self.name,
                                         font=dict(color="black", size=16,
                                                   family="Arial Black"),
                                         automargin=True))

        fig.show(renderer="browser")

    def save_simulation_data(self, output_folder):
        if self.dg_data:
            out_file = output_folder + self.name + ".pickle"
            object_file = open(out_file, 'wb')
            pickle.dump(self, object_file)
            object_file.close()
            print(f" Saving simulation data in: {out_file}")

