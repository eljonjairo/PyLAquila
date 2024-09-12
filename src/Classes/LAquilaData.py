"""
 Fault Class
"""
import pickle
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from scipy.signal import butter, lfilter
from scipy.io import loadmat

class Data:
    def __init__(self, name: str, nstats: int, nt: int, dt: float,
                 freqmin: float, freqmax: float):
        self.name = name  # Simulation dir name
        self.dt = dt  # Simulation traces dt
        self.nt = nt  # Simulation traces number of samples
        self.time = {}  # Time array for traces
        self.freq = {}  # Freq array for traces amplitude spectrum
        self.freqmin = freqmin  # Min freq of bandpass filter
        self.freqmax = freqmax  # Max freq of bandpass filter
        self.nstats = nstats  # Simulation number of stations
        self.dhF = None  # Size of square subfaults
        self.slip = None  # Slip distribution
        self.rise_t = None  # Rise time distribution
        self.rupt_t = None  # Rupture time distribution
        self.rupt_v = None  # Rupture velocity distribution
        self.data = {}  # DG Traces

    def set_velo_acc(self, velox, veloy, veloz, stations):
        self.time = np.linspace(0, self.nt * self.dt, num=self.nt)
        # Design a Butterworth bandpass filter with order 4
        nyq = 0.5 / self.dt
        [b, a] = butter(4, [self.freqmin / nyq, self.freqmax / nyq], 'band')
        for istat, stat in enumerate(stations):
            self.data[stat] = {}
            # Save station name
            self.data[stat]['name'] = stat

            # Save filtered velocities
            self.data[stat]['vx'] = lfilter(b, a, velox[istat, :])
            self.data[stat]['vy'] = lfilter(b, a, veloy[istat, :])
            self.data[stat]['vz'] = lfilter(b, a, veloz[istat, :])

            # Calculate Peak Ground Velocity
            self.data[stat]['PGVx'] = np.max(np.abs(self.data[stat]['vx']))
            self.data[stat]['PGVy'] = np.max(np.abs(self.data[stat]['vy']))
            self.data[stat]['PGVz'] = np.max(np.abs(self.data[stat]['vz']))

            # Calculate accelerations
            ax = np.gradient(velox[istat, :], self.dt)
            ay = np.gradient(veloy[istat, :], self.dt)
            az = np.gradient(veloz[istat, :], self.dt)

            # Save filtered velocities
            self.data[stat]['ax'] = lfilter(b, a, ax)
            self.data[stat]['ay'] = lfilter(b, a, ay)
            self.data[stat]['az'] = lfilter(b, a, az)

            # Calculate Peak Ground Velocity
            self.data[stat]['PGAx'] = np.max(np.abs(self.data[stat]['ax']))
            self.data[stat]['PGAy'] = np.max(np.abs(self.data[stat]['ay']))
            self.data[stat]['PGAz'] = np.max(np.abs(self.data[stat]['az']))

    def set_arias(self):
        g = 981  # Gravity in cm/s^2
        # Differentiate velocity to obtain acceleration (a = dv/dt)
        for key in self.data.keys():
            ax = self.data[key]['ax']
            ay = self.data[key]['ay']
            az = self.data[key]['az']

            # Calculate Arias Intensity
            arias_int_x = np.cumsum(ax ** 2) * (np.pi / (2 * g)) * self.dt
            arias_int_y = np.cumsum(ay ** 2) * (np.pi / (2 * g)) * self.dt
            arias_int_z = np.cumsum(az ** 2) * (np.pi / (2 * g)) * self.dt

            # Calculate Arias Duration
            # Find the t_i where the Arias Intensity reaches 5% and 95%
            arias_5_perc_x = 0.05 * arias_int_x[-1]
            arias_95_perc_x = 0.95 * arias_int_x[-1]
            arias_5_perc_y = 0.05 * arias_int_y[-1]
            arias_95_perc_y = 0.95 * arias_int_y[-1]
            arias_5_perc_z = 0.05 * arias_int_z[-1]
            arias_95_perc_z = 0.95 * arias_int_z[-1]

            # Find the indices where the cumulative Arias Intensity crosses 5% and 95%
            t_5_perc_x = ax[np.where(arias_int_x >= arias_5_perc_x)][0]
            t_95_perc_x = ax[np.where(arias_int_x >= arias_95_perc_x)][0]
            t_5_perc_y = ay[np.where(arias_int_y >= arias_5_perc_y)][0]
            t_95_perc_y = ay[np.where(arias_int_y >= arias_95_perc_y)][0]
            t_5_perc_z = az[np.where(arias_int_z >= arias_5_perc_z)][0]
            t_95_perc_z = az[np.where(arias_int_z >= arias_95_perc_z)][0]

            # Arias Duration is the time difference between 5% and 95% Arias Intensity
            arias_dur_x = t_95_perc_x - t_5_perc_x
            arias_dur_y = t_95_perc_y - t_5_perc_y
            arias_dur_z = t_95_perc_z - t_5_perc_z

            # Save results
            self.data[key]['ax'] = ax
            self.data[key]['ay'] = ay
            self.data[key]['az'] = az
            self.data[key]['arias_int_x_5perc'] = t_5_perc_x
            self.data[key]['arias_int_y_5perc'] = t_5_perc_y
            self.data[key]['arias_int_z_5perc'] = t_5_perc_z
            self.data[key]['arias_int_x_95perc'] = t_95_perc_x
            self.data[key]['arias_int_y_95perc'] = t_95_perc_y
            self.data[key]['arias_int_z_95perc'] = t_95_perc_z
            self.data[key]['arias_int_x'] = arias_int_x
            self.data[key]['arias_int_y'] = arias_int_y
            self.data[key]['arias_int_z'] = arias_int_z
            self.data[key]['arias_dur_x'] = arias_dur_x
            self.data[key]['arias_dur_y'] = arias_dur_y
            self.data[key]['arias_dur_z'] = arias_dur_z

    def set_velo_spec(self):
        self.freq = np.fft.fftfreq(self.nt, d=self.dt)[:self.nt // 2]
        for key in self.data.keys():
            # Perform FFT on the seismic trace
            fft_vx = np.fft.fft(self.data[key]['vx'])
            fft_vy = np.fft.fft(self.data[key]['vy'])
            fft_vz = np.fft.fft(self.data[key]['vz'])

            # Amplitude spectrum (take the magnitude of the FFT and normalize)
            self.data[key]['spec_vx'] = np.abs(fft_vx)[:self.nt // 2] * 2 / self.nt
            self.data[key]['spec_vy'] = np.abs(fft_vy)[:self.nt // 2] * 2 / self.nt
            self.data[key]['spec_vz'] = np.abs(fft_vz)[:self.nt // 2] * 2 / self.nt

    def set_src(self, src_folder: str):
        src_file = src_folder + self.name.split('LAquila_')[1].split('_1.0')[0] + '.mat'
        print(f" Loading src file: {src_file}")
        try:
            input_source = loadmat(src_file)
            self.dhF = input_source['Fault3D']['dhF'][-1][0]
            self.slip = input_source['Fault3D']['slip'][-1][0]
            self.rise_t = input_source['Fault3D']['riseT'][-1][0]
            self.rupt_t = input_source['Fault3D']['rupTime'][-1][0]
            self.rupt_v = input_source['Fault3D']['vrupt'][-1][0]

        except FileNotFoundError:
            print()
            print(f"The file {src_file} does not exist.")
            print()

    def plot_data(self):
        pcolors = ['blue', 'green', 'red', 'mediumblue', 'orange', 'magenta', 'darkviolet', 'brown',
                   'purple', 'orangered', 'navy', 'beige', 'grey']
        tickfont = dict(color="black", size=24, family="Arial Black")
        fig = (make_subplots(rows=4, cols=6, vertical_spacing=0.1, horizontal_spacing=0.05,
                             subplot_titles=("slip", "rise time", "rupture velocity"),
                             specs=[[{}, {}, {}, {"rowspan": 2}, {"rowspan": 2}, {"rowspan": 2}],
                                    [{}, {}, {}, None, None, None],
                                    [{}, {}, {}, {"rowspan": 2}, {"rowspan": 2}, {"rowspan": 2}],
                                    [{}, {}, {}, None, None, None]], print_grid=True))

        colorscale = [[0, 'rgb(250, 250, 250)'], [1.0, 'rgb(255, 255, 255)']]
        contours = dict(coloring='lines', showlabels=True)
        colorbar_slip = dict(lenmode='fraction', len=0.13, thickness=15,
                             bordercolor="black", tickvals=[0, 0.4, 0.8, 1.2, 1.6],
                             x=0.06, y=0.74, orientation='h')
        colorbar_rise = dict(lenmode='fraction', len=0.13, thickness=15,
                             tickvals=[0, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0],
                             bordercolor="black", x=0.24, y=0.74, orientation='h')
        colorbar_ruptv = dict(lenmode='fraction', len=0.13, thickness=15,
                              tickvals=[1000, 2000, 3000, 3000],
                              bordercolor="black", x=0.41, y=0.74, orientation='h')

        fig.add_trace(go.Heatmap(z=self.slip, hoverongaps=False, colorscale='viridis',
                                 colorbar=colorbar_slip, zmin=0, zmax=1.6),
                      row=1, col=1)
        fig.add_trace(go.Contour(z=self.rupt_t, contours=contours, line_width=2,
                                 showscale=False, colorscale=colorscale),
                      row=1, col=1)
        fig.add_trace(go.Heatmap(z=self.rupt_v, hoverongaps=False, colorscale='hot',
                                 zmin=1000, zmax=3000, colorbar=colorbar_ruptv),
                      row=1, col=3)
        fig.add_trace(go.Heatmap(z=self.rise_t, hoverongaps=False,
                                 zmin=0, zmax=4.0, colorbar=colorbar_rise),
                      row=1, col=2)

        time = self.time
        # freq = np.log(self.freq)
        freq = self.freq
        lev_vx = 0
        lev_vy = 0
        lev_vz = 0
        lev_ax = 0
        lev_ay = 0
        lev_az = 0
        for k, key in enumerate(self.data.keys()):
            name = self.data[key]['name']
            # Plot velocities
            vx = self.data[key]['vx']
            vy = self.data[key]['vy']
            vz = self.data[key]['vz']
            lev_vx += np.max(abs(vx)) * 1.2
            lev_vy += np.max(abs(vy)) * 1.2
            lev_vz += np.max(abs(vz)) * 1.2
            fig.add_trace(go.Scatter(y=vx - lev_vx, x=time, mode='lines', name=name,
                                     line=dict(color=pcolors[k], width=1)), row=1, col=4)
            fig.add_trace(go.Scatter(y=vy - lev_vy, x=time, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=1, col=5)
            fig.add_trace(go.Scatter(y=vz - lev_vz, x=time, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=1, col=6)

            # Plot velocities amplitude spectrum
            spec_vx = self.data[key]['spec_vx']
            spec_vy = self.data[key]['spec_vy']
            spec_vz = self.data[key]['spec_vz']
            fig.add_trace(go.Scatter(y=spec_vx, x=freq, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=2, col=1)
            fig.add_trace(go.Scatter(y=spec_vy, x=freq, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=2, col=2)
            fig.add_trace(go.Scatter(y=spec_vz, x=freq, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=2, col=3)
            # Plot accelerations
            ax = self.data[key]['ax'] * 1.2
            ay = self.data[key]['ay'] * 1.2
            az = self.data[key]['az'] * 1.2
            lev_ax += np.max(abs(ax))
            lev_ay += np.max(abs(ay))
            lev_az += np.max(abs(az))
            fig.add_trace(go.Scatter(y=ax - lev_ax, x=time, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=3, col=4)
            fig.add_trace(go.Scatter(y=ay - lev_ay, x=time, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=3, col=5)
            fig.add_trace(go.Scatter(y=az - lev_az, x=time, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=3, col=6)

            # Plot Arias Intensities
            arias_x = self.data[key]['arias_int_x']
            arias_y = self.data[key]['arias_int_y']
            arias_z = self.data[key]['arias_int_z']
            fig.add_trace(go.Scatter(y=arias_x, x=time, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=3, col=1)
            fig.add_trace(go.Scatter(y=arias_y, x=time, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=3, col=2)
            fig.add_trace(go.Scatter(y=arias_z, x=time, mode='lines', showlegend=False,
                                     line=dict(color=pcolors[k], width=1)), row=3, col=3)

        ndip, nstk = np.shape(self.slip)
        fig.update_xaxes(title='time (s)', row=1, col=4)
        fig.update_xaxes(title='time (s)', row=1, col=5)
        fig.update_xaxes(title='time (s)', row=1, col=6)
        fig.update_yaxes(title='vx (cm/s)', row=1, col=4)
        fig.update_yaxes(title='vy (cm/s)', row=1, col=5)
        fig.update_yaxes(title='vz (cm/s)', row=1, col=6)
        fig.update_xaxes(title='freq (Hz)', range=[0, self.freqmax * 1.2], row=2, col=1)
        fig.update_xaxes(title='freq (Hz)', range=[0, self.freqmax * 1.2], row=2, col=2)
        fig.update_xaxes(title='freq (Hz)', range=[0, self.freqmax * 1.2], row=2, col=3)
        fig.update_yaxes(title='Amplitude spec vx', row=2, col=1)
        fig.update_yaxes(title='Amplitude spec vy', row=2, col=2)
        fig.update_yaxes(title='Amplitude spec vz', row=2, col=3)
        fig.update_xaxes(title='time (s)', row=3, col=4)
        fig.update_xaxes(title='time (s)', row=3, col=5)
        fig.update_xaxes(title='time (s)', row=3, col=6)
        fig.update_yaxes(title='ax (cm/s2)', row=3, col=4)
        fig.update_yaxes(title='ay (cm/s2)', row=3, col=5)
        fig.update_yaxes(title='az (cm/s2)', row=3, col=6)
        fig.update_xaxes(title='time (s)', row=3, col=1)
        fig.update_xaxes(title='time (s)', row=3, col=2)
        fig.update_xaxes(title='time (s)', row=3, col=3)
        fig.update_yaxes(title='Arias Intensity x', row=3, col=1)
        fig.update_yaxes(title='Arias Intensity y', row=3, col=2)
        fig.update_yaxes(title='Arias Intensity z', row=3, col=3)
        fig.update_yaxes(autorange="reversed", scaleanchor="x",
                         range=[ndip, 0], row=1, col=1)
        fig.update_yaxes(autorange="reversed", scaleanchor="x",
                         range=[ndip, 0], row=1, col=2)
        fig.update_yaxes(autorange="reversed", scaleanchor="x",
                         range=[ndip, 0], row=1, col=3)
        fig.update_layout(title_text=self.name, titlefont=tickfont)
        fig.show(renderer="browser")

    def save_simulation_data(self, output_folder):
        out_file = output_folder + self.name + ".pickle"
        object_file = open(out_file, 'wb')
        pickle.dump(self, object_file)
        object_file.close()
        print(f" Saving simulation data in: {out_file}")
