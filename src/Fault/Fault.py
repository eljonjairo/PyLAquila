"""
   Fault Class

   """
import scipy
import utm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


class Fault:
    def __init__(self, name: str, dh_f: float, in_fault: dict):
        # Check the input dict
        assert bool(in_fault), f" input dict is empty. "

        # Assign object attributes
        self.name = name
        self.dh_f = dh_f

        for key, value in in_fault.items():
            exec("".join(["self.in_", key, " = value"]))

        self.in_ndip, self.in_nstk = self.in_slip.shape

    def __repr__(self):
        return f'''
                    name: {self.name} 
                    dh_f: {self.dh_f} 
                    nstk: {self.in_nstk} 
                    ndip: {self.in_ndip} '''

    def interpolate(self):
        """ Convert input fault coordinates to X, y in km """
        X_in, Y_in, t1, t2 = utm.from_latlon(self.in_lat, self.in_lon,
                                             33, 'T')/1000

    def plot_fault_inputs(self):
        fig = plt.figure(constrained_layout=True)
        ax1 = fig.add_subplot(211)
        ax1.set_xlabel(" Lon (°)")
        ax1.set_ylabel(" Lat (°)")
        ax1.set_aspect('equal', adjustable='box')
        ax1.set_title(" Input Slip ")
        levels = np.linspace(0, np.max(self.in_rupt_time), 9)
        cs = ax1.contour(self.in_lon, self.in_lat, self.in_rupt_time, levels,
                         colors=('w',), linewidths=(0.3,), origin='lower')
        ax1.clabel(cs, fmt='%2.1f', colors='w', fontsize=10)
        fp = ax1.pcolormesh(self.in_lon, self.in_lat, self.in_slip, cmap=cm.viridis)
        plt.colorbar(fp, location='bottom', shrink=.4)

        ax2 = fig.add_subplot(212)
        ax2.set_xlabel(" Lon (°)")
        ax2.set_ylabel(" Lat (°)")
        ax2.set_aspect('equal', adjustable='box')
        ax2.set_title("Rise time (s)")
        fp = ax2.pcolormesh(self.in_lon, self.in_lat, self.in_rise_time, cmap=cm.viridis)
        plt.colorbar(fp, location='bottom', shrink=.4)
        plt.show()


def load_mat_file(infile):
    """
    Load matlab structure file
    :param infile:
        infile: string =  matlab structure's name
    : return:
        in_fault: dict = input coords, slip, rupture time, rise time of the file

    """
    print()
    print(f" Loading matlab file from {"../Inputs/" + infile} ")
    print()

    # Load Matlab Finite Fault input file
    in_fault = scipy.io.loadmat("".join(["Inputs/", infile]))
    fault = dict(z=in_fault[infile]['geoZ'][-1][0],
                 lat=in_fault[infile]['geoLAT'][-1][0],
                 lon=in_fault[infile]['geoLON'][-1][0],
                 slip=in_fault[infile]['slipSPL'][-1][0],
                 rise_time=in_fault[infile]['riseSPL'][-1][0],
                 rupt_time=in_fault[infile]['timeSPL'][-1][0],
                 hypo_lon=in_fault[infile]['evLON'][-1][0][0][0],
                 hypo_lat=in_fault[infile]['evLAT'][-1][0][0][0],
                 hypo_z=in_fault[infile]['evDPT'][-1][0][0][0])

    return fault
