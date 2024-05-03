"""
    Fault Generator
    Params:
        dh_f: float = Size of output square subfaults
        .mat = Structure matlab file
        name: string = Name of structure matlab file
    Return:
        Fault object with slip, rupture time, and rise time distributions
        over the interpolated fault
"""

from src.Fault import Fault

if __name__ == '__main__':

    name = 'LAquilaCirella03'  # Output Fault object name

    dh_f: float = 0.5  # Output subfaults size in Km

    # Input matlab structure name
    file_name: str = "s2009LAQUIL03CIRE"
    in_fault = Fault.load_mat_file(file_name)

    # Creates a new instance of the Fault class
    LAquila_fault = Fault.Fault(name, dh_f, in_fault)
    print(LAquila_fault)
    LAquila_fault.plot_fault_inputs()

    LAquila_fault.interpolate()