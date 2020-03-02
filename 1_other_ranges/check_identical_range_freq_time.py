import numpy as np
import matplotlib.pyplot as pl
from matplotlib import patches
from scipy.optimize import curve_fit, root, fsolve
from scipy.signal import argrelextrema
from numpy import dot, outer
from os import listdir
from os.path import isfile, join

def main():
    
    # Set parameters
    n_minimax = 40                     # Number of minimax points
    R_min     = int(30000)
    R_max     = int(80000)

    path = "alpha_beta_of_N_"+str(n_minimax)

    path_time = "/home/wij17778/loctmp/01_minimax_python_time/1_other_ranges/"+path
    path_freq = "/home/wij17778/loctmp/02_minimax_python_frequency/1_other_ranges/"+path

    ranges_time = get_ranges(path_time, R_min, R_max)
    ranges_freq = get_ranges(path_freq, R_min, R_max)

    print("Common ranges:", list(set(ranges_time).intersection(ranges_freq)))

def get_ranges(path, R_min, R_max):

    files = [f for f in listdir(path) if isfile(join(path, f))]

    ranges = []

    for f in files:
        if(not f.startswith("alpha")): continue
        R_file = int(f[21:34])
        if(R_file > R_min and R_file < R_max):
           ranges.append(R_file)

    return ranges

if __name__ == "__main__":
    main()

