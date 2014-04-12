"""Plot results of exported grids, e.g. in cross-section with matplotlib"""

import sys, os
from matplotlib import use
use("Agg")
import matplotlib.pyplot as plt
import numpy as np

"""Settings"""
create_single_plots = True

# find all exported_grid*txt files in directory, and remove trailing commas

for filename in os.listdir("."):
    if "exported" in filename and "fixed" not in filename and "original" not in filename:
        lines = open(filename).readlines()
        new_name = os.path.splitext(filename)[0] + "_fixed.txt"
        file_new = open(new_name, "w")
        for line in lines:
            file_new.write(line[:-2] + "\n")
    if "delxyz" in filename:
        # determine nx, ny, nz
        lines = open(filename).readlines()
        nx = int(lines[0].split("*")[0])
        ny = int(lines[1].split("*")[0])
        nz = len(lines[2].split(","))-1
        
print nx, ny, nz

# now: open fixed files and store results in a list of numpy arrays

for filename in os.listdir("."):
    if "fixed" in filename:
        print filename
        try:
            grid = np.loadtxt(filename, delimiter = ',') 
        except ValueError:
            print "Could not load exported grid " + filename
            continue
        # reshape to match propert nx, ny, nz structure
        print grid.shape
        grid = grid.reshape(nz, ny, nx)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(grid[:,ny/2,:], interpolation='nearest', origin='lower_left')
        ax.set_title(filename)
        plt.savefig(os.path.splitext(filename)[0] + ".png")
        
        




