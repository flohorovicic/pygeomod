"""remove trailing comma from exported_grid files"""

import sys, os

filename = sys.argv[1]

lines = open(filename).readlines()

new_name = os.path.splitext(filename)[0] + "_fixed.csv"

file_new = open(new_name, "w")

for line in lines:
    file_new.write(line[:-2] + "\n")


