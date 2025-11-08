#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse
import glob

parser = argparse.ArgumentParser(description="Track maximum temperature from binary outputs.")
parser.add_argument("pattern", help="Glob for binary files (e.g. 'out/rank0_step*.bin')")
args = parser.parse_args()

files = sorted(glob.glob(args.pattern))
ts, Tmax = [], []

for f in files:
    with open(f, "rb") as fh:
        nx = int.from_bytes(fh.read(4), 'little')
        ny = int.from_bytes(fh.read(4), 'little')
        nz = int.from_bytes(fh.read(4), 'little')

        arr = np.frombuffer(fh.read(), dtype=np.float64)
        arr = arr.reshape((nx*ny*nz,))
        
        Tmax.append(arr.max())

    # crude time from step number in filename
    step = int(f.split("step")[-1].split(".")[0])
    ts.append(step)

plt.plot(ts, Tmax, marker="o")
plt.xlabel("Step")
plt.ylabel("Max T [K]")
plt.title("Max temperature evolution")
plt.grid()
plt.show()