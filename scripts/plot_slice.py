#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Plot a CSV slice produced by mpi-battery-heat.")
parser.add_argument("file", help="CSV file (out/rank*_xy_mid*.csv)")
parser.add_argument("--out", help="Output PNG filename (optional)")

args = parser.parse_args()

data = np.loadtxt(args.file, delimiter=",")

plt.imshow(data, origin="lower", cmap="inferno")
plt.colorbar(label="Temperature [K]")
plt.title(args.file)

if args.out:
    plt.savefig(args.out, dpi=150)
else:
    plt.show()
