#!/usr/bin/env python

"""Convergence check
"""

import sys
import os
import pickle
import numpy as np


pdb = sys.argv[1]
patience = int(sys.argv[2])

if os.path.isfile(f"../results/{pdb}_rmsd.pkl"):
    rmsds = pickle.load(open(f"../results/{pdb}_rmsd.pkl", "rb"))
    max_val_idx = np.argmax(rmsds)
    if len(rmsds) - 1 - max_val_idx >= patience:
        sys.exit("converged!")
