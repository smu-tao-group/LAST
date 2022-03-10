#!/usr/bin/env python

"""Seed structure selection
"""

import sys
import os
import glob
import pickle
import numpy as np
from tensorflow import keras
from sklearn.preprocessing import MinMaxScaler
import mdtraj as md
import statsmodels.api as sm


pdb = sys.argv[1]
iter_round = int(sys.argv[2])

# heavy atom indices
trajs = md.load(
    f"../trajs/{pdb}_aligned.dcd", top=f"../inputs/{pdb}.prmtop"
)
trajs = trajs.atom_slice(
    trajs.topology.select_atom_indices("heavy")
)

frames = trajs.n_frames // (iter_round + 1)
coors = trajs.xyz
coors = coors.reshape(trajs.n_frames, -1)

scaler = MinMaxScaler()
coors_scaled = scaler.fit_transform(coors)

# load latents and decoder
latents = pickle.load(
    open(f"../results/{pdb}_r{iter_round}_latents.pkl", "rb")
)
latent_dim = latents.shape[1]
decoder = keras.models.load_model(
    f"../models/{pdb}_r{iter_round}_decoder"
)

decoded_structure = decoder(latents).numpy()
decoded_structure = scaler.inverse_transform(decoded_structure) * 10
real_structure = scaler.inverse_transform(coors_scaled) * 10


def rmsd(p1, p2):
    """RMSD calculation
    """
    p1 = p1.reshape(-1, 3)
    p2 = p2.reshape(-1, 3)
    return np.sqrt(np.sum(np.square(p1 - p2)) / trajs.n_atoms)


# nonparametric fit
dens_u = sm.nonparametric.KDEMultivariate(
    data=latents, var_type='c' * latent_dim, bw='normal_reference'
)
cdf = dens_u.cdf()

# sort based on CDF
idxs = np.arange(0, latents.shape[0])
cdf_idx = np.stack((cdf, idxs), axis=1)
cdf_idx = sorted(cdf_idx, key=lambda x: x[0])

# detect outliers
# traverse from beginning
if os.path.isfile(f"../results/{pdb}_structure.pkl"):
    outliers = []
    outliers_structure = pickle.load(
        open(f"../results/{pdb}_structure.pkl", "rb")
    )
    idx = 0
else:
    outliers = [int(cdf_idx[0][1])]
    outliers_structure = [real_structure[outliers[-1]]]
    idx = 1

while idx < len(cdf_idx) and len(outliers) < 5:
    current_structure = real_structure[int(cdf_idx[idx][1])]
    current_latent = latents[int(cdf_idx[idx][1])]
    # whether this point is near within 1A
    near = False

    for out_structure in outliers_structure:
        if rmsd(current_structure, out_structure) < 1:
            near = True
            break

    if not near:
        outliers.append(int(cdf_idx[idx][1]))
        outliers_structure.append(real_structure[outliers[-1]])

    idx += 1


# traverse from the other side
if os.path.isfile(f"../results/{pdb}_structure.pkl"):
    idx = len(cdf_idx) - 1
else:
    outliers.append(int(cdf_idx[-1][1]))
    outliers_structure.append(real_structure[outliers[-1]])
    idx = len(cdf_idx) - 2

while idx >= 0 and len(outliers) < 10:
    current_structure = real_structure[int(cdf_idx[idx][1])]
    current_latent = latents[int(cdf_idx[idx][1])]
    near = False

    for out_structure in outliers_structure:
        if rmsd(current_structure, out_structure) < 1:
            near = True
            break

    if not near:
        outliers.append(int(cdf_idx[idx][1]))
        outliers_structure.append(real_structure[outliers[-1]])

    idx -= 1


# delete data to spare some space
del trajs
del coors


for i in range(len(outliers)):
    num = outliers[i]
    round_num = num // frames
    seed_num = (num - round_num * frames) // 100 + 1

    trajs = md.load(
        f'../trajs/{pdb}_r{round_num}_' + 's{0:02d}.dcd'.format(seed_num),
        top=f"../inputs/{pdb}.prmtop"
    )

    out_idx = num - round_num * frames - 100 * (seed_num - 1)
    trajs[out_idx].save_amberrst7(
        f'../rst/{pdb}_r{iter_round + 1}_rst.' + '{0:02d}'.format(i + 1)
    )


# save structure
pickle.dump(
    outliers_structure,
    open(f"../results/{pdb}_structure.pkl", "wb")
)

# save outliers idx
pickle.dump(
    outliers,
    open(f"../results/{pdb}_out_{iter_round}.pkl", "wb")
)

# save rmsd
trajs_dir = sorted(
    glob.glob(f"../trajs/{pdb}_r*.dcd"),
    key=lambda x: (int(x.split("_")[1][1:]), int(x.split("_")[2][1:-4])))

# get the latest iteration
trajs_dir = trajs_dir[-10:]

ref = md.load(f"../inputs/{pdb}.pdb")
ref = ref.atom_slice(
    ref.topology.select_atom_indices("alpha")
)

trajs = md.load(trajs_dir, top=f"../inputs/{pdb}.prmtop")
trajs = trajs.atom_slice(
    trajs.topology.select_atom_indices("alpha")
)

rmsd = md.rmsd(trajs, ref) * 10
rmsd_mean = np.mean(rmsd)

# save mean rmsd
if os.path.isfile(f"../results/{pdb}_rmsd.pkl"):
    rmsds = pickle.load(open(f"../results/{pdb}_rmsd.pkl", "rb"))
else:
    rmsds = []

rmsds.append(rmsd_mean)
pickle.dump(rmsds, open(f"../results/{pdb}_rmsd.pkl", "wb"))
