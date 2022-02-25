#!/usr/bin/env python

"""Train VAE model
"""

import sys
import pickle
import mdtraj as md
from sklearn.preprocessing import MinMaxScaler
from vae import build_vae


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

# scale
scaler = MinMaxScaler()
coors_scaled = scaler.fit_transform(coors)

encoder, decoder, vae = build_vae(
    input_dim=coors_scaled.shape[1:],
    encoder_neuron_nums=[512, 128, 32],
    latent_dim=2
)

history = vae.fit(
    x=coors_scaled, y=coors_scaled,
    shuffle=True, epochs=400, batch_size=16
)

# get latents
latents = encoder.predict(coors_scaled)[0]

# save latents, models
pickle.dump(latents, open(f"../results/{pdb}_r{iter_round}_latents.pkl", "wb"))
encoder.save(f"../models/{pdb}_r{iter_round}_encoder")
decoder.save(f"../models/{pdb}_r{iter_round}_decoder")
