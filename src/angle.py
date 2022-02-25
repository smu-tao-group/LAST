#!/usr/bin/env python

"""Calculate ADK angles
"""

import mdtraj as md
import numpy as np


def get_atoms(traj):
    """
    get atom slice
    """
    return traj.atom_slice(
        traj.topology.select_atom_indices("heavy")
    )


def vector_angle(vector1, vector2):
    """
    calculate angle of two vectors
    """
    unit_vector1 = vector1 / np.linalg.norm(vector1)
    unit_vector2 = vector2 / np.linalg.norm(vector2)
    dot_product = np.dot(unit_vector1, unit_vector2)
    angle_radian = np.arccos(dot_product)
    return np.degrees(angle_radian)


def get_traj_vector_angles(traj, selection):
    """
    get 2 vector angles
    """
    # selection : v1_tail, v1_head, v2_head, v3_tai, v3_head, v4_head
    v1_tail = traj.atom_slice(traj.topology.select(selection[0]))
    v1_head = traj.atom_slice(traj.topology.select(selection[1]))
    v2_head = traj.atom_slice(traj.topology.select(selection[2]))
    v3_tail = traj.atom_slice(traj.topology.select(selection[3]))
    v3_head = traj.atom_slice(traj.topology.select(selection[4]))
    v4_head = traj.atom_slice(traj.topology.select(selection[5]))

    # center of mass
    v1_tail_cor = md.compute_center_of_mass(v1_tail)
    v1_head_cor = md.compute_center_of_mass(v1_head)
    v2_head_cor = md.compute_center_of_mass(v2_head)
    v3_tail_cor = md.compute_center_of_mass(v3_tail)
    v3_head_cor = md.compute_center_of_mass(v3_head)
    v4_head_cor = md.compute_center_of_mass(v4_head)

    # calculate vectors
    v1 = v1_head_cor - v1_tail_cor
    v2 = v2_head_cor - v1_tail_cor
    v3 = v3_head_cor - v3_tail_cor
    v4 = v4_head_cor - v3_tail_cor

    # calculate angles
    v1_v2_angles = list(map(vector_angle, v1, v2))
    v3_v4_angles = list(map(vector_angle, v3, v4))

    return v1_v2_angles, v3_v4_angles


selection_1ake = [
    "resid 90 to 99", "resid 35 to 55", "resid 115 to 125",
    "resid 115 to 125", "resid 179 to 185", "resid 125 to 153"
]

selection_1dvr = [
    "resid 95 to 101 or resid 106 to 108", "resid 39 to 59",
    "resid 124 to 134", "resid 124 to 134", "resid 188 to 194",
    "resid 134 to 162"
]

selection_2ak3 = [
    "resid 95 to 104", "resid 41 to 61", "resid 119 to 129",
    "resid 119 to 129", "resid 183 to 189", "resid 129 to 157"
]


def get_angle(trajs_dir, top, stride=None, pdb="2ak3"):
    """
    get angle based on pdb
    """
    close_state = md.load(trajs_dir, top=top, stride=stride)
    close_state_alpha = get_atoms(close_state)

    if pdb == "2ak3":
        selection = [
            "resid 95 to 104", "resid 41 to 61", "resid 119 to 129",
            "resid 119 to 129", "resid 183 to 189", "resid 129 to 157"
        ]
    elif pdb == "1dvr":
        selection = [
            "resid 95 to 101 or resid 106 to 108", "resid 39 to 59",
            "resid 124 to 134", "resid 124 to 134", "resid 188 to 194",
            "resid 134 to 162"
        ]
    else:
        selection = [
            "resid 90 to 99", "resid 35 to 55", "resid 115 to 125",
            "resid 115 to 125", "resid 179 to 185", "resid 125 to 153"
        ]

    coors = get_traj_vector_angles(close_state_alpha, selection)
    return coors
