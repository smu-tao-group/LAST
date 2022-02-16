#!/usr/bin/env python
"""
preliminary MD simulations
"""

import sys
import mdtraj as md
from amber import Amber


pdb = sys.argv[1]

# do 100ps nvt
sim_amber = Amber()
sim_amber.set_top_file(f"../inputs/{pdb}.prmtop")
sim_amber.set_cor_file(f"../inputs/{pdb}.inpcrd")
sim_amber.set_simulation_time(total_steps=50000, report_interval=1000)
sim_amber.set_up_simulation(type='nvt', minimize_energy=True)
sim_amber.do_simulation(file_name=f'../trajs/{pdb}_100ps_nvt')

trajs = md.load(f'../trajs/{pdb}_100ps_nvt.dcd', top=f'../inputs/{pdb}.prmtop')
trajs[-1].save_amberrst7(f'../rst/{pdb}_100ps_nvt_rts')

# do 100ps npt
sim_amber.set_cor_file(f"../rst/{pdb}_100ps_nvt_rts")
sim_amber.set_simulation_time(total_steps=50000, report_interval=1000)
sim_amber.set_up_simulation(type='npt', minimize_energy=True)
sim_amber.do_simulation(file_name=f'../trajs/{pdb}_100ps_npt')

trajs = md.load(f'../trajs/{pdb}_100ps_npt.dcd', top=f'../inputs/{pdb}.prmtop')
trajs[-1].save_amberrst7(f'../rst/{pdb}_100ps_npt_rts')
