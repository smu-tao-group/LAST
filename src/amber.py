#!/usr/bin/env python

"""Amber simulation
"""

from openmm.app import (
    HBonds, PME, Simulation, DCDReporter, StateDataReporter
)
from openmm import LangevinMiddleIntegrator, MonteCarloBarostat
from openmm.unit import nanometer, kelvin, picosecond, femtosecond, bar
from openmm.app import AmberPrmtopFile, AmberInpcrdFile


class Amber():
    """
    Amber simulation in OpenMM
    """
    def __init__(self, total_steps=None, report_interval=None):
        self.top = None
        self.cor = None
        self.system = None
        self.total_steps = total_steps
        self.report_interval = report_interval
        self.integrator = None
        self.simulation = None

    def set_top_file(self, top):
        self.top = AmberPrmtopFile(top)

    def set_cor_file(self, cor):
        self.cor = AmberInpcrdFile(cor)

    def set_simulation_time(self, total_steps, report_interval):
        """Set simulation time params
        """
        # total simulation time = tstep * 2fs
        self.total_steps = total_steps
        # num of DCD frames = tstep / report_interval
        self.report_interval = report_interval

    def set_up_simulation(
        self, type='nvt', minimize_energy=True
    ):
        """Set up simulation details.
        """
        self.integrator = LangevinMiddleIntegrator(
            300 * kelvin, 1.0 / picosecond, 2.0 * femtosecond
        )

        # set up system
        self.system = self.top.createSystem(
            nonbondedMethod=PME, nonbondedCutoff=1.0 * nanometer,
            constraints=HBonds
        )

        if type == 'npt':
            # add force to keep NPT
            self.system.addForce(MonteCarloBarostat(1 * bar, 300 * kelvin))

        self.simulation = Simulation(
            self.top.topology, self.system, self.integrator
        )

        # load history if given check point file
        self.simulation.context.setPositions(self.cor.positions)

        if not minimize_energy:
            self.simulation.context.setVelocitiesToTemperature(300)

        if self.cor.boxVectors is not None:
            self.simulation.context.setPeriodicBoxVectors(
                *self.cor.boxVectors
            )

        if minimize_energy:
            self.simulation.minimizeEnergy()

    def do_simulation(self, file_name):
        """Conduct simulation.
        """
        self.simulation.reporters.append(DCDReporter(
            f"{file_name}.dcd", self.report_interval)
        )

        self.simulation.reporters.append(StateDataReporter(
            f"{file_name}.log", self.report_interval, step=True,
            potentialEnergy=True, temperature=True, volume=True,
            progress=True, remainingTime=True, speed=True,
            totalSteps=self.total_steps)
        )

        self.simulation.step(self.total_steps)
