import gromacs
import mdpow.equil, mdpow.fep
mdpow.equil.Simulation?
sim = mdpow.equil.Simulation(filename="./benzene/water.simulation")
G = mdpow.fep.Ghyd(simulation=sim, runtime=25)
G.setup(maxwarn=1)
G.save()

