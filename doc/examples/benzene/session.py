import mdpow.equil
S = mdpow.equil.WaterSimulation(molecule="BNZ")
S.topology("benzene.itp")
S.solvate(struct="benzene.pdb")
S.energy_minimize()
S.MD_relaxed(runtime=5)   # should be at least 1e3 ps for production not just 5 ps

# run simulation externally or use MDrunner
# (see docs for using mpi etc)
import gromacs
r = gromacs.run.MDrunner(dirname=S.dirs['MD_relaxed'], deffnm="md", c="md.pdb", cpi=True, append=True, v=True)
r.run()   # runs mdrun in the python shell


S.MD(runtime=10, qscript=['local.sh'])  # should be at least 10e3 ps for production, not just 10 ps
# run simulation
r = gromacs.run.MDrunner(dirname=S.dirs['MD_NPT'], deffnm="md", c="md.pdb", cpi=True, append=True, v=True)
r.run()   # runs mdrun in the python shell


import mdpow.fep
gwat = mdpow.fep.Ghyd(simulation=S, runtime=10)
gwat.setup()

# run multiple simulations on cluster



O = mdpow.equil.OctanolSimulation(molecule="BNZ")
O.topology("benzene.itp")
O.solvate(struct="benzene.pdb")
O.energy_minimize()
O.MD_relaxed(runtime=0.5)
