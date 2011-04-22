import mdpow.equil
S = mdpow.equil.WaterSimulation(molecule="BNZ")
S.topology("benzene.itp")
S.solvate(struct="benzene.pdb")
S.energy_minimize()
S.MD_relaxed(runtime=0.5)

# run simulation externally or use
# gromacs.run.mdrun() ...

S.MD(runtime=5, qscript=['local.sh'])

# run simulation

import mdpow.fep
gwat = mdpow.fep.Ghyd(simulation=S, runtime=10)
gwat.setup()

# run multiple simulations on cluster
