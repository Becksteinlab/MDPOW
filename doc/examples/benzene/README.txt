==================
 Example: Benzene
==================

Follow the commands in the file :file:`session.py`; you cannot simply
run it because at some point you have to run the individual
calculations, either on a cluster in parallel or sequentially on a
multi-core workstation.




Hydration free energy
=====================

Parameters
----------
- MD_Relax: 5 ps
- MD_NPT: 10 ps
- FEP windows: 250 ps 

Result
------
- Delta A_hyd = -3.37 +- 0.78 kJ/mol = -0.81 +/- 0.19
- exp: -0.87 kcal/mol

So, not too bad for very short simulations. 

With 5 ns per window we got -2.97 +/- 0.21 kJ/mol = -0.70 +/- 0.05
kca/mol, indicating that the better value with the shorter simulation
time is not converged and fortuituously close to the experimental
value.

How to do it (after having run all the simulations)::
  >>> G = mdpow.fep.Ghyd(filename="./FEP/water/Ghyd.fep")
  >>> G.analyze()
  -3.37211 (0.78367)
  >>> print G.results.DeltaA.Gibbs
  -3.37211 (0.78367)
  >>> mdpow.fep.kJ_to_kcal(G.results.DeltaA.Gibbs)
  -0.805954 (0.187302)

Plot dV/dlambda ::
 
  >>> import pylab
  >>> G.plot()
  >>> pylab.savefig("dVdl.pdf")

.. image:: ../dVdl/pdf
    
   Derivative of the potential energy with the coupling parameter
   lambda for the two parts of the free energy perturbation
   calculation. Left: Switching off of Coulomb interactions. Right:
   Switching off of the Lennard-Jones (van der Waals-type)
   interactions of the uncharged solute.

