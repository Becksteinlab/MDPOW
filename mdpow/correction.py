"""
Volume correction for running NVT instead of NPT

Compute the correction for Delta G_hyd for Delta F_hyd from NVT sims.

Required quantities for the correction DeltaW
=============================================

- volume of the solute
  * in water
  * in octanol

- volume of the simulation box of the *NVT* FEP simulation

- isothermal compressibility of
  * TIP4P water
  * 1-octanol (OPLS-AA)

Solute volume
-------------

1. average volume of the system (solute + water) in the NPT
   simulations with N_w water molecules: V_NPT

2. volume of a box with N_w water molecules at the same temperature
   and pressure is V_wat = N/N_A * M/rho (see :func:`volume`). The
   water density for TIP4P as a function of T and P (equation of
   state) only needs to be computed once.

   Right now we're actually just storing the value at T=300 K and 1
   bar: 993. kg/m^3 (computed from 20 ns MD with 296 water molecules).

3. V_solute = V_NPT - V_water

Similarly for 1-octanol as a solvent.


Volume of the FEP simulation box
--------------------------------

The FEP simulations are carried out at constant volume because
pressure coupling has been shown by Mobley, Shirts and co-workers
(ref?) to converge much more slowly and to be less accurate than NVT
unless one uses special barostats (such as Anderson).

The box volume is simply obtained from the last frame of the NPT
simulation.


Isothermal compressibility
--------------------------

kappa_T needs to be calculated for a given T and P. We compute it from
a long NPT simulation via the fluctuation formula and store the value.

"""
import tempfile
import numpy

import gromacs
from gromacs import GromacsError, MissingDataError
from gromacs.utilities import unlink_gmx

from recsql.convert import besttype
from numkit.observables import QuantityWithError

from mdpow.tables import N_Avogadro, molecular_weight, solvent_selections, solvent_density

import logging
logger = logging.getLogger("mdpow.correction")


testdata = """
Statistics over 10064096 steps [ 0.0000 through 20128.1900 ps ], 3 data sets
All statistics are over 2012820 points

Energy                      Average   Err.Est.       RMSD  Tot-Drift
-------------------------------------------------------------------------------
Pressure                   0.973056       0.08    616.702   0.358289  (bar)
Volume                      8.91529     0.0015  0.0856889 -0.000560363  (nm^3)
Density                     993.321       0.17    9.54136    0.04978  (kg/m^3)
"""


def parse_energy(data):
    """Put output from :program:`g_energy` into a recarray.

    Input *data* needs to look like ::

        Statistics over 10064096 steps [ 0.0000 through 20128.1900 ps ], 3 data sets
        All statistics are over 2012820 points

        Energy                      Average   Err.Est.       RMSD  Tot-Drift
        -------------------------------------------------------------------------------
        Pressure                   0.973056       0.08    616.702   0.358289  (bar)
        Volume                      8.91529     0.0015  0.0856889 -0.000560363  (nm^3)
        Density                     993.321       0.17    9.54136    0.04978  (kg/m^3)

    Only the lines starting with the keyword *Energy* up to the next
    blank line are considered. Data are considered white-space separated.
    """

    in_datablock = False
    records = []
    columns = None
    for line in data.split('\n'):
        line = line.strip()
        if line.startswith("----------") or len(line) == 0:
            continue
        if line.startswith("Energy"):
            columns = [s.replace('.','').replace('-','') for s in line.split()] + ["unit"]
            in_datablock = True
            continue
        if in_datablock:
            fields = line.split()
            if len(fields) == len(columns):
                records.append([besttype(s.replace('(','').replace(')',''))
                                for s in fields])
            else:
                in_datablock = False
    return numpy.rec.fromrecords(records, names=columns)

def get_observables(edr, observables=None, **kwargs):
    """Run :program:`g_energy` to get thermodynamic averages.

    :Arguments:
       *edr*
             Gromacs energy file
       *observables*
             list of energy (thermodynamic) averages to be computed;
             the default (for ``None``) is ::

               ["Pressure", "Volume", "Density", "Temperature"]
       *kwargs*
             other keywords are passed to
             :class:`gromacs.tools.G_energy`; for instance, *b* sets
             the start time (in ps) and *o* sets the output file name

    :Returns: :class:`~numpy.core.records.recarray` with the
              observable name in the first column ("*Energy*") and
              subsequent columns containing the exact averages
              produced by :program:`g_energy`.
    """
    if observables is None:
        observables = ["Pressure", "Volume", "Density", "Temperature"]
    logger.info("Analyzing observables %(observables)r in energy file %(edr)r...", vars())

    kwargs['input'] = observables
    kwargs['f'] = edr
    kwargs['stdout'] = False  # required for processing
    kwargs['stderr'] = False  # suppress stderr
    fd, kwargs['o'] = tempfile.mkstemp(suffix=".xvg")
    try:
        rc, data, stderr = gromacs.g_energy(**kwargs)
        if rc != 0:
            logger.error("g_energy failed")
            raise GromacsError("g_energy failed --- see output below:\n\n" + stderr + data)
    finally:
        unlink_gmx(kwargs['o'])
    return parse_energy(data)

def recwhere(r, value, name='Energy'):
    """Return record(s) in recarray *r* where named column *name* == *value*."""
    return r[r.field(name) == value]


def Observable(r, name, column='Energy', index=0):
    """Return a :class:`numkit.observables.QuantityWithError` for observable *name*.

    .. Note:: If for some unknown reason there is more than one entry
              in the recarray *r* that matches *name* in the field
              *column* then only the match at index *index* is
              returned. Other matches are silently ignored.
    """
    return QuantityWithError(recwhere(r, name, column).Average[index],
                             recwhere(r, name, column).ErrEst[index])

def v_mol(N, V):
    """Returns the molecular volume V/N."""
    return V/N

def volume(N, rho, M):
    """Volume from number of particles *N*, density *rho* and molecular weight *M*.

    Units:
     - *N* as number of particles (unit 1)
     - *rho* in kg m^-3 = g/L
     - *M* in g/mol

    :Returns: V = N/N_Avogadro * M / rho (in nm<sup>3</sup>)

    derivation
    # n = N/N_Avogadro
    # M = m/n = (rho*V)/n
    # V = n*M/rho = N*M/(N_A*rho)

    With input units:
      g/mol * mol * m^3 * (10^3 g)^-1 = 10^-3 m^3 = 10^-3 (10^9 nm)^3
      = 10^24 nm^3
    """
    return N*M/(rho*N_Avogadro) * 1e+24

def volume2(N, v_w):
    """Volume from number of particles *N* and molecular volume *v_w*.

    The molecular volume *v_w* is the volume per particle for given *T* and *P*.

    :Returns: V = N * v_w
    """
    return N * v_w

def calculate_solute_volume(solute_data, solvent_data=None,
                            N_solvent=None, rho_solvent=None, v_solvent=None,
                            solvent="water"):
    """Calculate the volume of the solute.

    Requires the formatted data from :def:`get_observables` for the
    the NPT simulation of the solute and either the corresponding data
    for the same simulation *without* the solute, *or* the number of
    solvent molecules, its density, and (if it's not water), its name
    (so that its molecular weight can be looked up in
    :data:`molecular_weight`).
    """
    solute_Volume = recwhere(solute_data, 'Volume')
    if len(solute_Volume) != 1:
        errmsg = "There should only be ONE Volume entry, not %r" % solute_Volume
        logger.error(errmsg)
        raise ValueError(errmsg)
    if not solute_Volume.unit[0] == "nm^3":
        errmsg = "Unit of Volume should be nm^3, not %r" % solute_Volume.unit[0]
        logger.error(errmsg)
        raise ValueError(errmsg)
    solute_Vbox = Observable(solute_data, "Volume")
    logger.info("Volume of the solute simulation box: %(solute_Vbox)r", vars())

    if solvent_data is not None:
        solvent_Vbox = Observable(solvent_data, 'Volume')
        # skipped sanity check for len...
        # Calling code should check that there are actually the same
        # number of solvent molecules in both simulations.
        logger.info("Calculating solvent box from simulation")
    elif not (N_solvent is None or v_solvent is None):
        solvent_Vbox = volume2(N_solvent, v_solvent)
        logger.info("Calculating solvent box from molecular volume")
        logger.info("Molecular volume: v=%(v_solvent)r nm^3", vars())
        logger.info("Number of solvent molecules: N=%(N_solvent)d", vars())
    elif not (N_solvent is None or rho_solvent is None):
        solvent_Vbox = volume(N_solvent, rho_solvent, molecular_weight[solvent])
        logger.info("Calculating solvent box from molecular weight and density")
        logger.info("Density: rho=%(rho_solvent)r kg/m^3", vars())
        logger.info("Molecular weight: M=%r g/mol", molecular_weight[solvent])
        logger.info("Number of solvent molecules: N=%(N_solvent)d", vars())
    else:
        errmsg = "No data for solvent"
        logger.error(errmsg)
        raise MissingDataError(errmsg)
    logger.info("Volume of the equivalent pure solvent simulation box: %(solvent_Vbox)r", vars())
    v_solute = solute_Vbox - solvent_Vbox
    logger.info("Volume of the solute: v_solute=%(v_solute)r nm^3", vars())
    return v_solute

def count_solvent_molecules(tpr, solvent="water"):
    """Count number of solvent  molecules from index file.

    .. Note:: Only solvent="water" or "octanol" supported. New
              selections can be added to :data:`solvent_selections`.
    """
    # This would be MUCH easier in MDAnalysis but I don't want to add
    # this as a dependency to mdpow. g_select would also work but that's
    # more recent and also requires reading of an xvg file.
    solventgroup = "solvent_%(solvent)s" % vars()   # should not start with a number for make_ndx!!
    try:
        cmd = ['keep 0', 'del 0', solvent_selections[solvent],
               'name 0 %(solventgroup)s' % vars(), 'q']
    except KeyError:
        raise NotImplementedError("solvent %(solvent)r not supported" % vars())

    fd, tmpndx = tempfile.mkstemp(suffix=".ndx")
    try:
        # silent, no screen output!!
        rc,output,junk = gromacs.cbook.make_ndx_captured(f=tpr, o=tmpndx,
                                                         input=cmd, stderr=False)
        if rc != 0:
            out = "\n".join(["GMX_FATAL: %s" % s for s in (junk + "\n\n" + output).split("\n")])
            errmsg = "Failed to build index. Look at the output below for hints:\n" + out
            logger.error(errmsg)
            raise GromacsError(errmsg)
    finally:
        unlink_gmx(tmpndx)
    try:
        ndx = gromacs.cbook.parse_ndxlist(output)
    except AttributeError:
        logger.error("Failed to find solvent %(solvent)r. Maybe specify solvent='octanol'?",
                     vars())
        raise ValueError("Failed to properly construct/parse make_ndx output.")
    N_solvent = None
    for rec in ndx:
        if rec['name'] == solventgroup:
            N_solvent = rec['natoms'] # == number of molecules as one atom per molecule
            break
    if N_solvent is None:
        errmsg = "Failed to find solvent %(solvent)r in %(tpr)r" % vars()
        logger.error(errmsg)
        raise ValueError(errmsg)
    logger.info("Found %(N_solvent)d %(solvent)s solvent molecules in %(tpr)r", vars())
    return N_solvent

def molecular_volume_analysis(tpr, edr, solvent="water", **kwargs):
    """Calculate molecular volume from NPT simulation data.

    v_mol = V/N

    with V as the average box volume (from the energy file *edr*) and
    the number of solvent molecules in the simulation (from the
    *tpr*). *kwargs* is passed to :func:`get_observables`.
    """
    kwargs.setdefault('stderr', False)
    obs = get_observables(edr, **kwargs)
    V = Observable(obs, "Volume")
    N = count_solvent_molecules(tpr, solvent=solvent)
    v = v_mol(N, V)
    logger.info("Molecular volume v_mol = %(v)r nm^3", vars())
    return v

class ThermodynamicAnalysis(object):
    """Analysis of thermodynamic properties of the equilibrium simulations.

    The main purpose of the class is to collect the edr and tpr file
    of a NPT simulation and calculate the solute volume.
    """
    def __init__(self, tpr, edr, solvent="water", **kwargs):
        self.tpr = tpr
        self.edr = edr
        self.solvent = solvent
        self.g_energy_kwargs = kwargs    # for g_energy

        self.N_solvent = self._count_solvent_molecules()
        self.data = self._get_observables(**kwargs)

        self.__v_solute = None

    def Observable(self, name):
        """Return a :class:`numkit.observables.QuantityWithError` for observable *name*."""
        column = "Energy"
        return QuantityWithError(recwhere(self.data, name, column).Average[0],
                                 recwhere(self.data, name, column).ErrEst[0])

    def _count_solvent_molecules(self):
        return count_solvent_molecules(self.tpr, solvent=self.solvent)

    def _get_observables(self, **kwargs):
        return get_observables(self.edr, **kwargs)

    def solute_volume(self, **kwargs):
        """Calculate the solute volume.

        Arguments are passed to :func:`calculate_solute_volume` unless they
        need to be set to values to match the system defined by this tpr (such
        as *N_solvent* and *solvent*). Useful keywords are *rho_solvent* or
        *v_solvent*.

        .. SeeAlso:: :attr:`~ThermodynamicAnalysis.v_solute` and :func:`calculate_solute_volume`
        """
        kwargs['N_solvent'] = self.N_solvent
        kwargs['solvent'] = self.solvent
        kwargs.setdefault('rho_solvent', solvent_density[self.solvent])
        self.__v_solute = calculate_solute_volume(self.data, **kwargs)  # update v_solute, too
        return self.__v_solute

    @property
    def v_solute(self):
        """Volume of the solute in nm^3.

        .. SeeAlso:: :meth:`~ThermodynamicAnalysis.solute_volume` (which
                     accepts *kwargs* and also updates this managed attribute).
        """
        if self.__v_solute is None:
            self.__v_solute = self.solute_volume()
        return self.__v_solute

    @property
    def v_solvent(self):
        """Volume of a solvent molecule V/N"""
        return self.volume/self.N_solvent

    @property
    def density(self):
        """Average total density of the simulation system in kg/m^3"""
        return self.Observable("Density")

    @property
    def volume(self):
        """Average total volume of the simulation system in nm^3"""
        return self.Observable("Volume")

    @property
    def temperature(self):
        """Average total temperature of the simulation system in K"""
        return self.Observable("Temperature")

    @property
    def pressure(self):
        """Average total pressure of the simulation system in bar"""
        return self.Observable("Pressure")

    def summary(self):
        """Log summary information"""
        logger.info("summary: tpr=%(tpr)r, edr=%(edr)r", vars(self))
        logger.info("summary: N_solvent = %(N_solvent)d", vars(self))
        logger.info("summary: T=%r K", self.temperature)
        logger.info("summary: P=%r bar", self.pressure)
        logger.info("summary: V=%r nm^3", self.volume)
        logger.info("summary: rho=%r kg/m^3 (total system density)", self.density)
        logger.info("summary: v_solvent = %r nm^3 (solvent particle volume)", self.v_solvent)
        logger.info("summary: v_solute = %r nm^3 (difference to ideal solvent box)", self.v_solute)

#------------------------------------------------------------
# from NVTcorrection.py
import numpy as np
from mdpow.tables import kappaT, bar

N_A = N_Avogadro ##6.0221415e23  # mol**-1

# base energy unit is kJ or kJ/mol

def kJ2kcal(x):
    return x/4.184

def energyUnit(x, **kwargs):
    if kwargs.get('unit', None) == 'kcal' or kwargs.get('kcal', False):
        return kJ2kcal(x)
    return x

def DeltaW(*args, **kwargs):
    """Correction for decoupling as used in POW

    (pretty much equivalent to the lowest order correction)
    """
    return DeltaW_A(*args, **kwargs)

def DeltaW_A(vs, V0, kappa=kappaT['DEFAULT'], **kwargs):
    """Exact correction for annihilation, see (Eq 13)"""
    x = vs/V0
    y = V0/(kappa/bar)*(x + np.log(1-x))
    return energyUnit(y, **kwargs)

def DeltaW_B(vs, V0, kappa=kappaT['DEFAULT'], **kwargs):
    """Exact correction for growing (Eq 12)"""
    x = vs/V0
    y = -V0/(kappa/bar)*(x - (x-1)*np.log(1-x))
    return energyUnit(y, **kwargs)

def pV(vs, P=1., **kwargs):
    y = vs*P*bar
    return energyUnit(y, **kwargs)

def DeltaW_lowestorder(vs, V0, kappa=kappaT['DEFAULT'], **kwargs):
    """Lowest order correction (Eq 16)"""
    # vs, V0: nm**3
    y = -0.5*vs**2/(kappa/bar * V0)
    return energyUnit(y, **kwargs)
