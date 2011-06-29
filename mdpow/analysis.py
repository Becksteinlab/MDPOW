# POW
# Copyright (c) 2010 Oliver Beckstein and Bogdan Iorga

"""
:mod:`mdpow.analysis` --- Collection of analysis and plotting functions
=======================================================================

Simple functions to quickly plot data. Typically, it works best if ran
interactively from the top level of the POW directory!

Experimental values are loaded from the targets list (``targets.csv``)
and computed values from the table in ``pow.txt``. See
:func:`plot_exp_vs_comp` for details.

.. Note: The header templates in :file:`lib/templates/_pow_header.txt` and
         :file:`lib/templates/_energies_header.txt` *define* the column names
         that are used in the code here, in particular in :class:`ComputeData`
         and :class:`ExpComp`.

Prepare data
------------

First copy the **computed results** , the ``pow.txt`` and ``energies.txt``
files that are produced by :program:`mdpow-pow`, into the ``data`` directories.

Then format them::

   ./lib/scripts/make_tables.sh data/*

The **experimental data** are taken from *targets.numbers*. In
:program:`numbers` export the table as *UTF-8* in *CSV* format to
``experimental/targets.csv``. This is only necessary if the experimental data
were changed. We only plot entries for which

 - a id number (first column *no*) is defined (should be unique)
 - a *logPow* value exists
 - a *itp_name* exists, which must correspond to the *molecule* name used in
   :mod:`mdpow.fep.Gsolv`


Making graphs
-------------

For the following, import the module::

  import mdpow.analysis


Octanol-water partition coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Plot results and save to a pdf file with :func:`plot_exp_vs_comp`::

  mdpow.analysis.plot_exp_vs_comp(figname="figs/logPow.pdf")

By default we also include the *SAMPL2* and reference (*Ref*) set
results. In practice some manual adjustments are required, e.g. ::

  mdpow.analysis.plot_exp_vs_comp()
  # resize window so that (huge!) legend fits
  ylim(-22,10)
  savefig("figs/logPow.pdf")

Using a file named ``exclusions.txt`` in the same directory as the
data file, one can exclude certain runs from appearing in the graph:
set the *exclusions* keyword to ``True``::

  pylab.clf()
  mdpow.analysis.plot_exp_vs_comp(exclusions=True, figname='figs/logPow_best.pdf')

In practice, manual fiddling is required such as resizing the graph::

  mdpow.analysis.plot_exp_vs_comp(exclusions=True)
  ylim(-14,10)
  savefig('figs/logPow_best.pdf')

The ``exclusions.txt`` files must contain a table such as ::

            Table[exclusions]: These sims are ignored.
            ======== ===========================================
            itp_name directory_regex
            ======== ===========================================
            AXX      .*
            5FH      benzylhyd$
            ======== ===========================================

then any simulation of a *molecule* equalling *itp_name* and which is
recorded with a *directory* matching the regular expression
*directory_regex* will be excluded from the analysis.


Solvation energies
~~~~~~~~~~~~~~~~~~

Plots that compare experimental hydration and octanol-solvation free
energies to computed values. DeltaG_hyd are only available for a few
compounds so we only plot a subset of all the compounds that we have
done.

Experimental octanol solvation free energies are computed from
experimental logPow and DeltaGhyd from

  logPow = -(DeltaGoct-DeltaGhyd)/kT * log10(e)

(see also :func:`gsolv2logpow`) as

  DeltaGoct = DeltaGhyd - kT*logPow / log10(e)

.. Warning:: In principle the temperature T of the logPow measurement and the
   DeltaGhyd measurement must be the same. For the time being we just *assume*
   that both were done at T=300K. Also, the error on DeltaGoct is not
   calculated properly at the moment, either, because the error on logPow is
   hard to quantify (based on the logKow_ data). We are estimating the error on
   kT*logPow/log10(e) error as 0.5 kcal/mol and combine it with the known
   experimental error for DeltaGhyd.

Plotting uses the :meth:`GsolvData.plot` method from :class:`GsolvData`::

   from pylab import *
   clf()
   G = mdpow.analysis.GsolvData()
   G.plot('hyd')
   # adjust things such as manually increasing window ...
   ylim(-180,20)
   savefig("figs/ghyd.pdf")

   clf()
   G.plot('oct')
   ylim(-150,20)
   savefig("figs/goct.pdf")

Right now, the plots are a bit messy but I opted to include the legend to make
it easier for us to understand the data. I had to manually increase the
plotting window to make things fit properly.

:class:`GsolvData` also honours the *exlusions* = ``True`` keyword argument.

.. _logKow: http://logkow.cisti.nrc.ca/logkow/

Functions
---------
.. autofunction:: plot_exp_vs_comp
.. autoclass:: GsolvData
   :members:
.. autoclass:: ExpComp
   :members:
.. autoclass:: ExpData
   :members:
.. autoclass:: ComputedData
   :members:
.. autofunction:: load_data
.. autofunction:: load_exp
.. autofunction:: gsolv2logpow
"""
import os.path
import numpy
import recsql
import logging
logger = logging.getLogger('mdpow.analysis')

from mdpow import kBOLTZ

#: Default paths to ``pow.txt`` for *experiments*, *run01*, *SAMPL2*, and *Ref*.
DEFAULTS_POW = {
    'experiments': "experimental/targets.csv",
    'run01': "data/run01/pow.txt",
    'SAMPL2': "data/SAMPL2/pow.txt",
    'Ref': "data/Ref/pow.txt",
    }

#: Default paths to ``energies.txt`` for *experiments*, *run01*, *SAMPL2*, and *Ref*.
DEFAULTS_E = {
    'experiments': "experimental/targets.csv",
    'run01': "data/run01/energies.txt",
    'SAMPL2': "data/SAMPL2/energies.txt",
    'Ref': "data/Ref/energies.txt",
    }

def gsolv2logpow(Gwat, Goct, unit='kcal/mol', temperature=300.):
    """Calculate logPow from the solvation free energies.

        logPow = -(Goct-Ghyd)/kT * log10(e)

    .. Note:: Default unit is kcal/mol, unlike the rest of mdpow, which
       uses kJ/mol. The reason is that most solvation free energies in
       the literature are quoted in kcal/mol.

    :Arguments:
       *Gwat*
           hydration free energy
       *Goct*
           octanol solvation free energy
       *temperature*
           temperature in K [300]
       *unit*
           unit of the energies, either "kcal/mol" or "kJ/mol";
           ["kcal/mol"]
    """
    if unit == 'kcal/mol':
        Gwat *= 4.184
        Goct *= 4.184
    return -(Goct-Gwat)/(kBOLTZ*temperature) * numpy.log10(numpy.e)

def load_data(filename=DEFAULTS_POW['run01'], **kwargs):
    """Load computed POW table and return :class:`recsql.SQLarray`.

    The data file is typically ``pow.txt``. It *must* contain a proper
    reST table. Use the ``_header`` and ``_footer`` files if you only
    have the raw output from :program:`mdpow-pow`.

    Furthermore, the column names (defined in the header and footer
    files!) are important because we use them here.
    """
    kwargs['filename'] = filename
    return recsql.SQLarray_fromfile(**kwargs)

def load_exp(filename=DEFAULTS_POW['experiments'], **kwargs):
    """Load experimental values table and return :class:`recsql.SQLarray`.

    To obtain ``targets.csv`` export ``targets.numbers`` in
    :program:`Numbers` as **UTF-8** (important!) in the **CSV* format
    (File->Export)
    """
    kwargs['filename'] = filename
    return recsql.SQLarray_fromfile(**kwargs)


def _plot_ideal(X, dy=0.5, dy2=None, padding=0.05, filled=True):
    """Plot lines perfect correlation (solid) and 'within *dy* units'.

    :param X:   array of *x* values to plot at
    :param dy:  plot dashed line a +/-dy of the ideal
    :param dy2: plot dotted line a +/-dy2 of the ideal
    :param filled:  ``True`` use fancy alpha-blended solid bars, ``False`` uses dashed and dotted lines
    :param padding: percentage of total extent in *X* that is added left and right
    """
    from pylab import plot, fill_between
    from numpy import array

    X = numpy.asarray(X)
    if len(X) < 2:
        raise ValueError("X must contain at least 2 values.")
    # add buffer
    if padding:
        dx = padding * numpy.abs(X[-1]-X[0])
        X[0] -= dx
        X[-1] += dx

    Y = X

    if filled:
        alpha = 0.1
        if dy2:
            fill_between(X, Y-dy2, Y+dy2, color="black", alpha=alpha)
        else:
            alpha = 0.2
        fill_between(X, Y-dy, Y+dy, color="black", alpha=alpha)
    else:
        plot(X, Y-dy, 'k--', X, Y+dy, 'k--')
        if dy2:
            plot(X, Y-dy2, 'k.', X, Y+dy2, 'k.')

    plot(X,Y,'k-')

def _finish(X, **kwargs):
    """Add axes labels and save to *figname*"""
    from pylab import xlabel, ylabel, savefig, xlim, ylim, gca
    import matplotlib as mpl

    figname = kwargs.pop('figname', None)

    _plot_ideal(X, dy=kwargs.pop('dy',0.5))
    xlim(kwargs.pop('xmin', -6), kwargs.pop('xmax',10))
    ylim(kwargs.pop('ymin', -10), kwargs.pop('ymax',10))

    ax = gca()
    locator = mpl.ticker.MultipleLocator(base=kwargs.pop('base',1.0))
    ax.xaxis.set_major_locator(locator)
    ax.yaxis.set_major_locator(locator)

    xlabel(r'experimental $\log P_{\rm OW}$')
    ylabel(r'computed $\log P_{\rm OW}$')
    if figname:
        savefig(figname)
        logger.info("Wrote plot to %r", figname)
    return figname


def plot_quick(a, **kwargs):
    """Quick hack to plot all exp_vs_comp data together."""
    from pylab import plot
    from numpy import array
    kwargs.setdefault('ylim',None) # ??
    kwargs.setdefault('ymax',None) # ??
    plot(a.recarray.exp, a.recarray.logP_OW, 'ro', **kwargs)  # data
    return _finish([-6,9], **kwargs)  # default extent for our logPow :-p

class ExpComp(object):
    """Database combining experimental with computed values."""

    def __init__(self, **kwargs):
        """
        :Keywords:
           *experiments*
               path to ``targets.csv``
           *data*
               list of ``pow.txt`` paths; default are the files for
               Ref, run01, SAMPL2 (stored in :data:`DEFAULTS_POW`)
           *exclusions*
               ``False``
                   Does nothing special.
               ``True``
                   Look for `exclusions.txt` in same directory as each data file.
                   If it contains a table such as::

                      Table[exclusions]: These sims are ignored.
                      ======== ===========================================
                      itp_name directory_regex
                      ======== ===========================================
                      AXX      .*
                      5FH      benzylhyd$
                      ======== ===========================================

                   then any simulation of a *molecule* equalling *itp_name*
                   and which is recorded with a *directory* matching the
                   regular expression *directory_regex* will be excluded from the analysis.
                   [``False``]
        """

        filename = kwargs.pop('experiments', DEFAULTS_POW['experiments'])
        self.exclusions = kwargs.pop('exclusions', False)
        data = kwargs.pop('data', [DEFAULTS_POW['Ref'], DEFAULTS_POW['run01'], DEFAULTS_POW['SAMPL2']])

        experimental = ExpData(filename=filename).data

        filename = data[0]
        # add to experimental db
        energies_name = "logPow_computed"
        exclusions_name = None
        computed = ComputedData(filename=filename,
                                name=energies_name,
                                connection=experimental.connection)

        for num,filename in enumerate(data[1:]):
            dbname = "logPow_computed_%d" % (num+1)
            compute2 = ComputedData(filename=filename,
                                    name=dbname,
                                    connection=experimental.connection)  # add to experimental db
            # merge compute2 into compute (will be dropped after init) via __del__
            computed.data.merge_table(dbname)

        if self.exclusions:
            # exclusion.txt files are not guaranteed to exist, so remember the first found:
            first_excl = None
            for num,filename in enumerate(data):
                exclusions = os.path.join(os.path.dirname(filename), 'exclusions.txt')
                if not os.path.exists(exclusions):
                    continue
                logger.info("Loading exclusions from %(exclusions)r.", vars())
                dbname = "exclusions_%d" % num
                excl = recsql.SQLarray_fromfile(exclusions, connection=experimental.connection, name=dbname)
                if first_excl:
                    first_excl.merge_table(dbname)
                else:
                    first_excl = excl
                    exclusions_name = dbname

        self.rawdb = experimental

        if self.exclusions and not exclusions_name is None:
            energies = self.rawdb.selection(
                """SELECT * FROM %(energies_name)s """
                """WHERE NOT directory IN """
                """(SELECT directory FROM %(energies_name)s """
                """LEFT JOIN %(exclusions_name)s USING (itp_name) """
                """WHERE MATCH(directory_regex, directory))""" % vars())
            energies_name = energies.name     # switch to the reduced table

        # combined (matched on the itp_name!)
        self.database = self.rawdb.selection(
            "SELECT "
            """CommonName AS name, directory AS comment, DeltaG0, """
            """mean,std,min,max,"""
            """E.logPow AS exp, C.logPow AS comp, C.errlogP AS errcomp """
            """FROM __self__ AS E """
            """LEFT JOIN %(energies_name)s AS C USING(itp_name) """
            """WHERE NOT (comp ISNULL OR C.itp_name ISNULL) """
            """GROUP BY comment ORDER BY no""" % vars())


    def plot(self, **kwargs):
        """Plot individual data points with legend.

        By default, the following should work:

        - Run from the top ``mdpow`` directory.
        - Prepare ``data/run01/pow.txt`` (must prepend header and append
          footer so that it is proper table). See :func:`load_data`.
        - Prepare ``experimental/targets.csv`` if it does not exist or if
          something changed. See :func:`load_exp` for details.

           *compoundnames*
               ``False`` puts the directory names in the legend, ``True`` uses the chemical
               compound names; "auto" chooses ``True`` if exclusions were applied ["auto"]
           *figname*
               write figure to *filename*; suffix determines file type
           *ymin*, *ymax*
               limits of the plot in the y direction (=computational results)
        """
        from pylab import figure, subplot, plot, legend, ylim, errorbar
        from matplotlib import colors, cm, rc
        import matplotlib

        # need large figure for the plot
        matplotlib.rc('figure', figsize=kwargs.pop('figsize', (8,12)))
        # default font
        matplotlib.rc('font', size=10)

        compoundnames = kwargs.pop("compoundnames", "auto")
        if compoundnames == "auto":
            if self.exclusions:
                compoundnames = True
            else:
                compoundnames = False
        if compoundnames:
            legendformat = r"%(name)s %(exp).1f/%(comp).1f$\pm$%(errcomp).1f"
        else:
            legendformat = r"%(comment)s %(exp).1f/%(comp).1f$\pm$%(errcomp).1f"

        # explicitly list columns so that we can safely expand the db with new columns
        c = self.database.SELECT("name,comment,DeltaG0,mean,std,min,max,exp,comp,errcomp")

        #subplot(121)
        norm = colors.normalize(0,len(c))
        for i, (name,comment,DeltaA0,xmean,xstd,xmin,xmax,exp,comp,errcomp) in enumerate(c):
            if exp is None or comp is None:
                continue

            # fix possibly missing values
            if xmean is None:
                xmean = comp
            if xstd is None:
                xstd = 0
            if xmin is None:
                xmin = xmean
            if xmax is None:
                xmax = xmean

            color = cm.jet(norm(i))
            label = legendformat % vars()
            #do not plot the aureola; pretty but does not contain information
            #plot(exp,comp, marker='o', markersize=14, color=color, markeredgewidth=0, alpha=0.1)
            plot(exp,comp, marker='o', markersize=5, color=color, label=label)
            xerr = numpy.abs(numpy.array([[xmin],[xmax]]) - exp)
            errorbar(exp,comp, xerr=xerr, yerr=errcomp, color=color, linewidth=1.5, capsize=0)

        # legend outside figure with bbox_to_anchor instead of  loc='lower right'
        #legend(ncol=3, numpoints=1, prop={'size':6}, bbox_to_anchor=(1.05,1), borderaxespad=0.)
        legend(ncol=3, numpoints=1, prop={'size':6}, mode="expand", loc='lower right')

        figname = _finish(self.database.limits('exp'), **kwargs)

        matplotlib.rcdefaults()  # restore defaults
        return figname


def plot_exp_vs_comp(**kwargs):
    """Plot computed logPow against experimental values from default files.

    Experimental values are stored in the reST table referenced byt
    the *experiments* keyword. *data* contains a list of ``pow.txt``
    tables for the calculated values.
    """
    logger.info("Plotting logPow.")
    expcomp = ExpComp(exclusions=kwargs.pop('exclusions',False))
    return expcomp.plot(**kwargs)

def unpackCSlist(s, convertor=float):
    """Unpack a comma-separated list in string form."""
    try:
        return map(convertor, s.split(','))
    except:
        return []


class ExpData(object):
    """Object that represents our experimental data.

    Access the raw data via :attr:`ExpData.rawdata` and a table enriched
    with statistics as :attr:`ExpData.data` (which is a
    :class:`recsql.SQLarray`).
    """
    def __init__(self, filename=DEFAULTS_POW['experiments'], **kwargs):
        """Load experimental values table and return :class:`recsql.SQLarray`.

        To obtain ``targets.csv`` export ``targets.numbers`` in
        :program:`Numbers` as **UTF-8** (important!) in the **CSV* format
        (File->Export)
        """
        logger.info("Loading experimental data from %(filename)r", vars())
        kwargs['filename'] = filename
        self.filename = filename
        # load original data from the targets list and set up database
        self.rawdata = recsql.SQLarray_fromfile(**kwargs)
        self.statistics = {}
        self._init_stats()  # generate stats and load statistics{}
        # add 'stats' table to the main database (via connection)
        # (must be kept around in self or the table gets gc/del'ed immediately again)
        self.__statsdb = recsql.SQLarray(
            name='__stats', records=list(self._stats_generator()),
            columns=self._stats_columns, connection=self.rawdata.connection)
        # generate another table '__experiments' in the database
        self.data = self.rawdata.selection(
            """SELECT no,name,CommonName,itp_name,CAS_RN,logPow,mean,std,min,max """
            """FROM __self__ NATURAL LEFT JOIN __stats""",
            name="__experiments")

        # access with
        #   self.data.SELECT('*')
        # or
        #   self.rawdata.selection('SELECT * FROM __experiments')

    def _init_stats(self):
        """Statistics calculated from the experimental values.

        Sets :attr:`ExpData._stats_columns` and
        :attr:`ExpData._stats_generator`.

        :Returns: dictionary; each entry is a numpy.array with the order
                  corresponding to the order of compounds listed n the key
                  'names'
         """
        def int_or_zero(n):
            try:
                return int(n)
            except TypeError:
                return 0

        r = {}  # results
        s = self.rawdata.SELECT("no,itp_name,logPow,other_logPow",
                                "WHERE NOT (itp_name ISNULL OR no ISNULL OR logPow ISNULL)")
        # explicit cast to a python int: work around broken sqlite... int64 not good?!
        # int_or_zero(): deal with incomplete (empty) data column
        #r['number'] = map(int_or_zero, s.no)
        r['number'] = map(int, s.no)
        r['itp_names'] = s.itp_name
        r['logPow'] = s.logPow.astype(float)
        # sort other values and add the main value to the list
        # (note: indices correspond to columns in SELECT!)
        all = [numpy.sort([row[2]] + unpackCSlist(row[3])) for row in s]
        r['mean'] = numpy.array([a.mean() for a in all])
        r['std'] =  numpy.array([a.std() for a in all])
        r['min'] = numpy.array([a.min() for a in all])
        r['max'] = numpy.array([a.max() for a in all])

        self.statistics.update(r)
        self.statistics['other_logPow'] = all

        # iterator for other things
        columns = r.keys()
        columns.sort()
        def records_generator():
            numrows = len(r['itp_names'])
            lengths = numpy.array([len(r[c]) for c in columns])
            if not numpy.all(lengths == numrows):
                raise TypeError("entries of dict r have different lengths")
            # could add a fix_type() wrapper-hack around each number that checks
            # if the type is sqlite-incompatible and then casts to a pure python type
            for i in xrange(numrows):
                yield tuple([r[c][i] for c in columns])
        self._stats_columns = columns
        self._stats_generator = records_generator
        return r

class ComputedData(object):
    """Object that represents computed data.

    Access via :attr:`ComputedData.data`.
    """

    def __init__(self, filename=DEFAULTS_POW['run01'], **kwargs):
        """Load computed POW table and return :class:`recsql.SQLarray`.

        The data file is typically ``pow.txt``. It *must* contain a proper
        reST table. Use the ``_header`` and ``_footer`` files if you only
        have the raw output from :program:`mdpow-pow`.

        Furthermore, the column names are important because we use them
        here.
        """
        logger.info("Loading computed data from %(filename)r.", vars())
        kwargs['filename'] = filename
        self.filename = filename
        self.rawdata = recsql.SQLarray_fromfile(**kwargs)
        self.data = self.rawdata


class GsolvData(object):
    """Solvation energies organized as a database.

    Plot either "hyd" or "oct" with :meth:`GsolvData.plot`.
    """

    def __init__(self, exp=DEFAULTS_E['experiments'], **kwargs):
        """Load experimental and simulation data.

        The defaults load all the data generated in the project so
        far. See :data:`DEFAULTS_E` in the python code.

        :Keywords:
         *exp*
            path to the experimental ``targets.csv`` file
         *data*
            list of simulation results (typically stored as reST ``energies.txt``).
         *temperature*
            temperature at which experimental measurements of logPow
            were presumed to be taken; needed for the calculations of
            the experimental DeltaG_oct from logPow and experimental [300.0]
            DeltaG_hyd.
         *exclusions*
            ``False`` does nothing special.
            ``True``: look for `exclusions.txt` in same directory as each data file.
            If it contains a table such as::
               Table[exclusions]: These sims are ignored.
               ======== ===========================================
               itp_name directory_regex
               ======== ===========================================
               AXX      .*
               5FH      benzylhyd$
               ======== ===========================================

            then any simulation of a *molecule* equalling *itp_name*
            and which is recorded with a *directory* matching the
            regular expression *directory_regex* will be excluded from the analysis.
            [``False``]
        """
        data = kwargs.pop('data', [DEFAULTS_E['Ref'], DEFAULTS_E['run01'], DEFAULTS_E['SAMPL2']])
        temperature = kwargs.pop('temperature', 300.0)  # in K, used to calculate exp G_oct
        self.exclusions = kwargs.pop('exclusions', False)

        logger.info("Loading experimental data from %(exp)r.", vars())
        experimental = recsql.SQLarray_fromfile(exp)

        # add to experimental db
        energies_name = "energies_computed"
        exclusions_name = None
        computed = ComputedData(filename=data[0],
                                name=energies_name,
                                connection=experimental.connection)

        self.computed = computed  # for testing ... needs to be kept against gc :-p

        for num,filename in enumerate(data[1:]):
            dbname = "energies_computed_%d" % (num+1)
            compute2 = ComputedData(filename=filename,
                                    name=dbname,
                                    connection=experimental.connection)  # add to experimental db
            # merge compute2 into compute (will be dropped after init) via __del__
            computed.data.merge_table(dbname)

        if self.exclusions:
            # exclusion.txt files are not guaranteed to exist, so remember the first found:
            first_excl = None
            for num,filename in enumerate(data):
                exclusions = os.path.join(os.path.dirname(filename), 'exclusions.txt')
                if not os.path.exists(exclusions):
                    continue
                logger.info("Loading exclusions from %(exclusions)r.", vars())
                dbname = "exclusions_%d" % num
                excl = recsql.SQLarray_fromfile(exclusions, connection=experimental.connection, name=dbname)
                if first_excl:
                    first_excl.merge_table(dbname)
                else:
                    first_excl = excl
                    exclusions_name = dbname

        self.rawdb = experimental  # all together

        if self.exclusions and not exclusions_name is None:
            # NOTE: change molecule --> itp_name so that I can use USING (itp_name)
            energies = self.rawdb.selection(
                """SELECT * FROM %(energies_name)s """
                """WHERE NOT directory IN """
                """(SELECT directory FROM %(energies_name)s """
                """LEFT JOIN %(exclusions_name)s ON itp_name = molecule """
                """WHERE MATCH(directory_regex, directory))""" % vars())
            energies_name = energies.name     # switch to the reduced table


        # Table for all hydration free energies with experimental values (unit: kJ/mol)
        # note: converting original Ghyd (in kcal/mol) into kJ/mol by multiplying with 4.184 kJ/kcal !!!
        # logPow = -(Goct-Ghyd)/kT * log10(e) ->   Goct = Ghyd - kT*logPow / log10(e)
        #
        # for error estimate I'd need the logPow error... but that requires _init_stats()...
        # so for a start I estimate the experimental logPow error as 0.5 log10 units.
        estimated_error_logPow = 0.5  # kcal/mol ESTIMATE!!!
        logger.warning("Errors on Delta G_oct use the estimated logPow error %(estimated_error_logPow)f (log10 units)!!! Needs to be done properly.", vars())
        kTlog10e = kBOLTZ*temperature/numpy.log10(numpy.e)  # in kJ/mol !!
        self.database = self.rawdb.selection(
            "SELECT no, CommonName, Ghyd * 4.184 AS ExpGhyd, IFNULL(error_Ghyd, 0.0) * 4.184 AS ExpGhydErr, "
            "       Ghyd*4.184 - ? * E.logPow AS ExpGoct, E.logPow AS ExpLogPow, "
            "       sqrt(pow(IFNULL(error_Ghyd, 0.0)*4.184,2) + pow(?*?,2))  AS ExpGoctErr, "
            "       W.DeltaG0 AS CompGhyd, W.errDG0 AS CompGhydErr, "
            "       O.DeltaG0 AS CompGoct, O.errDG0 AS CompGoctErr, "
            "       W.directory AS comment "
            "FROM __self__ AS E LEFT JOIN %(energies_name)s AS W ON itp_name = W.molecule "
            "                   LEFT JOIN %(energies_name)s AS O ON itp_name = O.molecule AND W.directory = O.directory "
            "WHERE W.solvent = 'water' AND O.solvent = 'octanol' AND NOT E.Ghyd ISNULL" % vars(),
            (kTlog10e, kTlog10e, estimated_error_logPow))


    def plot(self, mode, **kwargs):
        """Plot individual data points with legend.

           *mode*
               "hyd" or "oct"
           *compoundnames*
               ``False`` puts the directory names in the legend, ``True`` uses the chemical
               compound names; "auto" chooses ``True`` if exclusions were applied ["auto"]
           *figname*
               write figure to *filename*; suffix determines file type
           *ymin*, *ymax*
               limits of the plot in the y direction (=computational results)
        """
        from pylab import figure, plot, subplot, legend, ylim, errorbar, xlabel, ylabel, savefig
        from matplotlib import colors, cm, rc
        import matplotlib

        figname = kwargs.pop('figname', None)
        compoundnames = kwargs.pop("compoundnames", "auto")
        if compoundnames == "auto":
            if self.exclusions:
                compoundnames = True
            else:
                compoundnames = False
        if compoundnames:
            legendformat = r"%(name)s %(exp).1f/%(comp).1f$\pm$%(errcomp).1f"
        else:
            legendformat = r"%(comment)s %(exp).1f/%(comp).1f$\pm$%(errcomp).1f"

        # need large figure for the plot
        matplotlib.rc('figure', figsize=kwargs.pop('figsize', (10.5,12)))
        # default font
        matplotlib.rc('font', size=10)
        matplotlib.rc('mathtext', fontset='stixsans') # does not work?

        mode = mode.lower()
        if not mode in ("hyd","oct"):
            raise ValueError("mode must be either 'hyd' or 'oct', not %r" % mode)
        fields = "CommonName, comment, ExpG%(mode)s AS ExpG, ExpG%(mode)sErr AS ExpErr, CompG%(mode)s AS CompG, CompG%(mode)sErr AS CompErr" % vars()
        c = self.database.SELECT(fields)
        ExpG = "ExpG%(mode)s" % vars()

        #subplot(121)
        norm = colors.normalize(0,len(c))
        for i, (name,comment,exp,errexp,comp,errcomp) in enumerate(c):
            if exp is None or comp is None:  # should not be necessary...
                continue

            color = cm.jet(norm(i))
            label = legendformat % vars()
            #do not plot the aureola; pretty but does not contain information
            #plot(exp,comp, marker='o', markersize=14, color=color, markeredgewidth=0, alpha=0.1)
            plot(exp,comp, marker='o', markersize=5, color=color, label=label)
            errorbar(exp,comp, xerr=errexp, yerr=errcomp, color=color, linewidth=1.5, capsize=0)

        # legend outside figure with bbox_to_anchor instead of  loc='lower right'
        #legend(ncol=3, numpoints=1, prop={'size':6}, bbox_to_anchor=(1.05,1), borderaxespad=0.)
        legend(ncol=3, numpoints=1, prop={'size':6}, mode="expand", loc='lower right')

        # 1 kcal/mol = 4.184 kJ/mol band
        kcalmol = 4.184
        dy = kwargs.pop('dy', kcalmol)
        _plot_ideal(X=self.database.limits(ExpG), dy=dy, dy2=1.5*dy)

        xlabel(r'experimental $\Delta G_{\rm %(mode)s}$ (kJ/mol)' % vars())
        ylabel(r'computed $\Delta G_{\rm %(mode)s}$ (kJ/mol)' % vars())
        if figname:
            savefig(figname)
            logger.info("Wrote plot to %r", figname)

        matplotlib.rcdefaults()  # restore defaults
        return figname
