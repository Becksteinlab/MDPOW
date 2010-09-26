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

First copy the **computed results** , the ``pow.txt`` and ``energies.txt`` files that
are produced by :program:`mdpow-pow`, into the ``data`` directories.

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

Plot results and save to a pdf file with :func:`plot_exp_vs_comp`::

  mdpow.analysis.plot_exp_vs_comp(figname="figs/exp_vs_comp.pdf")

By default we also include the *SAMPL2* and reference (*Ref*) set results.

Remove the bad runs from ``pow.txt`` and save as ``pow_best.txt``. Then plot
again (this time excluding the SAMPL2 results)::

   pylab.clf()
   mdpow.analysis.plot_exp_vs_comp(data="data/run01/pow_best.txt", data2=None, figname='figs/run01/exp_vs_comp_best.pdf')

Note that the SAMPL2 results are excluded by setting ``data2=None``.


Functions
---------
.. autofunction:: plot_exp_vs_comp
.. autoclass:: ExpComp
.. autoclass:: ExpData
.. autoclass:: ComputedData
.. autofunction:: load_data
.. autofunction:: load_exp

"""

import numpy
import recsql
import logging
logger = logging.getLogger('mdpow.analysis')

#: Defaults paths to ``pow.txt`` for *experiments*, *run01*, *SAMPL2*, and *Ref*.
DEFAULTS_POW = {
    'experiments': "experimental/targets.csv",
    'run01': "data/run01/pow.txt",
    'SAMPL2': "data/SAMPL2/pow.txt",
    'Ref': "data/Ref/pow.txt",
    }

DEFAULTS_E = {
    'experiments': "experimental/targets.csv",
    'run01': "data/run01/energies.txt",
    'SAMPL2': "data/SAMPL2/energies.txt",
    'Ref': "data/Ref/energies.txt",
    }


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


def _plot_ideal(X=None, dy=0.5):
    from pylab import plot
    from numpy import array
    # plot perfect agreement (solid) and 'within 0.5 units'
    if X is None:
        X = array([-6,9])
    else:
        X = numpy.asarray(X)
    Y = X
    plot(X, Y-dy, 'k--', X, Y+dy, 'k--', X,Y,'k-')

def _finish(**kwargs):
    """Add axes labels and save to *figname*"""
    from pylab import xlabel, ylabel, savefig, xlim, ylim, gca
    import matplotlib as mpl

    figname = kwargs.pop('figname', None)

    _plot_ideal(dy=kwargs.pop('dy',0.5))
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
    return _finish(**kwargs)

class ExpComp(object):
    def __init__(self, **kwargs):
        """
        :Keywords:
           *experiments*
               path to ``targets.csv``
           *data*
               list of ``pow.txt`` paths; default are the files for 
               Ref, run01, SAMPL2 (stored in :data:`DEFAULTS_POW`)
        """

        filename = kwargs.pop('experiments', DEFAULTS_POW['experiments'])
        experimental = ExpData(filename=filename)

        filenames = kwargs.pop('data', [DEFAULTS_POW['Ref'], DEFAULTS_POW['run01'], DEFAULTS_POW['SAMPL2']])
        filename = filenames[0]
        # add to experimental db
        computed = ComputedData(filename=filename,
                                name="logPow_computed", 
                                connection=experimental.data.connection)

        for num,filename in enumerate(filenames[1:]):
            dbname = "logPow_compute_%d" % (num+1)
            compute2 = ComputedData(filename=filename,
                                    name=dbname, 
                                    connection=experimental.data.connection)  # add to experimental db
            # merge compute2 into compute (will be dropped after init) via __del__
            computed.data.merge_table(dbname)

        self.database = experimental
        
        # for debugging
        self._experimental = experimental
        self._computed = computed

    def plot(self, **kwargs):
        """Plot individual data points with legend.

        By default, the following should work:

        - Run from the top ``mdpow`` directory.
        - Prepare ``data/run01/pow.txt`` (must prepend header and append
          footer so that it is proper table). See :func:`load_data`.
        - Prepare ``experimental/targets.csv`` if it does not exist or if
          something changed. See :func:`load_exp` for details.

           *figname*
               write figure to *filename*; suffix determines file type
           *ymin*, *ymax*
               limits of the plot in the y direction (=computational results)
        """
        from pylab import figure, plot, legend, ylim, errorbar
        from matplotlib import colors, cm, rc
        import matplotlib

        # combined (matched on the itp_name!)
        c = self.database.data.SELECT(
            """CommonName AS name, directory AS comment, DeltaG0, """
            """mean,std,min,max,"""
            """__self__.logPow AS exp, C.logPow AS comp, C.errlogP AS errcomp""", 
            """LEFT JOIN logPow_computed AS C USING(itp_name) """
            """WHERE NOT (comp ISNULL OR C.itp_name ISNULL) """
            """GROUP BY comment ORDER BY no""")

        # need large figure for the plot
        matplotlib.rc('figure', figsize=kwargs.pop('figsize', (8,12)))
        # default font
        matplotlib.rc('font', size=10)

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
            label = r"%(comment)s %(exp).1f/%(comp).1f$\pm$%(errcomp).1f" % vars()
            plot(exp,comp, marker='o', markersize=14, color=color, markeredgewidth=0, alpha=0.1)
            plot(exp,comp, marker='o', markersize=5, color=color, label=label)
            xerr = numpy.abs(numpy.array([[xmin],[xmax]]) - exp)
            errorbar(exp,comp, xerr=xerr, yerr=errcomp, color=color, linewidth=1.5, capsize=0)

        legend(ncol=3, numpoints=1, loc='lower right', prop={'size':8})
        figname = _finish(**kwargs)

        matplotlib.rcdefaults()  # restore defaults
        return figname


def plot_exp_vs_comp(**kwargs):
    """Plot computed logPow against experimental values.

    Experimental values are stored in the reST table referenced byt
    the *experiments* keyword. *data* contains a list of ``pow.txt``
    tables for the calculated values.
    """
    expcomp = ExpComp(experiments=kwargs.pop("experiments",DEFAULTS_POW['experiments']),
                      data=kwargs.pop("data",[DEFAULTS_POW['Ref'], DEFAULTS_POW['run01'], DEFAULTS_POW['SAMPL2']]))
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
        kwargs['filename'] = filename
        self.filename = filename
        self.rawdata = recsql.SQLarray_fromfile(**kwargs)
        self.data = self.rawdata


# TODO: compute experimental Goct from Ghyd and logPow and compare to computed
#       Will show clearly for which compounds we should try simulating longer.

class GhydData(object):
    def __init__(self, exp=DEFAULTS_E['experiments'], 
                 data=[DEFAULTS_E['Ref'], DEFAULTS_E['run01'], DEFAULTS_E['SAMPL2']],
                 **kwargs):
        from mdpow import kBOLTZ

        temperature = kwargs.pop('temperature', 300.0)  # in K, used to calculate exp G_oct

        experimental = recsql.SQLarray_fromfile(exp)

        # add to experimental db
        computed = ComputedData(filename=data[0],
                                name="energies_computed", 
                                connection=experimental.connection)

        self.computed = computed  # for testing ... needs o be kept against gc :-p

        for num,filename in enumerate(data[1:]):
            dbname = "energies_computed_%d" % (num+1)
            compute2 = ComputedData(filename=filename,
                                    name=dbname, 
                                    connection=experimental.connection)  # add to experimental db
            # merge compute2 into compute (will be dropped after init) via __del__
            computed.data.merge_table(dbname)

        self.rawdb = experimental  # all together

        # Table for all hydration free energies with experimental values (unit: kJ/mol)
        # note: converting Ghyd(kcal/mol) into kJ/mol !!!
        # logPow = -(Goct-Ghyd)/kT * log10(e) ->   Goct = Ghyd - kT*logPow / log10(e)
        #
        # for error estimate I'd need the logPow error... but that requires _init_stats()...
        # so for a start I estimate it as 1 kcal/mol...
        kTlog10e = kBOLTZ*temperature/numpy.log10(numpy.e)
        self.database = self.rawdb.selection(
            "SELECT no, CommonName, Ghyd * 4.184 AS ExpGhyd, IFNULL(error_Ghyd, 0.0) * 4.184 AS ExpGhydErr, "
            "       Ghyd*4.184 - ? * E.logPow AS ExpGoct, E.logPow AS ExpLogPow, "
            "       sqrt(pow(IFNULL(error_Ghyd, 0.0)*4.184,2) + pow(4.184,2))  AS ExpGoctErr, "
            "       W.DeltaG0 AS CompGhyd, W.errDG0 AS CompGhydErr, "
            "       O.DeltaG0 AS CompGoct, O.errDG0 AS CompGoctErr, "
            "       W.directory AS comment "
            "FROM __self__ AS E LEFT JOIN energies_computed AS W ON itp_name = W.molecule "
            "                   LEFT JOIN energies_computed AS O ON itp_name = O.molecule "
            "WHERE W.solvent = 'water' AND O.solvent = 'octanol' AND NOT E.Ghyd ISNULL",
            (kTlog10e,))


    def plot(self, mode, **kwargs):
        """Plot individual data points with legend.

           *mode*
               "hyd" or "oct"
           *figname*
               write figure to *filename*; suffix determines file type
           *ymin*, *ymax*
               limits of the plot in the y direction (=computational results)
        """
        from pylab import figure, plot, legend, ylim, errorbar, xlabel, ylabel, savefig
        from matplotlib import colors, cm, rc
        import matplotlib

        figname = kwargs.pop('figname', None)

        # need large figure for the plot
        matplotlib.rc('figure', figsize=kwargs.pop('figsize', (8,12)))
        # default font
        matplotlib.rc('font', size=10)

        mode = mode.lower()
        if not mode in ("hyd","oct"):
            raise ValueError("mode must be either 'hyd' or 'oct', not %r" % mode)
        fields = "CommonName, comment, ExpG%(mode)s AS ExpG, ExpG%(mode)sErr AS ExpErr, CompG%(mode)s AS CompG, CompG%(mode)sErr AS CompErr" % vars()
        c = self.database.SELECT(fields)
        ExpG = "ExpG%(mode)s" % vars()

        norm = colors.normalize(0,len(c))
        for i, (name,comment,exp,errexp,comp,errcomp) in enumerate(c):
            if exp is None or comp is None:  # should not be necessary...
                continue

            color = cm.jet(norm(i))
            label = r"%(comment)s %(exp).1f/%(comp).1f$\pm$%(errcomp).1f" % vars()
            plot(exp,comp, marker='o', markersize=14, color=color, markeredgewidth=0, alpha=0.1)
            plot(exp,comp, marker='o', markersize=5, color=color, label=label)
            errorbar(exp,comp, xerr=errexp, yerr=errcomp, color=color, linewidth=1.5, capsize=0)

        legend(ncol=3, numpoints=1, loc='lower right', prop={'size':8})

        # 1 kcal/mol = 4.184 kJ/mol band
        _plot_ideal(X=self.database.limits(ExpG), dy=kwargs.pop('dy', 4.184))
        
        xlabel(r'experimental $\Delta G_{\rm %(mode)s}$ (kJ/mol)' % vars())
        ylabel(r'computed $\Delta G_{\rm %(mode)s}$ (kJ/mol)' % vars())
        if figname:
            savefig(figname)
            logger.info("Wrote plot to %r", figname)

        matplotlib.rcdefaults()  # restore defaults
        return figname
