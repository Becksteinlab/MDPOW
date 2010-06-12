# POW
# Copyright (c) 2010 Oliver Beckstein and Bogdan Iorga

"""
:mod:`analysis` --- Collection of analysis and plotting functions
=================================================================

Simple functions to quickly plot data. Typically, it works best if ran
interactively from the top level of the POW directory!

Experimental values are loaded from the targets list (``targets.csv``)
and computed values from the table in ``pow.txt``. See
:func:`plot_exp_vs_comp` for details.

.. Note: The header templates in :file:`lib/templates/_pow_header.txt` and
         :file:`lib/templates/_energies_header.txt` *define* the column names
         that are used in the code here, in particular in :class:`ComputeData`
         and :class:`ExpComp`.

Usage
-----

Plot results and save to a pdf file::

  mdpow.analysis.plot_exp_vs_comp(figname="figs/exp_vs_comp.pdf")

Remove the bad runs from  ``pow.txt`` and save as ``pow_best.txt``. Then plot again::

   pylab.clf()
   mdpow.analysis.plot_exp_vs_comp(data="data/run01/pow_best.txt", figname='figs/run01/exp_vs_comp_best.pdf')

By default we also include the SAMPL2 results,

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

DEFAULTS = {
    'experiments': "experimental/targets.csv",
    'run01': "data/run01/pow.txt",
    'SAMPL2': "data/SAMPL2/pow.txt",
    }

def load_data(filename=DEFAULTS['run01'], **kwargs):
    """Load computed POW table and return :class:`recsql.SQLarray`.

    The data file is typically ``pow.txt``. It *must* contain a proper
    reST table. Use the ``_header`` and ``_footer`` files if you only
    have the raw output from :program:`mdpow-pow`.

    Furthermore, the column names are important because we use them
    here.
    """
    kwargs['filename'] = filename
    return recsql.SQLarray_fromfile(**kwargs)

def load_exp(filename=DEFAULTS['experiments'], **kwargs):
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
           *data*, *data2*
               path to ``pow.txt``
        """
        # data and data2 are quick hacks to load both SAMPL2/pow.txt and run01/pow.txt
        experimental = ExpData(filename=kwargs.pop('experiments', DEFAULTS['experiments']))
        computed = ComputedData(filename=kwargs.pop('data', DEFAULTS['run01']),
                                name="logPow_computed", 
                                connection=experimental.data.connection)  # add to experimental db
        compute2 = ComputedData(filename=kwargs.pop('data2', DEFAULTS['SAMPL2']),
                                name="logPow_SAMPL2", 
                                connection=experimental.data.connection)  # add to experimental db
        # merge compute2 into compute (will be dropped after init) via __del__
        computed.data.merge_table("logPow_SAMPL2")

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
            """LEFT JOIN logPow_computed AS C using(itp_name) """
            """WHERE NOT (comp ISNULL OR C.itp_name ISNULL) """
            """GROUP BY comment ORDER BY no""")

        # need large figure for the plot
        matplotlib.rc('figure', figsize=kwargs.pop('figsize', (8,10)))
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
            label = "%(comment)s %(exp).1f/%(comp).1f" % vars()
            plot(exp,comp, marker='o', markersize=14, color=color, markeredgewidth=0, alpha=0.3)
            plot(exp,comp, marker='o', markersize=5, color=color, label=label)
            xerr = numpy.abs(numpy.array([[xmin],[xmax]]) - exp)
            errorbar(exp,comp, xerr=xerr, yerr=errcomp, color=color, linewidth=2, capsize=0)

        legend(ncol=3, numpoints=1, loc='lower right', prop={'size':8})
        figname = _finish(**kwargs)

        matplotlib.rcdefaults()  # restore defaults
        return figname


def plot_exp_vs_comp(**kwargs):
    expcomp = ExpComp(experiments=kwargs.pop("experiments",DEFAULTS['experiments']),
                      data=kwargs.pop("data",DEFAULTS['run01']),
                      data2=kwargs.pop("data2",DEFAULTS['SAMPL2']))
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
    def __init__(self, filename=DEFAULTS['experiments'], **kwargs):
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
        self.stats()  # generate stats and load statistics{}
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

    def stats(self):
        """Statistics calculated from the experimental values.

        :Arguments: *exp* is the database loaded with :func:`load_exp`.

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

    def __init__(self, filename=DEFAULTS['run01'], **kwargs):
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
