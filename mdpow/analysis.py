# POW
# Copyright (c) 2010 Oliver Beckstein and Bogdan Iorga

"""
:mod:`analysis` --- Collection of analysis and plotting functions
=================================================================

Simple functions to quickly plot data. Typically, it works best if ran
interactively from the top level of the POW directory!

Experimental values are loaded from the targets list (``targets.csv``)
and computed values from tthe table in ``pow.txt``. See
:func:`plot_exp_vs_comp` for details.

Usage
-----

Plot results and save to a pdf file::
   mdpow.analysis.plot_exp_vs_comp(figname="figs/run01/exp_vs_comp.pdf")

Remove the bad runs from  ``pow.txt`` and save as ``pow_best.txt``. Then plot again::
   pylab.clf()
   mdpow.analysis.plot_exp_vs_comp(data="data/run01/pow_best.txt", figname='figs/run01/exp_vs_comp_best.pdf')

Functions
---------

.. autofunction:: load_data
.. autofunction:: load_exp
.. autofunction:: plot_exp_vs_comp
"""

import recsql
import logging
logger = logging.getLogger('mdpow.analysis')

DEFAULTS = {
    'experiments': "experimental/targets.csv",
    'data': "data/run01/pow.txt",
    }

def load_data(filename=DEFAULTS['data'], **kwargs):
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

def plot_exp_vs_comp(**kwargs):
    """Plot individual data points with legend.

    By default, the following should work:
    - Run from the top mdpow directory.
    - Prepare ``data/run01/pow.txt`` (must prepend header and append
      footer so that it is proper table). See :func:`load_data`.
    - Prepare ``experimental/targets.csv`` if it does not exist or if
      something changed. See :func:`load_exp`

    :Keywords:
       *experiments*
           path to ``targets.csv``
       *data*
           path to ``pow.txt``
       *figname*
           write figure to *filename*; suffix determines file type
       *ymin*, *ymax*
           limits of the plot in the y direction (=computational results)
    """
    from pylab import figure, plot, legend, ylim
    from matplotlib import colors, cm, rc
    import matplotlib

    experimental = load_exp(filename=kwargs.pop('experiments', DEFAULTS['experiments']))
    computed = load_data(filename=kwargs.pop('data', DEFAULTS['data']),
                         name="logPow_computed", connection=experimental.connection)  # add to experimental db

    # combined (matched on the itp_name!)
    c = experimental.SELECT("CommonName AS name, directory AS comment, DeltaG0, "
                            "__self__.logPow AS exp, C.logPow AS comp", 
                            "LEFT JOIN logPow_computed AS C using(itp_name)")

    # need large figure for the plot
    matplotlib.rc('figure', figsize=kwargs.pop('figsize', (8,10)))
    # default font
    matplotlib.rc('font', size=10)

    norm = colors.normalize(0,len(c))
    for i, (name,comment,DeltaA0,exp,comp) in enumerate(c):
        if exp is None or comp is None:
            continue
        color = cm.jet(norm(i))
        label = "%(comment)s %(exp).1f/%(comp).1f" % vars()
        plot(exp,comp, marker='o', markersize=10, color=color, label=label)

    legend(ncol=3, numpoints=1, loc='lower right', prop={'size':8})
    figname = _finish(**kwargs)

    matplotlib.rcdefaults()  # restore defaults
    return figname
