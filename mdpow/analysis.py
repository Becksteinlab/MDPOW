# POW
# Copyright (c) 2010 Oliver Beckstein and Bogdan Iorga

"""
:mod:`analysis` --- Collection of analysis and plotting functions
=================================================================

Simple functions to quickly plot data. Typically, it works best if ran
interactively from the top level of the POW directory!

Right now we only look at the data in :doc:`doc/results/exp_vs_comp.txt`.

Usage
-----

Load the data file::

   results = mdpow.analysis.load_data("doc/results/exp_vs_comp.txt")

Plot results and save to a pdf file::

   mdpow.analysis.plot_exp_vs_comp(results, figname="figs/run01/exp_vs_comp.pdf")


Functions
---------

.. autofunction:: load_data
.. autofunction:: plot_exp_vs_comp
.. autofunction:: plot_quick
"""

import recsql
import logging
logger = logging.getLogger('mdpow.analysis')

def load_data(filename="doc/results/exp_vs_comp.txt"):
    """Load exp_vs_comp table and return :class:`recsql.SQLarray`."""
    return recsql.SQLarray(filename=filename)

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

def plot_exp_vs_comp(a, **kwargs):
    """Plot individual data points with legend.

    :Keywords:
       *figname*
           write figure to *filename*; suffix determines file type
       *ymin*, *ymax*
           limits of the plot in the y direction (=computational results)
    """
    from pylab import figure, plot, legend, ylim
    from matplotlib import colors, cm, rc
    import matplotlib

    # need large figure for the plot
    matplotlib.rc('figure', figsize=kwargs.pop('figsize', (8,10)))
    # default font
    matplotlib.rc('font', size=10)

    norm = colors.normalize(0,len(a))    
    for i, (mol,DeltaA0,comp,exp,name,comment) in enumerate(a.recarray):
        color = cm.jet(norm(i))
        label = "%(name)s %(exp).1f/%(comp).1f" % vars()
        plot(exp,comp, marker='o', markersize=10, color=color, label=label)

    legend(ncol=3, numpoints=1, loc='lower right', prop={'size':8})
    figname = _finish(**kwargs)

    matplotlib.rcdefaults()  # restore defaults
    return figname
