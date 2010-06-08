# mdpow.maths -- helper functions for analysis
# Copyright (c) 2010 Oliver Beckstein <orbeckst@gmail.com>
# Released under the GNU Public License version 3 (or higher)
"""
:mod:`mdpow.numkit` --- Mathematical analysis helper functions
==============================================================

.. autofunction:: Pearson
.. autofunction:: linfit
.. autofunction:: autocorrelation_fft
.. autofunction:: averaged_autocorrelation
.. autofunction:: tcorrel

Examples
--------

Autocorrelation time (time when ACF becomes 0 for the first time)::
  R = gromacs.formats.XVG("./md.xvg")
  acf = mdpow.numkit.autocorrelation_fft(R.array[1])
  where(acf <= 0)[0][0]
 
Alternatively, fit an exponential to the ACF and extract the time constant.

"""

import numpy
import scipy.integrate

# functions copied from hop.utilities

def Pearson_r(x,y):
    """Pearson's r (correlation coefficient)

      Pearson(x,y) --> correlation coefficient

    *x* and *y* are arrays of same length.
    
    Historical note:
    Naive implementation of Pearson's r:

      Ex = scipy.stats.mean(x)
      Ey = scipy.stats.mean(y)
      covxy = numpy.sum((x-Ex)*(y-Ey))
      r = covxy/math.sqrt(numpy.sum((x-Ex)**2)*numpy.sum((y-Ey)**2))
    """
    return numpy.corrcoef(x,y)[1,0]

def linfit(x,y,dy=[]):
    """Fit a straight line y = a + bx to the data in x and y; errors
    on y should be provided in dy in order to assess the goodness of
    the fit and derive errors on the parameters.

      linfit(x,y[,dy]) --> result_dict

    Fit y = a + bx to the data in x and y by analytically minimizing
    chi^2. dy holds the standard deviations of the individual y_i. If
    dy is not given, they are assumed to be constant (note that in
    this case Q is set to 1 and it is meaningless and chi2 is
    normalised to unit standard deviation on all points!).

    Returns the parameters a and b, their uncertainties sigma_a and
    sigma_b, and their correlation coefficient r_ab; it also returns
    the chi-squared statistic and the goodness-of-fit probability Q
    (that the fit would have chi^2 this large or larger; Q < 10^-2
    indicates that the model is bad --- Q is the probability that a
    value of chi-square as _poor_ as the calculated statistic chi2
    should occur by chance.)
  
    result_dict::
       intercept, sigma_intercept    a +/- sigma_a
       slope, sigma_slope            b +/- sigma_b
       parameter_correlation         correlation coefficient r_ab
                                     between a and b
       chi_square                    chi^2 test statistic
       Q_fit                         goodness-of-fit probability

    Based on 'Numerical Recipes in C', Ch 15.2.
    """
    import scipy.stats

    n = len(x)
    m = len(y)
    if n != m:
        raise ValueError("lengths of x and y must match: %s != %s" % (n, m))
    
    try:
        have_dy = (len(dy) > 0)
    except TypeError:
        have_dy = False

    if not have_dy:
        dy = numpy.ones((n),numpy.float)

    x  = numpy.asarray(x)
    y  = numpy.asarray(y)
    dy = numpy.asarray(dy)

    s2  = dy*dy
    S   = numpy.add.reduce(1/s2)
    Sx  = numpy.add.reduce(x/s2)
    Sy  = numpy.add.reduce(y/s2)
    Sxx = numpy.add.reduce(x*x/s2)
    Sxy = numpy.add.reduce(x*y/s2)

    t   = (x - Sx/S)/dy
    Stt = numpy.add.reduce(t*t)

    b = numpy.add.reduce(t*y/dy)/Stt
    a = (Sy - Sx*b)/S

    sa = numpy.sqrt((1 + (Sx*Sx)/(S*Stt))/S)
    sb = numpy.sqrt(1/Stt)

    covab = -Sx/(S*Stt)
    r = covab/(sa*sb)

    chi2 = numpy.add.reduce(((y-a-b*x)/dy)**2)
    if not have_dy:
        # estimate error if none were provided
        sigmadata = numpy.sqrt(chi2/(n-2))
        sa *= sigmadata
        sb *= sigmadata
        Q = 1.0
    else:
        Q = scipy.stats.chisqprob(chi2,n-2)

    return {"intercept":a,"slope":b,
            "sigma_intercept":sa,"sigma_slope":sb,
            "parameter_correlation":r, "chi_square":chi2, "Q":Q}

def autocorrelation_fft(series,include_mean=False,periodic=False,
                        normalize=True,start=None,stop=None,**kwargs):
    """Calculate the auto correlation function.
    
       autocorrelation_fft(series,include_mean=False,**kwargs) --> acf

    The time series is correlated with itself across its whole length. It is
    0-padded and the ACF is corrected for the 0-padding (the values for larger
    lags are increased) unless mode='valid' (see below).  Only the
    [0,len(series)[ interval is returned. The series is normalized to its 0-th
    element.

    Note that the series for mode='same'|'full' is inaccurate for long times
    and should probably be truncated at 1/2*len(series). Alternatively, only
    sample a subseries with the stop keyword.

    :Arguments:
      *series*
        (time) series, a 1D numpy array
      *include_mean*
        ``False``: subtract mean(series) from series [``False``]
      *periodic*
        ``False``: corrected for 0-padding; ``True``: return as is
        [``False``]
      *normalize*
        ``True`` divides by acf[0] so that the first element is 1;
        ``False`` leaves un-normalized [``True``]
      *start*,*stop*
        If set, calculate the ACF of series[start:stop] with series;
        in this case mode='valid' is enforced
      *mode*
        "full" | "same" | "valid": see :func:`scipy.signal.fftconvolve`
      *kwargs*
        other keyword arguments for :func:`scipy.signal.fftconvolve`
    """
    import scipy.signal
    kwargs.setdefault('mode','full')
    series = numpy.squeeze(series.astype(float)).copy()   # must copy because de-meaning modifies it
    if not include_mean:
        mean = series.mean()
        series -= mean
    if start or stop:
        kwargs['mode'] = 'valid'   # makes sense to only return the true unpolluted ACF
        start = start or 0
        stop = stop or len(series)
        if start >= stop:
            raise ValueError('Must be start < stop but start = %(start)d >= stop = %(stop)d.' 
                             % locals())

    ac = scipy.signal.fftconvolve(series,series[stop:start:-1,...],**kwargs)

    if kwargs['mode'] == 'valid':
        # origin at start+1
        if normalize:
            norm = ac[start+1] or 1.0   # to guard against ACFs of zero arrays
        else:
            norm = 1.0
        # Note that this is periodic (and symmetric) over result[0,stop-start+1] and so we 
        # only return one half:
        ##return numpy.concatenate( (ac[start+1:], ac[:start+1]) )[:len(ac)/2] / norm
        # ac is symmetric around start+1 so we average the two halves (just in case):
        ac[:] = numpy.concatenate( (ac[start+1:], ac[:start+1]) ) / norm
        ac = numpy.resize(ac,len(ac)+1)   # make space for replicated 0-th element
        ac[-1] = ac[0]
        if len(ac) % 2 == 1:   
            # orig ac was even
            return 0.5*(ac[:len(ac)/2] + ac[:len(ac)/2:-1])
        else:
            # orig ac was odd: replicate the least important datapoint for second half
            return 0.5*(ac[:len(ac)/2] + ac[:len(ac)/2-1:-1])
    else:
        origin = ac.shape[0]/2        # should work for both odd and even len(series)
        ac = ac[origin:]
        assert len(ac) <= len(series), "Oops: len(ac)=%d  len(series)=%d" % (len(ac),len(series))
        if not periodic:
            ac *= len(series)/(len(series) - 1.0*scipy.arange(len(ac)))   # correct for 0 padding
        if normalize:
            norm = ac[0] or 1.0  # to guard against ACFs of zero arrays
        else:
            norm = 1.0
        return ac/norm

def averaged_autocorrelation(series,length=None,sliding_window=None,**kwargs):
    """Calculates the averaged ACF of a series.

      averaged_autocorrelation(series,length=None,sliding_window=None) --> mean(acf), std(acf)

    Calculate the ACF of a series for only a fraction of the total length,
    *length* but repeat the calculation by setting the origin progressively
    every *sliding_window* steps and average over all the ACFs.

    :Arguments:
      *series*
          time series (by default, mean will be removed)
      *length*
          length (in frames) of the ACF [``1/2*len(series)``]
      *sliding_window*
          repeat ACF calculation every *sliding_window* frames
          [``len(series)/100``]
      *kwargs*
          additional arguments to :func:`autocorrelation_fft`
    """
    kwargs.pop('start',None) # must filter those kwargs as we set them ourselves
    kwargs.pop('stop',None)
    nframes = len(series)
    length = length or nframes/2
    _length = nframes - length   # _length is the length of the comparison series
    sliding_window = sliding_window or nframes/100
    # note: do NOT be tempted to change nframes-_length to nframes-_length+1
    #       (this will make the last acf 1 step longer, see series[stop:start:-1] !)
    acfs = numpy.array([autocorrelation_fft(series,start=start,stop=start+_length,**kwargs) 
                 for start in xrange(0,nframes-_length,sliding_window)])
    return acfs.mean(axis=0), acfs.std(axis=0)


def tcorrel(x,y,nstep=100):
    """Calculate the correlation time and an estimate of the error.

    The autocorrelation function is calculated via FFT on every *nstep* of the
    data. It is assumed to decay exponentially, f(t) = exp(-t/tau) and the
    decay constant is estimated as the integral of the ACF from the start up to
    its first root.

    See Frenkel and Smit, Academic Press, San Diego 2002, p526.

    .. Note:: *nstep* should be set sufficiently large so that there are less
              than ~50,000 entries in the input.
    """
    if x.shape != y.shape:
        raise TypeError("x and y must be y(x), i.e. same shape")
    _x = x[::nstep]  # do not run acf on all data: takes too long
    _y = y[::nstep]  # and does not improbe accuracy
    acf = autocorrelation_fft(_y, normalize=False)
    try:
        i0 = numpy.where(acf <= 0)[0][0]  # first root of acf
    except IndexError:
        i0 = -1   # use last value as best estimate
    t0 = _x[i0]
    # integral of the _normalized_ acf
    tc = scipy.integrate.simps(acf[:i0]/acf[0], x=_x[:i0])
    # error estimate [Frenkel & Smit, p526]
    sigma = numpy.sqrt(2*tc*acf[0]/(x[-1] - x[0]))

    return {'tc':tc, 't0':t0, 'sigma':sigma, 
            't':_x[:i0], 'acf':acf[:i0]}     # for debugging




class FitFunc(object):
    """Fit a function f to data (x,y) using the method of least squares.

    The function is fitted when the object is created, using
    :func:`scipy.optimize.leastsq`. One must derive from the base class
    :class:`FitFunc` and override the :meth:`FitFunc.f_factory` (including
    the definition of an appropriate local :func:`fitfunc` function) and
    :meth:`FitFunc.initial_values` appropriately. See the examples for a
    linear fit :class:`FitLin`, a 1-parameter exponential fit :class:`FitExp`,
    or a 3-parameter double exponential fit :class:`FitExp2`.

    The object provides two attributes
     :attr:`FitFunc.parameters`
           list of parameters of the fit
     :attr:`FitFunc.message`
           message from :func:`scipy.optimize.leastsq`

    After a successful fit, the fitted function can be applied to any data (a
    1D-numpy array) with :meth:`FitFunc.fit`.
     
    """
    def __init__(self,x,y):
        import scipy.optimize
        _x = numpy.asarray(x)
        _y = numpy.asarray(y)
        p0 = self.initial_values()
        fitfunc = self.f_factory()
        def errfunc(p,x,y):
            return  fitfunc(p,x) - y     # residuals        
        p,msg = scipy.optimize.leastsq(errfunc,p0[:],args=(_x,_y))
        try:
            p[0]
            self.parameters = p
        except (TypeError,IndexError,):
            # TypeError for int p, IndexError for numpy scalar (new scipy)
            self.parameters = [p]
        self.message = msg

    def f_factory(self):
        """Stub for fit function factory, which returns the fit function.
        Override for derived classes.
        """
        def fitfunc(p,x):
            # return f(p,x); should be a numpy ufunc
            raise NotImplementedError("base class must be extended for each fit function")
        return fitfunc

    def initial_values(self):
        """List of initital guesses for all parameters p[]"""
        # return [1.0, 2.0, 0.5]
        raise NotImplementedError("base class must be extended for each fit function")    

    def fit(self,x):
        """Applies the fit to all *x* values"""
        fitfunc = self.f_factory()
        return fitfunc(self.parameters,numpy.asarray(x))

class FitExp(FitFunc):
    """y = f(x) = exp(-p[0]*x)"""
    def f_factory(self):
        def fitfunc(p,x):
            return numpy.exp(-p[0]*x)   # exp(-B*x)
        return fitfunc
    def initial_values(self):
        return [1.0]
    def __repr__(self):
        return "<FitExp "+str(self.parameters)+">"

class FitExp2(FitFunc):
    """y = f(x) = p[0]*exp(-p[1]*x) + (1-p[0])*exp(-p[2]*x)"""
    def f_factory(self):
        def fitfunc(p,x):
            return p[0]*numpy.exp(-p[1]*x) + (1-p[0])*numpy.exp(-p[2]*x)
        return fitfunc
    def initial_values(self):
        return [0.5,0.1,1e-4]
    def __repr__(self):
        return "<FitExp2"+str(self.parameters)+">"

class FitLin(FitFunc):
    """y = f(x) = p[0]*x + p[1]"""
    def f_factory(self):
        def fitfunc(p,x):
            return p[0]*x + p[1]
        return fitfunc
    def initial_values(self):
        return [1.0,0.0]
    def __repr__(self):
        return "<FitLin"+str(self.parameters)+">"
