# -*- coding: utf-8 -*-
from numpy import arange, asarray, concatenate, convolve, exp, mat, zeros
from numpy.linalg import pinv
from scipy.misc import factorial
from matplotlib.pyplot import figure, plot


def deconvolve(x, y, sigma, noise, wsize=3):
    """ Deconvolves function y(x) by factor sigma.
        Noisy data  => large wsize (so derivatives don't fluctuate wildly)
        Large wsize => curve overshoots at base of peaks and peaks cannot
        become so narrow
        => We want smallest wsize possible while still getting smooth results
    """
    sigma = float(sigma)
    ret = zeros(y.size)
    dx = x[1] - x[0]

    # Find optimal window size for smoothing of y
    w0 = 3
    dw = durbin_watson(y, savitzky_golay(y, w0, 2, dx))
    while True:
        dw0 = durbin_watson(y, savitzky_golay(y, w0+2, 2, dx))
        if (abs(dw0-1) < abs(dw-1)):
            dw = dw0
            w0 = w0 + 2
        else:
            break

    i = 0
    while True:
        w = wsize*2*i + w0
        y0 = concatenate((zeros(w), y, zeros(w)))
        temp = savitzky_golay(y0, w, 2*i, dx, deriv=2*i)[w:-w]
        c = (-sigma**2/4)**i / factorial(i)
        temp *= c

        if max(abs(temp)) > max(abs(y)):    # Process probably won't converge
            print('Deconvolution Failed')
            return None

        elif max(abs(temp)) > noise:
            ret += temp
            i += 1

        else:
            break

        # Remove noise along x=0 at each iteration
        # (otherwise noise may be amplified)
        ret[ret < noise] = 0

    return ret


def durbin_watson(y1, y2):
    """ Returns the Durbin-Watson parameter for a pair of curves.
        The closer DW is to 1, the better y2 approximates y1
    """
    n = asarray(y1).size
    num = 0.
    den = (y1[0]-y2[0])**2

    for i in range(1, n):
        temp = (y1[i]-y2[i]) - (y1[i-1]-y2[i-1])
        num += temp**2
        den += (y1[i]-y2[i])**2

    n = float(n)

    return (n/(n-1)) * (num / (4*den))


def G(x, n, p):
    """ Returns the evaluation at x of the sum of n Gaussian curves of the form
        A * exp(-(x-m)^2/s).
        For the kth curve (k=0...n-1), A=p[3*k], m=p[3*k+1], s=p[3*k+2].
    """
    if n == 0:
        return 0
    else:
        return p[0] * exp(-(x-p[1])**2 / p[2]) + G(x, n-1, p[3:])


def multi_deconvolve(x, y, sigs, noise, wsize=3):
    """ Apply multiple recursive deconvolutions, with deconvolution parameters
        given by sigs. Applying many small recursive deconvolutions is often
        more stable than apply one large deconvolution.
    """
    n = len(sigs)

    if n > 1:
        yd = deconvolve(x, y, sigs[-1], noise, wsize)
        if yd is None:
            return None
        else:
            return multi_deconvolve(x, yd, sigs[:-1], noise, wsize)
    else:
        return deconvolve(x, y, sigs[0], noise, wsize)


def savitzky_golay(y, window_size, order, dx, deriv=0):
    """ Uses the Savitsky-Golay algorithm to return a smoothed version of the
        specified derivative of y
    """
    half_window = (window_size - 1) // 2
    b = mat([[k**i for i in range(order+1)]
             for k in range(-half_window, half_window+1)])
    m = pinv(b).A[deriv] * (1/dx)**deriv * factorial(deriv)

    firstvals = y[0] + abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] - abs(y[-half_window-1:-1][::-1] - y[-1])

    y = concatenate((firstvals, y, lastvals))

    return convolve(m[::-1], y, mode="valid")


if __name__ == "__main__":
    """ A demonstration of the deconvolution of a superposition of two
        overlapping Gaussian curves
    """
    x = arange(-10, 10, 0.025)
    y = G(x, 2, [100, 0, 1, 70, 1.5, 1])
    figure(1)
    plot(x, y, c='b')                           # Plot sum of curves
    plot(x, G(x, 1, [100, 0, 1]), c='r')        # Plot curve 1
    plot(x, G(x, 1, [70, 1.5, 1]), c='r')       # Plot curve 2
    figure(2)
    plot(x, deconvolve(x, y, 0.5, 0.1), c='b')  # Plot deconvolution
    plot(x, G(x, 1, [100, 0, 1]), c='r')        # Plot curve 1
    plot(x, G(x, 1, [70, 1.5, 1]), c='r')       # Plot curve 2
