GaussianDeconvolution
=====================

The problem of fitting Gaussian curves to noisy data is common in signal processing. This problem can be difficult if multiple overlapping curves need to be fitted, such as in Figure 1 below.

![Figure 1](/images/figure1.png)

The blue curve is the superposition of the two red curves. In this case, it's difficult to identify and locate the peak on the right, as it's obscured by the right shoulder of the peak on the left.

### Gaussian Convolution

It's [a well known result](http://mathworld.wolfram.com/Convolution.html) that the convolution of a Gaussian function *f* with a Gaussian kernel *g* results in another Gaussian function whose mean is the sum of the means of *f* and *g*, and whose variance is the sum of the variances of *f* and *g*.

Importantly, if

![Equation 1](/images/equation1.png)

then

![Equation 2](/images/equation2.png)

Convolution is a linear operation, so if we convolve the superposition of two Gaussian curves, *f1+f2*, with a Gaussian kernel *g* with mean 0, the result is again the superposition of two Gaussian curves, with the same means as *f1* and *f2*, but larger variances and smaller heights. This operation is demonstrated in Figure 2 below.

![Figure 2](/images/figure2.png)

Clearly, it is easier to identify, locate, and fit the peaks in the graph on the left-hand side, than in the graph on the right-hand side. Given the graph on the right-hand side, if it were possible to perform the reverse operation to obtain the graph on the left-hand side, we could fit Gaussian curves to the result, and we would then be able to derive the means, widths and heights of the original curves from the means, widths and heights of the processed curves.

### Gaussian Deconvolution

It's [a well known result](http://en.wikipedia.org/wiki/Weierstrass_transform#The_inverse) that if

![Equation 3](/images/equation3.png)

then

![Equation 4](/images/equation4.png)

Thus, by a change of variables, if *H(x)* is the convolution of a function *h(x)* with the Gaussian kernel *g(x)* given in (1), then

![Equation 5](/images/equation5.png)

This "deconvolution" operation is implemented here, using [Savitzky–Golay filters](http://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter) to find and smooth the derivatives in the sum, and [Durbin–Watson statistics](http://en.wikipedia.org/wiki/Durbin%E2%80%93Watson_statistic) to find the optimal window sizes for the Savitzky–Golay filters.