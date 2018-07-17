<img src="README_files/banner.png" width="4032" />

The `integrated` R package estimates demographic vital rates from diverse and distinct data types. The integrated modelling process is designed to be flexible and customisable, and supports several different process models and likelihoods. Users can combine any number of data types in a single analysis.

Inference in the `integrated` package is enabled by [greta](https://greta-dev.github.io/greta/), which supports a range a approximate or fully Bayesian inference tools. `greta` uses the TensorFlow software library to support fast computation of complex probabilistic models. 

You can install the `integrated` package directly from GitHub using the devtools package:

``` r
devtools::install_github('jdyen/integrated.pkg')
```

Please leave feedback, bug reports or feature requests at the GitHub [issues page](https://github.com/jdyen/integrated.pkg/issues). 

[![build status](https://travis-ci.org/jdyen/integrated.svg?branch=master)](https://travis-ci.org/jdyen/integrated.pkg) [![codecov.io](https://codecov.io/github/jdyen/integrated.pkg/coverage.svg?branch=master)](https://codecov.io/github/jdyen/integrated.pkg?branch=master) [![license](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
 
