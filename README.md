# Uncertainty propagation in C++
[![CircleCI](https://dl.circleci.com/status-badge/img/gh/martukas/uncertain/tree/master.svg?style=svg)](https://dl.circleci.com/status-badge/redirect/gh/martukas/uncertain/tree/master)
[![Coverage Status](https://coveralls.io/repos/github/martukas/uncertain/badge.svg)](https://coveralls.io/github/martukas/uncertain)
[![codecov](https://codecov.io/gh/martukas/uncertain/branch/master/graph/badge.svg?token=0FGFZJ7QEF)](https://codecov.io/gh/martukas/uncertain)

A set of drop-in classes that keep track of uncertainties along with calculations.
Most math functions are reimplemented, allowing one to use these classes in place of the
primitive arithmetic `double` type without extensive rewriting of code.

## Word of warning

:warning: Reported code coverage is optimistic. No proper unit tests have been implemented. :warning:

## Background
This code is an adaptation of code produced by Evan Manning for NASA/JPL. The original project was described in an article
for C/C++ Users Journal in March 1996, available [here](http://www.pennelynn.com/Documents/CUJ/HTML/14.03/MANNING/MANNING.HTM).
The code has been released by Caltech under the following [license](LICENSE).