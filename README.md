# Uncertainty propagation in C++
[![Build Status](https://travis-ci.com/martukas/uncertain.svg?branch=master)](https://travis-ci.com/martukas/uncertain)
[![Coverage Status](https://coveralls.io/repos/github/martukas/uncertain/badge.svg?branch=master)](https://coveralls.io/github/martukas/uncertain?branch=master)

A set of drop-in classes that keep track of uncertainties along with calculations.
Most math functions are reimplemented, allowing one to use these classes in place of the
primitive arithmetic `double` type without extensive rewriting of code.

## Examples

Add some examples here...

## Background
This code is an adaptation of code produced by Evan Manning for NASA/JPL. The original project was described in an article
for C/C++ Users Journal in March 1996, available [here](http://www.pennelynn.com/Documents/CUJ/HTML/14.03/MANNING/MANNING.HTM).
The code has been released by Caltech under the following [license](LICENSE).