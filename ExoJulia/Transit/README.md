15 April 2016 - Eric Agol

This directory contains routines for computing transit light
curves in Julia.

The workhorse routine is occultquad.jl which implements the
computations from Mandel & Agol (2002) (correcting typos in
the paper).

The code has been tested from p=0.01 to p=100 with a range
of values and special cases, and works to high precision in
all test cases.

To give it a whirl, type Julia> include("test_quad.jl")
which will show plots of flux & the derivative of the
flux (added to one to display along with flux) for the
fiducial parameters given at the end of the test_quad.jl
routine.
