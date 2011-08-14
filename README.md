##### Chuckie !["Rugrat + BoxLib = Chuckie"](https://github.com/mynameisfiber/chuckie/raw/master/chuckie.png)

# What is this?

This is a code written in fortran, meant to be compiled with [f2py](http://www.scipy.org/F2py) and used in python (or any other language that can interface with fortran ;).  The code itself is meant to solve the [einstein field equations](http://en.wikipedia.org/wiki/Einstein_field_equations) using the [BSSN 3+1 formulation](http://en.wikipedia.org/wiki/BSSN_formalism)

# Why Chuckie?

This code is based off of some work done by [Paul Duffel](http://duffell.org/) called Rugrat.  Rugrat was written in C, so moving it to fortran would require a name change, but nothing too drastic.

# How do I use it?

Easy!  First you must have the basic requirements:

* python (>2.4)
* numpy
* gfortran (or other fortran compiler)

To compile, simply type:

    $ make evolve.so

And now you will have a module that is loadable in python!  How do you use that module?  Well, look in the examples folder to see some simple usage and the test folder for some actual physics!
