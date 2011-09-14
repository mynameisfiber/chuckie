#!/bin/env python

import numpy
from evolve import evolve

data = np.ones( (64, 64, 64, 24) )

while True:
  data = evolve.advance(U, ng=2, dx=1.0, dt=1.0)
  subroutine advance (U, Unew, n, m, o, ncomp, ng, dx, dt)
