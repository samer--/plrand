# plrand
## Pseudorandom generator with various distributions

This module provides for the sampling of pseudorandom numbers from a
variety of distributions, while maintaining complete and explicit control
of the random generator state. The generator is a high quality, long
period combined multiple recursive generator by Pierre L'Ecuyer,
with a period of about 2^191.

The state can be initialised from the truly random /dev/random if this 
device exists.

The state can be efficiently advanced a known number of steps - this
supports stream splitting, by partitioning the sequence of 2^191
values into subsequences.

Predicates for sampling from and computing probabilities from several
distributions are included, but not exported -- there is a nicer interface
to them in progress -- but in the mean time, see examples/test_crp.pl
for usage.

See plrand.pl and crp.pl for more information.

### PREREQUISITES

GSL (GNU Scientific Library) for Zeta distribution and probability computations [optional]
Pack genutils to run examples/test_crp.pl


### ACKNOWLEDGEMENTS

This distribution contains a modified version of the RngStreams library
by Pierre L'Ecuyer and Richard Simard. See cpp/RngStream.{h,c} in this
distribution.
