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

### PREREQUISITES

SWI Prolog
GSL (GNU Scientific Library) for Zeta distribution [optional]


### INSTALLATION

Edit the variables in the top half of the root Makefile if necessary.

If you don't have GSL, you can set HAVE_GSL to 0 to disable this feature.
In fact, all that happens in this case is that calls to two GSL functions
are replaced with NAN so that attempts to sample from the Zeta distribution
will not work.

The build process installs a Prolog module file and a foreign library
to ~/lib/prolog by default. If you wish to change this, edit the root Makefile
accordingly and be sure that the referenced directories are in your
Prolog file_search_path.

In the root directory of this package, type

	$ make install


### ACKNOWLEDGEMENTS

This distribution contains a modified version of the RngStreams library
by Pierre L'Ecuyer and Richard Simard. See cpp/RngStream.{h,c} in this
distribution.
