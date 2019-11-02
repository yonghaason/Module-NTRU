Identity-Based Encryption over Module-NTRU Lattices
===========

This software is a computer program which purpose is to provide 
a proof-of-concept implementation of PKC 2020 submission
"A New Trapdoor over Module-NTRU Lattices and its Application to ID-based Encryption".

This software is based on Thomas Prest's original implementation
governed by the CeCILL license under French law and abiding by the rules of distribution of free software.
The original implementation can be found in "https://github.com/tprest/Lattice-IBE/"


Warning
=======
This code is not to be considered secure, efficient or fully portable. Its purpose is not to be used for actual encryption, but to provide the research community a tool to verify, analyze and reproduce the statements made in our paper.

How to use?
===========

To modify the parameters, edit the values N0 and q0 in params.h.

To run on an Unix machine with g++:
```
$ make
$ ./IBE
```

If GMP and NTL are not in a standard directory, you have to modify the CCFLAGS and LDFLAGS in the Makefile to indicate where they are.
