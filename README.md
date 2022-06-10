# EfficientBowedString
Companion page for the paper *Efficient Simulation of the Bowed String in Modal Form*.

This repository stores sample audios synthesised by solving the stiff string equation with the modal algorithm detailed in the paper. The prefix notes indicate which string was simulated. The suffix `Free` means that the string was bowed and then left free to vibrate: this way the decaying sound can be percieved. `Stop` specifies that the string was bowed and then stopped. `5th` indicates that the string base lentgh was reduced of 1/3, thus the output sound corresponds to the pythagorean 5th. The strings parameters were taken from: C. Desvages, Physical modelling of the bowed string and applications to sound synthesis, Ph.D. thesis, The University of Edinburgh, 2018.

In all cases the bow velocity started from zero and reached a maximum value through a linear ramp, was kept constant for a while and then decreased linearly until zero.
The stopping of vibration was implemented by keeping the bowing force constant after the bow velocity reached zero. The Free case, on the other hand, was obtained by setting the force to zero as well.
