# EfficientBowedString
Companion page for the paper *Efficient Simulation of the Bowed String in Modal Form*.

This repository stores Matlab code (folder `MatlabCode`) and sample audios (folder `Sounds`) synthesised by solving the stiff string equation with the modal algorithm detailed in the paper.  

###Paper Errata: The curve in Figure 1(c) of the paper, along with the corresponding formula in the caption, is
incorrect. The correct version is stored here under the name `frictionCurves`. The correct formula is: $`\phi_\mathrm{C}(\eta)=\sign(\eta)(0.4 e^{-|\eta|/0.7}+0.35)`$.

The folder `RealTimeImpl` stores a real-time working string made with JUCE

The names of the scripts indicate:  
`bowedStiffStringModal`: stiff string algorithm solved with the modal non-iterative solver\
`bowedStiffStringFDTD`: stiff string algorithm solved with the non-iterative solver applied to the first-order system (FDTD), with Sherman-Morrison\
`SOIT`: ideal string solved with an iterative solver applied to the second-order system (FDTD)\
`FOIT`: ideal string solved with an iterative solver applied to the first-order system (FDTD)\
`FONIT`: ideal string solved with the non-iterative solver applied to the first-order system (FDTD), No Sherman-Morrison\
`MOD`: ideal string solved with the non-iterative solver applied to the modal system  


The prefix notes indicate which string was simulated. The suffix `Free` means that the string was bowed and then left free to vibrate: this way the decaying sound can be percieved. `Stop` specifies that the string was bowed and then stopped. `5th` indicates that the string base lentgh was reduced of 1/3, thus the output sound corresponds to the pythagorean 5th. The strings physical parameters were taken from: C. Desvages, Physical modelling of the bowed string and applications to sound synthesis, Ph.D. thesis, The University of Edinburgh, 2018.

In all cases the bow velocity started from zero and reached a maximum value through a linear ramp, was kept constant for a while and then decreased linearly until zero.
The stopping of vibration was implemented by keeping the bowing force constant (and not null) after the bow velocity reached zero. The Free case, on the other hand, was obtained by setting the force to zero as well. The amplitude envelope of the samples reflects the forcing type: in the "Stop" case, the string is put into vibration, and has a fast decay after the bow is stopped. In the "Free" case, the vibration continues after the bow is stopped, and the decay is slower.
