# waveEq1DMatlab
This is a small collection of Matlab scripts and functions that can illustrate some of the basic ideas/results of accurate interface treatment in numerically solving 1D wave equation problems.

I'd suggest users start by downloading the several scripts and functions included here, and run the runDriver1DWE.m driver script.  This sets up and solves a 1-D wave equation problem in the periodic interval [-1,1).  A Gaussian-shaped wave front starts in the left side of the domain and encounters a heterogeneous region (new modeled material) at x = 0.  Part of the wave reflects off the interface where the model parameters jump (at x = 0), and part of the wave continues through to the right.  The extent to which the wave is reflected/transmitted depends on the contrast in wave speed and density across the interface (and, to a certain degree - especially in the case of a very thin spacing between the two interfaces here - the width of the heterogeneous region).

After running the driver script described above, the user is free to visualize data in the form of individual snapshots at a designated point in (unitless) time through the plot1DWEresults.m script.  Alternatively, the waveMovieMaker1D.m script is set up to create an .avi movie via Matlab's videowriter function.

After you confirm that an initial run using the default settings works out just fine, feel free to explore different changes to parameters such as c2, rho2, and interfaceWidth in the driver script (runDriver1DWE.m).  This will change the wave speed, density, and width of the heterogeneous region that starts at x = 0.

What sets this code apart from a simple, traditional FD solution is its ability to accurately model the changes that occur to wave shape as that wave crosses interfaces, even when those interfaces are close enough together that stencils cross more than one of them.
