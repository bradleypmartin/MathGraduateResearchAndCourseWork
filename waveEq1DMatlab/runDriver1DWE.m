% This script will execute a 1D wave equation solution/simulation in the
% periodic interval [-1,1), using the FD4wave1DAC code referenced below.
% With the default parameter settings below, the simulation will store data
% between 0 and 1.5 time units, during which an initial wave front present
% at x = -0.5 travels in the positive x-direction and encounters a new
% modeled material between x = 0 and x = 0.5.  Part of the incident wave
% will be reflected, and part will travel through the heterogeneous region.

% In the default setting, the heterogeneous region features 2x the wave
% speed of the rest of the domain, and density is uniform everywhere (value
% 1.0).  The user is free to explore alternative options as described
% below.  One of the advantages of this particular solution is that it can
% accurately handle not only the immediate changes in model parameters, but
% also a diminishing thickness (even below grid resolution) of the
% heterogeneous layer.  Have fun and try out examples with very thin layers
% and high parameter contrasts!

% After running this script, you may view snapshots of the result with the
% plot1DWEResults.m script or create an .avi video with the
% waveMovieMaker1D.m script.

% USER PARAMETERS
endtime  = 1.5;       % Current simulation end (unitless time)
slowMult = 5;         % Slowdown factor of results for video making
interfaceWidth = 0.5; % Width of heterogenous layer starting at x = 0
                      % (should be positive)

c2 = 2.0; 
% Wave speed within heterogeneous region (elsewhere, c = 1); you
% may change this parameter, but be careful about how much you do
% so; a significant increase in "c" may demand a smaller time
% step than the default setting (for stability).
          
rho2 = 1; 
% Density within the heterogeneous region

N = 400;             
% Number of equispaced nodes used for the discrete FD
% solution (spaced between x = -1 and x = 1)

FDorder = 4;         % Order of FD method used for the solution
doubleNaiveFlag = 0; % By default (0), the method is set up to handle an
                     % interface "double-cross," and as stated above, very
                     % thin heterogeneous layers may be treated accurately.
                     % Setting this flag to "1" can show what effect
                     % removing this feature will have (in the case of a
                     % subgrid- or near-grid-resolution heterogeneous
                     % region.

tstep = 0.001*(400/N); % length of discrete time step (delta-t)

% BLACK BOX COMPUTATION (modify at your own risk! ;) )

% scaling time step conservatively to prepare for smooth 60 FPS video
invTstep = 1/tstep;
frameDivide = ceil(invTstep/(60*slowMult));
newinvTstep = 60*slowMult*frameDivide;
tstep = 1/newinvTstep;

% setting up equispaced arrangement of nodes along the problem (periodic)
% interval [-1,1)
prexgrid = linspace(-1,1,N+1);
xgrid = prexgrid(1,1:N)'+(prexgrid(1,2)-prexgrid(1,1))/2;

% calling the function FD4wave1DAC for solution computation
[uxt,fxt,timeVec] = ...
    FD4wave1DAC(xgrid,tstep,endtime,c2,rho2,FDorder,interfaceWidth,...
    doubleNaiveFlag);