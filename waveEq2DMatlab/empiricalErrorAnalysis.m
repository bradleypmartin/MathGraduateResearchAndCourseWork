% This script can be used to compare a given experimental simulation's
% results to that of a reference solution, and thereby get an estimate for
% solution error (most likely on a coarser node set than the reference).

% This analysis compares y-oriented particle velocities of the two
% solutions.

% Data from refined node set (right now set as an example on ~313,000 nodes
load('trial8HO.mat')
N = size(xnodes,1);
REFwave = datastore(1+1*N:2*N,4);
xnodesREF = xnodes;
ynodesREF = ynodes;

% Data from "query" node set (right now set as an example on 2500 nodes)
load('trial7HO.mat')
N = size(xnodes,1);
QUERYwave = datastore(1+1*N:2*N,4);
xnodesQUERY = xnodes;
ynodesQUERY = ynodes;

% creation of interpolation operator
tic;
interpmatrix = createRbfInterpolantStencils(xnodesREF,ynodesREF,...
    xnodesQUERY,ynodesQUERY,0,1,0,1,0.2,49,6);
toc;

% restricting error analysis to top part of domain (away from interface)
nullVec = (ynodesQUERY < 0.95).*(ynodesQUERY > 0.60);
error = (interpmatrix*REFwave-QUERYwave).*nullVec;

L2errorEst = (sum(error.^2)/size(xnodesQUERY,1)).^0.5

% plotting (restricted) error estimate vs x,y
figure(1)
resultspcolor(xnodesQUERY,ynodesQUERY,error,0,1,0,1,0.2,200)
colormap('jet')
colorbar
%view([-85,15])
%caxis([-0.002 0.002]); axis([0 1 0 1 -0.002 0.002])
xlabel('x'); ylabel('y')



