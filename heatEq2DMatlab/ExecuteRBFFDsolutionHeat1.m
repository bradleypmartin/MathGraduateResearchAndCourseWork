% EXECUTERBFFDSOLUTION is a script written to actually run a wave model
% once dx/dy and hyperviscosity RBF-FD operators have been made and domain
% properties have been set.  These tasks can be carried out in a simple way
% with the EXEPREPRBF6 script.

exLocX = 0.5;
exLocY = 0.75;
timeSigma = 1;
spaceSigma = 23;
exTime = 0.1;

const = max(2*round(N/2500),1);

endtime = 0.11;               % end time of model
k = 0.0001*1/const;                   % time step
storesteps = (1:100*const:1000*const+1)';       % vector of time steps to save (for visualization)
%storesteps = (1:150:1801)';

% Running the simulation with RBFEWAVE6MSEXE

datastore = RBFHeat1exe(bigLmatrix,endtime,k,...
    initdatavec,storesteps,exLocX,exLocY,timeSigma,spaceSigma,...
    xnodes,ynodes,exTime,N,numIntNodes,rho1,rho2,yLI,yUI,closedFlag);
