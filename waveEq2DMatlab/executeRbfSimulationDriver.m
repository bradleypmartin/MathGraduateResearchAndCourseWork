
spaceSigma = 23;  % "Tightness" of initial condition in problem domain
                  % (higher = more high-frequency components; lower =
                  % smoother)

FDorder = 4;  % order of hybrid FD stencils away from interface (I'd leave
              % this at 6 or lower for now; I forget exactly how I declared
              % the changeover between RBF-FD and traditional FD stencils)

relH = (2500/N)^0.5; % average node spacing

gamma = 2.4E-11*relH^(2*dampLk-1); %LK = 3 (you'll have to change this
                                 % hyperviscosity scaling based on Lk

% "back-of-the-envelope" style calculation to figure out dt scaling
const = 2*floor((N/2500)^0.5);
const2 = round(N/2500);
const3 = round(log(const2)/log(2));

endtime = 0.32;               % end time of model

if mod(const3,2) == 0
    k = 0.01/const;
    storesteps = 1+const*[0 10 20 30]';
else
    k = 0.01/const/1.4;
    storesteps = 1+const*[0 14 28 42]';
end

explFlag = 0;

datastore = EWERBF1hybHO(bigLmatrix,bigDampMatrix,gamma,endtime,k,...
    initdatavec,storesteps,spaceSigma,...
    xnodes,ynodes,N,explFlag,FDflagMat,lambdaMuVec,rhoVec,FDorder);
