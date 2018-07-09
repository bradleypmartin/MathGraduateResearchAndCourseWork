naiveFlag = 0;
CGflag = 0; CGiterations = 1000; 
FOflag = 1; FOalpha = 1000; FOloc = 0;
ZOflag = 0; ZOalpha = 1; ZOradius = 0.025;
RhoKnownFlag = 1; Rhoalpha = 1000; knownRhoVec = 1./[ones(100,1); 1.5*ones(100,1)];
intKnownFlag = 1;

%intLocs = [0.31001; 0.545001; 0.792001];
%intLocs = [0.31001; 0.53001; 0.77001];
%intLocs = [0.31001; 0.52001; 0.75001];
%intLocs = [0.31001; 0.435001; 0.525001];
intLocs = [0.50001];

numNewtons = 6;

for k = 1:numNewtons

    intStack = intStackSetup(intLocs,xgrid,ufxtK,bigSize,N);
    
    [xgrid ufxtKP1 bigOperator bigJacobian F] = ...
        FD4STInv9(N,tstep,endtime,dampConst,naiveFlag,...
        knownu,knownf,knownuInd,knownfInd,knownuAlpha,knownfAlpha,ufxtK,...
        CGflag,CGiterations,spaceSigma,timeSigma,exLoc,exTime,...
        FOflag,FOalpha,FOloc,ZOflag,ZOalpha,ZOradius,...
        intKnownFlag,intStack,RhoKnownFlag,Rhoalpha,...
        ZOswitchu,ZOswitchf,knownRhoVec);

    ufxtK = ufxtKP1;

end