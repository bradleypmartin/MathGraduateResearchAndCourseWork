function datastore = ...
    RBFHeat2exe(bigLmatrix,endtime,k,initdatavec,storesteps,...
    exLocX,exLocY,timeSigma,spaceSigma,xnodes,ynodes,...
    exTime,N,numIntNodes,k1,k2,yLI,yUI,closedFlag)

timesteps = ceil(endtime/k);  % total number of time steps through which
                              % we integrate data

% below, we initialize the storage array DATASTORE to store our data.

datavecCurrent = zeros(size(initdatavec,1),1);
datastore = zeros(size(initdatavec,1),size(storesteps,1));

datastore(:,1) = datavecCurrent;
sConst = 2;

% we combine the block differential operator with a block hyperviscosity
% operator to make one unified operator with which we can time step.

boundVec = [zeros(N-numIntNodes,1); (1-cos(sConst*pi*xnodes(N-numIntNodes+1:N,1)))/2];
boundNullVec = [ones(N-numIntNodes,1); zeros(numIntNodes,1)];

tic;

% In the time stepping FOR loop below, we integrate via RK4.

for time_step = 2:timesteps
    
    time_step
    expTime = (time_step-2)*k-0.05;
    
    %%% RK4 step 1
    
    d1 = bigLmatrix*datavecCurrent;
    step2data = datavecCurrent+k/2*d1;
    locT = timeSigma*(expTime+k/2);
    step2data = boundNullVec.*step2data + (1/(1+exp(-locT)))*boundVec;
    
    d2 = bigLmatrix*step2data;
    step3data = datavecCurrent+k/2*d2;
    locT = timeSigma*(expTime+k/2);
    step3data = boundNullVec.*step3data + (1/(1+exp(-locT)))*boundVec;
    
    d3 = bigLmatrix*step3data;
    step4data = datavecCurrent+k*d3;
    locT = timeSigma*(expTime+k);
    step4data = boundNullVec.*step4data + (1/(1+exp(-locT)))*boundVec;
    
    d4 = bigLmatrix*step4data;
    
    datavecCurrent = datavecCurrent+k/6*(d1+2*d2+2*d3+d4);
    datavecCurrent = boundNullVec.*datavecCurrent + (1/(1+exp(-locT)))*boundVec;
    
    for m = 1:size(storesteps,1)
        if storesteps(m,1) == time_step
            datastore(:,m) = datavecCurrent;
        end
    end
    
end

disp('For time stepping:')
toc;
disp(' ')

%%% End function RBFewave6msexe %%%

end

    