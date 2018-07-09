function datastore = ...
    RBFHeat1exe(bigLmatrix,endtime,k,initdatavec,storesteps,...
    exLocX,exLocY,timeSigma,spaceSigma,xnodes,ynodes,...
    exTime,N,numIntNodes,k1,k2,yLI,yUI,closedFlag)

timesteps = ceil(endtime/k);  % total number of time steps through which
                              % we integrate data

% below, we initialize the storage array DATASTORE to store our data.

if closedFlag > 0
    C1 = (pi^2+timeSigma/k1)^0.5;
    C2 = (pi^2+timeSigma/k2)^0.5;
    sConst = 1;
else
    C1 = (4*pi^2+timeSigma/k1)^0.5;
    C2 = (4*pi^2+timeSigma/k2)^0.5;
    sConst = 2;
end

A = [1 1 0 0 0 0;
    exp(C1*yLI) exp(-C1*yLI) -exp(C2*yLI) -exp(-C2*yLI) 0 0;
    C1*k1*exp(C1*yLI) -C1*k1*exp(-C1*yLI) -C2*k2*exp(C2*yLI) C2*k2*exp(-C2*yLI) 0 0;
    0 0 exp(C2*yUI) exp(-C2*yUI) -exp(C1*yUI) -exp(-C1*yUI);
    0 0 C2*k2*exp(C2*yUI) -C2*k2*exp(-C2*yUI) -C1*k1*exp(C1*yUI) C1*k1*exp(-C1*yUI);
    0 0 0 0 exp(C1) exp(-C1)];

cVec2 = A\([0 0 0 0 0 1]');

Z1 = sin(sConst*pi*xnodes).*(cVec2(1,1)*exp(ynodes*C1)+cVec2(2,1)*exp(-ynodes*C1));
Z2 = sin(sConst*pi*xnodes).*(cVec2(3,1)*exp(ynodes*C2)+cVec2(4,1)*exp(-ynodes*C2));
Z3 = sin(sConst*pi*xnodes).*(cVec2(5,1)*exp(ynodes*C1)+cVec2(6,1)*exp(-ynodes*C1));

datavecCurrent = (ynodes < yLI).*Z1 + (ynodes >= yLI).*(ynodes < yUI).*Z2 + ...
    (ynodes >= yUI).*Z3;

datastore = zeros(size(initdatavec,1),size(storesteps,1));

datastore(:,1) = datavecCurrent;

% we combine the block differential operator with a block hyperviscosity
% operator to make one unified operator with which we can time step.

boundVec = [zeros(N-numIntNodes,1); sin(sConst*pi*xnodes(N-numIntNodes+1:N,1))];
boundNullVec = [ones(N-numIntNodes,1); zeros(numIntNodes,1)];

tic;

% In the time stepping FOR loop below, we integrate via RK4.

for time_step = 2:timesteps
    
    time_step
    expTime = (time_step-2)*k;
    
    %%% RK4 step 1
    
    d1 = bigLmatrix*datavecCurrent;
    step2data = datavecCurrent+k/2*d1;
    step2data = boundNullVec.*step2data + exp(timeSigma*(expTime+k/2))*boundVec;
    
    d2 = bigLmatrix*step2data;
    step3data = datavecCurrent+k/2*d2;
    step3data = boundNullVec.*step3data + exp(timeSigma*(expTime+k/2))*boundVec;
    
    d3 = bigLmatrix*step3data;
    step4data = datavecCurrent+k*d3;
    step4data = boundNullVec.*step4data + exp(timeSigma*(expTime+k))*boundVec;
    
    d4 = bigLmatrix*step4data;
    
    datavecCurrent = datavecCurrent+k/6*(d1+2*d2+2*d3+d4);
    datavecCurrent = boundNullVec.*datavecCurrent + exp(timeSigma*(expTime+k))*boundVec;
    
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

    