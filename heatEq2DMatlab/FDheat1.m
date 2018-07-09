function [datastore xnodes ynodes] = ...
    FDheat1(N1D,endtime,k,HOflag,storesteps)

timesteps = ceil(endtime/k);  % total number of time steps through which
% we integrate data

% below, we initialize the storage array DATASTORE to store our data.

N = N1D^2;
xgrid = linspace(0,1,N1D);
[X Y] = meshgrid(xgrid,xgrid);
h = 1/(N1D-1);

Dymat = zeros(N1D,N1D);

if HOflag == 1
    preweightsLL = weights(-2,-2:2,1);
    preweightsRR = weights(2,-2:2,1);
    preweightsL = weights(-1,-2:2,1);
    preweightsR = weights(1,-2:2,1);
    preweightsC = weights(0,-2:2,1);
    Dymat(2,1:5) = preweightsL(2,:)/h;
    Dymat(N1D-1,N1D-4:N1D) = preweightsR(2,:)/h;
    Dymat(1,1:5) = preweightsLL(2,:)/h;
    Dymat(N1D,N1D-4:N1D) = preweightsRR(2,:)/h;
    for counter1 = 3:N1D-2
        Dymat(counter1,counter1-2:counter1+2) = preweightsC(2,:)/h;
    end
    Dxmat = Dymat';
else
    preweightsL = weights(-1,-1:1,1);
    preweightsR = weights(1,-1:1,1);
    preweights = weights(0,-1:1,1);
    Dymat(1,1:3) = preweightsL(2,:)/h;
    Dymat(N1D,N1D-2:N1D) = preweightsR(2,:)/h;
    for counter1 = 2:N1D-1
        Dymat(counter1,counter1-1:counter1+1) = preweights(2,:)/h;
    end
    Dxmat = Dymat';
end

yLI = 0.6;
yUI = 0.8;

k1 = 1;
k2 = 1/5;

kMat = ones(size(X));
kMat = kMat+(Y > 0.6).*(Y < 0.8)*(k2-k1);

timeSigma = 1;

C1 = (pi^2+timeSigma/k1)^0.5;
C2 = (pi^2+timeSigma/k2)^0.5;
sConst = 1;

A = [1 1 0 0 0 0;
    exp(C1*yLI) exp(-C1*yLI) -exp(C2*yLI) -exp(-C2*yLI) 0 0;
    C1*k1*exp(C1*yLI) -C1*k1*exp(-C1*yLI) -C2*k2*exp(C2*yLI) C2*k2*exp(-C2*yLI) 0 0;
    0 0 exp(C2*yUI) exp(-C2*yUI) -exp(C1*yUI) -exp(-C1*yUI);
    0 0 C2*k2*exp(C2*yUI) -C2*k2*exp(-C2*yUI) -C1*k1*exp(C1*yUI) C1*k1*exp(-C1*yUI);
    0 0 0 0 exp(C1) exp(-C1)];

cVec2 = A\([0 0 0 0 0 1]');

Z1 = sin(sConst*pi*X).*(cVec2(1,1)*exp(Y*C1)+cVec2(2,1)*exp(-Y*C1));
Z2 = sin(sConst*pi*X).*(cVec2(3,1)*exp(Y*C2)+cVec2(4,1)*exp(-Y*C2));
Z3 = sin(sConst*pi*X).*(cVec2(5,1)*exp(Y*C1)+cVec2(6,1)*exp(-Y*C1));

dataCurrent = (Y < yLI).*Z1 + (Y >= yLI).*(Y < yUI).*Z2 + ...
    (Y >= yUI).*Z3;

datastore = zeros(N,size(storesteps,1));
datastore(:,1) = reshape(dataCurrent,N,1);

xnodes = reshape(X,N,1);
ynodes = reshape(Y,N,1);

boundNull = ones(size(X));
boundNull(1:N1D,1) = zeros(N1D,1);
boundNull(1:N1D,N1D) = zeros(N1D,1);
boundNull(1,1:N1D) = zeros(1,N1D);
boundNull(N1D,1:N1D) = zeros(1,N1D);

boundFunc = zeros(size(X));
boundFunc(N1D,1:N1D) = sin(pi*X(N1D,1:N1D));

Dymat = sparse(Dymat);
Dxmat = sparse(Dxmat);

for time_step = 2:timesteps
    
    time_step
    expTime = (time_step-2)*k;
    
    %%% RK4 step 1
    
    d1 = Dymat*(kMat.*(Dymat*dataCurrent))+(kMat.*(dataCurrent*Dxmat))*Dxmat;
    step2data = dataCurrent+k/2*d1;
    step2data = boundNull.*step2data + exp(timeSigma*(expTime+k/2))*boundFunc;
    
    d2 = Dymat*(kMat.*(Dymat*step2data))+(kMat.*(step2data*Dxmat))*Dxmat;
    step3data = dataCurrent+k/2*d2;
    step3data = boundNull.*step3data + exp(timeSigma*(expTime+k/2))*boundFunc;
    
    d3 = Dymat*(kMat.*(Dymat*step3data))+(kMat.*(step3data*Dxmat))*Dxmat;
    step4data = dataCurrent+k*d3;
    step4data = boundNull.*step4data + exp(timeSigma*(expTime+k))*boundFunc;
    
    d4 = Dymat*(kMat.*(Dymat*step4data))+(kMat.*(step4data*Dxmat))*Dxmat;
    
    dataCurrent = dataCurrent+k/6*(d1+2*d2+2*d3+d4);
    dataCurrent = boundNull.*dataCurrent + exp(timeSigma*(expTime+k))*boundFunc;
    
    for m = 1:size(storesteps,1)
        if storesteps(m,1) == time_step
            datastore(:,m) = reshape(dataCurrent,N,1);
        end
    end
    
end

end

