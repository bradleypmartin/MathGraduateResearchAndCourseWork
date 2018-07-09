% Script RESULTSDISPLAY helps display data from RBF-FD modeling of wave
% transport (2D elastic wave equation).

timeindex = 11;
%
timeSigma = 1;
rho1 = 1;
rho2 = 1/5;
closedFlag = 0;
yLI = 0.6;
yUI = 0.8;
%}
%k = 0.001;
% Index of the time point at which you want to examine
% data (index in the data storage vector)

N = size(xnodes,1);

dataexp = datastore(:,timeindex);

figure(1)
resultssurf(xnodes,ynodes,...
    real(dataexp),0,1,0,1,0.2,200);
colormap('jet')
colorbar
xlabel('x'); ylabel('y')
%}

if closedFlag > 0
    C1 = (pi^2+timeSigma/rho1)^0.5;
    C2 = (pi^2+timeSigma/rho2)^0.5;
    sConst = 1;
else
    C1 = (4*pi^2+timeSigma/rho1)^0.5;
    C2 = (4*pi^2+timeSigma/rho2)^0.5;
    sConst = 2;
end

A = [1 1 0 0 0 0;
    exp(C1*yLI) exp(-C1*yLI) -exp(C2*yLI) -exp(-C2*yLI) 0 0;
    C1*rho1*exp(C1*yLI) -C1*rho1*exp(-C1*yLI) -C2*rho2*exp(C2*yLI) C2*rho2*exp(-C2*yLI) 0 0;
    0 0 exp(C2*yUI) exp(-C2*yUI) -exp(C1*yUI) -exp(-C1*yUI);
    0 0 C2*rho2*exp(C2*yUI) -C2*rho2*exp(-C2*yUI) -C1*rho1*exp(C1*yUI) C1*rho1*exp(-C1*yUI);
    0 0 0 0 exp(C1) exp(-C1)];

cVec2 = A\([0 0 0 0 0 1]');

Z1 = exp(0.1*timeSigma)*sin(sConst*pi*xnodes).*(cVec2(1,1)*exp(ynodes*C1)+cVec2(2,1)*exp(-ynodes*C1));
Z2 = exp(0.1*timeSigma)*sin(sConst*pi*xnodes).*(cVec2(3,1)*exp(ynodes*C2)+cVec2(4,1)*exp(-ynodes*C2));
Z3 = exp(0.1*timeSigma)*sin(sConst*pi*xnodes).*(cVec2(5,1)*exp(ynodes*C1)+cVec2(6,1)*exp(-ynodes*C1));

dataAn = (ynodes < yLI).*Z1 + (ynodes >= yLI).*(ynodes < yUI).*Z2 + ...
    (ynodes >= yUI).*Z3;

error = dataAn-dataexp;

figure(2)
resultssurf(xnodes,ynodes,...
    real(dataAn),0,1,0,1,0.2,200);
colormap('jet')
colorbar
xlabel('x'); ylabel('y')

figure(3)
resultssurf(xnodes,ynodes,...
    real(error),0,1,0,1,0.2,200);
colormap('jet')
colorbar
xlabel('x'); ylabel('y')
%caxis([-0.01 0.01])

L2error = sum(error.^2/N).^0.5
maxerror = max(abs(error))
