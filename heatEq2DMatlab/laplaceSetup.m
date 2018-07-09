k1 = rho1;
k2 = rho2;
trialSigma = 0;

if closedFlag > 0
    C1 = (pi^2+trialSigma/k1)^0.5;
    C2 = (pi^2+trialSigma/k2)^0.5;
    sConst = 1;
else
    C1 = (4*pi^2+trialSigma/k1)^0.5;
    C2 = (4*pi^2+trialSigma/k2)^0.5;
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

dataExact = (ynodes < yLI).*Z1 + (ynodes >= yLI).*(ynodes < yUI).*Z2 + ...
    (ynodes >= yUI).*Z3;

R1 = (C1^2-(sConst*pi)^2)*sin(sConst*pi*xnodes).*(cVec2(1,1)*exp(ynodes*C1)+cVec2(2,1)*exp(-ynodes*C1))*k1;
R2 = (C2^2-(sConst*pi)^2)*sin(sConst*pi*xnodes).*(cVec2(3,1)*exp(ynodes*C2)+cVec2(4,1)*exp(-ynodes*C2))*k2;
R3 = (C1^2-(sConst*pi)^2)*sin(sConst*pi*xnodes).*(cVec2(5,1)*exp(ynodes*C1)+cVec2(6,1)*exp(-ynodes*C1))*k1;

RHS = (ynodes < yLI).*R1 + (ynodes >= yLI).*(ynodes < yUI).*R2 + ...
    (ynodes >= yUI).*R3;

RHSa = (ynodes < yLI).*R1/k1^2 + (ynodes >= yLI).*(ynodes < yUI).*R2/k2^2 + ...
    (ynodes >= yUI).*R3/k1^2;

boundVec = [zeros(N-numIntNodes,1); sin(sConst*pi*xnodes(N-numIntNodes+1:N,1))];
boundNullVec = [ones(N-numIntNodes,1); zeros(numIntNodes,1)];

Ltrunc = bigLmatrix(1:N-4*numIntNodes,1:N-4*numIntNodes);
RHStruncTemp = RHS-bigLmatrix*boundVec;
RHStrunc = RHStruncTemp(1:N-4*numIntNodes,1); 

intExp = Ltrunc\RHStrunc;
dataexp = [intExp; zeros(3*numIntNodes,1); sin(sConst*pi*xnodes(N-numIntNodes+1:N,1))];

% we combine the block differential operator with a block hyperviscosity
% operator to make one unified operator with which we can time step.

error = dataExact-dataexp;

L2error = sum(error.^2/N).^0.5
maxerror = max(abs(error))

figure(1);
resultssurf(xnodes,ynodes,dataexp,0,1,0,1,0.1,200);

figure(2);
resultssurf(xnodes,ynodes,error,0,1,0,1,0.1,200);

%{
figure(1);
resultssurf(xnodes,ynodes,RHS,0,1,0,1,0.1,200);

figure(2);
resultssurf(xnodes,ynodes,RHSa,0,1,0,1,0.1,200);
%}