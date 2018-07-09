function [xgrid ufxtKP1 bigOperator bigJacobian F] = ...
    FD4STInv9(N,tstep,endtime,dampConst,naiveFlag,...
    knownu,knownf,knownuInd,knownfInd,knownuAlpha,knownfAlpha,...
    ufxtK,CGflag,CGiterations,spaceSigma,timeSigma,exLoc,exTime,...
    FOflag,FOalpha,FOloc,ZOflag,ZOalpha,ZOradius,...
    intKnownFlag,intStack,RhoKnownflag,Rhoalpha,ZOswitchu,ZOswitchf,...
        knownRhoVec)

% Edit 4 - 140404 - changing code to 1-D simplification of 2-D EWE;
% Edit 5 - 140405 - enabling toggle between naive and AC solution
% 140408 - change to first iteration of inversion
% 140410 - inclusion of CG solver for Newton step

% Edit 6 - 140410 - change to Ricker wavelet forcing cond. on f
% Edit 7 - 140413 - first trial of 1st order regularization
% Edit 8 - 140416 - restricts 1st order reg. to end of domain
% Edit 9 - 140416 - adds zeroth order reg. at knownInd and exLoc

h = 1/N;
numInterfaces = size(intStack,1);

uxM1 = zeros(numInterfaces,5);
uxM2 = zeros(numInterfaces,5);
uxM3 = zeros(numInterfaces,5);
uxM4 = zeros(numInterfaces,5);

fxM1 = zeros(numInterfaces,5);
fxM2 = zeros(numInterfaces,5);
fxM3 = zeros(numInterfaces,5);
fxM4 = zeros(numInterfaces,5);

for interfaceCounter = 1:numInterfaces

    c1 = intStack(interfaceCounter,4);
    c2 = intStack(interfaceCounter,5);

    p1 = intStack(interfaceCounter,2);
    p2 = intStack(interfaceCounter,3);

    superStackMid = zeros(5,6,4,2);  % u = 1; v = 2;

    for l = 1:2
        for k = 1:4
            for i = 1:5
                for j = 1:5
                    constMult = 1;
                    if j <= 5-k

                        if l == 1 % u conditions
                            if i == 2
                                constMult = p2/p1*(c2/c1)^2;
                            end
                            if i == 3
                                constMult = (c2/c1)^2;
                            end
                            if i == 4
                                constMult = p2/p1*(c2/c1)^4;
                            end
                            if i == 5
                                constMult = (c2/c1)^4;
                            end
                        end

                        if l == 2 % v condtions
                            if i == 2
                                constMult = p1/p2;
                            end
                            if i == 3
                                constMult = (c2/c1)^2;
                            end
                            if i == 4
                                constMult = p1/p2*(c2/c1)^2;
                            end
                            if i == 5
                                constMult = (c2/c1)^4;
                            end
                        end
                    end
                    superStackMid(i,j,k,l)= constMult*(-1.5-4+j+k)^(i-1);
                    
                end

                % rhs alteration
                constMult = 1;
                if k < 3
                    if l == 1 % u conditions
                        if i == 2
                            constMult = p2/p1*(c2/c1)^2;
                        end 
                        if i == 3
                            constMult = (c2/c1)^2;
                        end
                        if i == 4
                            constMult = p2/p1*(c2/c1)^4;
                        end
                        if i == 5
                            constMult = (c2/c1)^4;
                        end
                    end

                    if l == 2 % v condtions
                        if i == 2
                            constMult = p1/p2;
                        end
                        if i == 3
                            constMult = (c2/c1)^2;
                        end
                        if i == 4
                            constMult = p1/p2*(c2/c1)^2;
                        end
                        if i == 5
                            constMult = (c2/c1)^4;
                        end
                    end
                end

                superStackMid(i,6,k,l) = constMult*(i-1)*(-2.5+k)^(i-2);

            end
        end
    end

    weightsMid = zeros(5,4,2);
    

    for k = 1:4
        for l = 1:2
            weightsMid(:,k,l) = superStackMid(:,1:5,k,l)\superStackMid(:,6,k,l);
        end
    end
    
    uxM1(interfaceCounter,:) = weightsMid(:,1,1)'/h;
    uxM2(interfaceCounter,:) = weightsMid(:,2,1)'/h;
    uxM3(interfaceCounter,:) = weightsMid(:,3,1)'/h;
    uxM4(interfaceCounter,:) = weightsMid(:,4,1)'/h;

    fxM1(interfaceCounter,:) = weightsMid(:,1,2)'/h;
    fxM2(interfaceCounter,:) = weightsMid(:,2,2)'/h;
    fxM3(interfaceCounter,:) = weightsMid(:,3,2)'/h;
    fxM4(interfaceCounter,:) = weightsMid(:,4,2)'/h;
    
    weightsC = weights(0,-2:2,1);
    ufxC = weightsC(2,1:5)/h;

    if naiveFlag == 1
        uxM1(interfaceCounter,:) = ufxC; 
        uxM2(interfaceCounter,:) = ufxC; 
        uxM3(interfaceCounter,:) = ufxC; 
        uxM4(interfaceCounter,:) = ufxC;
        fxM1(interfaceCounter,:) = ufxC; 
        fxM2(interfaceCounter,:) = ufxC; 
        fxM3(interfaceCounter,:) = ufxC; 
        fxM4(interfaceCounter,:) = ufxC;
    end
    
end

xgrid = linspace(0,1,N+1)';
xgrid = xgrid(1:N,1)+h/2;

dampMult = zeros(size(xgrid));

intIndices = zeros(numInterfaces,1);
tempCounter = 1;

for m = 1:N
    for intCounter = 1:numInterfaces
        if xgrid(m,1) > intStack(intCounter,1)
            if m > 1
                if xgrid(m-1,1) < intStack(intCounter,1)
                    intIndices(tempCounter,1) = m-1;
                    tempCounter = tempCounter+1;
                end
            end
        end
    end
    if xgrid(m,1)>0.9
        dampMult(m,1) = (10*(xgrid(m,1)-0.9))^2*dampConst;
    end
end

timesteps = ceil(endtime/tstep);

%ux0=exp(-1200*(xgrid-1/4).^2);
%fx0=-exp(-1200*(xgrid-1/4).^2);

ux0 = zeros(size(xgrid));
fx0 = zeros(size(xgrid));

%[1/12 -2/3 0 2/3 -1/12]/h

weightsC = weights(0,-2:2,1);
weightsL = weights(-1,-2:2,1);
weightsR = weights(1,-2:2,1);
weights2L = weights(-1.5,-2:2,1);
weights2R = weights(1.5,-2:2,1);

ux2L = weightsL(2,2:5);
ux2L = (ux2L-weightsL(2,1)*weights2L(2,2:5)/weights2L(2,1))/h;

uxL = weightsC(2,2:5);
uxL = (uxL-weightsC(2,1)*weights2L(2,2:5)/weights2L(2,1))/h;

uxR = weightsC(2,1:4);
uxR = (uxR-weightsC(2,5)*weights2R(2,1:4)/weights2R(2,5))/h;

ux2R = weightsR(2,1:4);
ux2R = (ux2R-weightsR(2,5)*weights2R(2,1:4)/weights2R(2,5))/h;

fx2L = weightsL(2,2:5);
fx2L = (fx2L-weightsL(2,1)*weights2L(1,2:5)/weights2L(1,1))/h;

fxL = weightsC(2,2:5);
fxL = (fxL-weightsC(2,1)*weights2L(1,2:5)/weights2L(1,1))/h;

fxR = weightsC(2,1:4);
fxR = (fxR-weightsC(2,5)*weights2R(1,1:4)/weights2R(1,5))/h;

fx2R = weightsR(2,1:4);
fx2R = (fx2R-weightsR(2,5)*weights2R(1,1:4)/weights2R(1,5))/h;

ufxC = weightsC(2,1:5)/h;

weightsT1 = weights(-1,-2:2,1)/tstep;
weightsTC = weights(0,-2:2,1)/tstep;
weightsTEm1 = weights(1,-2:2,1)/tstep;
weightsTE = weights(2,-2:2,1)/tstep;

uftT1 = weightsT1(2,1:5);
uftC = weightsTC(2,1:5);
uftEm1 = weightsTEm1(2,1:5);
uftE = weightsTE(2,1:5);

bigSize = N*timesteps;

RHS = zeros(2*bigSize+size(knownu,1)*size(knownu,2)+...
    size(knownf,1)*size(knownf,2),1);
sparsei = [];
sparsej = [];
sparseval = [];
counter1 = 1;

sparseiJac = [];
sparsejJac = [];
sparsevalJac = [];
counter2 = 1;

for m = 1:bigSize
    xindex = mod(m-1,N)+1;
    tindex = floor((m-1)/N)+1;
    if tindex == 1  % first timestep
        if xindex < 3 % left side
            if xindex == 1 % leftmost side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:4                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-2);
                    sparsevalJac(counter2,1) = uftT1(1,k+1);
                    counter2 = counter2+1;
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparseval(counter1,1) = -fx2L(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparsevalJac(counter2,1) = -fx2L(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fx2L(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-1),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);
                
                % ft-c2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-2);
                    sparsevalJac(counter2,1) = uftT1(1,k+1);
                    counter2 = counter2+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparseval(counter1,1) = -ux2L(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparsevalJac(counter2,1) = -ux2L(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -ux2L(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-1),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m+bigSize,1) = -uftT1(1,1)*fx0(xindex,1);
                
            end
            
            if xindex == 2 % next-to-left side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:4                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-2);
                    sparsevalJac(counter2,1) = uftT1(1,k+1);
                    counter2 = counter2+1;
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-2),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);
                
                % ft-c2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-2);
                    sparsevalJac(counter2,1) = uftT1(1,k+1);
                    counter2 = counter2+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-2),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m+bigSize,1) = -uftT1(1,1)*fx0(xindex,1);
                
            end
        end
        
        if xindex > N-2 % right side
            if xindex == N % right-most side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:4                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-2);
                    sparsevalJac(counter2,1) = uftT1(1,k+1);
                    counter2 = counter2+1;
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparseval(counter1,1) = -fx2R(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparsevalJac(counter2,1) = -fx2R(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fx2R(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-4),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-2);
                    sparsevalJac(counter2,1) = uftT1(1,k+1);
                    counter2 = counter2+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparseval(counter1,1) = -ux2R(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparsevalJac(counter2,1) = -ux2R(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -ux2R(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-4),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m+bigSize,1) = -uftT1(1,1)*fx0(xindex,1);
                
            end
            
            if xindex == N-1 % next-to-left side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:4                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-2);
                    sparsevalJac(counter2,1) = uftT1(1,k+1);
                    counter2 = counter2+1;
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxR(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparsevalJac(counter2,1) = -fxR(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fxR(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);
                
                % ft-c2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-2);
                    sparsevalJac(counter2,1) = uftT1(1,k+1);
                    counter2 = counter2+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxR(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparsevalJac(counter2,1) = -uxR(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -uxR(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m+bigSize,1) = -uftT1(1,1)*fx0(xindex,1);
                
            end
        end
        
        if xindex > 2 && xindex < N-1 % center
            tempIntFlag = 0;
            for intCounter = 1:numInterfaces
                if xindex > intIndices(intCounter,1)-2 && ...
                        xindex < intIndices(intCounter,1)+3
                    tempIntFlag = intCounter;
                end
            end
                if tempIntFlag > 0
                    if xindex == intIndices(tempIntFlag,1)-1
                        uxTemp = uxM1(tempIntFlag,:);
                        fxTemp = fxM1(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)
                        uxTemp = uxM2(tempIntFlag,:);
                        fxTemp = fxM2(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)+1
                        uxTemp = uxM3(tempIntFlag,:);
                        fxTemp = fxM3(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)+2
                        uxTemp = uxM4(tempIntFlag,:);
                        fxTemp = fxM4(tempIntFlag,:);
                    end
                    for k = 1:4                   
                        % ut
                        sparsei(counter1,1) = m;      
                        sparsej(counter1,1) = xindex+2*N*(tindex+k-2);
                        sparseval(counter1,1) = uftT1(1,k+1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;      
                        sparsejJac(counter2,1) = xindex+2*N*(tindex+k-2);
                        sparsevalJac(counter2,1) = uftT1(1,k+1);
                        counter2 = counter2+1;
                    end
                    for k = 1:5
                        % -fx/p
                        sparsei(counter1,1) = m;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparseval(counter1,1) = -fxTemp(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparsevalJac(counter2,1) = -fxTemp(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*bigSize;
                        sparsevalJac(counter2,1) = -fxTemp(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*u
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1);
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);

                    % ft-pc2ux+dampMult*f = RHS        
                    for k = 1:4                   
                        % ft
                        sparsei(counter1,1) = m+bigSize;      
                        sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                        sparseval(counter1,1) = uftT1(1,k+1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;      
                        sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-2);
                        sparsevalJac(counter2,1) = uftT1(1,k+1);
                        counter2 = counter2+1;
                    end
                    for k = 1:5
                        % -pc2ux
                        sparsei(counter1,1) = m+bigSize;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparseval(counter1,1) = -uxTemp(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparsevalJac(counter2,1) = -uxTemp(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+N+2*bigSize;
                        sparsevalJac(counter2,1) = -uxTemp(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*f
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    RHS(m+bigSize,1) = -uftT1(1,1)*fx0(xindex,1);
                else
                % ut-fx/p+dampMult*u = RHS        
                    for k = 1:4                   
                        % ut
                        sparsei(counter1,1) = m;      
                        sparsej(counter1,1) = xindex+2*N*(tindex+k-2);
                        sparseval(counter1,1) = uftT1(1,k+1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;      
                        sparsejJac(counter2,1) = xindex+2*N*(tindex+k-2);
                        sparsevalJac(counter2,1) = uftT1(1,k+1);
                        counter2 = counter2+1;
                    end
                    for k = 1:5
                        % -fx/p
                        sparsei(counter1,1) = m;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparseval(counter1,1) = -ufxC(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*bigSize;
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*u
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1);
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);

                    % ft-pc2ux+dampMult*f = RHS        
                    for k = 1:4                   
                        % ft
                        sparsei(counter1,1) = m+bigSize;      
                        sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                        sparseval(counter1,1) = uftT1(1,k+1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;      
                        sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-2);
                        sparsevalJac(counter2,1) = uftT1(1,k+1);
                        counter2 = counter2+1;
                    end
                    for k = 1:5
                        % -pc2ux
                        sparsei(counter1,1) = m+bigSize;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparseval(counter1,1) = -ufxC(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+N+2*bigSize;
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*f
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    RHS(m+bigSize,1) = -uftT1(1,1)*fx0(xindex,1);
                end
            
        end
    end
    
    if tindex == 2  % second timestep
        if xindex < 3 % left side
            if xindex == 1 % leftmost side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:4                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-3);
                    sparsevalJac(counter2,1) = uftC(1,k+1);
                    counter2 = counter2+1;
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparseval(counter1,1) = -fx2L(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparsevalJac(counter2,1) = -fx2L(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fx2L(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-1),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-3);
                    sparsevalJac(counter2,1) = uftC(1,k+1);
                    counter2 = counter2+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparseval(counter1,1) = -ux2L(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparsevalJac(counter2,1) = -ux2L(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -ux2L(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-1),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
            
            if xindex == 2 % next-to-left side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:4                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-3);
                    sparsevalJac(counter2,1) = uftC(1,k+1);
                    counter2 = counter2+1;
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-2),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-3);
                    sparsevalJac(counter2,1) = uftC(1,k+1);
                    counter2 = counter2+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparseval(counter1,1) = -uxL(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparsevalJac(counter2,1) = -uxL(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -uxL(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-2),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
        end
        
        if xindex > N-2 % right side
            if xindex == N % right-most side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:4                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-3);
                    sparsevalJac(counter2,1) = uftC(1,k+1);
                    counter2 = counter2+1;
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparseval(counter1,1) = -fx2R(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparsevalJac(counter2,1) = -fx2R(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fx2R(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-4),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-3);
                    sparsevalJac(counter2,1) = uftC(1,k+1);
                    counter2 = counter2+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparseval(counter1,1) = -ux2R(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparsevalJac(counter2,1) = -ux2R(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -ux2R(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-4),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
            
            if xindex == N-1 % next-to-left side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:4                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-3);
                    sparsevalJac(counter2,1) = uftC(1,k+1);
                    counter2 = counter2+1;
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxR(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparsevalJac(counter2,1) = -fxR(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fxR(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-3);
                    sparsevalJac(counter2,1) = uftC(1,k+1);
                    counter2 = counter2+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxR(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparsevalJac(counter2,1) = -uxR(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -uxR(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
        end
        
        if xindex > 2 && xindex < N-1 % center
            tempIntFlag = 0;
            for intCounter = 1:numInterfaces
                if xindex > intIndices(intCounter,1)-2 && ...
                        xindex < intIndices(intCounter,1)+3
                    tempIntFlag = intCounter;
                end
            end
                if tempIntFlag > 0
                    if xindex == intIndices(tempIntFlag,1)-1
                        uxTemp = uxM1(tempIntFlag,:);
                        fxTemp = fxM1(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)
                        uxTemp = uxM2(tempIntFlag,:);
                        fxTemp = fxM2(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)+1
                        uxTemp = uxM3(tempIntFlag,:);
                        fxTemp = fxM3(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)+2
                        uxTemp = uxM4(tempIntFlag,:);
                        fxTemp = fxM4(tempIntFlag,:);
                    end
                    % ut-fx/p+dampMult*u = RHS        
                    for k = 1:4                   
                        % ut
                        sparsei(counter1,1) = m;      
                        sparsej(counter1,1) = xindex+2*N*(tindex+k-3);
                        sparseval(counter1,1) = uftC(1,k+1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;      
                        sparsejJac(counter2,1) = xindex+2*N*(tindex+k-3);
                        sparsevalJac(counter2,1) = uftC(1,k+1);
                        counter2 = counter2+1;
                    end
                    for k = 1:5
                        % -fx/p
                        sparsei(counter1,1) = m;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparseval(counter1,1) = -fxTemp(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparsevalJac(counter2,1) = -fxTemp(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*bigSize;
                        sparsevalJac(counter2,1) = -fxTemp(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*u
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1);
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                    % ft-pc2ux+dampMult*f = RHS        
                    for k = 1:4                   
                        % ft
                        sparsei(counter1,1) = m+bigSize;      
                        sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                        sparseval(counter1,1) = uftC(1,k+1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;      
                        sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-3);
                        sparsevalJac(counter2,1) = uftC(1,k+1);
                        counter2 = counter2+1;
                    end
                    for k = 1:5
                        % -pc2ux
                        sparsei(counter1,1) = m+bigSize;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparseval(counter1,1) = -uxTemp(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparsevalJac(counter2,1) = -uxTemp(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+N+2*bigSize;
                        sparsevalJac(counter2,1) = -uxTemp(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*f
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                else
                    % ut-fx/p+dampMult*u = RHS        
                    for k = 1:4                   
                        % ut
                        sparsei(counter1,1) = m;      
                        sparsej(counter1,1) = xindex+2*N*(tindex+k-3);
                        sparseval(counter1,1) = uftC(1,k+1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;      
                        sparsejJac(counter2,1) = xindex+2*N*(tindex+k-3);
                        sparsevalJac(counter2,1) = uftC(1,k+1);
                        counter2 = counter2+1;
                    end
                    for k = 1:5
                        % -fx/p
                        sparsei(counter1,1) = m;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparseval(counter1,1) = -ufxC(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*bigSize;
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*u
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1);
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                    % ft-pc2ux+dampMult*f = RHS        
                    for k = 1:4                   
                        % ft
                        sparsei(counter1,1) = m+bigSize;      
                        sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                        sparseval(counter1,1) = uftC(1,k+1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;      
                        sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-3);
                        sparsevalJac(counter2,1) = uftC(1,k+1);
                        counter2 = counter2+1;
                    end
                    for k = 1:5
                        % -pc2ux
                        sparsei(counter1,1) = m+bigSize;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparseval(counter1,1) = -ufxC(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+N+2*bigSize;
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*f
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                end
            
        end
    end
    
    if tindex > 2 && tindex < timesteps-1 % middle timesteps
        if xindex < 3 % left side
            if xindex == 1 % leftmost side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-4);
                    sparsevalJac(counter2,1) = uftC(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparseval(counter1,1) = -fx2L(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparsevalJac(counter2,1) = -fx2L(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fx2L(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-1),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-4);
                    sparsevalJac(counter2,1) = uftC(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparseval(counter1,1) = -ux2L(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparsevalJac(counter2,1) = -ux2L(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -ux2L(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-1),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
            
            if xindex == 2 % next-to-left side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-4);
                    sparsevalJac(counter2,1) = uftC(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-2),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-4);
                    sparsevalJac(counter2,1) = uftC(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparseval(counter1,1) = -uxL(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparsevalJac(counter2,1) = -uxL(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -uxL(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-2),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
        end
        
        if xindex > N-2 % right side
            if xindex == N % right-most side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-4);
                    sparsevalJac(counter2,1) = uftC(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparseval(counter1,1) = -fx2R(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparsevalJac(counter2,1) = -fx2R(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fx2R(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-4),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-4);
                    sparsevalJac(counter2,1) = uftC(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparseval(counter1,1) = -ux2R(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparsevalJac(counter2,1) = -ux2R(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -ux2R(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-4),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
            
            if xindex == N-1 % next-to-left side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-4);
                    sparsevalJac(counter2,1) = uftC(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxR(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparsevalJac(counter2,1) = -fxR(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fxR(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-4);
                    sparsevalJac(counter2,1) = uftC(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxR(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparsevalJac(counter2,1) = -uxR(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -uxR(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
        end
        
        if xindex > 2 && xindex < N-1 % center
            tempIntFlag = 0;
            for intCounter = 1:numInterfaces
                if xindex > intIndices(intCounter,1)-2 && ...
                            xindex < intIndices(intCounter,1)+3
                    tempIntFlag = intCounter;    
                end
            end
                if tempIntFlag > 0
                    if xindex == intIndices(tempIntFlag,1)-1
                        uxTemp = uxM1(tempIntFlag,:);
                        fxTemp = fxM1(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)
                        uxTemp = uxM2(tempIntFlag,:);
                        fxTemp = fxM2(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)+1
                        uxTemp = uxM3(tempIntFlag,:);
                        fxTemp = fxM3(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)+2
                        uxTemp = uxM4(tempIntFlag,:);
                        fxTemp = fxM4(tempIntFlag,:);
                    end
                    % ut-fx/p+dampMult*u = RHS        
                    for k = 1:5                   
                        % ut
                        sparsei(counter1,1) = m;      
                        sparsej(counter1,1) = xindex+2*N*(tindex+k-4);
                        sparseval(counter1,1) = uftC(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;      
                        sparsejJac(counter2,1) = xindex+2*N*(tindex+k-4);
                        sparsevalJac(counter2,1) = uftC(1,k);
                        counter2 = counter2+1;

                        % -fx/p
                        sparsei(counter1,1) = m;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparseval(counter1,1) = -fxTemp(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparsevalJac(counter2,1) = -fxTemp(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*bigSize;
                        sparsevalJac(counter2,1) = -fxTemp(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*u
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1);
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                    % ft-pc2ux+dampMult*f = RHS        
                    for k = 1:5                   
                        % ft
                        sparsei(counter1,1) = m+bigSize;      
                        sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                        sparseval(counter1,1) = uftC(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;      
                        sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-4);
                        sparsevalJac(counter2,1) = uftC(1,k);
                        counter2 = counter2+1;

                        % -pc2ux
                        sparsei(counter1,1) = m+bigSize;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparseval(counter1,1) = -uxTemp(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparsevalJac(counter2,1) = -uxTemp(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+N+2*bigSize;
                        sparsevalJac(counter2,1) = -uxTemp(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*f
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                else
                    % ut-fx/p+dampMult*u = RHS        
                    for k = 1:5                   
                        % ut
                        sparsei(counter1,1) = m;      
                        sparsej(counter1,1) = xindex+2*N*(tindex+k-4);
                        sparseval(counter1,1) = uftC(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;      
                        sparsejJac(counter2,1) = xindex+2*N*(tindex+k-4);
                        sparsevalJac(counter2,1) = uftC(1,k);
                        counter2 = counter2+1;

                        % -fx/p
                        sparsei(counter1,1) = m;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparseval(counter1,1) = -ufxC(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*bigSize;
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*u
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1);
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                    % ft-pc2ux+dampMult*f = RHS        
                    for k = 1:5                   
                        % ft
                        sparsei(counter1,1) = m+bigSize;      
                        sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                        sparseval(counter1,1) = uftC(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;      
                        sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-4);
                        sparsevalJac(counter2,1) = uftC(1,k);
                        counter2 = counter2+1;

                        % -pc2ux
                        sparsei(counter1,1) = m+bigSize;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparseval(counter1,1) = -ufxC(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+N+2*bigSize;
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*f
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                end
            
        end
    end
    
    if tindex == timesteps-1  % second-to-last timestep
        if xindex < 3 % left side
            if xindex == 1 % leftmost side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-5);
                    sparsevalJac(counter2,1) = uftEm1(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparseval(counter1,1) = -fx2L(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparsevalJac(counter2,1) = -fx2L(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fx2L(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-1),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-5);
                    sparsevalJac(counter2,1) = uftEm1(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparseval(counter1,1) = -ux2L(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparsevalJac(counter2,1) = -ux2L(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -ux2L(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-1),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
            
            if xindex == 2 % next-to-left side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-5);
                    sparsevalJac(counter2,1) = uftEm1(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-2),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-5);
                    sparsevalJac(counter2,1) = uftEm1(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparseval(counter1,1) = -uxL(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparsevalJac(counter2,1) = -uxL(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -uxL(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-2),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
        end
        
        if xindex > N-2 % right side
            if xindex == N % right-most side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-5);
                    sparsevalJac(counter2,1) = uftEm1(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparseval(counter1,1) = -fx2R(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparsevalJac(counter2,1) = -fx2R(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fx2R(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-4),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-5);
                    sparsevalJac(counter2,1) = uftEm1(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparseval(counter1,1) = -ux2R(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparsevalJac(counter2,1) = -ux2R(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -ux2R(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-4),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
            
            if xindex == N-1 % next-to-left side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-5);
                    sparsevalJac(counter2,1) = uftEm1(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxR(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparsevalJac(counter2,1) = -fxR(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fxR(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-5);
                    sparsevalJac(counter2,1) = uftEm1(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxR(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparsevalJac(counter2,1) = -uxR(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -uxR(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
        end
        
        if xindex > 2 && xindex < N-1 % center
            tempIntFlag = 0;
            for intCounter = 1:numInterfaces
                if xindex > intIndices(intCounter,1)-2 && ...
                        xindex < intIndices(intCounter,1)+3 
                    tempIntFlag = intCounter;
                end
            end
                if tempIntFlag > 0
                    if xindex == intIndices(tempIntFlag,1)-1
                        uxTemp = uxM1(tempIntFlag,:);
                        fxTemp = fxM1(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)
                        uxTemp = uxM2(tempIntFlag,:);
                        fxTemp = fxM2(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)+1
                        uxTemp = uxM3(tempIntFlag,:);
                        fxTemp = fxM3(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)+2
                        uxTemp = uxM4(tempIntFlag,:);
                        fxTemp = fxM4(tempIntFlag,:);
                    end
                    % ut-fx/p+dampMult*u = RHS        
                    for k = 1:5                   
                        % ut
                        sparsei(counter1,1) = m;      
                        sparsej(counter1,1) = xindex+2*N*(tindex+k-5);
                        sparseval(counter1,1) = uftEm1(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;      
                        sparsejJac(counter2,1) = xindex+2*N*(tindex+k-5);
                        sparsevalJac(counter2,1) = uftEm1(1,k);
                        counter2 = counter2+1;

                        % -fx/p
                        sparsei(counter1,1) = m;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparseval(counter1,1) = -fxTemp(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparsevalJac(counter2,1) = -fxTemp(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*bigSize;
                        sparsevalJac(counter2,1) = -fxTemp(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*u
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1);
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                    % ft-pc2ux+dampMult*f = RHS        
                    for k = 1:5                   
                        % ft
                        sparsei(counter1,1) = m+bigSize;      
                        sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                        sparseval(counter1,1) = uftEm1(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;      
                        sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-5);
                        sparsevalJac(counter2,1) = uftEm1(1,k);
                        counter2 = counter2+1;

                        % -pc2ux
                        sparsei(counter1,1) = m+bigSize;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparseval(counter1,1) = -uxTemp(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparsevalJac(counter2,1) = -uxTemp(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+N+2*bigSize;
                        sparsevalJac(counter2,1) = -uxTemp(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*f
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                else  
                    % ut-fx/p+dampMult*u = RHS        
                    for k = 1:5                   
                        % ut
                        sparsei(counter1,1) = m;      
                        sparsej(counter1,1) = xindex+2*N*(tindex+k-5);
                        sparseval(counter1,1) = uftEm1(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;      
                        sparsejJac(counter2,1) = xindex+2*N*(tindex+k-5);
                        sparsevalJac(counter2,1) = uftEm1(1,k);
                        counter2 = counter2+1;

                        % -fx/p
                        sparsei(counter1,1) = m;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparseval(counter1,1) = -ufxC(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*bigSize;
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*u
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1);
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                    % ft-pc2ux+dampMult*f = RHS        
                    for k = 1:5                   
                        % ft
                        sparsei(counter1,1) = m+bigSize;      
                        sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                        sparseval(counter1,1) = uftEm1(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;      
                        sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-5);
                        sparsevalJac(counter2,1) = uftEm1(1,k);
                        counter2 = counter2+1;

                        % -pc2ux
                        sparsei(counter1,1) = m+bigSize;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparseval(counter1,1) = -ufxC(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+N+2*bigSize;
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*f
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                end
            
        end
    end
    
    if tindex == timesteps  % last timestep
        if xindex < 3 % left side
            if xindex == 1 % leftmost side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-6);
                    sparsevalJac(counter2,1) = uftE(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparseval(counter1,1) = -fx2L(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparsevalJac(counter2,1) = -fx2L(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fx2L(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-1),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-6);
                    sparsevalJac(counter2,1) = uftE(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparseval(counter1,1) = -ux2L(1,k)*ufxtK(xindex+N+2*bigSize);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparsevalJac(counter2,1) = -ux2L(1,k)*ufxtK(xindex+N+2*bigSize);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -ux2L(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-1));
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
            
            if xindex == 2 % next-to-left side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-6);
                    sparsevalJac(counter2,1) = uftE(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fxL(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-2),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-6);
                    sparsevalJac(counter2,1) = uftE(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparseval(counter1,1) = -uxL(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparsevalJac(counter2,1) = -uxL(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -uxL(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-2),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
        end
        
        if xindex > N-2 % right side
            if xindex == N % right-most side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-6);
                    sparsevalJac(counter2,1) = uftE(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparseval(counter1,1) = -fx2R(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparsevalJac(counter2,1) = -fx2R(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fx2R(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-4),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-6);
                    sparsevalJac(counter2,1) = uftE(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparseval(counter1,1) = -ux2R(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparsevalJac(counter2,1) = -ux2R(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -ux2R(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-4),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
            
            if xindex == N-1 % next-to-left side
                
                % ut-fx/p+dampMult*u = RHS        
                for k = 1:5                   
                    % ut
                    sparsei(counter1,1) = m;      
                    sparsej(counter1,1) = xindex+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;      
                    sparsejJac(counter2,1) = xindex+2*N*(tindex+k-6);
                    sparsevalJac(counter2,1) = uftE(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxR(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparsevalJac(counter2,1) = -fxR(1,k)*ufxtK(xindex+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*bigSize;
                    sparsevalJac(counter2,1) = -fxR(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                    counter2 = counter2+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;      
                    sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-6);
                    sparsevalJac(counter2,1) = uftE(1,k);
                    counter2 = counter2+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxR(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter1 = counter1+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparsevalJac(counter2,1) = -uxR(1,k)*ufxtK(xindex+N+2*bigSize,1);
                    counter2 = counter2+1;
                    
                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+N+2*bigSize;
                    sparsevalJac(counter2,1) = -uxR(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                    counter2 = counter2+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                
                sparseiJac(counter2,1) = m+bigSize;
                sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                sparsevalJac(counter2,1) = dampMult(xindex,1);
                counter2 = counter2+1;
                % = RHS
                %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                
            end
        end
        
        if xindex > 2 && xindex < N-1 % center
            tempIntFlag = 0;
            for intCounter = 1:numInterfaces
                if xindex > intIndices(intCounter,1)-2 && ...
                        xindex < intIndices(intCounter,1)+3                
                    tempIntFlag = intCounter;
                end
            end
                if tempIntFlag > 0
                    if xindex == intIndices(tempIntFlag,1)-1
                        uxTemp = uxM1(tempIntFlag,:);
                        fxTemp = fxM1(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)
                        uxTemp = uxM2(tempIntFlag,:);
                        fxTemp = fxM2(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)+1
                        uxTemp = uxM3(tempIntFlag,:);
                        fxTemp = fxM3(tempIntFlag,:);
                    end
                    if xindex == intIndices(tempIntFlag,1)+2
                        uxTemp = uxM4(tempIntFlag,:);
                        fxTemp = fxM4(tempIntFlag,:);
                    end
                    % ut-fx/p+dampMult*u = RHS        
                    for k = 1:5                   
                        % ut
                        sparsei(counter1,1) = m;      
                        sparsej(counter1,1) = xindex+2*N*(tindex+k-6);
                        sparseval(counter1,1) = uftE(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;      
                        sparsejJac(counter2,1) = xindex+2*N*(tindex+k-6);
                        sparsevalJac(counter2,1) = uftE(1,k);
                        counter2 = counter2+1;

                        % -fx/p
                        sparsei(counter1,1) = m;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparseval(counter1,1) = -fxTemp(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparsevalJac(counter2,1) = -fxTemp(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*bigSize;
                        sparsevalJac(counter2,1) = -fxTemp(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*u
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1);
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                    % ft-pc2ux+dampMult*f = RHS        
                    for k = 1:5                   
                        % ft
                        sparsei(counter1,1) = m+bigSize;      
                        sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                        sparseval(counter1,1) = uftE(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;      
                        sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-6);
                        sparsevalJac(counter2,1) = uftE(1,k);
                        counter2 = counter2+1;

                        % -pc2ux
                        sparsei(counter1,1) = m+bigSize;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparseval(counter1,1) = -uxTemp(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparsevalJac(counter2,1) = -uxTemp(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+N+2*bigSize;
                        sparsevalJac(counter2,1) = -uxTemp(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*f
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                else
                    % ut-fx/p+dampMult*u = RHS        
                    for k = 1:5                   
                        % ut
                        sparsei(counter1,1) = m;      
                        sparsej(counter1,1) = xindex+2*N*(tindex+k-6);
                        sparseval(counter1,1) = uftE(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;      
                        sparsejJac(counter2,1) = xindex+2*N*(tindex+k-6);
                        sparsevalJac(counter2,1) = uftE(1,k);
                        counter2 = counter2+1;

                        % -fx/p
                        sparsei(counter1,1) = m;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparseval(counter1,1) = -ufxC(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N+(k-3);
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m;
                        sparsejJac(counter2,1) = xindex+2*bigSize;
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*N*(tindex-1)+N+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*u
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1);
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1);
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                    % ft-pc2ux+dampMult*f = RHS        
                    for k = 1:5                   
                        % ft
                        sparsei(counter1,1) = m+bigSize;      
                        sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                        sparseval(counter1,1) = uftE(1,k);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;      
                        sparsejJac(counter2,1) = xindex+N+2*N*(tindex+k-6);
                        sparsevalJac(counter2,1) = uftE(1,k);
                        counter2 = counter2+1;

                        % -pc2ux
                        sparsei(counter1,1) = m+bigSize;
                        sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparseval(counter1,1) = -ufxC(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter1 = counter1+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+(k-3);
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+N+2*bigSize,1);
                        counter2 = counter2+1;

                        sparseiJac(counter2,1) = m+bigSize;
                        sparsejJac(counter2,1) = xindex+N+2*bigSize;
                        sparsevalJac(counter2,1) = -ufxC(1,k)*ufxtK(xindex+2*N*(tindex-1)+(k-3),1);
                        counter2 = counter2+1;
                    end
                    % + dampMult*f
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                    sparseval(counter1,1) = dampMult(xindex,1);
                    counter1 = counter1+1;

                    sparseiJac(counter2,1) = m+bigSize;
                    sparsejJac(counter2,1) = xindex+2*N*(tindex-1)+N;
                    sparsevalJac(counter2,1) = dampMult(xindex,1);
                    counter2 = counter2+1;
                    % = RHS
                    %RHS(m+bigSize,1) = -uftC(1,1)*fx0(xindex,1);
                end
            
        end
    end
    
    RHS(m+bigSize,1) = RHS(m+bigSize,1)-...
        2/((3*timeSigma)^0.5*pi^0.25)*...
        (1-(tindex*tstep-exTime)^2/timeSigma^2)*...
        exp(-(tindex*tstep-exTime)^2/(2*timeSigma^2))*...
        (1/(spaceSigma*(2*pi)^0.5))*...
        exp(-0.5*((xgrid(xindex,1)-exLoc)/spaceSigma)^2);
    
end

knownRowsCounter = 0;
for m = 1:size(knownu,1)
    % known u
    for uCounter = 1:size(knownuInd,1)
        knownRowsCounter = knownRowsCounter+1;
        sparsei(counter1,1) = knownRowsCounter+2*bigSize;
        sparsej(counter1,1) = knownuInd(uCounter,1)+2*N*(m-1);
        sparseval(counter1,1) = knownuAlpha;
        counter1 = counter1+1;

        sparseiJac(counter2,1) = knownRowsCounter+2*bigSize;
        sparsejJac(counter2,1) = knownuInd(uCounter,1)+2*N*(m-1);
        sparsevalJac(counter2,1) = knownuAlpha;
        counter2 = counter2+1;

        RHS(knownRowsCounter+2*bigSize,1) = knownuAlpha*knownu(m,uCounter);
    end
    % known f
    for fCounter = 1:size(knownfInd,1)
        knownRowsCounter = knownRowsCounter+1;
        sparsei(counter1,1) = knownRowsCounter+2*bigSize;
        sparsej(counter1,1) = knownfInd(fCounter,1)+N+2*N*(m-1);
        sparseval(counter1,1) = knownfAlpha;
        counter1 = counter1+1;

        sparseiJac(counter2,1) = knownRowsCounter+2*bigSize;
        sparsejJac(counter2,1) = knownfInd(fCounter,1)+N+2*N*(m-1);
        sparsevalJac(counter2,1) = knownfAlpha;
        counter2 = counter2+1;

        RHS(knownRowsCounter+2*bigSize,1) = knownfAlpha*knownf(m,fCounter);
    end
end

augRowsCounter = 0;
if FOflag == 1
   if intKnownFlag == 0
       for m = 1:N
           if xgrid(m,1) > FOloc
               if m ~= N
               %if m ~= N/2 && m ~= N
                   % known FO mult - 1/p
                   sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsej(counter1,1) = 2*bigSize+m;
                   sparseval(counter1,1) = -FOalpha;
                   counter1 = counter1+1;

                   sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsej(counter1,1) = 2*bigSize+m+1;
                   sparseval(counter1,1) = FOalpha;
                   counter1 = counter1+1;

                   sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsejJac(counter2,1) = 2*bigSize+m;
                   sparsevalJac(counter2,1) = -FOalpha;
                   counter2 = counter2+1;

                   sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsejJac(counter2,1) = 2*bigSize+m+1;
                   sparsevalJac(counter2,1) = FOalpha;
                   counter2 = counter2+1;

                   augRowsCounter = augRowsCounter+1;

                   % known FO mult - pc2
                   sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsej(counter1,1) = 2*bigSize+m+N;
                   sparseval(counter1,1) = -FOalpha;
                   counter1 = counter1+1;

                   sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsej(counter1,1) = 2*bigSize+m+1+N;
                   sparseval(counter1,1) = FOalpha;
                   counter1 = counter1+1;

                   sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsejJac(counter2,1) = 2*bigSize+m+N;
                   sparsevalJac(counter2,1) = -FOalpha;
                   counter2 = counter2+1;

                   sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsejJac(counter2,1) = 2*bigSize+m+1+N;
                   sparsevalJac(counter2,1) = FOalpha;
                   counter2 = counter2+1;

                   augRowsCounter = augRowsCounter+1;
               end
           end    
       end
   else
       for m = 1:N
           FOexception = 0;
           if xgrid(m,1) > FOloc
               for k4 = 1:numInterfaces
                   if m == intIndices(k4,1)
                       FOexception = 1;
                   end
               end
               if FOexception ~= 1 && m ~= N
                   % known FO mult - 1/p
                   sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsej(counter1,1) = 2*bigSize+m;
                   sparseval(counter1,1) = -FOalpha;
                   counter1 = counter1+1;

                   sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsej(counter1,1) = 2*bigSize+m+1;
                   sparseval(counter1,1) = FOalpha;
                   counter1 = counter1+1;

                   sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsejJac(counter2,1) = 2*bigSize+m;
                   sparsevalJac(counter2,1) = -FOalpha;
                   counter2 = counter2+1;

                   sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsejJac(counter2,1) = 2*bigSize+m+1;
                   sparsevalJac(counter2,1) = FOalpha;
                   counter2 = counter2+1;

                   augRowsCounter = augRowsCounter+1;

                   % known FO mult - pc2
                   sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsej(counter1,1) = 2*bigSize+m+N;
                   sparseval(counter1,1) = -FOalpha;
                   counter1 = counter1+1;

                   sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsej(counter1,1) = 2*bigSize+m+1+N;
                   sparseval(counter1,1) = FOalpha;
                   counter1 = counter1+1;

                   sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsejJac(counter2,1) = 2*bigSize+m+N;
                   sparsevalJac(counter2,1) = -FOalpha;
                   counter2 = counter2+1;

                   sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+1;
                   sparsejJac(counter2,1) = 2*bigSize+m+1+N;
                   sparsevalJac(counter2,1) = FOalpha;
                   counter2 = counter2+1;

                   augRowsCounter = augRowsCounter+1;
               end
           end    
       end
   end
end

ZORowsCounter = 0;
if ZOflag == 1
   
   for m = 1:N
       tempFlag = 0;
       if abs(xgrid(m,1) - exLoc) < ZOradius
           tempFlag = tempFlag+1;
       end
       for counter4 = 1:ZOswitchu
           if abs(xgrid(knownuInd(counter4,1)) - xgrid(m,1)) < ZOradius
               tempFlag = tempFlag+1;
           end
       end
       for counter4 = 1:ZOswitchf
           if abs(xgrid(knownfInd(counter4,1)) - xgrid(m,1)) < ZOradius
               tempFlag = tempFlag+1;
           end
       end
       if tempFlag > 0  

           % known ZO mult - 1/p
           sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+ZORowsCounter+1;
           sparsej(counter1,1) = 2*bigSize+m;
           sparseval(counter1,1) = ZOalpha;
           counter1 = counter1+1;

           sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+ZORowsCounter+1;
           sparsejJac(counter2,1) = 2*bigSize+m;
           sparsevalJac(counter2,1) = ZOalpha;
           counter2 = counter2+1;

           ZORowsCounter = ZORowsCounter+1;

           % known FO mult - pc2
           sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+ZORowsCounter+1;
           sparsej(counter1,1) = 2*bigSize+m+N;
           sparseval(counter1,1) = ZOalpha;
           counter1 = counter1+1;

           sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+ZORowsCounter+1;
           sparsejJac(counter2,1) = 2*bigSize+m+N;
           sparsevalJac(counter2,1) = ZOalpha;
           counter2 = counter2+1;

           ZORowsCounter = ZORowsCounter+1; 
       end
         
   end
end

RhoRowsCounter = 0;
if RhoKnownflag == 1
   
   for m = 1:N  

           % known ZO mult - 1/p
           sparsei(counter1,1) = 2*bigSize+knownRowsCounter+augRowsCounter+RhoRowsCounter+ZORowsCounter+1;
           sparsej(counter1,1) = 2*bigSize+m;
           sparseval(counter1,1) = Rhoalpha;
           counter1 = counter1+1;

           sparseiJac(counter2,1) = 2*bigSize+knownRowsCounter+augRowsCounter+RhoRowsCounter+ZORowsCounter+1;
           sparsejJac(counter2,1) = 2*bigSize+m;
           sparsevalJac(counter2,1) = Rhoalpha;
           counter2 = counter2+1;

           RhoRowsCounter = RhoRowsCounter+1; 
         
   end
end

if FOflag == 0 && ZOflag == 0 && RhoKnownflag == 0
    bigOperator = sparse(sparsei,sparsej,sparseval,2*bigSize+knownRowsCounter,2*bigSize+2*N);
    bigJacobian = sparse(sparseiJac,sparsejJac,sparsevalJac,2*bigSize+knownRowsCounter,2*bigSize+2*N);
else
    bigOperator = sparse(sparsei,sparsej,sparseval,2*bigSize+knownRowsCounter+augRowsCounter+ZORowsCounter+RhoRowsCounter,2*bigSize+2*N);
    bigJacobian = sparse(sparseiJac,sparsejJac,sparsevalJac,2*bigSize+knownRowsCounter+augRowsCounter+ZORowsCounter+RhoRowsCounter,2*bigSize+2*N);
    RHS = [RHS; zeros(augRowsCounter,1); ones(ZORowsCounter,1)*ZOalpha; knownRhoVec*Rhoalpha];
end
tic;
F = -(bigOperator*ufxtK-RHS);
if CGflag == 0
    ufxtKP1 = bigJacobian\F+ufxtK;
else
    deltaPrev = zeros(2*bigSize+2*N,1); pPrev = zeros(2*bigSize+2*N,1);
    betaPrev = 0; sPrev = -F; rPrev = bigJacobian'*sPrev;
    
    for k = 1:CGiterations
        pNext = -rPrev+betaPrev*pPrev;
        alpha = (rPrev'*rPrev)/((bigJacobian*pNext)'*(bigJacobian*pNext));
        deltaNext = deltaPrev+alpha*pNext;
        sNext = sPrev+alpha*bigJacobian*pNext;
        rNext = bigJacobian'*sNext;
        betaNext = (rNext'*rNext)/(rPrev'*rPrev);
        
        deltaPrev = deltaNext; pPrev = pNext; betaPrev = betaNext;
        sPrev = sNext; rPrev = rNext;
    end
    
    ufxtKP1 = deltaNext+ufxtK;
end
toc;
%

