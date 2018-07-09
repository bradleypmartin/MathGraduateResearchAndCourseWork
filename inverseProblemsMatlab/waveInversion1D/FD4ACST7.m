function [xgrid ufxt bigOperator bigSize] = ...
    FD4ACST7(N,tstep,endtime,dampConst,naiveFlag,...
    spaceSigma,timeSigma,exLoc,exTime,intStack)

% Edit 4 - 140404 - changing code to 1-D simplification of 2-D EWE;
% Edit 5 - 140405 - enabling toggle between naive and AC solution
% Edit 6 - 140410 - addition of Ricker wavelet forcing condition
% Edit 7 - 140418 - adds "intstack" as an input - allows placement of an
% arbitrary number of interfaces in 1-D unit interval.  Instack is a 3 x k
% matrix, where the rows indicate individual interfaces.  The first column
% holds the x-location of the interface.  The second column holds value of
% c to the right of that interface (until the next interface) and the third
% entry holds the value of rho to the right of that interface (until the
% next)
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

    c2 = intStack(interfaceCounter,2);
    if interfaceCounter == 1
        c1 = 1;
    else
        c1 = intStack(interfaceCounter-1,2);
    end

    p2 = intStack(interfaceCounter,3);
    if interfaceCounter == 1
        p1 = 1;
    else
        p1 = intStack(interfaceCounter-1,3);
    end

    superStackMid = zeros(5,6,4,2); % u = 1; f = 2

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

cmult = ones(size(xgrid));
pmult = ones(size(xgrid));
dampMult = zeros(size(xgrid));

intIndices = zeros(numInterfaces,1); % intIndices will store the index
                                     % of the left-most data node that
                                     % hasn't crossed the given int. (row)
tempCounter = 1;
                                     
for m = 1:N
    for intCounter = 1:numInterfaces
        if xgrid(m,1) > intStack(intCounter,1)
            cmult(m,1) = intStack(intCounter,2)^2;
            pmult(m,1) = intStack(intCounter,3);
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

ux0 = zeros(size(xgrid));
fx0 = zeros(size(xgrid));

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

weightsT1 = weights(-1,-2:2,1)/tstep;
weightsTC = weights(0,-2:2,1)/tstep;
weightsTEm1 = weights(1,-2:2,1)/tstep;
weightsTE = weights(2,-2:2,1)/tstep;

uftT1 = weightsT1(2,1:5);
uftC = weightsTC(2,1:5);
uftEm1 = weightsTEm1(2,1:5);
uftE = weightsTE(2,1:5);

bigSize = N*timesteps;

RHS = zeros(2*bigSize,1);
sparsei = [];
sparsej = [];
sparseval = [];
counter1 = 1;

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
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparseval(counter1,1) = -fx2L(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);
                
                % ft-c2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparseval(counter1,1) = -ux2L(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);
                
                % ft-c2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparseval(counter1,1) = -fx2R(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparseval(counter1,1) = -ux2R(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxR(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);
                
                % ft-c2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxR(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:5
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxTemp(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);

                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                end
                for k = 1:5
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxTemp(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:5
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -ufxC(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftT1(1,1)*ux0(xindex,1);

                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-2);
                    sparseval(counter1,1) = uftT1(1,k+1);
                    counter1 = counter1+1;
                end
                for k = 1:5
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -ufxC(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparseval(counter1,1) = -fx2L(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparseval(counter1,1) = -ux2L(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparseval(counter1,1) = -uxL(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparseval(counter1,1) = -fx2R(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparseval(counter1,1) = -ux2R(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxR(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxR(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:5
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxTemp(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                end
                for k = 1:5
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxTemp(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:5
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -ufxC(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:4                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-3);
                    sparseval(counter1,1) = uftC(1,k+1);
                    counter1 = counter1+1;
                end
                for k = 1:5
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -ufxC(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparseval(counter1,1) = -fx2L(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparseval(counter1,1) = -ux2L(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparseval(counter1,1) = -uxL(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparseval(counter1,1) = -fx2R(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparseval(counter1,1) = -ux2R(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxR(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxR(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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

                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxTemp(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;

                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxTemp(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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

                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -ufxC(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-4);
                    sparseval(counter1,1) = uftC(1,k);
                    counter1 = counter1+1;

                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -ufxC(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparseval(counter1,1) = -fx2L(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparseval(counter1,1) = -ux2L(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparseval(counter1,1) = -uxL(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparseval(counter1,1) = -fx2R(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparseval(counter1,1) = -ux2R(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxR(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxR(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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

                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxTemp(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;

                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxTemp(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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

                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -ufxC(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-5);
                    sparseval(counter1,1) = uftEm1(1,k);
                    counter1 = counter1+1;

                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -ufxC(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-1);
                    sparseval(counter1,1) = -fx2L(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-1);
                    sparseval(counter1,1) = -ux2L(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-2);
                    sparseval(counter1,1) = -fxL(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-2);
                    sparseval(counter1,1) = -uxL(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-4);
                    sparseval(counter1,1) = -fx2R(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-4);
                    sparseval(counter1,1) = -ux2R(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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
                end
                for k = 1:4
                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxR(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);
                
                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;
                end
                for k = 1:4
                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxR(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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

                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -fxTemp(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;

                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -uxTemp(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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

                    % -fx/p
                    sparsei(counter1,1) = m;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+N+(k-3);
                    sparseval(counter1,1) = -ufxC(1,k)/pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*u
                sparsei(counter1,1) = m;
                sparsej(counter1,1) = xindex+2*N*(tindex-1);
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
                % = RHS
                %RHS(m,1) = -uftC(1,1)*ux0(xindex,1);

                % ft-pc2ux+dampMult*f = RHS        
                for k = 1:5                   
                    % ft
                    sparsei(counter1,1) = m+bigSize;      
                    sparsej(counter1,1) = xindex+N+2*N*(tindex+k-6);
                    sparseval(counter1,1) = uftE(1,k);
                    counter1 = counter1+1;

                    % -pc2ux
                    sparsei(counter1,1) = m+bigSize;
                    sparsej(counter1,1) = xindex+2*N*(tindex-1)+(k-3);
                    sparseval(counter1,1) = -ufxC(1,k)*cmult(xindex,1)*pmult(xindex,1);
                    counter1 = counter1+1;
                end
                % + dampMult*f
                sparsei(counter1,1) = m+bigSize;
                sparsej(counter1,1) = xindex+2*N*(tindex-1)+N;
                sparseval(counter1,1) = dampMult(xindex,1);
                counter1 = counter1+1;
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

bigOperator = sparse(sparsei,sparsej,sparseval,2*bigSize,2*bigSize);
tic;
ufxt = bigOperator\RHS;
toc;
%

