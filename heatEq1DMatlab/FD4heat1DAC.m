function [xgrid uxt timeVec] = FD4heat1DAC(N,tstep,endtime,centerChange,rampFactor)

h = 2/N;

xgrid = linspace(-1,1,N+1)';

k1 = 1/9;
k2 = 1;

kmult = k1*ones(size(xgrid));

for m = 1:N+1
    if xgrid(m,1)>0
        kmult(m,1) = k2;
    end
end

superStackMid = zeros(5,6,4);

for k = 1:4
    for i = 1:5
        for j = 1:5
            constMult = 1;
            if j <= 5-k
                if i > 1
                    constMult = (k2/k1);
                end
                if i > 3
                    constMult = (k2/k1)^2;
                end
            end
            superStackMid(i,j,k)= constMult*(-1.5-4+j+k)^(i-1);
        end
        
        % rhs alteration
        constMult = 1;
        if k < 3
            if i > 1
                constMult = (k2/k1);
            end
            if i > 3
                constMult = (k2/k1)^2;
            end
        end
        
        superStackMid(i,6,k) = constMult*(i-1)*(i-2)*(-2.5+k)^(i-3);
        
    end
end

weightsMid = zeros(5,4);

for k = 1:4
    weightsMid(:,k) = superStackMid(:,1:5,k)\superStackMid(:,6,k)/h^2;
end

timesteps = ceil(endtime/tstep);

uxt = zeros(N+1,timesteps);
uxtTemp = zeros(N+1,1);
timeVec = zeros(timesteps,1);

preWeightsCentered = weights(0,-2:2,2)/h^2;
preWeightsLeft = weights(-1,-2:2,2)/h^2;
preWeightsRight = weights(1,-2:2,2)/h^2;

uxxmatrix = diag(preWeightsCentered(3,3)*ones(N+1,1),0)+...
    diag(preWeightsCentered(3,2)*ones(N,1),-1)+...
    diag(preWeightsCentered(3,4)*ones(N,1),1)+...
    diag(preWeightsCentered(3,1)*ones(N-1,1),-2)+...
    diag(preWeightsCentered(3,5)*ones(N-1,1),2);

uxxmatrix(2,:) = zeros(1,N+1);
uxxmatrix(N,:) = zeros(1,N+1);

uxxmatrix(2,1:5) = preWeightsLeft(3,:);
uxxmatrix(N,N-3:N+1) = preWeightsRight(3,:);

uxxmatrix(N/2,:) = zeros(1,N+1);
uxxmatrix(N/2,N/2-2:N/2+2) = weightsMid(:,1)';

uxxmatrix(N/2+1,:) = zeros(1,N+1);
uxxmatrix(N/2+1,N/2-1:N/2+3) = weightsMid(:,2)';

uxxmatrix(N/2+2,:) = zeros(1,N+1);
uxxmatrix(N/2+2,N/2:N/2+4) = weightsMid(:,3)';

uxxmatrix(N/2+3,:) = zeros(1,N+1);
uxxmatrix(N/2+3,N/2+1:N/2+5) = weightsMid(:,4)';

uxxmatrix = sparse(uxxmatrix);

%RK4 timestep

for n = 2:timesteps

    prevTime = (n-2)*tstep;
    timeVec(n,1) = (n-1)*tstep;
    
    uxtTemp(1,1) = 1/(1+exp(-rampFactor*(prevTime-centerChange)));
    uxtTemp(N+1,1) = 0;
    
    d1u = uxxmatrix*uxtTemp.*kmult;
    step2u = uxtTemp+tstep/2*d1u;
    
    step2u(1,1) = 1/(1+exp(-rampFactor*(prevTime+tstep/2-centerChange)));
    step2u(N+1,1) = 0;
    
    d2u = uxxmatrix*step2u.*kmult;
    step3u = uxtTemp+tstep/2*d2u;
    
    step3u(1,1) = 1/(1+exp(-rampFactor*(prevTime+tstep/2-centerChange)));
    step3u(N+1,1) = 0;
    
    d3u = uxxmatrix*step3u.*kmult;
    step4u = uxtTemp+tstep*d3u;
    
    step4u(1,1) = 1/(1+exp(-rampFactor*(prevTime+tstep-centerChange)));
    step4u(N+1,1) = 0;
    
    d4u = uxxmatrix*step4u.*kmult;
    
    uxtTemp = uxtTemp+tstep/6*(d1u+2*d2u+2*d3u+d4u);
    
    uxtTemp(1,1) = 1/(1+exp(-rampFactor*(prevTime+tstep-centerChange)));
    uxtTemp(N+1,1) = 0;
    
    uxt(:,n) = uxtTemp;

end
