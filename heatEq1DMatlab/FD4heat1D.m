function [xgrid uxt timeVec] = FD4heat1D(N,tstep,endtime,centerChange,rampFactor)

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
