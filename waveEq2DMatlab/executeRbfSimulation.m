function datastore = ...
    executeRbfSimulation(bigLmatrix,bigDampMatrix,gamma,endtime,k,initdatavec,storesteps,...
    spaceSigma,xnodes,ynodes,N,explFlag,FDflagMat,lambdaMuVec,rhoVec,FDorder)

% This function carries out RK4 time stepping of an RBF-FD/traditional FD
% hybrid solution in our doubly-periodic unit square.

timesteps = ceil(endtime/k);  % total number of time steps through which
                              % we integrate data
                            
N1D = round(size(xnodes,1)^0.5);
h = 1/N1D;

lambdamu = reshape(lambdaMuVec,N1D,N1D);
rho = reshape(rhoVec,N1D,N1D);

datastore = zeros(size(initdatavec,1),size(storesteps,1));

if explFlag == 0
    initdatavec(1+N:2*N,1) = exp(-spaceSigma^2*(ynodes-0.75).^2);
    initdatavec(1+2*N:3*N,1) = 1/3^0.5*exp(-spaceSigma^2*(ynodes-0.75).^2);
    initdatavec(1+4*N:5*N,1) = 3^0.5*exp(-spaceSigma^2*(ynodes-0.75).^2);
else
    initdatavec(1+N:2*N,1) = (ynodes-0.75).*exp(-spaceSigma^2*((ynodes-0.75).^2+(xnodes-0.5).^2));
    initdatavec(1+0*N:1*N,1) = (xnodes-0.5).*exp(-spaceSigma^2*((ynodes-0.5).^2+(xnodes-0.5).^2));
end

uxt = reshape(initdatavec(1+0*N:1*N,1),N1D,N1D);
vxt = reshape(initdatavec(1+1*N:2*N,1),N1D,N1D);
s1xt = reshape(initdatavec(1+2*N:3*N,1),N1D,N1D);
s2xt = reshape(initdatavec(1+3*N:4*N,1),N1D,N1D);
s3xt = reshape(initdatavec(1+4*N:5*N,1),N1D,N1D);

datastore(:,1) = initdatavec;
datavecCurrent = initdatavec;

Dxmatrix = zeros(N1D,N1D);
Dymatrix = zeros(N1D,N1D);

preygrid = linspace(0,1,N1D+1);
ygrid = preygrid+h/2;

for rowCounter = 1:N1D
    defaultFDorder = 6;
    maxFDorder = 0;
    if ygrid(1,rowCounter) < 0.23
        maxFDorder = round((0.23-ygrid(1,rowCounter))/h)-1;
    end
    if ygrid(1,rowCounter) > 0.27 && ygrid(1,rowCounter) < 0.48
        oneFDorder = -(round((0.27-ygrid(1,rowCounter))/h)+1);
        twoFDorder = round((0.48-ygrid(1,rowCounter))/h)-1;
        maxFDorder = min(oneFDorder,twoFDorder);
    end
    if ygrid(1,rowCounter) > 0.52
        maxFDorder = -(round((0.52-ygrid(1,rowCounter))/h)+1);
    end
    tempFDorder = max(maxFDorder,defaultFDorder);
    locFDorder = min(tempFDorder,FDorder);
    
    if mod(locFDorder,2) == 1
        locFDorder = locFDorder-1;
    end
    
    preWeights = BengtsWeights(0,-locFDorder/2:locFDorder/2,1);
    preColIdx = (-locFDorder/2:locFDorder/2)+(rowCounter-1);
    colIdx = mod(preColIdx,N1D)+1;
    
    for idxCounter = 1:size(colIdx,2)
        Dymatrix(rowCounter,colIdx(1,idxCounter)) = ...
            preWeights(2,idxCounter)/h;
    end
    
    preWeights = BengtsWeights(0,-FDorder/2:FDorder/2,1);
    preColIdx = (-FDorder/2:FDorder/2)+(rowCounter-1);
    colIdx = mod(preColIdx,N1D)+1;
    for idxCounter = 1:size(colIdx,2)
        Dxmatrix(rowCounter,colIdx(1,idxCounter)) = ...
            preWeights(2,idxCounter)/h;
    end
end

Dxmatrix = sparse(Dxmatrix);

% we combine the block differential operator with a block hyperviscosity
% operator to make one unified operator with which we can time step.

tic;

% In the time stepping FOR loop below, we integrate via RK4.

for time_step = 2:timesteps
    
    time_step
    
    %%% RK4 Step 1 %%%
    d1 = bigLmatrix*datavecCurrent + gamma*(bigDampMatrix*datavecCurrent);
    %step2data = datavecCurrent+k/2*d1;
    
    ux = -uxt*Dxmatrix;
    uy = Dymatrix*uxt;
    
    %vx = (Dxmatrix*(vxt(:,:)'))';
    vx = -vxt*Dxmatrix;
    vy = Dymatrix*vxt;
    
    s1x = -s1xt*Dxmatrix;
    %s1y = Dxmatrix*s1xt(:,:,z-1);
    
    s2x = -s2xt*Dxmatrix;
    s2y = Dymatrix*s2xt;

    %s3x = (Dxmatrix*(s3xt(:,:,z-1)'))';
    s3y = Dymatrix*s3xt;

    d1u = FDflagMat.*(s1x+s2y)./rho;
    d1v = FDflagMat.*(s2x+s3y)./rho;
    d1s1 = FDflagMat.*((3*lambdamu).*ux+lambdamu.*vy);
    d1s2 = FDflagMat.*(lambdamu.*uy+lambdamu.*vx);
    d1s3 = FDflagMat.*(lambdamu.*ux+(3*lambdamu).*vy);

    step2u = uxt(:,:)+k/2*(d1u + reshape(d1(1+0*N:1*N,1),N1D,N1D));
    step2v = vxt(:,:)+k/2*(d1v + reshape(d1(1+1*N:2*N,1),N1D,N1D));
    step2s1 = s1xt(:,:)+k/2*(d1s1 + reshape(d1(1+2*N:3*N,1),N1D,N1D));
    step2s2 = s2xt(:,:)+k/2*(d1s2 + reshape(d1(1+3*N:4*N,1),N1D,N1D));
    step2s3 = s3xt(:,:)+k/2*(d1s3 + reshape(d1(1+4*N:5*N,1),N1D,N1D));
    
    step2data = [reshape(step2u,N1D*N1D,1); reshape(step2v,N1D*N1D,1); ...
        reshape(step2s1,N1D*N1D,1); reshape(step2s2,N1D*N1D,1); ...
        reshape(step2s3,N1D*N1D,1)];
    
    %%% RK4 Step 2 %%%
    d2 = bigLmatrix*step2data + gamma*(bigDampMatrix*step2data);
    %step3data = datavecCurrent+k/2*d2;
    
    ux = -step2u*Dxmatrix;
    uy = Dymatrix*step2u;
    
    vx = -step2v*Dxmatrix;
    vy = Dymatrix*step2v;
    
    s1x = -step2s1*Dxmatrix;
    %s1y = Dxmatrix*step2s1;
    
    s2x = -step2s2*Dxmatrix;
    s2y = Dymatrix*step2s2;

    %s3x = (Dxmatrix*(step2s3'))';
    s3y = Dymatrix*step2s3;

    d2u = FDflagMat.*(s1x+s2y)./rho;
    d2v = FDflagMat.*(s2x+s3y)./rho;
    d2s1 = FDflagMat.*((3*lambdamu).*ux+lambdamu.*vy);
    d2s2 = FDflagMat.*(lambdamu.*uy+lambdamu.*vx);
    d2s3 = FDflagMat.*(lambdamu.*ux+(3*lambdamu).*vy);
    
    step3u = uxt(:,:)+k/2*(d2u + reshape(d2(1+0*N:1*N,1),N1D,N1D));
    step3v = vxt(:,:)+k/2*(d2v + reshape(d2(1+1*N:2*N,1),N1D,N1D));
    step3s1 = s1xt(:,:)+k/2*(d2s1 + reshape(d2(1+2*N:3*N,1),N1D,N1D));
    step3s2 = s2xt(:,:)+k/2*(d2s2 + reshape(d2(1+3*N:4*N,1),N1D,N1D));
    step3s3 = s3xt(:,:)+k/2*(d2s3 + reshape(d2(1+4*N:5*N,1),N1D,N1D));
    
    step3data = [reshape(step3u,N1D*N1D,1); reshape(step3v,N1D*N1D,1); ...
        reshape(step3s1,N1D*N1D,1); reshape(step3s2,N1D*N1D,1); ...
        reshape(step3s3,N1D*N1D,1)];
    
    %%% RK4 Step 3 %%%
    d3 = bigLmatrix*step3data + gamma*(bigDampMatrix*step3data);
    %step4data = datavecCurrent+k*d3;
    
    ux = -step3u*Dxmatrix;
    uy = Dymatrix*step3u;
    
    vx = -step3v*Dxmatrix;
    vy = Dymatrix*step3v;
    
    s1x = -step3s1*Dxmatrix;
    %s1y = Dxmatrix*step3s1;
    
    s2x = -step3s2*Dxmatrix;
    s2y = Dymatrix*step3s2;

    %s3x = (Dxmatrix*(step3s3'))';
    s3y = Dymatrix*step3s3;

    d3u = FDflagMat.*(s1x+s2y)./rho;
    d3v = FDflagMat.*(s2x+s3y)./rho;
    d3s1 = FDflagMat.*((3*lambdamu).*ux+lambdamu.*vy);
    d3s2 = FDflagMat.*(lambdamu.*uy+lambdamu.*vx);
    d3s3 = FDflagMat.*(lambdamu.*ux+(3*lambdamu).*vy);

    step4u = uxt(:,:)+k*(d3u + reshape(d3(1+0*N:1*N,1),N1D,N1D));
    step4v = vxt(:,:)+k*(d3v + reshape(d3(1+1*N:2*N,1),N1D,N1D));
    step4s1 = s1xt(:,:)+k*(d3s1 + reshape(d3(1+2*N:3*N,1),N1D,N1D));
    step4s2 = s2xt(:,:)+k*(d3s2 + reshape(d3(1+3*N:4*N,1),N1D,N1D));
    step4s3 = s3xt(:,:)+k*(d3s3 + reshape(d3(1+4*N:5*N,1),N1D,N1D));
    
    step4data = [reshape(step4u,N1D*N1D,1); reshape(step4v,N1D*N1D,1); ...
        reshape(step4s1,N1D*N1D,1); reshape(step4s2,N1D*N1D,1); ...
        reshape(step4s3,N1D*N1D,1)];
    
    %%% RK4 Step 4 %%%
    d4 = bigLmatrix*step4data + gamma*(bigDampMatrix*step4data);
    
    ux = -step4u*Dxmatrix;
    uy = Dymatrix*step4u;
    
    vx = -step4v*Dxmatrix;
    vy = Dymatrix*step4v;
    
    s1x = -step4s1*Dxmatrix;
    %s1y = Dxmatrix*step4s1;
    
    s2x = -step4s2*Dxmatrix;
    s2y = Dymatrix*step4s2;

    %s3x = (Dxmatrix*(step4s3'))';
    s3y = Dymatrix*step4s3;

    d4u = FDflagMat.*(s1x+s2y)./rho;
    d4v = FDflagMat.*(s2x+s3y)./rho;
    d4s1 = ((3*lambdamu).*ux+lambdamu.*vy).*FDflagMat;
    d4s2 = (lambdamu.*uy+lambdamu.*vx).*FDflagMat;
    d4s3 = (lambdamu.*ux+(3*lambdamu).*vy).*FDflagMat;
    
    %%% RK4 final step %%%
    
    d1Vec = [reshape(d1u,N1D*N1D,1); reshape(d1v,N1D*N1D,1); ...
        reshape(d1s1,N1D*N1D,1); reshape(d1s2,N1D*N1D,1); ...
        reshape(d1s3,N1D*N1D,1)];
    
    d2Vec = [reshape(d2u,N1D*N1D,1); reshape(d2v,N1D*N1D,1); ...
        reshape(d2s1,N1D*N1D,1); reshape(d2s2,N1D*N1D,1); ...
        reshape(d2s3,N1D*N1D,1)];
    
    d3Vec = [reshape(d3u,N1D*N1D,1); reshape(d3v,N1D*N1D,1); ...
        reshape(d3s1,N1D*N1D,1); reshape(d3s2,N1D*N1D,1); ...
        reshape(d3s3,N1D*N1D,1)];
    
    d4Vec = [reshape(d4u,N1D*N1D,1); reshape(d4v,N1D*N1D,1); ...
        reshape(d4s1,N1D*N1D,1); reshape(d4s2,N1D*N1D,1); ...
        reshape(d4s3,N1D*N1D,1)];
    
    datavecCurrent = datavecCurrent+k/6*(d1+2*d2+2*d3+d4 + ...
        d1Vec + 2*d2Vec + 2*d3Vec + d4Vec);
    
    uxt = reshape(datavecCurrent(1+0*N:1*N,1),N1D,N1D);
    vxt = reshape(datavecCurrent(1+1*N:2*N,1),N1D,N1D);
    s1xt = reshape(datavecCurrent(1+2*N:3*N,1),N1D,N1D);
    s2xt = reshape(datavecCurrent(1+3*N:4*N,1),N1D,N1D);
    s3xt = reshape(datavecCurrent(1+4*N:5*N,1),N1D,N1D);
    
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

    