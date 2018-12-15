function [bigLmatrix, bigDampMatrix, initdatavec, xnodes, ynodes, lambdaMuVec, rhoVec, tripleVec, FDflagMat] = ...
EWE2DRbfPrep(N,GAshp,stencilsize,polydegree,...
    xlb,xub,ylb,yub,lambdaMu1,lambdaMu2,rho1,rho2,xnodes,ynodes,nearestneighbors,...
    totalmovements,plotflag,numIntNodes,intWrapFlag,...
    polydegreeInt,stencilsizeInt,naiveFlag,dampLk,...
    const1,const2,const3,yLI,yUI,queryFactor,allRBFflag)

% This function takes care of operator creation for 2D EWE solutions in the
% doubly-periodic unit square.

lambdaMuVec = zeros(N,1);
rhoVec = zeros(N,1);

% Below: creating a new node set of (xnodes,ynodes) (if those are needed)

if size(xnodes,1) == 1
    
    N1D = round(N^0.5);
    
    prexgrid = linspace(xlb,xub,N1D+1);
    xgrid = prexgrid(1,1:N1D)+(prexgrid(1,2)-prexgrid(1,1))/2;
    
    [X Y] = meshgrid(xgrid,xgrid);
    
    hCart = (xub-xlb)/N1D;
    
    moveFlag = ((abs(Y-curvedinterface1(X,yLI)) < (8*hCart)) + (abs(Y-curvedinterface2(X,yUI)) < (8*hCart))) > 0;
    antiFDflagMat = ((abs(Y-curvedinterface1(X,yLI)) < (12*hCart)) + (abs(Y-curvedinterface2(X,yUI)) < (12*hCart))) > 0;
    FDflagMat = 1-antiFDflagMat;
    RBFflagVec = reshape(antiFDflagMat,N,1);
    
    xnodesInt = zeros(12*numIntNodes,1);
    ynodesInt = zeros(12*numIntNodes,1);
    
    if allRBFflag > 0
        moveFlag = ones(N1D,N1D);
        RBFflagVec = ones(N,1);
    end
    
    if intWrapFlag == 1
       
       hInt = 1/numIntNodes;
       prexgrid = linspace(0,1,numIntNodes+1);
       xgrid = (prexgrid(1,1:numIntNodes))';
       [y1 yprime1 theta1] = curvedinterface1(xgrid,yLI);
       [y2 yprime2 theta2] = curvedinterface2(xgrid,yUI);
       xnodesInt(1:numIntNodes,1) = xgrid + sin(theta1)*hInt*0.5;
       ynodesInt(1:numIntNodes,1) = (curvedinterface1(xgrid,yLI)-0.5*cos(theta1)*hInt);
       xnodesInt(1+numIntNodes:2*numIntNodes,1) = xgrid - sin(theta1)*hInt*0.5;
       ynodesInt(1+numIntNodes:2*numIntNodes,1) = (curvedinterface1(xgrid,yLI)+0.5*hInt*cos(theta1));
       xnodesInt(1+2*numIntNodes:3*numIntNodes,1) = xgrid + sin(theta2)*hInt*0.5;
       ynodesInt(1+2*numIntNodes:3*numIntNodes,1) = (curvedinterface2(xgrid,yUI)-0.5*hInt*cos(theta2));
       xnodesInt(1+3*numIntNodes:4*numIntNodes,1) = xgrid - sin(theta2)*hInt*0.5;
       ynodesInt(1+3*numIntNodes:4*numIntNodes,1) = (curvedinterface2(xgrid,yUI)+0.5*hInt*cos(theta2));
       
       xnodesInt(1+4*numIntNodes:5*numIntNodes,1) = xgrid + sin(theta1)*hInt*(0.5+3^0.5/2)+hInt*0.5;
       ynodesInt(1+4*numIntNodes:5*numIntNodes,1) = (curvedinterface1(xgrid,yLI)-(0.5+3^0.5/2)*cos(theta1)*hInt);
       xnodesInt(1+5*numIntNodes:6*numIntNodes,1) = xgrid - sin(theta1)*hInt*(0.5+3^0.5/2)+hInt*0.5;
       ynodesInt(1+5*numIntNodes:6*numIntNodes,1) = (curvedinterface1(xgrid,yLI)+(0.5+3^0.5/2)*hInt*cos(theta1));
       xnodesInt(1+6*numIntNodes:7*numIntNodes,1) = xgrid + sin(theta2)*hInt*(0.5+3^0.5/2)+hInt*0.5;
       ynodesInt(1+6*numIntNodes:7*numIntNodes,1) = (curvedinterface2(xgrid,yUI)-(0.5+3^0.5/2)*hInt*cos(theta2));
       xnodesInt(1+7*numIntNodes:8*numIntNodes,1) = xgrid - sin(theta2)*hInt*(0.5+3^0.5/2)+hInt*0.5;
       ynodesInt(1+7*numIntNodes:8*numIntNodes,1) = (curvedinterface2(xgrid,yUI)+(0.5+3^0.5/2)*hInt*cos(theta2));
       
       xnodesInt(1+8*numIntNodes:9*numIntNodes,1) = xgrid + sin(theta1)*hInt*(0.5+3^0.5);
       ynodesInt(1+8*numIntNodes:9*numIntNodes,1) = (curvedinterface1(xgrid,yLI)-(0.5+3^0.5)*cos(theta1)*hInt);
       xnodesInt(1+9*numIntNodes:10*numIntNodes,1) = xgrid - sin(theta1)*hInt*(0.5+3^0.5);
       ynodesInt(1+9*numIntNodes:10*numIntNodes,1) = (curvedinterface1(xgrid,yLI)+(0.5+3^0.5)*hInt*cos(theta1));
       xnodesInt(1+10*numIntNodes:11*numIntNodes,1) = xgrid + sin(theta2)*hInt*(0.5+3^0.5);
       ynodesInt(1+10*numIntNodes:11*numIntNodes,1) = (curvedinterface2(xgrid,yUI)-(0.5+3^0.5)*hInt*cos(theta2));
       xnodesInt(1+11*numIntNodes:12*numIntNodes,1) = xgrid - sin(theta2)*hInt*(0.5+3^0.5);
       ynodesInt(1+11*numIntNodes:12*numIntNodes,1) = (curvedinterface2(xgrid,yUI)+(0.5+3^0.5)*hInt*cos(theta2));
       %{
       xnodes(1+12*numIntNodes:13*numIntNodes,1) = xgrid + sin(theta1)*hInt*(0.5+1.5*3^0.5)+hInt*0.5;
       ynodes(1+12*numIntNodes:13*numIntNodes,1) = (curvedinterface1(xgrid)-(0.5+1.5*3^0.5)*cos(theta1)*hInt);
       xnodes(1+13*numIntNodes:14*numIntNodes,1) = xgrid - sin(theta1)*hInt*(0.5+1.5*3^0.5)+hInt*0.5;
       ynodes(1+13*numIntNodes:14*numIntNodes,1) = (curvedinterface1(xgrid)+(0.5+1.5*3^0.5)*hInt*cos(theta1));
       xnodes(1+14*numIntNodes:15*numIntNodes,1) = xgrid + sin(theta2)*hInt*(0.5+1.5*3^0.5)+hInt*0.5;
       ynodes(1+14*numIntNodes:15*numIntNodes,1) = (curvedinterface2(xgrid)-(0.5+1.5*3^0.5)*hInt*cos(theta2));
       xnodes(1+15*numIntNodes:16*numIntNodes,1) = xgrid - sin(theta2)*hInt*(0.5+1.5*3^0.5)+hInt*0.5;
       ynodes(1+15*numIntNodes:16*numIntNodes,1) = (curvedinterface2(xgrid)+(0.5+1.5*3^0.5)*hInt*cos(theta2));
       %}
       
       intWrapCounter = 1;
       
       
       scatterVec = randperm(N)';
       for k = 1:N
           if intWrapCounter <= (12*numIntNodes)
               tempMatIdx = scatterVec(k,1);
               m = mod(tempMatIdx-1,N1D)+1;
               n = floor((tempMatIdx-1)/N1D)+1;
               if moveFlag(m,n) == 1
                   moveFlag(m,n) = 0;
                   X(m,n) = xnodesInt(intWrapCounter,1);
                   Y(m,n) = ynodesInt(intWrapCounter,1);
                   intWrapCounter = intWrapCounter + 1;
               end
           end
       end
       
    end
    
    xnodes = reshape(X,N,1);
    ynodes = reshape(Y,N,1);
    moveFlagVec = reshape(moveFlag,N,1);
    
    xnodes = xnodes + ((rand(N,1)-0.5)*0.05*hCart).*moveFlagVec;
    ynodes = ynodes + ((rand(N,1)-0.5)*0.05*hCart).*moveFlagVec;
    
    tic;

    for movement_iteration = 1:totalmovements
        movement_iteration
        [xnodes ynodes] = mos2dsqperiodic7(xnodes,ynodes,...
        (0.05*(2500/N)^0.5)/movement_iteration,nearestneighbors,xlb,xub,ylb,yub,...
        numIntNodes,intWrapFlag,yLI,yUI,moveFlagVec,allRBFflag);
    
        % here, we may enable or disable plotting the nodes while the node set
        % is being created.
    
        if plotflag == 1
            figure(1)
    
            plot(xnodes,ynodes,'.k')
            axis([xlb xub ylb yub])
            hold on;
            xgridPlot = linspace(0,1,100);
            ygridPlot1 = curvedinterface1(xgridPlot,yLI);
            ygridPlot2 = curvedinterface2(xgridPlot,yUI);
            plot(xgridPlot,ygridPlot1,'--k',xgridPlot,ygridPlot2,'--k');
            hold off;
        end
    end
    disp('For node set creation:')
    toc;
    disp(' ')
end

% Below: assigning the values of lambda and mu at each node

seqVec1 = zeros(N,1);
seqVec2 = zeros(N,1);

numPolysIntPlus = (polydegreeInt + 2)*(polydegreeInt+3)/2;
expStack = zeros(numPolysIntPlus,N);

xCoeff = [0];
yCoeff = [0];

for k = 1:(polydegreeInt+1)
    xCoeff = [xCoeff; (k:-1:0)'];
    yCoeff = [yCoeff; (0:k)'];
end

for k = 1:N
    if ynodes(k,1) <= curvedinterface2(xnodes(k,1),yUI) && ynodes(k,1) >= curvedinterface1(xnodes(k,1),yLI) %curvedinterface(xnodes(k,1))
        lambdaMuVec(k,1) = lambdaMu2+const1*sin(const2*pi*xnodes(k,1))*sin(const3*pi*ynodes(k,1));
        seqVec2(k,1) = 1;
        %expStack(1,k) = lambdaMu2+const1*sin(const2*pi*xnodes(k,1))*sin(const3*pi*ynodes(k,1));
        
    else
        lambdaMuVec(k,1) = lambdaMu1;
        seqVec1(k,1) = 1;
        %expStack(1,k) = lambdaMu1;
    end
    
    if (abs(ynodes(k,1)-curvedinterface1(xnodes(k,1),yLI)) <= (queryFactor*hInt)) || ...
            (abs(ynodes(k,1)-curvedinterface2(xnodes(k,1),yUI)) <= (queryFactor*hInt))
        
        if ynodes(k,1) < (curvedinterface1(xnodes(k,1),yLI)+curvedinterface2(xnodes(k,1),yUI))/2
            [xClosest theta] = pointFinder1(xnodes(k,1),ynodes(k,1),yLI);
            yClosest = curvedinterface1(xClosest,yLI);
        else
            [xClosest theta] = pointFinder2(xnodes(k,1),ynodes(k,1),yUI);
            yClosest = curvedinterface2(xClosest,yUI);
        end
        
        xRedef = [cos(theta); -sin(theta)];
        yRedef = [sin(theta); cos(theta)];
        
        multByXprime = zeros(numPolysIntPlus);
        multByYprime = zeros(numPolysIntPlus);
        
        rowCounterX = 2;
        skipCounterX = 1;
        skipLimitX = 1;
        
        rowCounterY = 3;
        skipCounterY = 1;
        skipLimitY = 1;
        
        for colCounter = 1:numPolysIntPlus
            if rowCounterX <= numPolysIntPlus
                multByXprime(rowCounterX,colCounter) = 1;
                rowCounterX = rowCounterX+1;
                if skipCounterX == skipLimitX
                    rowCounterX = rowCounterX+1;
                    skipCounterX = 1;
                    skipLimitX = skipLimitX+1;
                else
                    skipCounterX = skipCounterX+1;
                end
            end
            if rowCounterY <= numPolysIntPlus
                multByYprime(rowCounterY,colCounter) = 1;
                rowCounterY = rowCounterY+1;
                if skipCounterY == skipLimitY
                    rowCounterY = rowCounterY+1;
                    skipCounterY = 1;
                    skipLimitY = skipLimitY+1;
                else
                    skipCounterY = skipCounterY+1;
                end
            end
        end
        
        multByX = xRedef(1,1)*multByXprime + xRedef(2,1)*multByYprime;
        multByY = yRedef(1,1)*multByXprime + yRedef(2,1)*multByYprime;
        
        expTransMatrix = zeros(numPolysIntPlus);
        
        for z = 1:numPolysIntPlus
            tempMatrix = (multByX^xCoeff(z,1))*(multByY^yCoeff(z,1));
            expTransMatrix(:,z) = tempMatrix(:,1);
        end
        
        expStack(1,k) = lambdaMu2+const1*sin(const2*pi*xClosest)*sin(const3*pi*yClosest);
        for expCounter = 2:numPolysIntPlus
            xSign = 1;
            ySign = 1;
            xCoeffTemp = xCoeff(expCounter,1);
            yCoeffTemp = yCoeff(expCounter,1);
            xConst = (pi*const2)^xCoeffTemp;
            yConst = (pi*const3)^yCoeffTemp;
            xVar = sin(const2*pi*xClosest);
            yVar = sin(const3*pi*yClosest);
            if (mod(xCoeffTemp,4) > 1)
                xSign = -1;
            end
            if (mod(yCoeffTemp,4) > 1)
                ySign = -1;
            end
            if mod(xCoeffTemp,2) == 1
                xVar = cos(const2*pi*xClosest);
            end
            if mod(yCoeffTemp,2) == 1
                yVar = cos(const3*pi*yClosest);
            end
            
            expStack(expCounter,k) = const1*xSign*ySign*xVar*yVar*xConst*yConst/(factorial(xCoeffTemp)*factorial(yCoeffTemp));
            
        end
        
        expStack(:,k) = expTransMatrix*expStack(:,k);
        
    else
        
        expStack(1,k) = lambdaMuVec(k,1);
        
    end
    
end

for k = 1:N
    if ynodes(k,1) <= curvedinterface2(xnodes(k,1),yUI) && ynodes(k,1) >= curvedinterface1(xnodes(k,1),yLI) %curvedinterface(xnodes(k,1))
        lambdaMuVec(k,1) = lambdaMu2+const1*sin(const2*pi*xnodes(k,1))*sin(const3*pi*ynodes(k,1));
        rhoVec(k,1) = rho2;
        seqVec2(k,1) = 1;
    else
        lambdaMuVec(k,1) = lambdaMu1;
        rhoVec(k,1) = rho1;
        seqVec1(k,1) = 1;
    end
end

% Next, we create sparse, diagonal lambda and mu matrices to aid us in
% creating the large differential operator that we will call to convect
% waves.

% IC creation

p = zeros(5*N,1);  

initdatavec = p;

% Below: setting up dx/dy and hyperviscosity matrix operators

tic;
%
[bigLmatrix tripleVec] = ...
    createRBFLoperator1(xnodes,ynodes,xlb,xub,ylb,yub,GAshp,...
    stencilsize,polydegree,...
    lambdaMu1,rho1,lambdaMu2,rho2,polydegreeInt,...
    stencilsizeInt,naiveFlag,1,yLI,yUI,queryFactor,0,RBFflagVec,...
    lambdaMuVec,expStack);

bigDampMatrix = ...
    createRBFLoperator1(xnodes,ynodes,xlb,xub,ylb,yub,GAshp,...
    stencilsize,polydegree,...
    lambdaMu1,rho1,lambdaMu2,rho2,polydegreeInt,...
    stencilsizeInt,naiveFlag,dampLk,yLI,yUI,queryFactor,1,RBFflagVec,...
    lambdaMuVec,expStack);

disp('For matrix operator creation:')
toc;
disp(' ')

end


function [y yprime theta] = curvedinterface1(x,yLI)

% function CURVEDINTERFACE1 creates a simple, slightly curved C-inf interface that
% is periodic on our square.

y = 0.02*sin(2*pi*x)+yLI*ones(size(x));
yprime = 0.04*pi*cos(2*pi*x);
theta = atan(yprime);

end

function [y yprime theta] = curvedinterface2(x,yUI)

% function CURVEDINTERFACE2 creates a simple, slightly curved C-inf interface that
% is periodic on our square.

y = 0.02*sin(2*pi*x)+yUI*ones(size(x));
yprime = 0.04*pi*cos(2*pi*x);
theta = atan(yprime);

end


function [xnew ynew] = mos2dsqperiodic7(x,y,delta,nnn,xlb,xub,ylb,yub,...
    numIntNodes,intWrapFlag,yLI,yUI,moveFlagVec,allRBFflag)

% function MOS2DSQPERIODIC7 moves RBF-FD nodes through one step of
% electrostatic repulsion (by distance DELTA).

% ver 7: one interface at y = 0.5; another at 0.25.

% IN: 

% x: x-coordinates of RBF-FD nodes
% y: y-coordinates "             "
% delta: distance to move nodes along force vector determined by electrostatic
% repulsion
% nnn: number of nearest neighbors to use in calculating electrostatic
% force vector
% xlb, etc. : lower and upper bounds of periodic domain in x and y
% numIntNodes: number of interface nodes to throw down
% intWrapFlag: if 0: don't wrap int's; if 1: do.

% OUT:

% xnew: new x-coordinates for RBF-FD nodes
% ynew: new y-coordinates "              "

xcorr = zeros(size(x));  % xcorr and ycorr will hold the corrections
ycorr = zeros(size(y));  % (displacements) to move nodes by one step.

nnn = nnn+1;  % when performing knnsearches, we'll throw out the returned
              % identity indices (in the first slot of index matrices)

% The first (small) knn search is for finding out how much of the domain we
% should "tile" to ensure correct identification of periodic nearest
% neighbors near the boundaries.
              
[idx dist] = knnsearch([x y],[x(1:20,1) y(1:20,1)],'k',nnn);

avgdist = sum(dist(:,size(dist,2)))/size(dist,1);
safedist = 2*avgdist;

% The first FOR loop (below) finds out how many nodes we have to "tile."

counter = 0;

for m = 1:size(x,1)
    
    if xub-x(m,1)<=safedist
        counter = counter+1;
    end
    
    if x(m,1)-xlb<=safedist
        counter = counter+1;
    end
    
    if yub-y(m,1)<=safedist
        counter = counter+1;
    end
    
    if y(m,1)-ylb<=safedist
        counter = counter+1;
    end
    
    if xub-x(m,1)<=safedist
        if yub-y(m,1)<=safedist
        counter = counter+1;
        end
    end
    
    if xub-x(m,1)<=safedist
        if y(m,1)-ylb<=safedist
        counter = counter+1;
        end
    end
    
    if x(m,1)-xlb<=safedist
        if yub-y(m,1)<=safedist
        counter = counter+1;
        end
    end
    
    if x(m,1)-xlb<=safedist
        if y(m,1)-ylb<=safedist
        counter = counter+1;
        end
    end
end

% Here, we initialize OVERLAP vectors to hold coordinates for tiled nodes.

xoverlap = zeros(counter,1);
yoverlap = zeros(counter,1);

xspan = xub-xlb;
yspan = yub-ylb;

% The second FOR loop (below) actually assigns appropriate coordinates to
% "tiled" nodes.

counter = 1;

for m = 1:size(x,1)
    
    if xub-x(m,1)<=safedist
        xoverlap(counter,1)=x(m,1)-xspan;
        yoverlap(counter,1)=y(m,1);
        counter = counter+1;
    end
    
    if x(m,1)-xlb<=safedist
        xoverlap(counter,1)=x(m,1)+xspan;
        yoverlap(counter,1)=y(m,1);
        counter = counter+1;
    end
    
    if yub-y(m,1)<=safedist
        yoverlap(counter,1)=y(m,1)-yspan;
        xoverlap(counter,1)=x(m,1);
        counter = counter+1;
    end
    
    if y(m,1)-ylb<=safedist
        yoverlap(counter,1)=y(m,1)+yspan;
        xoverlap(counter,1)=x(m,1);
        counter = counter+1;
    end
    
    if xub-x(m,1)<=safedist
        if yub-y(m,1)<=safedist
        xoverlap(counter,1)=x(m,1)-xspan;
        yoverlap(counter,1)=y(m,1)-yspan;
        counter = counter+1;
        end
    end
    
    if xub-x(m,1)<=safedist
        if y(m,1)-ylb<=safedist
        xoverlap(counter,1)=x(m,1)-xspan;
        yoverlap(counter,1)=y(m,1)+yspan;
        counter = counter+1;
        end
    end
    
    if x(m,1)-xlb<=safedist
        if y(m,1)-ylb<=safedist
        xoverlap(counter,1)=x(m,1)+xspan;
        yoverlap(counter,1)=y(m,1)+yspan;
        counter = counter+1;
        end
    end
    
    if x(m,1)-xlb<=safedist
        if yub-y(m,1)<=safedist
        xoverlap(counter,1)=x(m,1)+xspan;
        yoverlap(counter,1)=y(m,1)-yspan;
        counter = counter+1;
        end
    end
    
end

xaug = [x; xoverlap];  % x- and y- coordinates augmented with tiled nodes
yaug = [y; yoverlap];
    
% Now, we actually find nearest neighbors in the augmented domain

[idx dist] = knnsearch([xaug yaug],[x y],'k',nnn);

idxmod = idx(:,2:nnn);
distmod = dist(:,2:nnn);

oneoverr5 = 1./(distmod.^5); % we use a relation of 1/r^4 in calculating
                             % an electrostatic repulsion vector.  One
                             % power of radius is already included in the
                             % numerator (effectively)      
    
for m = 1:size(x,1)
    
    if moveFlagVec(m,1) ~= 0
        
        for n = 1:nnn-1
            
            % we sum up contributions to the force vector from NNN nearest
            % neighbors...
            
            xcorr(m,1) = xcorr(m,1)+(x(m,1)-xaug(idxmod(m,n),1)) ...
                *oneoverr5(m,n);
            ycorr(m,1) = ycorr(m,1)+(y(m,1)-yaug(idxmod(m,n),1)) ...
                *oneoverr5(m,n);
            
        end
        
        % ... and then normalize those force vectors.
        
        xycorrlength = (xcorr(m,1)^2+ycorr(m,1)^2)^0.5;
        xcorr(m,1) = xcorr(m,1)/xycorrlength;
        ycorr(m,1) = ycorr(m,1)/xycorrlength;
        
    end
    
end

% we move each node by a distance DELTA.

xnew = x+xcorr*delta;
ynew = y+ycorr*delta;

% And finally, we make sure to "wrap" the nodes to the other side of the
% domain if they exit it through the movement step.

hCart = (xub-xlb)/size(x,1)^0.5;
hInt = 1/numIntNodes;

for k = 1:size(xnew,1)
    
    if moveFlagVec(k,1) ~= 0
        if allRBFflag == 0
            if ynew(k,1)-curvedinterface2(xnew(k,1),yUI) >= 8.3*hCart
                xnew(k,1) = xlb+rand*(xlb-xub);
                temp = rand;
                if temp < 0.25
                    ynew(k,1) = curvedinterface1(xnew(k,1),yLI)-4*hCart;
                else
                    if temp < 0.5
                        ynew(k,1) = curvedinterface1(xnew(k,1),yLI)+4*hCart;
                    else
                        if temp < 0.75
                            ynew(k,1) = curvedinterface2(xnew(k,1),yUI)-4*hCart;
                        else
                            ynew(k,1) = curvedinterface2(xnew(k,1),yUI)+4*hCart;
                        end
                    end
                end
            end
            if ynew(k,1)-curvedinterface1(xnew(k,1),yLI) <= -8.3*hCart
                xnew(k,1) = xlb+rand*(xlb-xub);
                temp = rand;
                if temp < 0.25
                    ynew(k,1) = curvedinterface1(xnew(k,1),yLI)-4*hCart;
                else
                    if temp < 0.5
                        ynew(k,1) = curvedinterface1(xnew(k,1),yLI)+4*hCart;
                    else
                        if temp < 0.75
                            ynew(k,1) = curvedinterface2(xnew(k,1),yUI)-4*hCart;
                        else
                            ynew(k,1) = curvedinterface2(xnew(k,1),yUI)+4*hCart;
                        end
                    end
                end
            end
            %
            if (ynew(k,1) < curvedinterface2(xnew(k,1),yUI)) && ...
                    (ynew(k,1) > curvedinterface1(xnew(k,1),yLI))
                if (ynew(k,1)-curvedinterface2(xnew(k,1),yUI) <= -8.3*hCart) && ...
                        ynew(k,1)-curvedinterface1(xnew(k,1),yLI) >= 8.3*hCart
                    xnew(k,1) = xlb+rand*(xlb-xub);
                    temp = rand;
                    if temp < 0.25
                        ynew(k,1) = curvedinterface1(xnew(k,1),yLI)-4*hCart;
                    else
                        if temp < 0.5
                            ynew(k,1) = curvedinterface1(xnew(k,1),yLI)+4*hCart;
                        else
                            if temp < 0.75
                                ynew(k,1) = curvedinterface2(xnew(k,1),yUI)-4*hCart;
                            else
                                ynew(k,1) = curvedinterface2(xnew(k,1),yUI)+4*hCart;
                            end
                        end
                    end
                end
            end
        end
        %}
        if intWrapFlag ~= 0
            
            if abs(ynew(k,1)-curvedinterface2(xnew(k,1),yUI)) <= (3^0.5+0.5)*hInt
                xnew(k,1) = xlb+rand*(xlb-xub);
                temp = rand;
                if temp < 0.25
                    ynew(k,1) = curvedinterface1(xnew(k,1),yLI)-4*hCart;
                else
                    if temp < 0.5
                        ynew(k,1) = curvedinterface1(xnew(k,1),yLI)+4*hCart;
                    else
                        if temp < 0.75
                            ynew(k,1) = curvedinterface2(xnew(k,1),yUI)-4*hCart;
                        else
                            ynew(k,1) = curvedinterface2(xnew(k,1),yUI)+4*hCart;
                        end
                    end
                end
            end
            if abs(ynew(k,1)-curvedinterface1(xnew(k,1),yLI)) <= (3^0.5+0.5)*hInt
                xnew(k,1) = xlb+rand*(xlb-xub);
                temp = rand;
                if temp < 0.25
                    ynew(k,1) = curvedinterface1(xnew(k,1),yLI)-4*hCart;
                else
                    if temp < 0.5
                        ynew(k,1) = curvedinterface1(xnew(k,1),yLI)+4*hCart;
                    else
                        if temp < 0.75
                            ynew(k,1) = curvedinterface2(xnew(k,1),yUI)-4*hCart;
                        else
                            ynew(k,1) = curvedinterface2(xnew(k,1),yUI)+4*hCart;
                        end
                    end
                end
            end

        end
    
    end
    
    if xnew(k,1) >= xub
        xnew(k,1) = xnew(k,1)-xspan;
    end
    if xnew(k,1) <= xlb
        xnew(k,1) = xnew(k,1)+xspan;
    end
    if ynew(k,1) >= yub
        ynew(k,1) = ynew(k,1)-yspan;
    end
    if ynew(k,1) <= ylb
        ynew(k,1) = ynew(k,1)+yspan;
    end
    
    
        
end

%%% End function mos2dsqperiodic %%%

end



function [sparseLmatrix tripleVec] = ...
    createRBFLoperator1(xnodes,ynodes,xlb,xub,ylb,yub,shp,stencilsize,...
    polydegree,lambdaMu1,rho1,lambdaMu2,rho2,polydegreeInt,...
    stencilsizeInt,naiveFlag,Lk,yLI,yUI,queryFactor,hyperFlag,RBFflagVec,...
    lambdaMuVec,expStack)

% First, we define the size of our operator.

N1 = size(xnodes,1);

% Next, we call FINDPERIODICNEIGHBORS6 to find periodic nearest neighbors for all 
% RBF-FD nodes in the domain.  ORIGIDXSTACK holds the original indices
% (locations in [xnodes ynodes] vector pair) for the
% STENCILSIZE nearest neighbors (columns) of each of the N1 RBF-FD evaluation nodes
% (rows).

origidxstack = findperiodicneighbors6(xnodes,ynodes,xnodes,ynodes,xlb,xub,ylb,yub,stencilsize);
origidxstackInt = findperiodicneighbors6(xnodes,ynodes,xnodes,ynodes,xlb,xub,ylb,yub,stencilsizeInt);

% the KNNSEARCH below outputs the squared euclidean distance between each
% RBF-FD node and its 3rd nearest neighbor in the stack MINDISTANCES2.  This
% will be used to normalize distances in each RBF-FD stencil (improving
% condition of the system we need to solve).

% Distances to the FURTHEST neighbor of a stencil
% will be used to normalize distances in polynomial evaluation.

[idx distances] = ...
    knnsearch([xnodes ynodes],[xnodes ynodes],'k',stencilsize);
mindistances2 = distances(:,4).^2;
normfactorint = distances(:,stencilsize);

% sparse storage initialization (apologies for the blatant hardcoding here!
% :) )

sparseiL1 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL1 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL1 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL2 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL2 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL2 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL3 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL3 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL3 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL4 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL4 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL4 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL5 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL5 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL5 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL6 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL6 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL6 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL7 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL7 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL7 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL8 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL8 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL8 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL9 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL9 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL9 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL10 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL10 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL10 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL11 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL11 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL11 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL12 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL12 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL12 = zeros(max(stencilsize,stencilsizeInt),N1);

sparseiL13 = ones(max(stencilsize,stencilsizeInt),N1);   % initializing sparse indices
sparsejL13 = ones(max(stencilsize,stencilsizeInt),N1);
sparsewL13 = zeros(max(stencilsize,stencilsizeInt),N1);
                                   
xspan = xub-xlb;                   
yspan = yub-ylb;                     

tripleVec = zeros(N1,1);

onesStencil = ones(stencilsize,1);
onesStencilInt = ones(stencilsizeInt,1);

stencilsizeTemp = stencilsizeInt;

onesPad = ones(stencilsize-stencilsizeInt,1);
zerosPad = zeros(stencilsize-stencilsizeInt,1);

parfor m = 1:N1
    
    lambdaMuTemp = lambdaMuVec(m,1);
    
    if RBFflagVec(m,1) ~= 0
        xpositionsInt = zeros(1,stencilsizeInt);
        ypositionsInt = zeros(1,stencilsizeInt);
        stencilZoneVec = zeros(1,stencilsizeInt);
        
        xpositions = zeros(1,stencilsize);
        ypositions = zeros(1,stencilsize);
        
        if ynodes(m,1) >= curvedinterface1(xnodes(m,1),yLI) && ynodes(m,1) <= curvedinterface2(xnodes(m,1),yUI)
            rhoTemp = rho2;
            %lambdaMuTemp = lambdaMu2;
            rhoEval = rho2;
            rhoAcross = rho1;
        else
            rhoTemp = rho1;
            %lambdaMuTemp = lambdaMu1;
            rhoEval = rho1;
            rhoAcross = rho2;
        end
        
        A = ones(stencilsize,stencilsize);
        RHS = zeros(stencilsize,1);
        
        RHSdx = RHS;
        RHSdy = RHS;
        
        Aint = ones(stencilsizeInt,stencilsizeInt);
        RHSint = zeros(stencilsizeInt,1);
        RHSdxint = RHSint;
        RHSdyint = RHSint;
        
        seqFlag = 0;
        %
        hApprox = 1/N1^0.5;
        if (abs(ynodes(m,1)-curvedinterface1(xnodes(m,1),yLI)) < queryFactor*hApprox) || (abs(ynodes(m,1)-curvedinterface2(xnodes(m,1),yUI)) < queryFactor*hApprox)
            seqFlag = 1;
        end
        %}
        
        
        
        if naiveFlag == 1
            %tempidxstack = origidxstack(m,:);
            seqFlag = 0;
        else
            %tempidxstack = origidxstackSeq1(m,:);
        end
        
        if seqFlag == 0
            
            tempidxstack = origidxstack(m,:);
            
            for n = 1:stencilsize
                
                xpositions(1,n)=xnodes(tempidxstack(1,n),1);
                ypositions(1,n)=ynodes(tempidxstack(1,n),1);
                
                % Below, each node in every stencil is checked for the need to
                % "tile" or "ghost" that node to another part of the domain.  In
                % the future, this step could probably be incorporated into
                % FINDPERIODICNEIGHBORS6 for speed.
                
                if abs(xpositions(1,1)-(xpositions(1,n)+xspan))<...
                        abs(xpositions(1,1)-xpositions(1,n))
                    xpositions(1,n)=xpositions(1,n)+xspan;
                end
                
                if abs(xpositions(1,1)-(xpositions(1,n)-xspan))<...
                        abs(xpositions(1,1)-xpositions(1,n))
                    xpositions(1,n)=xpositions(1,n)-xspan;
                end
                
                if abs(ypositions(1,1)-(ypositions(1,n)+yspan))<...
                        abs(ypositions(1,1)-ypositions(1,n))
                    ypositions(1,n)=ypositions(1,n)+yspan;
                end
                
                if abs(ypositions(1,1)-(ypositions(1,n)-yspan))<...
                        abs(ypositions(1,1)-ypositions(1,n))
                    ypositions(1,n)=ypositions(1,n)-yspan;
                end
            end
            
            distmin2 = mindistances2(m,1);  % DISTMIN2 normalizes distances for
            % each RBF-FD stencil so we don't have
            % to choose optimal shape parameters for each stencil separately - with
            % this done and with IMQ RBFs, a uniform SHP of 0.2 works quite well.
            
            % the Laguerre polynomial, coefficient, and dimension of the space (2)
            % below will help evaluate the Laplacian of Gaussian RBFs, as described
            % in Fornberg and Lehto (2011).
            
            pk = zeros(Lk+1,1);
            coeff = shp^(2*Lk)/(distmin2)^Lk;
            d = 2;
            
            for i = 1:stencilsize        % creating the A matrix.
                for j = 1:stencilsize
                    A(i,j) = ...
                        exp(-shp^2*...
                        ((xpositions(1,j)-xpositions(1,i))^2+...
                        (ypositions(1,j)-ypositions(1,i))^2)/distmin2);
                end
                
                if hyperFlag > 0
                    % preparing to create the RHS, for hypervisc. weights
                    
                    er2dbrmin2 = shp^2*((xpositions(1,1)-xpositions(1,i))^2+...
                        (ypositions(1,1)-ypositions(1,i))^2)/distmin2;
                    pk(1,1) = 1;
                    pk(2,1) = 4*er2dbrmin2-2*d;
                    
                    for k = 2:Lk
                        pk(k+1,1) = 4*(er2dbrmin2-2*(k-1)-d/2)*pk(k,1)...
                            -8*(k-1)*(2*(k-1)-2+d)*pk(k-1,1);
                    end
                    
                    % creating the RHS via Laguerre polynomial (Fornberg and Lehto,
                    % 2011)
                    
                    RHS(i,1) = coeff*pk(Lk+1)*A(i,1);
                else
                    RHSdx(i,1) = -2*shp^2/distmin2*(xpositions(1,1)-xpositions(1,i))*A(i,1);
                    RHSdy(i,1) = -2*shp^2/distmin2*(ypositions(1,1)-ypositions(1,i))*A(i,1);
                end
            end
            
            % appending the A matrix with polynomials via AUGHYPERWITHPOLYS6
            if hyperFlag > 0
                hyperweightstemp = augHypwithpolys1(A,RHS,xpositions,ypositions,...
                    xpositions(1,1),ypositions(1,1),polydegree,1/normfactorint(m,1),Lk);
                
                % setting up indices and weights for sparse hyp. operator creation
                
                % u hyp
                sparseiL1(:,m) = m*onesStencil;
                sparsejL1(:,m) = tempidxstack';
                sparsewL1(:,m) = (-1)^(Lk+1)*hyperweightstemp';
                
                % v hyp
                sparseiL2(:,m) = m*onesStencil+N1;
                sparsejL2(:,m) = tempidxstack'+N1;
                sparsewL2(:,m) = (-1)^(Lk+1)*hyperweightstemp';
                
                % f hyp
                sparseiL3(:,m) = m*onesStencil+2*N1;
                sparsejL3(:,m) = tempidxstack'+2*N1;
                sparsewL3(:,m) = (-1)^(Lk+1)*hyperweightstemp';
                
                % g hyp
                sparseiL4(:,m) = m*onesStencil+3*N1;
                sparsejL4(:,m) = tempidxstack'+3*N1;
                sparsewL4(:,m) = (-1)^(Lk+1)*hyperweightstemp';
                
                % h hyp
                sparseiL5(:,m) = m*onesStencil+4*N1;
                sparsejL5(:,m) = tempidxstack'+4*N1;
                sparsewL5(:,m) = (-1)^(Lk+1)*hyperweightstemp';
                
                
            else
                [dxweightstemp dyweightstemp] = augDxDywithpolys1(A,RHSdx,RHSdy,xpositions,ypositions,...
                    xpositions(1,1),ypositions(1,1),polydegree,1/normfactorint(m,1),Lk);
                
                % setting up indices and weights for sparse hyp. operator creation
                
                
                
                % ROW 1: ut = (fx + gy)/rho
                sparseiL1(:,m) = m*onesStencil+0*N1;
                sparsejL1(:,m) = tempidxstack'+2*N1;
                sparsewL1(:,m) = dxweightstemp'/rhoTemp;
                
                sparseiL2(:,m) = m*onesStencil+0*N1;
                sparsejL2(:,m) = tempidxstack'+3*N1;
                sparsewL2(:,m) = dyweightstemp'/rhoTemp;
                
                % ROW 2: vt = (gx + fy)/rho
                sparseiL3(:,m) = m*onesStencil+1*N1;
                sparsejL3(:,m) = tempidxstack'+3*N1;
                sparsewL3(:,m) = dxweightstemp'/rhoTemp;
                
                sparseiL4(:,m) = m*onesStencil+1*N1;
                sparsejL4(:,m) = tempidxstack'+4*N1;
                sparsewL4(:,m) = dyweightstemp'/rhoTemp;
                
                % ROW 3: ft = (lambda+2*mu)*ux + lambda*vy
                sparseiL5(:,m) = m*onesStencil+2*N1;
                sparsejL5(:,m) = tempidxstack'+0*N1;
                sparsewL5(:,m) = dxweightstemp'*3*lambdaMuTemp;
                
                sparseiL6(:,m) = m*onesStencil+2*N1;
                sparsejL6(:,m) = tempidxstack'+1*N1;
                sparsewL6(:,m) = dyweightstemp'*lambdaMuTemp;
                
                % ROW 4: gt = mu*uy + mu*vx
                sparseiL7(:,m) = m*onesStencil+3*N1;
                sparsejL7(:,m) = tempidxstack'+0*N1;
                sparsewL7(:,m) = dyweightstemp'*lambdaMuTemp;
                
                sparseiL8(:,m) = m*onesStencil+3*N1;
                sparsejL8(:,m) = tempidxstack'+1*N1;
                sparsewL8(:,m) = dxweightstemp'*lambdaMuTemp;
                
                % ROW 5: ht = lambda*ux + (lambda+2*mu)*vy
                sparseiL9(:,m) = m*onesStencil+4*N1;
                sparsejL9(:,m) = tempidxstack'+0*N1;
                sparsewL9(:,m) = dxweightstemp'*lambdaMuTemp;
                
                sparseiL10(:,m) = m*onesStencil+4*N1;
                sparsejL10(:,m) = tempidxstack'+1*N1;
                sparsewL10(:,m) = dyweightstemp'*3*lambdaMuTemp;
                
            end
            
            %%% START HERE AFTER LUNCH - recall new functions for aughyp, augdx/y %%%
            
        else
            
            tempidxstack = origidxstackInt(m,:);
            
            tripleCheck = [0 0 0];
            
            %tempidxstack = origidxstackInt(m,:);
            
            for n = 1:stencilsizeInt
                
                xpositionsInt(1,n)=xnodes(origidxstackInt(m,n),1);
                ypositionsInt(1,n)=ynodes(origidxstackInt(m,n),1);
                
                stencilZoneVec(1,n) = 1;
                
                if ypositionsInt(1,n) > curvedinterface1(xpositionsInt(1,n),yLI)
                    stencilZoneVec(1,n) = 2;
                    if ypositionsInt(1,n) > curvedinterface2(xpositionsInt(1,n),yUI)
                        stencilZoneVec(1,n) = 3;
                    end
                end
                
                if stencilZoneVec(1,n) == 1
                    tripleCheck(1,1) = 1;
                end
                if stencilZoneVec(1,n) == 2
                    tripleCheck(1,2) = 1;
                end
                if stencilZoneVec(1,n) == 3
                    tripleCheck(1,3) = 1;
                end
                
                % Below, each node in every stencil is checked for the need to
                % "tile" or "ghost" that node to another part of the domain.  In
                % the future, this step could probably be incorporated into
                % FINDPERIODICNEIGHBORS6 for speed.
                
                if abs(xpositionsInt(1,1)-(xpositionsInt(1,n)+xspan))<...
                        abs(xpositionsInt(1,1)-xpositionsInt(1,n))
                    xpositionsInt(1,n)=xpositionsInt(1,n)+xspan;
                end
                
                if abs(xpositionsInt(1,1)-(xpositionsInt(1,n)-xspan))<...
                        abs(xpositionsInt(1,1)-xpositionsInt(1,n))
                    xpositionsInt(1,n)=xpositionsInt(1,n)-xspan;
                end
                
                if abs(ypositionsInt(1,1)-(ypositionsInt(1,n)+yspan))<...
                        abs(ypositionsInt(1,1)-ypositionsInt(1,n))
                    ypositionsInt(1,n)=ypositionsInt(1,n)+yspan;
                end
                
                if abs(ypositionsInt(1,1)-(ypositionsInt(1,n)-yspan))<...
                        abs(ypositionsInt(1,1)-ypositionsInt(1,n))
                    ypositionsInt(1,n)=ypositionsInt(1,n)-yspan;
                end
                
            end
            
            tripleFlag = tripleCheck(1,1)*tripleCheck(1,3);
            
            tripleVec(m,1) = tripleFlag;
            
            evalZone = 1;
            
            if ypositionsInt(1,1) < (curvedinterface1(xpositionsInt(1,1),yLI)+curvedinterface2(xpositionsInt(1,1),yUI))/2
                [xClosest theta] = pointFinder1(xpositionsInt(1,1),ypositionsInt(1,1),yLI);
                yClosest = curvedinterface1(xClosest,yLI);
                %yIntLoc = 0.25;
                if ypositionsInt(1,1) > curvedinterface1(xpositionsInt(1,1),yLI)
                    evalZone = 2;
                end
            else
                [xClosest theta] = pointFinder2(xpositionsInt(1,1),ypositionsInt(1,1),yUI);
                yClosest = curvedinterface2(xClosest,yUI);
                %yIntLoc = 0.5;
                
                evalZone = 3;
                
                if ypositionsInt(1,1) > curvedinterface2(xpositionsInt(1,1),yUI)
                    evalZone = 4;
                end
            end
            
            xeval = xpositionsInt(1,1);
            yeval = ypositionsInt(1,1);
            
            Tpoints = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            
            %xpositionsInt = xpositionsInt-xClosest;
            %ypositionsInt = ypositionsInt-yClosest;
            ypositionsIntWarp = ypositionsInt;
            
            
        
        zoneWidth = yUI-yLI;
        %
        if tripleFlag == 0
            for n = 1:stencilsizeInt
                tempVec = Tpoints*[(xpositionsInt(1,n)-xClosest); (ypositionsInt(1,n)-yClosest)];
                xpositionsInt(1,n) = tempVec(1,1);
                ypositionsInt(1,n) = tempVec(2,1);
                ypositionsIntWarp(1,n) = ypositionsInt(1,n);
                %{
                if ((ypositionsInt(1,1)*ypositionsInt(1,n)) >= 0)
                    
                    ypositionsIntWarp(1,n) = ypositionsInt(1,n);
                else
                    ypositionsIntWarp(1,n) = rhoAcross/rhoEval*ypositionsInt(1,n);
                end
                %}
            end
        else
            for n = 1:stencilsizeInt
                tempVec = Tpoints*[(xpositionsInt(1,n)-xClosest); (ypositionsInt(1,n)-yClosest)];
                xpositionsInt(1,n) = tempVec(1,1);
                ypositionsInt(1,n) = tempVec(2,1);
                ypositionsIntWarp(1,n) = ypositionsInt(1,n);
                %{
                if evalZone == 1
                    if stencilZoneVec(1,n) == 2
                        ypositionsIntWarp(1,n) = rho2/rho1*ypositionsInt(1,n);
                    end
                    if stencilZoneVec(1,n) == 3
                        ypositionsIntWarp(1,n) = rho2/rho1*zoneWidth + ...
                            ypositionsInt(1,n)-zoneWidth;
                    end
                end
                if evalZone == 2
                    if stencilZoneVec(1,n) == 1
                        ypositionsIntWarp(1,n) = rho1/rho2*ypositionsInt(1,n);
                    end
                    if stencilZoneVec(1,n) == 3
                        ypositionsIntWarp(1,n) = zoneWidth + ...
                            (ypositionsInt(1,n)-zoneWidth)*rho1/rho2;
                    end
                end
                if evalZone == 3
                    if stencilZoneVec(1,n) == 1
                        ypositionsIntWarp(1,n) = -zoneWidth + ...
                            (ypositionsInt(1,n)+zoneWidth)*rho1/rho2;
                    end
                    if stencilZoneVec(1,n) == 3
                        ypositionsIntWarp(1,n) = rho1/rho2*ypositionsInt(1,n);
                    end
                end
                if evalZone == 4
                    if stencilZoneVec(1,n) == 2
                        ypositionsIntWarp(1,n) = rho2/rho1*ypositionsInt(1,n);
                    end
                    if stencilZoneVec(1,n) == 1
                        ypositionsIntWarp(1,n) = -rho2/rho1*zoneWidth + ...
                            ypositionsInt(1,n)+zoneWidth;
                    end
                end
                %}
            end
        end
        
            
            
            distmin2 = mindistances2(m,1);  % DISTMIN2 normalizes distances for
            % each RBF-FD stencil so we don't have
            % to choose optimal shape parameters for each stencil separately - with
            % this done and with IMQ RBFs, a uniform SHP of 0.2 works quite well.
            
            % the Laguerre polynomial, coefficient, and dimension of the space (2)
            % below will help evaluate the Laplacian of Gaussian RBFs, as described
            % in Fornberg and Lehto (2011).
            
            [c1Matrix c2Matrix c3Matrix c4Matrix pFuncsUV pFuncsFGH nullMatrix1 nullMatrix2] = ...
                continuityCreator(xpositionsInt,ypositionsInt,...
                xeval,yeval,polydegreeInt,lambdaMu1,rho1,lambdaMu2,rho2,xClosest,...
                Lk,yLI,yUI,tripleFlag,1/normfactorint(m,1),evalZone,...
                expStack(:,m),theta);
            
            pk = zeros(Lk+1,1);
            
            coeff = shp^(2*Lk)/(distmin2)^Lk;
            d = 2;
            
            for i = 1:stencilsizeInt        % creating the A matrix.
                for j = 1:stencilsizeInt
                    Aint(i,j) = ...
                        exp(-shp^2*...
                        ((xpositionsInt(1,j)-xpositionsInt(1,i))^2+...
                        (ypositionsIntWarp(1,j)-ypositionsIntWarp(1,i))^2)/distmin2);
                    
                end
                
                % preparing to create the RHS.
                
                er2dbrmin2 = shp^2*((xpositionsInt(1,1)-xpositionsInt(1,i))^2+...
                    (ypositionsIntWarp(1,1)-ypositionsIntWarp(1,i))^2)/distmin2;
                
                if hyperFlag > 0
                    
                    pk(1,1) = 1;
                    pk(2,1) = 4*er2dbrmin2-2*d;
                    
                    for k = 2:Lk
                        pk(k+1,1) = 4*(er2dbrmin2-2*(k-1)-d/2)*pk(k,1)...
                            -8*(k-1)*(2*(k-1)-2+d)*pk(k-1,1);
                        
                    end
                    
                    % creating the RHS via Laguerre polynomial (Fornberg and Lehto,
                    % 2011)
                    
                    % row 1
                    
                    RHSint(i,1) = coeff*pk(Lk+1)*Aint(i,1);
                else
                    RHSdxint(i,1) = -2*shp^2/distmin2*(xpositionsInt(1,1)-xpositionsInt(1,i))*Aint(i,1);
                    RHSdyint(i,1) = -2*shp^2/distmin2*(ypositionsIntWarp(1,1)-ypositionsIntWarp(1,i))*Aint(i,1);
                end
            end
            
            if hyperFlag > 0
                [row1w,row2w,row3w,row4w,row5w] = ...
                    augHypwithpolysInt1(Aint,RHSint,xpositionsInt,ypositionsInt,...
                    xeval,yeval,polydegreeInt,1/normfactorint(m,1),...
                    lambdaMu1,rho1,lambdaMu2,rho2,xClosest,Lk,c1Matrix,c2Matrix,c3Matrix,c4Matrix,...
                    yLI,yUI,evalZone,stencilZoneVec,tripleFlag,pFuncsUV,pFuncsFGH);
                
                
                Tbig = [(cos(theta))^2 2*sin(theta)*cos(theta) (sin(theta))^2; ...
                    -sin(theta)*cos(theta) ((cos(theta))^2-(sin(theta))^2) sin(theta)*cos(theta); ...
                    (sin(theta))^2 -2*sin(theta)*cos(theta) (cos(theta))^2];
                
                Tsmall = [cos(theta) sin(theta); -sin(theta) cos(theta)];
                
                for n = 1:stencilsizeInt
                    temp12 = [row1w(1,n) row1w(1,n+stencilsizeInt); ...
                        row2w(1,n) row2w(1,n+stencilsizeInt)];
                    
                    temp345 = [row3w(1,n) row3w(1,n+stencilsizeInt) row3w(1,n+2*stencilsizeInt); ...
                        row4w(1,n) row4w(1,n+stencilsizeInt) row4w(1,n+2*stencilsizeInt); ...
                        row5w(1,n) row5w(1,n+stencilsizeInt) row5w(1,n+2*stencilsizeInt)];
                    
                    temp12trans = Tsmall\temp12*Tsmall;
                    temp345trans = Tbig\temp345*Tbig;
                    
                    row1w(1,n) = temp12trans(1,1);
                    row1w(1,n+stencilsizeInt) = temp12trans(1,2);
                    
                    row2w(1,n) = temp12trans(2,1);
                    row2w(1,n+stencilsizeInt) = temp12trans(2,2);
                    
                    row3w(1,n) = temp345trans(1,1);
                    row3w(1,n+stencilsizeInt) = temp345trans(1,2);
                    row3w(1,n+2*stencilsizeInt) = temp345trans(1,3);
                    
                    row4w(1,n) = temp345trans(2,1);
                    row4w(1,n+stencilsizeInt) = temp345trans(2,2);
                    row4w(1,n+2*stencilsizeInt) = temp345trans(2,3);
                    
                    row5w(1,n) = temp345trans(3,1);
                    row5w(1,n+stencilsizeInt) = temp345trans(3,2);
                    row5w(1,n+2*stencilsizeInt) = temp345trans(3,3);
                    
                end
                
                % setting up indices and weights for sparse hyp. operator creation
                
                
                
                % u hyp
                sparseiL1(:,m) = [m*onesStencilInt+0*N1; onesPad];
                sparsejL1(:,m) = [tempidxstack'+0*N1; onesPad];
                sparsewL1(:,m) = [(-1)^(Lk+1)*row1w(1,1:stencilsizeTemp)'; zerosPad];
                
                sparseiL2(:,m) = [m*onesStencilInt+0*N1; onesPad];
                sparsejL2(:,m) = [tempidxstack'+1*N1; onesPad];
                sparsewL2(:,m) = [(-1)^(Lk+1)*row1w(1,1+stencilsizeTemp:2*stencilsizeTemp)'; zerosPad];
                
                % v hyp
                sparseiL3(:,m) = [m*onesStencilInt+1*N1; onesPad];
                sparsejL3(:,m) = [tempidxstack'+0*N1; onesPad];
                sparsewL3(:,m) = [(-1)^(Lk+1)*row2w(1,1:stencilsizeTemp)'; zerosPad];
                
                sparseiL4(:,m) = [m*onesStencilInt+1*N1; onesPad];
                sparsejL4(:,m) = [tempidxstack'+1*N1; onesPad];
                sparsewL4(:,m) = [(-1)^(Lk+1)*row2w(1,1+stencilsizeTemp:2*stencilsizeTemp)'; zerosPad];
                
                % f hyp
                sparseiL5(:,m) = [m*onesStencilInt+2*N1; onesPad];
                sparsejL5(:,m) = [tempidxstack'+2*N1; onesPad];
                sparsewL5(:,m) = [(-1)^(Lk+1)*row3w(1,1:stencilsizeTemp)'; zerosPad];
                
                sparseiL6(:,m) = [m*onesStencilInt+2*N1; onesPad];
                sparsejL6(:,m) = [tempidxstack'+3*N1; onesPad];
                sparsewL6(:,m) = [(-1)^(Lk+1)*row3w(1,1+stencilsizeTemp:2*stencilsizeTemp)'; zerosPad];
                
                sparseiL7(:,m) = [m*onesStencilInt+2*N1; onesPad];
                sparsejL7(:,m) = [tempidxstack'+4*N1; onesPad];
                sparsewL7(:,m) = [(-1)^(Lk+1)*row3w(1,1+2*stencilsizeTemp:3*stencilsizeTemp)'; zerosPad];
                
                % g hyp
                sparseiL8(:,m) = [m*onesStencilInt+3*N1; onesPad];
                sparsejL8(:,m) = [tempidxstack'+2*N1; onesPad];
                sparsewL8(:,m) = [(-1)^(Lk+1)*row4w(1,1:stencilsizeTemp)'; zerosPad];
                
                sparseiL9(:,m) = [m*onesStencilInt+3*N1; onesPad];
                sparsejL9(:,m) = [tempidxstack'+3*N1; onesPad];
                sparsewL9(:,m) = [(-1)^(Lk+1)*row4w(1,1+stencilsizeTemp:2*stencilsizeTemp)'; zerosPad];
                
                sparseiL10(:,m) = [m*onesStencilInt+3*N1; onesPad];
                sparsejL10(:,m) = [tempidxstack'+4*N1; onesPad];
                sparsewL10(:,m) = [(-1)^(Lk+1)*row4w(1,1+2*stencilsizeTemp:3*stencilsizeTemp)'; zerosPad];
                
                % h hyp
                sparseiL11(:,m) = [m*onesStencilInt+4*N1; onesPad];
                sparsejL11(:,m) = [tempidxstack'+2*N1; onesPad];
                sparsewL11(:,m) = [(-1)^(Lk+1)*row5w(1,1:stencilsizeTemp)'; zerosPad];
                
                sparseiL12(:,m) = [m*onesStencilInt+4*N1; onesPad];
                sparsejL12(:,m) = [tempidxstack'+3*N1; onesPad];
                sparsewL12(:,m) = [(-1)^(Lk+1)*row5w(1,1+stencilsizeTemp:2*stencilsizeTemp)'; zerosPad];
                
                sparseiL13(:,m) = [m*onesStencilInt+4*N1; onesPad];
                sparsejL13(:,m) = [tempidxstack'+4*N1; onesPad];
                sparsewL13(:,m) = [(-1)^(Lk+1)*row5w(1,1+2*stencilsizeTemp:3*stencilsizeTemp)'; zerosPad];
                
            else
                [row1w,row2w,row3w,row4w,row5w] = ...
                    augDxDywithpolysInt1(Aint,RHSdxint,RHSdyint,xpositionsInt,ypositionsInt,...
                    xeval,yeval,polydegreeInt,1/normfactorint(m,1),...
                    lambdaMu1,rho1,lambdaMu2,rho2,xClosest,Lk,c1Matrix,c2Matrix,c3Matrix,c4Matrix,...
                    yLI,yUI,evalZone,stencilZoneVec,tripleFlag,pFuncsUV,pFuncsFGH,lambdaMuVec(m,1));
                
                Tbig = [(cos(theta))^2 2*sin(theta)*cos(theta) (sin(theta))^2; ...
                    -sin(theta)*cos(theta) ((cos(theta))^2-(sin(theta))^2) sin(theta)*cos(theta); ...
                    (sin(theta))^2 -2*sin(theta)*cos(theta) (cos(theta))^2];
                
                Tsmall = [cos(theta) sin(theta); -sin(theta) cos(theta)];
                
                
                for n = 1:stencilsizeInt
                    temp12 = [row1w(1,n) row1w(1,n+stencilsizeInt) row1w(1,n+2*stencilsizeInt); ...
                        row2w(1,n) row2w(1,n+stencilsizeInt) row2w(1,n+2*stencilsizeInt)];
                    
                    temp345 = [row3w(1,n) row3w(1,n+stencilsizeInt); ...
                        row4w(1,n) row4w(1,n+stencilsizeInt); ...
                        row5w(1,n) row5w(1,n+stencilsizeInt)];
                    
                    temp12trans = Tsmall\temp12*Tbig;
                    temp345trans = Tbig\temp345*Tsmall;
                    
                    row1w(1,n) = temp12trans(1,1);
                    row1w(1,n+stencilsizeInt) = temp12trans(1,2);
                    row1w(1,n+2*stencilsizeInt) = temp12trans(1,3);
                    
                    row2w(1,n) = temp12trans(2,1);
                    row2w(1,n+stencilsizeInt) = temp12trans(2,2);
                    row2w(1,n+2*stencilsizeInt) = temp12trans(2,3);
                    
                    row3w(1,n) = temp345trans(1,1);
                    row3w(1,n+stencilsizeInt) = temp345trans(1,2);
                    
                    row4w(1,n) = temp345trans(2,1);
                    row4w(1,n+stencilsizeInt) = temp345trans(2,2);
                    
                    row5w(1,n) = temp345trans(3,1);
                    row5w(1,n+stencilsizeInt) = temp345trans(3,2);
 
                end
                % Row1 (ut = fx/rho + gy/rho)
                sparseiL1(:,m) = [m*onesStencilInt+0*N1; onesPad];
                sparsejL1(:,m) = [tempidxstack'+2*N1; onesPad];
                sparsewL1(:,m) = [row1w(1,1:stencilsizeTemp)'; zerosPad];
                
                sparseiL2(:,m) = [m*onesStencilInt+0*N1; onesPad];
                sparsejL2(:,m) = [tempidxstack'+3*N1; onesPad];
                sparsewL2(:,m) = [row1w(1,1+stencilsizeTemp:2*stencilsizeTemp)'; zerosPad];
                
                sparseiL3(:,m) = [m*onesStencilInt+0*N1; onesPad];
                sparsejL3(:,m) = [tempidxstack'+4*N1; onesPad];
                sparsewL3(:,m) = [row1w(1,1+2*stencilsizeTemp:3*stencilsizeTemp)'; zerosPad];
                
                % Row2 (vt = gx/rho + hy/rho)
                sparseiL4(:,m) = [m*onesStencilInt+1*N1; onesPad];
                sparsejL4(:,m) = [tempidxstack'+2*N1; onesPad];
                sparsewL4(:,m) = [row2w(1,1:stencilsizeTemp)'; zerosPad];
                
                sparseiL5(:,m) = [m*onesStencilInt+1*N1; onesPad];
                sparsejL5(:,m) = [tempidxstack'+3*N1; onesPad];
                sparsewL5(:,m) = [row2w(1,1+stencilsizeTemp:2*stencilsizeTemp)'; zerosPad];
                
                sparseiL6(:,m) = [m*onesStencilInt+1*N1; onesPad];
                sparsejL6(:,m) = [tempidxstack'+4*N1; onesPad];
                sparsewL6(:,m) = [row2w(1,1+2*stencilsizeTemp:3*stencilsizeTemp)'; zerosPad];
                
                % Row3 (ft = 3*lambdamu*ux + lambdamu*vy)
                sparseiL7(:,m) = [m*onesStencilInt+2*N1; onesPad];
                sparsejL7(:,m) = [tempidxstack'+0*N1; onesPad];
                sparsewL7(:,m) = [row3w(1,1:stencilsizeTemp)'; zerosPad];
                
                sparseiL8(:,m) = [m*onesStencilInt+2*N1; onesPad];
                sparsejL8(:,m) = [tempidxstack'+1*N1; onesPad];
                sparsewL8(:,m) = [row3w(1,1+stencilsizeTemp:2*stencilsizeTemp)'; zerosPad];
                
                % Row4 (gt = lambdamu*uy + lambdamu*vx)
                sparseiL9(:,m) = [m*onesStencilInt+3*N1; onesPad];
                sparsejL9(:,m) = [tempidxstack'+0*N1; onesPad];
                sparsewL9(:,m) = [row4w(1,1:stencilsizeTemp)'; zerosPad];
                
                sparseiL10(:,m) = [m*onesStencilInt+3*N1; onesPad];
                sparsejL10(:,m) = [tempidxstack'+1*N1; onesPad];
                sparsewL10(:,m) = [row4w(1,1+stencilsizeTemp:2*stencilsizeTemp)'; zerosPad];
                
                % Row5 (ht = lambdamu*ux + 3*lambdamu*vy)
                sparseiL11(:,m) = [m*onesStencilInt+4*N1; onesPad];
                sparsejL11(:,m) = [tempidxstack'+0*N1; onesPad];
                sparsewL11(:,m) = [row5w(1,1:stencilsizeTemp)'; zerosPad];
                
                sparseiL12(:,m) = [m*onesStencilInt+4*N1; onesPad];
                sparsejL12(:,m) = [tempidxstack'+1*N1; onesPad];
                sparsewL12(:,m) = [row5w(1,1+stencilsizeTemp:2*stencilsizeTemp)'; zerosPad];
                
            end
        end
    end
end

% and finally, creating the sparse RBF-FD hyperviscosity operator.

sparseiL1 = reshape(sparseiL1,N1*stencilsize,1);
sparseiL2 = reshape(sparseiL2,N1*stencilsize,1);
sparseiL3 = reshape(sparseiL3,N1*stencilsize,1);
sparseiL4 = reshape(sparseiL4,N1*stencilsize,1);
sparseiL5 = reshape(sparseiL5,N1*stencilsize,1);
sparseiL6 = reshape(sparseiL6,N1*stencilsize,1);
sparseiL7 = reshape(sparseiL7,N1*stencilsize,1);
sparseiL8 = reshape(sparseiL8,N1*stencilsize,1);
sparseiL9 = reshape(sparseiL9,N1*stencilsize,1);
sparseiL10 = reshape(sparseiL10,N1*stencilsize,1);
sparseiL11 = reshape(sparseiL11,N1*stencilsize,1);
sparseiL12 = reshape(sparseiL12,N1*stencilsize,1);
sparseiL13 = reshape(sparseiL13,N1*stencilsize,1);

sparsejL1 = reshape(sparsejL1,N1*stencilsize,1);
sparsejL2 = reshape(sparsejL2,N1*stencilsize,1);
sparsejL3 = reshape(sparsejL3,N1*stencilsize,1);
sparsejL4 = reshape(sparsejL4,N1*stencilsize,1);
sparsejL5 = reshape(sparsejL5,N1*stencilsize,1);
sparsejL6 = reshape(sparsejL6,N1*stencilsize,1);
sparsejL7 = reshape(sparsejL7,N1*stencilsize,1);
sparsejL8 = reshape(sparsejL8,N1*stencilsize,1);
sparsejL9 = reshape(sparsejL9,N1*stencilsize,1);
sparsejL10 = reshape(sparsejL10,N1*stencilsize,1);
sparsejL11 = reshape(sparsejL11,N1*stencilsize,1);
sparsejL12 = reshape(sparsejL12,N1*stencilsize,1);
sparsejL13 = reshape(sparsejL13,N1*stencilsize,1);

sparsewL1 = reshape(sparsewL1,N1*stencilsize,1);
sparsewL2 = reshape(sparsewL2,N1*stencilsize,1);
sparsewL3 = reshape(sparsewL3,N1*stencilsize,1);
sparsewL4 = reshape(sparsewL4,N1*stencilsize,1);
sparsewL5 = reshape(sparsewL5,N1*stencilsize,1);
sparsewL6 = reshape(sparsewL6,N1*stencilsize,1);
sparsewL7 = reshape(sparsewL7,N1*stencilsize,1);
sparsewL8 = reshape(sparsewL8,N1*stencilsize,1);
sparsewL9 = reshape(sparsewL9,N1*stencilsize,1);
sparsewL10 = reshape(sparsewL10,N1*stencilsize,1);
sparsewL11 = reshape(sparsewL11,N1*stencilsize,1);
sparsewL12 = reshape(sparsewL12,N1*stencilsize,1);
sparsewL13 = reshape(sparsewL13,N1*stencilsize,1);

sparseiL = [sparseiL1; sparseiL2; sparseiL3; sparseiL4; sparseiL5; ...
    sparseiL6; sparseiL7; sparseiL8; sparseiL9; sparseiL10; ...
    sparseiL11; sparseiL12; sparseiL13];

sparsejL = [sparsejL1; sparsejL2; sparsejL3; sparsejL4; sparsejL5; ...
    sparsejL6; sparsejL7; sparsejL8; sparsejL9; sparsejL10; ...
    sparsejL11; sparsejL12; sparsejL13];

sparsewL = [sparsewL1; sparsewL2; sparsewL3; sparsewL4; sparsewL5; ...
    sparsewL6; sparsewL7; sparsewL8; sparsewL9; sparsewL10; ...
    sparsewL11; sparsewL12; sparsewL13];

sparseLmatrix = sparse(sparseiL,sparsejL,sparsewL,5*N1,5*N1);

%%% End function CREATERBFHYPOPERATORS6 %%%

end


function RBFFDweights = ...
    augHypwithpolys1(A,RHS,xpositions,ypositions,xeval,yeval,...
    polydegree,normfactor,Lk)


numpolys = (polydegree+1)*(polydegree+2)/2;  % total number of terms added
n = size(A,1);                               % stencil size
augmatrix = zeros(n + numpolys);             % augmented (plus polys) matrix

augPblock = ones(n,numpolys);     % upper right and lower left blocks of LHS
augRHSblock = zeros(numpolys,1);   % additional vector for RHS (poly evals)

normfactor1 = normfactor;

xpositionsmod = normfactor1*(xpositions-xeval)'; % normalized x-coordinates
ypositionsmod = normfactor1*(ypositions-yeval)'; % normalized y-coordinates

counter1 = 0;   % counters help with poly evaluation
counter2 = 1;

for j = 2:numpolys    % here are the poly evaluations (at stencil nodes)
    
    augPblock(:,j)=xpositionsmod.^(counter2-counter1).*...
        ypositionsmod.^(counter1);

    counter1 = counter1+1;

    if counter1 > counter2
        counter1 = 0;
        counter2 = counter2+1;
    end

end

if Lk == 1
augRHSblock(4,1) = 2*normfactor1^2;
augRHSblock(6,1) = 2*normfactor1^2;
end

% building the augmented matrix...

augmatrix(1:n,1:n) = A;
augmatrix(1:n,(n+1):(n+numpolys)) = augPblock;
augmatrix((n+1):(n+numpolys),1:n) = augPblock';

newRHS = [RHS; augRHSblock];   % ...and augmented RHS

RBFFDweightstemp = (augmatrix\newRHS)';   % inverting...

RBFFDweights = RBFFDweightstemp(1,1:n);   % ...and truncating the weights.

%%% END function AUGHYPERWITHPOLYS6 %%%

end

function [dxWeights dyWeights] = ...
    augDxDywithpolys1(A,RHSdx,RHSdy,xpositions,ypositions,xeval,yeval,...
    polydegree,normfactor,Lk)


numpolys = (polydegree+1)*(polydegree+2)/2;  % total number of terms added
n = size(A,1);                               % stencil size
augmatrix = zeros(n + numpolys);             % augmented (plus polys) matrix

augPblock = ones(n,numpolys);     % upper right and lower left blocks of LHS
augRHSblockX = zeros(numpolys,1);   % additional vector for RHS (poly evals)
augRHSblockY = zeros(numpolys,1);   % additional vector for RHS (poly evals)

normfactor1 = normfactor;

xpositionsmod = normfactor1*(xpositions-xeval)'; % normalized x-coordinates
ypositionsmod = normfactor1*(ypositions-yeval)'; % normalized y-coordinates

counter1 = 0;   % counters help with poly evaluation
counter2 = 1;

for j = 2:numpolys    % here are the poly evaluations (at stencil nodes)
    
    augPblock(:,j)=xpositionsmod.^(counter2-counter1).*...
        ypositionsmod.^(counter1);

    counter1 = counter1+1;

    if counter1 > counter2
        counter1 = 0;
        counter2 = counter2+1;
    end

end

augRHSblockX(2,1) = normfactor1;
augRHSblockY(3,1) = normfactor1;

% building the augmented matrix...

augmatrix(1:n,1:n) = A;
augmatrix(1:n,(n+1):(n+numpolys)) = augPblock;
augmatrix((n+1):(n+numpolys),1:n) = augPblock';

newRHSX = [RHSdx; augRHSblockX];   % ...and augmented RHS
newRHSY = [RHSdy; augRHSblockY];   % ...and augmented RHS

RBFFDweightstempX = (augmatrix\newRHSX)';   % inverting...
RBFFDweightstempY = (augmatrix\newRHSY)';   % inverting...

dxWeights = RBFFDweightstempX(1,1:n);   % ...and truncating the weights.
dyWeights = RBFFDweightstempY(1,1:n);   % ...and truncating the weights.

%%% END function AUGHYPERWITHPOLYS6 %%%

end


function [c1Matrix c2Matrix c3Matrix c4Matrix pFuncsUV pFuncsFGH nullMatrix1 nullMatrix2] = ...
    continuityCreator(xpositionsInt,ypositionsInt,xeval,yeval,...
    polydegreeInt,lambdaMu1,rho1,lambdaMu2,rho2,xClosest,Lk,yLI,yUI,tripleFlag,normfactor,evalZone,...
    expVec,theta)

numpolysPlus = (polydegreeInt+2)*(polydegreeInt+3)/2;  % total number of terms added
numpolys = (polydegreeInt+1)*(polydegreeInt+2)/2;

multByY = zeros(numpolysPlus);
multByX = zeros(numpolysPlus);

rowCounterX = 2;
rowCounterY = 3;
skipCounterX = 1;
skipCounterY = 1;
skipLimitX = 1;
skipLimitY = 1;

for colCounter = 1:numpolysPlus
    if rowCounterX <= numpolysPlus
        multByX(rowCounterX,colCounter) = 1;
        rowCounterX = rowCounterX+1;
        if skipCounterX == skipLimitX
            rowCounterX = rowCounterX+1;
            skipCounterX = 1;
            skipLimitX = skipLimitX+1;
        else
            skipCounterX = skipCounterX+1;
        end
    end
    if rowCounterY <= numpolysPlus
        multByY(rowCounterY,colCounter) = 1;
        rowCounterY = rowCounterY+1;
        if skipCounterY == skipLimitY
            rowCounterY = rowCounterY+1;
            skipCounterY = 1;
            skipLimitY = skipLimitY+1;
        else
            skipCounterY = skipCounterY+1;
        end
    end
end

% posInd gives us which side of which interface we're on.

dxBlock = zeros(numpolysPlus);
dyBlock = dxBlock;

dxBlock(2,1) = 1;
dyBlock(3,1) = 1;

intLoc1 = 0;
intLoc2 = 0;

if tripleFlag > 0
    if (evalZone == 1) || (evalZone == 2)
        intLoc2 = normfactor*(yUI-yLI);
    end
    if (evalZone == 3) || (evalZone == 4)
        intLoc1 = normfactor*(yLI-yUI);
    end
end

xPower = []; yPower = [];

for k = 0:polydegreeInt+1
   yPower = [yPower; (0:k)'];
   xPower = [xPower; (k:-1:0)'];
end

overallPower = xPower+yPower;

LMblock = zeros(numpolysPlus);

for k = 1:numpolysPlus
    LMblock = LMblock + expVec(k,1)*(multByX^xPower(k,1))*(multByY^yPower(k,1))/(normfactor^overallPower(k,1));
end

for k = 2:polydegreeInt+1
    startIndexi = (k*(k+1))/2+1;
    finalIndexi = ((k+1)*(k+2))/2;
    
    startIndexj = ((k-1)*k)/2+1;
    finalIndexj = (k*(k+1))/2;
    
    dxBlock(startIndexi:startIndexi+(k-1),startIndexj:finalIndexj) = ...
        diag(fliplr(1:k));
    dyBlock(startIndexi+1:finalIndexi,startIndexj:finalIndexj) = ...
        diag(1:k);
end
%

nullMatrix1 = [];
nullMatrix2 = [];

blockSize = size(dxBlock,1);
zeroBlock = zeros(blockSize);

dxBlock = dxBlock';
dyBlock = dyBlock';

lambdaMuLowerTemp = lambdaMu1;
lambdaMuUpperTemp = LMblock;
rhoLowerTemp = rho1;
rhoUpperTemp = rho2;

bigOp1 = [zeroBlock zeroBlock dxBlock/rhoLowerTemp dyBlock/rhoLowerTemp zeroBlock; ...
    zeroBlock zeroBlock zeroBlock dxBlock/rhoLowerTemp dyBlock/rhoLowerTemp; ...
    3*lambdaMuLowerTemp*dxBlock lambdaMuLowerTemp*dyBlock zeroBlock zeroBlock zeroBlock; ...
    lambdaMuLowerTemp*dyBlock lambdaMuLowerTemp*dxBlock zeroBlock zeroBlock zeroBlock; ...
    lambdaMuLowerTemp*dxBlock 3*lambdaMuLowerTemp*dyBlock zeroBlock zeroBlock zeroBlock];

bigOp2 = [zeroBlock zeroBlock dxBlock/rhoUpperTemp dyBlock/rhoUpperTemp zeroBlock; ...
    zeroBlock zeroBlock zeroBlock dxBlock/rhoUpperTemp dyBlock/rhoUpperTemp; ...
    3*lambdaMuUpperTemp*dxBlock lambdaMuUpperTemp*dyBlock zeroBlock zeroBlock zeroBlock; ...
    lambdaMuUpperTemp*dyBlock lambdaMuUpperTemp*dxBlock zeroBlock zeroBlock zeroBlock; ...
    lambdaMuUpperTemp*dxBlock 3*lambdaMuUpperTemp*dyBlock zeroBlock zeroBlock zeroBlock];

rowCounter = 1;

% lower equivalence

for opPowerCounter = 0:polydegreeInt+1
    
    polyLimit = polydegreeInt-opPowerCounter;
    polyLimNo = (polyLimit+1)*(polyLimit+2)/2;
    
    tempOp1 = bigOp1^opPowerCounter;
    tempOp2 = -bigOp2^opPowerCounter;
    if mod(opPowerCounter,2) == 0
        tempStack = [tempOp1(1+0*blockSize:2*blockSize,1+0*blockSize:2*blockSize) ...
            tempOp2(1+0*blockSize:2*blockSize,1+0*blockSize:2*blockSize)];
        
        for k1 = 1:numpolysPlus
            if yPower(k1,1) > 0
                for k2 = 1:numpolysPlus
                    if (xPower(k2,1) == xPower(k1,1)) && (yPower(k2,1) == 0)
                        tempStack(k2,:) = tempStack(k2,:)+tempStack(k1,:)*intLoc1^yPower(k1,1);
                        tempStack(k2+numpolysPlus,:) = tempStack(k2+numpolysPlus,:)+tempStack(k1+numpolysPlus,:)*intLoc1^yPower(k1,1);
                    end
                end
            end
        end
        
        colCounter = 1;
        advCounter = 1;
        for k = 0:polydegreeInt+1-opPowerCounter
            nullMatrix1(rowCounter,:) = tempStack(colCounter,:);
            rowCounter = rowCounter+1;
            nullMatrix1(rowCounter,:) = tempStack(colCounter+numpolysPlus,:);
            colCounter = colCounter+advCounter;
            advCounter = advCounter+1;
            rowCounter = rowCounter+1;
        end
    else
        %{
        tempStack = [tempOp1(1+0*blockSize:1*blockSize,1+2*blockSize:3*blockSize); ...
            tempOp2(1+0*blockSize:1*blockSize,1+2*blockSize:3*blockSize)];
        colCounter = 1;
        advCounter = 1;
        for k = 0:polydegreeInt-opPowerCounter
            nullMatrix(rowCounter,:) = tempStack(:,colCounter)';
            colCounter = colCounter+advCounter;
            advCounter = advCounter+1;
            rowCounter = rowCounter+1;
        end
        %}
        %
        tempStack = [tempOp1(1+3*blockSize:5*blockSize,1+0*blockSize:2*blockSize) ...
            tempOp2(1+3*blockSize:5*blockSize,1+0*blockSize:2*blockSize)];
        
        for k1 = 1:numpolysPlus
            if yPower(k1,1) > 0
                for k2 = 1:numpolysPlus
                    if (xPower(k2,1) == xPower(k1,1)) && (yPower(k2,1) == 0)
                        tempStack(k2,:) = tempStack(k2,:)+tempStack(k1,:)*intLoc1^yPower(k1,1);
                        tempStack(k2+numpolysPlus,:) = tempStack(k2+numpolysPlus,:)+tempStack(k1+numpolysPlus,:)*intLoc1^yPower(k1,1);
                    end
                end
            end
        end
        
        colCounter = 1;
        advCounter = 1;
        for k = 0:polydegreeInt+1-opPowerCounter
            nullMatrix1(rowCounter,:) = tempStack(colCounter,:);
            rowCounter = rowCounter+1;
            nullMatrix1(rowCounter,:) = tempStack(colCounter+numpolysPlus,:);
            colCounter = colCounter+advCounter;
            advCounter = advCounter+1;
            rowCounter = rowCounter+1;
        end
        %}
    end
end

lambdaMuLowerTemp = LMblock;
lambdaMuUpperTemp = lambdaMu1;
rhoLowerTemp = rho2;
rhoUpperTemp = rho1;

bigOp1 = [zeroBlock zeroBlock dxBlock/rhoLowerTemp dyBlock/rhoLowerTemp zeroBlock; ...
    zeroBlock zeroBlock zeroBlock dxBlock/rhoLowerTemp dyBlock/rhoLowerTemp; ...
    3*lambdaMuLowerTemp*dxBlock lambdaMuLowerTemp*dyBlock zeroBlock zeroBlock zeroBlock; ...
    lambdaMuLowerTemp*dyBlock lambdaMuLowerTemp*dxBlock zeroBlock zeroBlock zeroBlock; ...
    lambdaMuLowerTemp*dxBlock 3*lambdaMuLowerTemp*dyBlock zeroBlock zeroBlock zeroBlock];

bigOp2 = [zeroBlock zeroBlock dxBlock/rhoUpperTemp dyBlock/rhoUpperTemp zeroBlock; ...
    zeroBlock zeroBlock zeroBlock dxBlock/rhoUpperTemp dyBlock/rhoUpperTemp; ...
    3*lambdaMuUpperTemp*dxBlock lambdaMuUpperTemp*dyBlock zeroBlock zeroBlock zeroBlock; ...
    lambdaMuUpperTemp*dyBlock lambdaMuUpperTemp*dxBlock zeroBlock zeroBlock zeroBlock; ...
    lambdaMuUpperTemp*dxBlock 3*lambdaMuUpperTemp*dyBlock zeroBlock zeroBlock zeroBlock];

rowCounter = 1;

% upper equivalence

for opPowerCounter = 0:polydegreeInt+1
    tempOp1 = bigOp1^opPowerCounter;
    tempOp2 = -bigOp2^opPowerCounter;
    if mod(opPowerCounter,2) == 0
        tempStack = [tempOp1(1+0*blockSize:2*blockSize,1+0*blockSize:2*blockSize) ...
            tempOp2(1+0*blockSize:2*blockSize,1+0*blockSize:2*blockSize)];
        
        for k1 = 1:numpolysPlus
            if yPower(k1,1) > 0
                for k2 = 1:numpolysPlus
                    if (xPower(k2,1) == xPower(k1,1)) && (yPower(k2,1) == 0)
                        tempStack(k2,:) = tempStack(k2,:)+tempStack(k1,:)*intLoc2^yPower(k1,1);
                        tempStack(k2+numpolysPlus,:) = tempStack(k2+numpolysPlus,:)+tempStack(k1+numpolysPlus,:)*intLoc2^yPower(k1,1);
                    end
                end
            end
        end
        
        colCounter = 1;
        advCounter = 1;
        for k = 0:polydegreeInt+1-opPowerCounter
            nullMatrix2(rowCounter,:) = tempStack(colCounter,:);
            rowCounter = rowCounter+1;
            nullMatrix2(rowCounter,:) = tempStack(colCounter+numpolysPlus,:);
            colCounter = colCounter+advCounter;
            advCounter = advCounter+1;
            rowCounter = rowCounter+1;
        end
    else
        %{
        tempStack = [tempOp1(1+0*blockSize:1*blockSize,1+2*blockSize:3*blockSize); ...
            tempOp2(1+0*blockSize:1*blockSize,1+2*blockSize:3*blockSize)];
        colCounter = 1;
        advCounter = 1;
        for k = 0:polydegreeInt-opPowerCounter
            nullMatrix(rowCounter,:) = tempStack(:,colCounter)';
            colCounter = colCounter+advCounter;
            advCounter = advCounter+1;
            rowCounter = rowCounter+1;
        end
        %}
        %
        tempStack = [tempOp1(1+3*blockSize:5*blockSize,1+0*blockSize:2*blockSize) ...
            tempOp2(1+3*blockSize:5*blockSize,1+0*blockSize:2*blockSize)];
        
        for k1 = 1:numpolysPlus
            if yPower(k1,1) > 0
                for k2 = 1:numpolysPlus
                    if (xPower(k2,1) == xPower(k1,1)) && (yPower(k2,1) == 0)
                        tempStack(k2,:) = tempStack(k2,:)+tempStack(k1,:)*intLoc2^yPower(k1,1);
                        tempStack(k2+numpolysPlus,:) = tempStack(k2+numpolysPlus,:)+tempStack(k1+numpolysPlus,:)*intLoc2^yPower(k1,1);
                    end
                end
            end
        end
        
        colCounter = 1;
        advCounter = 1;
        for k = 0:polydegreeInt+1-opPowerCounter
            nullMatrix2(rowCounter,:) = tempStack(colCounter,:);
            rowCounter = rowCounter+1;
            nullMatrix2(rowCounter,:) = tempStack(colCounter+numpolysPlus,:);
            colCounter = colCounter+advCounter;
            advCounter = advCounter+1;
            rowCounter = rowCounter+1;
        end
        %}
    end
end

for k = 1:2*numpolysPlus
    nullMatrix1(k,:) = nullMatrix1(k,:)/max(abs(nullMatrix1(k,:)));
    nullMatrix2(k,:) = nullMatrix2(k,:)/max(abs(nullMatrix2(k,:)));
end

c1Matrix = nullMatrix1(1:2*numpolysPlus,1:2*numpolysPlus);
c2Matrix = -nullMatrix1(1:2*numpolysPlus,1+2*numpolysPlus:4*numpolysPlus);

c3Matrix = nullMatrix2(1:2*numpolysPlus,1:2*numpolysPlus);
c4Matrix = -nullMatrix2(1:2*numpolysPlus,1+2*numpolysPlus:4*numpolysPlus);

numpolysUV = 2*numpolys;
numpolysFGH = 2*numpolysPlus-3;

pFuncsUV = zeros(6*numpolys,numpolysUV);
pFuncsFGH = zeros(9*numpolys,numpolysFGH);

pFuncsUV(1+2*numpolys:4*numpolys,1:2*numpolys) = eye(2*numpolys);
fghCenterTemp = bigOp1(1+2*numpolysPlus:5*numpolysPlus,1+0*numpolysPlus:2*numpolysPlus);
fghCenterTemp = [fghCenterTemp(:,2:numpolysPlus) fghCenterTemp(:,3+numpolysPlus:2*numpolysPlus)];
fghCenterTemp = [fghCenterTemp(1:numpolys,:); ...
    fghCenterTemp(1+numpolysPlus:numpolys+numpolysPlus,:); ...
    fghCenterTemp(1+2*numpolysPlus:numpolys+2*numpolysPlus,:)];
pFuncsFGH(1+3*numpolys:6*numpolys,:) = fghCenterTemp;

uvLowerTemp = c1Matrix\c2Matrix;
fghLowerTemp = bigOp2(1+2*numpolysPlus:5*numpolysPlus,1+0*numpolysPlus:2*numpolysPlus)*uvLowerTemp;
uvLowerTemp = [uvLowerTemp(:,1:numpolys) uvLowerTemp(:,1+numpolysPlus:numpolys+numpolysPlus)];
uvLowerTemp = [uvLowerTemp(1:numpolys,:); uvLowerTemp(1+numpolysPlus:numpolys+numpolysPlus,:)];
fghLowerTemp = [fghLowerTemp(:,2:numpolysPlus) fghLowerTemp(:,3+numpolysPlus:2*numpolysPlus)];
fghLowerTemp = [fghLowerTemp(1:numpolys,:); ...
    fghLowerTemp(1+numpolysPlus:numpolys+numpolysPlus,:); ...
    fghLowerTemp(1+2*numpolysPlus:numpolys+2*numpolysPlus,:)];

pFuncsUV(1+0*numpolys:2*numpolys,:) = uvLowerTemp;
pFuncsFGH(1+0*numpolys:3*numpolys,:) = fghLowerTemp;

uvUpperTemp = c4Matrix\c3Matrix;
fghUpperTemp = bigOp2(1+2*numpolysPlus:5*numpolysPlus,1+0*numpolysPlus:2*numpolysPlus)*uvUpperTemp;
uvUpperTemp = [uvUpperTemp(:,1:numpolys) uvUpperTemp(:,1+numpolysPlus:numpolys+numpolysPlus)];
uvUpperTemp = [uvUpperTemp(1:numpolys,:); uvUpperTemp(1+numpolysPlus:numpolys+numpolysPlus,:)];
fghUpperTemp = [fghUpperTemp(:,2:numpolysPlus) fghUpperTemp(:,3+numpolysPlus:2*numpolysPlus)];
fghUpperTemp = [fghUpperTemp(1:numpolys,:); ...
    fghUpperTemp(1+numpolysPlus:numpolys+numpolysPlus,:); ...
    fghUpperTemp(1+2*numpolysPlus:numpolys+2*numpolysPlus,:)];

pFuncsUV(1+4*numpolys:6*numpolys,:) = uvUpperTemp;
pFuncsFGH(1+6*numpolys:9*numpolys,:) = fghUpperTemp;

%%% END function continuityCreator %%%

end

function [row1w,row2w,row3w,row4w,row5w] = ...
    augHypwithpolysInt1(A,RHS,xpositionsInt,ypositionsInt,xeval,yeval,...
    polydegreeInt,normfactor,lambdaMu1,rho1,lambdaMu2,rho2,xClosest,Lk,c1Matrix,c2Matrix,c3Matrix,c4Matrix,...
    yLI,yUI,evalZone,stencilZoneVec,tripleFlag,pFuncsUV,pFuncsFGH)

normfactor1 = normfactor;

numpolys = (polydegreeInt+1)*(polydegreeInt+2)/2;  % total number of terms added
numpolysPlus = (polydegreeInt+2)*(polydegreeInt+3)/2;

numpolysUV = 2*numpolys;
numpolysFGH = 2*numpolysPlus-3;

n = size(A,1);                               % stencil size (total both sides)
zeroBlock = zeros(size(A));

preAugPblockUV = zeros(2*n,numpolysUV);
preAugPblockFGH = zeros(3*n,numpolysFGH);

% posInd gives us which side of which interface we're on.

rholambdaMu2LowerTemp = rho1*lambdaMu1^2;
rholambdaMu2UpperTemp = rho2*lambdaMu2^2;
rhoLowerTemp = rho1;
rhoUpperTemp = rho2;

if yeval > (curvedinterface1(xeval,yLI)+curvedinterface2(xeval,yUI))/2
    rholambdaMu2LowerTemp = rho2*lambdaMu2^2;
    rholambdaMu2UpperTemp = rho1*lambdaMu1^2;
    rhoLowerTemp = rho2;
    rhoUpperTemp = rho1;
end

xpositionsmod = normfactor1*(xpositionsInt)'; % normalized x-coordinates
upperEvalFlag = 0;
%
if yeval < (curvedinterface1(xeval,yLI)+curvedinterface2(xeval,yUI))/2
    ypositionsmod = normfactor1*(ypositionsInt)';
    if yeval > curvedinterface1(xeval,yLI)
        upperEvalFlag = 1;
    end
                                        % normalized y-coordinates
else                                               % (to nearest int.)
    ypositionsmod = normfactor1*(ypositionsInt)';
    if yeval > curvedinterface2(xeval,yUI)
        upperEvalFlag = 1;
    end
end

%%% MODIFY HEAVILY for u,v,f,g,h %%%
% moved to ContinuityCreator 
%{
pFuncsUV = zeros(6*numpolys,numpolysUV);
pFuncsFGH = zeros(9*numpolys,numpolysFGH);
pFuncs(1+numpolys:2*numpolys,1:numpolys) = eye(numpolys);

pFuncs(1:numpolys,1:numpolys) = c1Matrix\c2Matrix;
pFuncs(1+2*numpolys:3*numpolys,1:numpolys) = c4Matrix\c3Matrix;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xCoefInd = [];
yCoefInd = [];

for k3 = 0:polydegreeInt
   yCoefInd = [yCoefInd; (0:k3)'];
   xCoefInd = [xCoefInd; (k3:-1:0)'];
end

%
for m = 1:n
    
    if stencilZoneVec(1,m) == 1
        
    % defining polys - below interface zone
        
        for z = 1:size(pFuncsUV,2)
     
            % u and v (lower)
            for k1 = 1:numpolys
                
                preAugPblockUV(m,z)=preAugPblockUV(m,z)+...
                    pFuncsUV(k1+0*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockUV(m+n,z)=preAugPblockUV(m+n,z)+...
                    pFuncsUV(k1+1*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
        for z = 1:size(pFuncsFGH,2)
     
            % f,g,h (lower)
            for k1 = 1:numpolys
                
                preAugPblockFGH(m,z)=preAugPblockFGH(m,z)+...
                    pFuncsFGH(k1+0*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+n,z)=preAugPblockFGH(m+n,z)+...
                    pFuncsFGH(k1+1*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+2*n,z)=preAugPblockFGH(m+2*n,z)+...
                    pFuncsFGH(k1+2*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
    end
    
    % defining polys - mid zone    
        
    if stencilZoneVec(1,m) == 2
        
        for z = 1:size(pFuncsUV,2)
     
            % u and v (lower)
            for k1 = 1:numpolys
                
                preAugPblockUV(m,z)=preAugPblockUV(m,z)+...
                    pFuncsUV(k1+2*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockUV(m+n,z)=preAugPblockUV(m+n,z)+...
                    pFuncsUV(k1+3*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
        for z = 1:size(pFuncsFGH,2)
     
            % f,g,h (lower)
            for k1 = 1:numpolys
                
                preAugPblockFGH(m,z)=preAugPblockFGH(m,z)+...
                    pFuncsFGH(k1+3*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+n,z)=preAugPblockFGH(m+n,z)+...
                    pFuncsFGH(k1+4*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+2*n,z)=preAugPblockFGH(m+2*n,z)+...
                    pFuncsFGH(k1+5*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
    end
    
    % defining polys - top zone    
        
    if stencilZoneVec(1,m) == 3
        
        for z = 1:size(pFuncsUV,2)
     
            % u and v (lower)
            for k1 = 1:numpolys
                
                preAugPblockUV(m,z)=preAugPblockUV(m,z)+...
                    pFuncsUV(k1+4*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockUV(m+n,z)=preAugPblockUV(m+n,z)+...
                    pFuncsUV(k1+5*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
        for z = 1:size(pFuncsFGH,2)
     
            % f,g,h (lower)
            for k1 = 1:numpolys
                
                preAugPblockFGH(m,z)=preAugPblockFGH(m,z)+...
                    pFuncsFGH(k1+6*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+n,z)=preAugPblockFGH(m+n,z)+...
                    pFuncsFGH(k1+7*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+2*n,z)=preAugPblockFGH(m+2*n,z)+...
                    pFuncsFGH(k1+8*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
    end
end

preAugRHSblockRow1 = zeros(size(preAugPblockUV,2),1);
preAugRHSblockRow2 = zeros(size(preAugPblockUV,2),1);
preAugRHSblockRow3 = zeros(size(preAugPblockFGH,2),1);
preAugRHSblockRow4 = zeros(size(preAugPblockFGH,2),1);
preAugRHSblockRow5 = zeros(size(preAugPblockFGH,2),1);

xCoef2Der = xCoefInd.*(xCoefInd-1);
yCoef2Der = yCoefInd.*(yCoefInd-1);
if Lk == 1
    if evalZone == 1

        for z = 1:size(pFuncsUV,2)

            % lower u and v
            for k1 = 1:numpolys
                if xCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsUV(k1+0*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                    preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsUV(k1+1*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                end
                if yCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsUV(k1+0*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                    preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsUV(k1+1*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                end
            end
        end
        
        for z = 1:size(pFuncsFGH,2)

            % lower f,g,h
            for k1 = 1:numpolys
                if xCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsFGH(k1+0*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                    preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsFGH(k1+1*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                    preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsFGH(k1+2*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                end
                if yCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsFGH(k1+0*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                    preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsFGH(k1+1*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                    preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsFGH(k1+2*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                end
            end
        end

    end
    
    if (evalZone == 2) || (evalZone == 3)
        for z = 1:size(pFuncsUV,2)

            % lower u and v
            for k1 = 1:numpolys
                if xCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsUV(k1+2*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                    preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsUV(k1+3*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                end
                if yCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsUV(k1+2*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                    preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsUV(k1+3*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                end
            end
        end
        
        for z = 1:size(pFuncsFGH,2)

            % lower f,g,h
            for k1 = 1:numpolys
                if xCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsFGH(k1+3*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                    preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsFGH(k1+4*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                    preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsFGH(k1+5*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                end
                if yCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsFGH(k1+3*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                    preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsFGH(k1+4*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                    preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsFGH(k1+5*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                end
            end
        end
    end
    
    if evalZone == 4
        for z = 1:size(pFuncsUV,2)

            % lower u and v
            for k1 = 1:numpolys
                if xCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsUV(k1+4*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                    preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsUV(k1+5*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                end
                if yCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsUV(k1+4*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                    preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsUV(k1+5*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                end
            end
        end
        
        for z = 1:size(pFuncsFGH,2)

            % lower f,g,h
            for k1 = 1:numpolys
                if xCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsFGH(k1+6*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                    preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsFGH(k1+7*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                    preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                        normfactor1^2*xCoef2Der(k1,1)*pFuncsFGH(k1+8*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-2)*ypositionsmod(1,1)^yCoefInd(k1,1);
                end
                if yCoef2Der(k1,1) ~= 0
                    preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsFGH(k1+6*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                    preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsFGH(k1+7*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                    preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                        normfactor1^2*yCoef2Der(k1,1)*pFuncsFGH(k1+8*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-2);
                end
            end
        end
    end
end
% building the augmented matrices...

Auv = [A zeroBlock; zeroBlock A];

augAmatrixUV = [Auv preAugPblockUV; ...
    preAugPblockUV' zeros(size(preAugPblockUV,2))];

Afgh = [A zeroBlock zeroBlock; zeroBlock A zeroBlock; zeroBlock zeroBlock A];

augAmatrixFGH = [Afgh preAugPblockFGH; ...
    preAugPblockFGH' zeros(size(preAugPblockFGH,2))];

% storing weights

zeroRblock = zeros(size(RHS,1),1);

RHS1 = [RHS; zeroRblock; preAugRHSblockRow1];
RHS2 = [zeroRblock; RHS; preAugRHSblockRow2];

RHS3 = [RHS; zeroRblock; zeroRblock; preAugRHSblockRow3];
RHS4 = [zeroRblock; RHS; zeroRblock; preAugRHSblockRow4];
RHS5 = [zeroRblock; zeroRblock; RHS; preAugRHSblockRow5];

row1wtemp = (augAmatrixUV\RHS1)'; 
row1w = row1wtemp(1,1:2*n);

row2wtemp = (augAmatrixUV\RHS2)'; 
row2w = row2wtemp(1,1:2*n);

row3wtemp = (augAmatrixFGH\RHS3)'; 
row3w = row3wtemp(1,1:3*n);

row4wtemp = (augAmatrixFGH\RHS4)'; 
row4w = row4wtemp(1,1:3*n);

row5wtemp = (augAmatrixFGH\RHS5)'; 
row5w = row5wtemp(1,1:3*n);

%%% END function AUGHYPERWITHPOLYSINT2 %%%

end

function [row1w,row2w,row3w,row4w,row5w] = ...
    augDxDywithpolysInt1(A,RHSdx,RHSdy,xpositionsInt,ypositionsInt,xeval,yeval,...
    polydegreeInt,normfactor,lambdaMu1,rho1,lambdaMu2,rho2,xClosest,Lk,c1Matrix,c2Matrix,c3Matrix,c4Matrix,...
    yLI,yUI,evalZone,stencilZoneVec,tripleFlag,pFuncsUV,pFuncsFGH,lambdaMuEval)

normfactor1 = normfactor;

numpolys = (polydegreeInt+1)*(polydegreeInt+2)/2;  % total number of terms added
numpolysPlus = (polydegreeInt+2)*(polydegreeInt+3)/2;

numpolysUV = 2*numpolys;
numpolysFGH = 2*numpolysPlus-3;

n = size(A,1);                               % stencil size (total both sides)
zeroBlock = zeros(size(A));

preAugPblockUV = zeros(2*n,numpolysUV);
preAugPblockFGH = zeros(3*n,numpolysFGH);

% posInd gives us which side of which interface we're on.

xpositionsmod = normfactor1*(xpositionsInt)'; % normalized x-coordinates
upperEvalFlag = 0;
%
%lambdaMuEval = lambdaMu1;
rhoEval = rho1;

if yeval < (curvedinterface1(xeval,yLI)+curvedinterface2(xeval,yUI))/2
    ypositionsmod = normfactor1*(ypositionsInt)';
    if yeval > curvedinterface1(xeval,yLI)
        upperEvalFlag = 1;
        %lambdaMuEval = lambdaMu2;
        rhoEval = rho2;
    end
                                        % normalized y-coordinates
else                                               % (to nearest int.)
    ypositionsmod = normfactor1*(ypositionsInt)';
    if yeval > curvedinterface2(xeval,yUI)
        upperEvalFlag = 1;
    else
        %lambdaMuEval = lambdaMu2;
        rhoEval = rho2;
    end
end

xCoefInd = [];
yCoefInd = [];

for k3 = 0:polydegreeInt
   yCoefInd = [yCoefInd; (0:k3)'];
   xCoefInd = [xCoefInd; (k3:-1:0)'];
end

%
for m = 1:n
    
    if stencilZoneVec(1,m) == 1
        
    % defining polys - below interface zone
        
        for z = 1:size(pFuncsUV,2)
     
            % u and v (lower)
            for k1 = 1:numpolys
                
                preAugPblockUV(m,z)=preAugPblockUV(m,z)+...
                    pFuncsUV(k1+0*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockUV(m+n,z)=preAugPblockUV(m+n,z)+...
                    pFuncsUV(k1+1*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
        for z = 1:size(pFuncsFGH,2)
     
            % f,g,h (lower)
            for k1 = 1:numpolys
                
                preAugPblockFGH(m,z)=preAugPblockFGH(m,z)+...
                    pFuncsFGH(k1+0*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+n,z)=preAugPblockFGH(m+n,z)+...
                    pFuncsFGH(k1+1*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+2*n,z)=preAugPblockFGH(m+2*n,z)+...
                    pFuncsFGH(k1+2*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
    end
    
    % defining polys - mid zone    
        
    if stencilZoneVec(1,m) == 2
        
        for z = 1:size(pFuncsUV,2)
     
            % u and v (lower)
            for k1 = 1:numpolys
                
                preAugPblockUV(m,z)=preAugPblockUV(m,z)+...
                    pFuncsUV(k1+2*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockUV(m+n,z)=preAugPblockUV(m+n,z)+...
                    pFuncsUV(k1+3*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
        for z = 1:size(pFuncsFGH,2)
     
            % f,g,h (lower)
            for k1 = 1:numpolys
                
                preAugPblockFGH(m,z)=preAugPblockFGH(m,z)+...
                    pFuncsFGH(k1+3*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+n,z)=preAugPblockFGH(m+n,z)+...
                    pFuncsFGH(k1+4*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+2*n,z)=preAugPblockFGH(m+2*n,z)+...
                    pFuncsFGH(k1+5*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
    end
    
    % defining polys - top zone    
        
    if stencilZoneVec(1,m) == 3
        
        for z = 1:size(pFuncsUV,2)
     
            % u and v (lower)
            for k1 = 1:numpolys
                
                preAugPblockUV(m,z)=preAugPblockUV(m,z)+...
                    pFuncsUV(k1+4*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockUV(m+n,z)=preAugPblockUV(m+n,z)+...
                    pFuncsUV(k1+5*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
        for z = 1:size(pFuncsFGH,2)
     
            % f,g,h (lower)
            for k1 = 1:numpolys
                
                preAugPblockFGH(m,z)=preAugPblockFGH(m,z)+...
                    pFuncsFGH(k1+6*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+n,z)=preAugPblockFGH(m+n,z)+...
                    pFuncsFGH(k1+7*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);
                preAugPblockFGH(m+2*n,z)=preAugPblockFGH(m+2*n,z)+...
                    pFuncsFGH(k1+8*numpolys,z)*xpositionsmod(m,1)^xCoefInd(k1,1)*ypositionsmod(m,1)^yCoefInd(k1,1);

            end
        end
        
    end
end

preAugRHSblockRow1 = zeros(size(preAugPblockFGH,2),1);
preAugRHSblockRow2 = zeros(size(preAugPblockFGH,2),1);
preAugRHSblockRow3 = zeros(size(preAugPblockUV,2),1);
preAugRHSblockRow4 = zeros(size(preAugPblockUV,2),1);
preAugRHSblockRow5 = zeros(size(preAugPblockUV,2),1);

xCoef2Der = xCoefInd.*(xCoefInd-1);
yCoef2Der = yCoefInd.*(yCoefInd-1);

if evalZone == 1
    
    for z = 1:size(pFuncsFGH,2)
        
        % lower u and v
        for k1 = 1:numpolys
            if xCoefInd(k1,1) ~= 0
                preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsFGH(k1+0*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)/rhoEval;
                preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsFGH(k1+1*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)/rhoEval;
            end
            if yCoefInd(k1,1) ~= 0
                preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsFGH(k1+1*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)/rhoEval;
                preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsFGH(k1+2*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)/rhoEval;
            end
        end
    end
    
    for z = 1:size(pFuncsUV,2)
        
        % lower f,g,h
        for k1 = 1:numpolys
            if xCoefInd(k1,1) ~= 0
                preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsUV(k1+0*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)*3*lambdaMuEval;
                preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsUV(k1+1*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)*lambdaMuEval;
                preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsUV(k1+0*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)*lambdaMuEval;
            end
            if yCoefInd(k1,1) ~= 0
                preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsUV(k1+1*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)*lambdaMuEval;
                preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsUV(k1+0*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)*lambdaMuEval;
                preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsUV(k1+1*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)*3*lambdaMuEval;
            end
        end
    end
    
end

if (evalZone == 2) || (evalZone == 3)
    for z = 1:size(pFuncsFGH,2)
        
        % lower u and v
        for k1 = 1:numpolys
            if xCoefInd(k1,1) ~= 0
                preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsFGH(k1+3*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)/rhoEval;
                preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsFGH(k1+4*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)/rhoEval;
            end
            if yCoefInd(k1,1) ~= 0
                preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsFGH(k1+4*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)/rhoEval;
                preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsFGH(k1+5*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)/rhoEval;
            end
        end
    end
    
    for z = 1:size(pFuncsUV,2)
        
        % lower f,g,h
        for k1 = 1:numpolys
            if xCoefInd(k1,1) ~= 0
                preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsUV(k1+2*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)*3*lambdaMuEval;
                preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsUV(k1+3*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)*lambdaMuEval;
                preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsUV(k1+2*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)*lambdaMuEval;
            end
            if yCoefInd(k1,1) ~= 0
                preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsUV(k1+3*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)*lambdaMuEval;
                preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsUV(k1+2*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)*lambdaMuEval;
                preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsUV(k1+3*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)*3*lambdaMuEval;
            end
        end
    end
end

if evalZone == 4
    for z = 1:size(pFuncsFGH,2)
        
        % lower u and v
        for k1 = 1:numpolys
            if xCoefInd(k1,1) ~= 0
                preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsFGH(k1+6*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)/rhoEval;
                preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsFGH(k1+7*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)/rhoEval;
            end
            if yCoefInd(k1,1) ~= 0
                preAugRHSblockRow1(z,1)=preAugRHSblockRow1(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsFGH(k1+7*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)/rhoEval;
                preAugRHSblockRow2(z,1)=preAugRHSblockRow2(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsFGH(k1+8*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)/rhoEval;
            end
        end
    end
    
    for z = 1:size(pFuncsUV,2)
        
        % lower f,g,h
        for k1 = 1:numpolys
            if xCoefInd(k1,1) ~= 0
                preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsUV(k1+4*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)*3*lambdaMuEval;
                preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsUV(k1+5*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)*lambdaMuEval;
                preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                    normfactor1*xCoefInd(k1,1)*pFuncsUV(k1+4*numpolys,z)*xpositionsmod(1,1)^(xCoefInd(k1,1)-1)*ypositionsmod(1,1)^yCoefInd(k1,1)*lambdaMuEval;
            end
            if yCoefInd(k1,1) ~= 0
                preAugRHSblockRow3(z,1)=preAugRHSblockRow3(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsUV(k1+5*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)*lambdaMuEval;
                preAugRHSblockRow4(z,1)=preAugRHSblockRow4(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsUV(k1+4*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)*lambdaMuEval;
                preAugRHSblockRow5(z,1)=preAugRHSblockRow5(z,1)+...
                    normfactor1*yCoefInd(k1,1)*pFuncsUV(k1+5*numpolys,z)*xpositionsmod(1,1)^xCoefInd(k1,1)*ypositionsmod(1,1)^(yCoefInd(k1,1)-1)*3*lambdaMuEval;
            end
        end
    end
end

% building the augmented matrices...

Auv = [A zeroBlock; zeroBlock A];

augAmatrixUV = [Auv preAugPblockUV; ...
    preAugPblockUV' zeros(size(preAugPblockUV,2))];

Afgh = [A zeroBlock zeroBlock; zeroBlock A zeroBlock; zeroBlock zeroBlock A];

augAmatrixFGH = [Afgh preAugPblockFGH; ...
    preAugPblockFGH' zeros(size(preAugPblockFGH,2))];

% storing weights

zeroRblock = zeros(size(RHSdx,1),1);

RHS1 = [RHSdx/rhoEval; RHSdy/rhoEval; zeroRblock; preAugRHSblockRow1];
RHS2 = [zeroRblock; RHSdx/rhoEval; RHSdy/rhoEval; preAugRHSblockRow2];

RHS3 = [RHSdx*3*lambdaMuEval; RHSdy*lambdaMuEval; preAugRHSblockRow3];
RHS4 = [RHSdy*lambdaMuEval; RHSdx*lambdaMuEval; preAugRHSblockRow4];
RHS5 = [RHSdx*lambdaMuEval; RHSdy*3*lambdaMuEval; preAugRHSblockRow5];

row1wtemp = (augAmatrixFGH\RHS1)'; 
row1w = row1wtemp(1,1:3*n);

row2wtemp = (augAmatrixFGH\RHS2)'; 
row2w = row2wtemp(1,1:3*n);

row3wtemp = (augAmatrixUV\RHS3)'; 
row3w = row3wtemp(1,1:2*n);

row4wtemp = (augAmatrixUV\RHS4)'; 
row4w = row4wtemp(1,1:2*n);

row5wtemp = (augAmatrixUV\RHS5)'; 
row5w = row5wtemp(1,1:2*n);

%%% END function AUGDXDYPOLYSINT2 %%%

end

function [origidxstack distances] = ...
    findperiodicneighbors6(oldxnodes, oldynodes, newxnodes, newynodes, ...
    xlb, xub, ylb, yub, stencilsize)

% FINDPERIODICNEIGHBORS6 is a "tiling" method for finding periodic nearest
% neighbors of each RBF node (NEWNODES) within 
% SOME set of RBF nodes (OLDNODES)  in a periodic square - which makes it useful for local RBF 
% interpolation during dynamic node refinement (or other applications).

% INPUTS:

% oldxnodes: (N1,1) vector of x-coordinates for RBF-FD nodes.  Neighbors
%            for all NEWNODES will be found from within THIS set (paired 
%            with OLDYNODES).
% oldynodes: Same as above, for OLD y-coordinates
% newxnodes: (N2,1) vector of NEW x-coordinates.  Neighbors for each of
%            these will be found from within OLDNODES.
% newynodes: Same as above, for NEW y-coordinates
% xlb, etc.: lower and upper bounds for the problem domain in x and y
% stencilsize: # of nearest neighbors to find for a node in the domain
%              (including that node itself, if it's in the set!)

% OUTPUTS:

% origidxstack: a matrix of indices of the STENCILSIZE nearest neighbors 
%               in OLDNODES (columns) to each of the N2 NEW RBF-FD nodes 
%               (rows).  These indices point to the nearest neighbors as 
%               seen in the ORIGINAL (OLD) set of RBF-nodes
% distances:    stack of Euclidean distances that corresponds to the
%               nearest-neighbor matchings in ORIGIDXSTACK above

%%% Begin function FINDPERIODICNEIGHBORS6 %%%

% First, average distance to the "STENCILSIZE"th nearest neighbor is found
% for 20 random RBF-FD nodes...

[testidx testdist] = ...
    knnsearch([oldxnodes oldynodes],...
    [oldxnodes(1:20,1) oldynodes(1:20,1)],'k',stencilsize);

avgdist = sum(testdist(:,size(testdist,2)))/size(testdist,1);

safedist = 2*avgdist;   % ... then, this distance is doubled...

xspan = xub-xlb;        % ... and below, if nodes are located within this
yspan = yub-ylb;        % "safedist" of the domain boundary, they are
                        % "ghosted" or "replicated on the other side of
xoverlap = [];          % the periodic domain so that the subsequent
yoverlap = [];          % KNNSEARCH produces the correct identity of
origidx = [];           % nearest neighbors for each node.
counter = 1;

for m = 1:size(oldxnodes,1)
    
    if xub-oldxnodes(m,1)<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)-xspan;
        yoverlap(counter,1)=oldynodes(m,1);
        origidx(counter,1) = m;
        counter = counter+1;
    end
    
    if oldxnodes(m,1)-xlb<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)+xspan;
        yoverlap(counter,1)=oldynodes(m,1);
        origidx(counter,1) = m;
        counter = counter+1;
    end
    
    if yub-oldynodes(m,1)<=safedist
        yoverlap(counter,1)=oldynodes(m,1)-yspan;
        xoverlap(counter,1)=oldxnodes(m,1);
        origidx(counter,1) = m;
        counter = counter+1;
    end
    
    if oldynodes(m,1)-ylb<=safedist
        yoverlap(counter,1)=oldynodes(m,1)+yspan;
        xoverlap(counter,1)=oldxnodes(m,1);
        origidx(counter,1) = m;
        counter = counter+1;
    end
    
    if xub-oldxnodes(m,1)<=safedist
        if yub-oldynodes(m,1)<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)-xspan;
        yoverlap(counter,1)=oldynodes(m,1)-yspan;
        origidx(counter,1) = m;
        counter = counter+1;
        end
    end
    
    if xub-oldxnodes(m,1)<=safedist
        if oldynodes(m,1)-ylb<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)-xspan;
        yoverlap(counter,1)=oldynodes(m,1)+yspan;
        origidx(counter,1) = m;
        counter = counter+1;
        end
    end
    
    if oldxnodes(m,1)-xlb<=safedist
        if yub-oldynodes(m,1)<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)+xspan;
        yoverlap(counter,1)=oldynodes(m,1)+yspan;
        origidx(counter,1) = m;
        counter = counter+1;
        end
    end
    
    if oldxnodes(m,1)-xlb<=safedist
        if oldynodes(m,1)-ylb<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)+xspan;
        yoverlap(counter,1)=oldynodes(m,1)-yspan;
        origidx(counter,1) = m;
        counter = counter+1;
        end
    end
    
end

xaug = [oldxnodes; xoverlap];                % xaug and yaug hold the x-
yaug = [oldynodes; yoverlap];                % and y-coordinates of both
origidx = [(1:size(oldxnodes,1))'; origidx]; % the original and "tiled"
                                             % OLD RBF-FD nodes.
                                             
% ...now we search within this augmented list of both original and "tiled"
% nodes to find the appropriate, periodic nearest neighbors...                                             

[idxtemp distances] = knnsearch([xaug yaug],[newxnodes newynodes],...
    'k',stencilsize);

origidxstack = zeros(size(idxtemp));

for m = 1:size(idxtemp,1)                      % ...and finally we create
    for n = 1:size(idxtemp,2)                       % a stack of original
        origidxstack(m,n)=origidx(idxtemp(m,n),1);  % indices for the
    end                                             % STENCILSIZE nearest
end                                                 % neighbors (columns)
                               % to each RBF-FD node in the domain (rows).

%%% End function FINDPERIODICNEIGHBORS6 %%%
                               
end


function [origidxstack origidxstack2 distances1 distances2] = ...
    findperiodicneighborsSeq1(oldxnodes, oldynodes, newxnodes, newynodes, ...
    xlb, xub, ylb, yub, stencilsize, seqVec1, seqVec2)

% this is a modification of the routine described below -
% meant for limiting stencils to one side of an interface.

% FINDPERIODICNEIGHBORS6 is a "tiling" method for finding periodic nearest
% neighbors of each RBF node (NEWNODES) within 
% SOME set of RBF nodes (OLDNODES)  in a periodic square - which makes it useful for local RBF 
% interpolation during dynamic node refinement (or other applications).

% INPUTS:

% oldxnodes: (N1,1) vector of x-coordinates for RBF-FD nodes.  Neighbors
%            for all NEWNODES will be found from within THIS set (paired 
%            with OLDYNODES).
% oldynodes: Same as above, for OLD y-coordinates
% newxnodes: (N2,1) vector of NEW x-coordinates.  Neighbors for each of
%            these will be found from within OLDNODES.
% newynodes: Same as above, for NEW y-coordinates
% xlb, etc.: lower and upper bounds for the problem domain in x and y
% stencilsize: # of nearest neighbors to find for a node in the domain
%              (including that node itself, if it's in the set!)

% OUTPUTS:

% origidxstack: a matrix of indices of the STENCILSIZE nearest neighbors 
%               in OLDNODES (columns) to each of the N2 NEW RBF-FD nodes 
%               (rows).  These indices point to the nearest neighbors as 
%               seen in the ORIGINAL (OLD) set of RBF-nodes
% distances:    stack of Euclidean distances that corresponds to the
%               nearest-neighbor matchings in ORIGIDXSTACK above

%%% Begin function FINDPERIODICNEIGHBORS6 %%%

% First, average distance to the "STENCILSIZE"th nearest neighbor is found
% for 20 random RBF-FD nodes...

[testidx testdist] = ...
    knnsearch([oldxnodes oldynodes],...
    [oldxnodes(1:20,1) oldynodes(1:20,1)],'k',stencilsize);

avgdist = sum(testdist(:,size(testdist,2)))/size(testdist,1);

safedist = 2*avgdist;   % ... then, this distance is doubled...

xspan = xub-xlb;        % ... and below, if nodes are located within this
yspan = yub-ylb;        % "safedist" of the domain boundary, they are
                        % "ghosted" or "replicated on the other side of
xoverlap = [];          % the periodic domain so that the subsequent
yoverlap = [];          % KNNSEARCH produces the correct identity of
origidx = [];           % nearest neighbors for each node.
seqVec1overlap = [];
seqVec2overlap = [];
counter = 1;

for m = 1:size(oldxnodes,1)
    
    if xub-oldxnodes(m,1)<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)-xspan;
        yoverlap(counter,1)=oldynodes(m,1);
        origidx(counter,1) = m;
        seqVec1overlap(counter,1) = seqVec1(m,1);
        seqVec2overlap(counter,1) = seqVec2(m,1);
        counter = counter+1;
    end
    
    if oldxnodes(m,1)-xlb<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)+xspan;
        yoverlap(counter,1)=oldynodes(m,1);
        origidx(counter,1) = m;
        seqVec1overlap(counter,1) = seqVec1(m,1);
        seqVec2overlap(counter,1) = seqVec2(m,1);
        counter = counter+1;
    end
    
    if yub-oldynodes(m,1)<=safedist
        yoverlap(counter,1)=oldynodes(m,1)-yspan;
        xoverlap(counter,1)=oldxnodes(m,1);
        origidx(counter,1) = m;
        seqVec1overlap(counter,1) = seqVec1(m,1);
        seqVec2overlap(counter,1) = seqVec2(m,1);
        counter = counter+1;
    end
    
    if oldynodes(m,1)-ylb<=safedist
        yoverlap(counter,1)=oldynodes(m,1)+yspan;
        xoverlap(counter,1)=oldxnodes(m,1);
        origidx(counter,1) = m;
        seqVec1overlap(counter,1) = seqVec1(m,1);
        seqVec2overlap(counter,1) = seqVec2(m,1);
        counter = counter+1;
    end
    
    if xub-oldxnodes(m,1)<=safedist
        if yub-oldynodes(m,1)<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)-xspan;
        yoverlap(counter,1)=oldynodes(m,1)-yspan;
        origidx(counter,1) = m;
        seqVec1overlap(counter,1) = seqVec1(m,1);
        seqVec2overlap(counter,1) = seqVec2(m,1);
        counter = counter+1;
        end
    end
    
    if xub-oldxnodes(m,1)<=safedist
        if oldynodes(m,1)-ylb<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)-xspan;
        yoverlap(counter,1)=oldynodes(m,1)+yspan;
        origidx(counter,1) = m;
        seqVec1overlap(counter,1) = seqVec1(m,1);
        seqVec2overlap(counter,1) = seqVec2(m,1);
        counter = counter+1;
        end
    end
    
    if oldxnodes(m,1)-xlb<=safedist
        if yub-oldynodes(m,1)<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)+xspan;
        yoverlap(counter,1)=oldynodes(m,1)+yspan;
        origidx(counter,1) = m;
        seqVec1overlap(counter,1) = seqVec1(m,1);
        seqVec2overlap(counter,1) = seqVec2(m,1);
        counter = counter+1;
        end
    end
    
    if oldxnodes(m,1)-xlb<=safedist
        if oldynodes(m,1)-ylb<=safedist
        xoverlap(counter,1)=oldxnodes(m,1)+xspan;
        yoverlap(counter,1)=oldynodes(m,1)-yspan;
        origidx(counter,1) = m;
        seqVec1overlap(counter,1) = seqVec1(m,1);
        seqVec2overlap(counter,1) = seqVec2(m,1);
        counter = counter+1;
        end
    end
    
end

xaug = [oldxnodes; xoverlap];                % xaug and yaug hold the x-
yaug = [oldynodes; yoverlap];                % and y-coordinates of both
origidx = [(1:size(oldxnodes,1))'; origidx]; % the original and "tiled"
                                             % OLD RBF-FD nodes.
seqVec1aug = [seqVec1; seqVec1overlap];
seqVec2aug = [seqVec2; seqVec2overlap];

isolateVec1 = 4*seqVec2aug;
isolateVec2 = 4*seqVec1aug;

% ...now we search within this augmented list of both original and "tiled"
% nodes to find the appropriate, periodic nearest neighbors...                                             

[idxtemp1 distances1] = knnsearch([xaug+isolateVec1 yaug],...
    [newxnodes newynodes],...
    'k',stencilsize);

[idxtemp2 distances2] = knnsearch([xaug+isolateVec2 yaug],...
    [newxnodes newynodes],...
    'k',stencilsize);

origidxstack = zeros(size(idxtemp1));
origidxstack2 = zeros(size(idxtemp1));

for m = 1:size(idxtemp1,1)                      % ...and finally we create
    for n = 1:size(idxtemp1,2)                  % a stack of original
        if seqVec1(m,1) == 1
            origidxstack(m,n)=origidx(idxtemp1(m,n),1);
            origidxstack2(m,n)=origidx(idxtemp2(m,n),1);
        else
            origidxstack(m,n)=origidx(idxtemp2(m,n),1);
            origidxstack2(m,n)=origidx(idxtemp1(m,n),1);
        end
        % indices for the
    end                                             % STENCILSIZE nearest
end                                                 % neighbors (columns)
                               % to each RBF-FD node in the domain (rows).

%%% End function FINDPERIODICNEIGHBORS6 %%%
                               
end


function sparseblockmatrix = sparseblockmaker(origmatrix,rowindex,colindex,blockrowscols)

% function SPARSEBLOCKMAKER creates a big square block matrix (sparse) with
% a smaller square block matrix (origmatrix) as its only nonzero block.
% One can call this function a number of times to efficiently create a
% big, sparse block matrix.

% IN:

% origmatrix: original sparse matrix
% rowindex: the block row index that ORIGMATRIX will occupy in the sparse block
% matrix
% colindex: the block col index that ORIGMATRIX will occupy in the sparse
% block matrix
% blockrowscols: number of block rows and columns in the sparse block
% matrix

% OUT:

% sparseblockmatrix: sparse block matrix including origmatrix as one block.

%%% Begin function sparseblockmatrix %%%

Norig = size(origmatrix,1);           % size of original matrix (square)
Nnew = Norig*blockrowscols;           % size of new big matrix
onesvec = ones(size(origmatrix,1),1); % used to place orig. matrix in the
indexvec = (1:Norig)';                % correct block

% below, leftmatrix and rightmatrix left- and right-multiply the original
% matrix to create the new sparse block matrix.

leftmatrix = sparse(indexvec+(rowindex-1)*Norig,indexvec,onesvec,Nnew,Norig);
rightmatrix = sparse(indexvec,indexvec+(colindex-1)*Norig,onesvec,Norig,Nnew);

sparseblockmatrix = leftmatrix*origmatrix*rightmatrix;

%%% End function sparseblockmatrix %%%

end

function [xClosest theta] = pointFinder1(x,y,yLI)
    
    %xClosest = x; theta = 0;

    %
    xClosest = []; yClosest = [];
    
    xClosest = [xClosest; x]; yClosest = [yClosest; curvedinterface1(x,yLI)];
    theta = []; theta = [theta; 0];
    
    for k = 2:20
        [yeval yprime] = curvedinterface1(xClosest(k-1,1),yLI);
        
        b = yeval-yprime*xClosest(k-1,1);
        
        theta = atan(yprime);
        
        transform = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        invTransf = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        
        newVec1 = transform*([x; y-b]);
        newVec2 = [newVec1(1,1); 0];
        newVec3 = invTransf*newVec2 + [0; b];
        
        xClosest = [xClosest; newVec3(1,1)];
        yClosest = [yClosest; curvedinterface1(newVec3(1,1),yLI)];  
    end
    xClosest = newVec3(1,1);
  %}
end

function [xClosest theta] = pointFinder2(x,y,yUI)
    %xClosest = x; theta = 0;
    %
    xClosest = []; yClosest = [];
    
    xClosest = [xClosest; x]; yClosest = [yClosest; curvedinterface2(x,yUI)];
    theta = []; theta = [theta; 0];
    
    for k = 2:20
        [yeval yprime] = curvedinterface2(xClosest(k-1,1),yUI);
        
        b = yeval-yprime*xClosest(k-1,1);
        
        theta = atan(yprime);
        
        transform = [cos(theta) sin(theta); -sin(theta) cos(theta)];
        invTransf = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        
        newVec1 = transform*([x; y-b]);
        newVec2 = [newVec1(1,1); 0];
        newVec3 = invTransf*newVec2 + [0; b];
        
        xClosest = [xClosest; newVec3(1,1)];
        yClosest = [yClosest; curvedinterface2(newVec3(1,1),yUI)];  
    end
    xClosest = newVec3(1,1);
    %}
end