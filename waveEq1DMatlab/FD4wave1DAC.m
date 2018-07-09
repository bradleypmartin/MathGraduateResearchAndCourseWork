function [uxt,fxt,timeVec] = ...
    FD4wave1DAC(xgrid,tstep,endtime,c2,rho2,FDorder,interfaceWidth,...
    doubleNaiveFlag)

% For easy execution/trial of this function, use the runDriver1DWE.m
% script.

% The FD4wave1DAC function can be used to solve a particular 1D wave
% transport PDE/IVP in the periodic interval [-1,1).  A Gaussian wave front
% begins at x = -0.5 in this interval, and travels with wave speed c = 1
% toward a heterogeneous region that begins at x = 0 and ends at the
% user-defined x-location interfaceWidth.

% INPUTS:

% xgrid:   column vector of x-locations for discretized (FD) solution
% tstep:   discrete time step (delta-t)
% endtime: end time value (unitless) for simulation (1.5 by default)
% c2:      wave speed within the heterogeneous region (1.0 elsewhere)
% rho2:    density    "                                             "
% FDorder: order (width-1) of FD stecils used for solution updating
% interfaceWidth:  width of the heterogeneous region beginning at x = 0
% doubleNaiveFlag: governs organized (0) vs. naive (~= 0) treatment of
%                  heterogeneous layers that are small relative to the
%                  spacing of discrete data nodes in xgrid

% OUTPUTS:

% uxt: N x (#timesteps) matrix of u data (particle velocity) within the 1D
%      domain (u(x,t))
% fxt: N x (#timesteps) matrix of f data (stress) within the 1D
%      domain (f(x,t))
% timeVec: (#timesteps) x 1 column vector of (unitless) time values at
%          which snapshots of u and f are stored in uxt, fxt (resp.)

% # of data nodes
N = size(xgrid,1);

% values of density and wave speed outside heterogeneous region
rho1 = 1;
c1 = 1;

% vectors for c^2, density, and c values at each data point
c2mult = ones(size(xgrid));
rhomult = ones(size(xgrid));
cmult = ones(size(xgrid));

% defining parameter values throughout domain (changing their value in
% heterogeneous region)
for m = 1:N
    if xgrid(m,1) >= 0 && xgrid(m,1) < interfaceWidth
        c2mult(m,1)= c2^2;
        cmult(m,1) = c2;
        rhomult(m,1) = rho2;
    end
end

% # of discrete RK4 time steps to complete in simulation
timesteps = ceil(endtime/tstep)+1;

% initialization of output data
uxt = zeros(N,timesteps);
fxt = zeros(N,timesteps);
timeVec = zeros(1,timesteps);

% t = 0 (initial) data in domain (Gaussian wave front centered at x = -1/2)
uxt(:,1)=-exp(-600*(xgrid+1/2).^2);
fxt(:,1)=exp(-600*(xgrid+1/2).^2);

% differential operator initialization
uxmatrix = zeros(size(xgrid,1));
fxmatrix = uxmatrix;

% in the loop below, counter1 will progress through each of the N data
% nodes, forming stencils (updating equations) for u and f at each point.
for counter1 = 1:N
    
    % defining (periodic) indexing of stencils
    preStencilIndices = (counter1-FDorder/2):(counter1+FDorder/2);
    trueStencilIndices = mod(preStencilIndices-1,N)+1;
    
    % defining c and density for the data point
    cTemp = cmult(counter1,1);
    rhoTemp = rhomult(counter1,1);
    
    % computing relative position of nodes in stencil
    stencilTempX = xgrid(trueStencilIndices,1);
    evalX = xgrid(counter1,1);
    stencilX = stencilTempX - 2*((evalX-stencilTempX) < -1) + ...
        2*((evalX-stencilTempX) > 1);
    
    % determining need for special treatment (for stencils that cross the
    % interface into the heterogeneous region)
    ACindicator = (stencilX(1,1) <= 0)*(stencilX(FDorder+1,1) > 0) + ...
        (stencilX(1,1) <= interfaceWidth)*(stencilX(FDorder+1,1) ...
        > interfaceWidth);
    
    % determining need for additional work in the case of a stencil that
    % crosses TWO interfaces (in the case of a very thin heterogeneous
    % region)
    DoubleIndicator = (stencilX(1,1) <= 0)*(stencilX(FDorder+1,1) ...
        > interfaceWidth);
    ACindicator = ACindicator + DoubleIndicator;
    
    % if the user wishes to naively treat a thin heterogeneous region,
    % these flags can shunt a solution back in that direction.
    if (DoubleIndicator > 0) && (doubleNaiveFlag > 0)
        ACindicator = 0;
        DoubleIndicator = 0;
    end
    
    if ACindicator == 0
        
        % No special treatment is needed for stencils that do not cross an
        % interface.  Standard finite difference stencil creation is
        % carried out here.
        
        stencilX = stencilX-evalX;
        stencilScale = 1/max(abs(stencilX));
        scaledX = stencilX*stencilScale;
        preWeightsX = weights(0,scaledX',1);
        weightsX = preWeightsX(2,:)*stencilScale;
        uxmatrix(counter1,trueStencilIndices) = weightsX;
        fxmatrix(counter1,trueStencilIndices) = weightsX;
    else
        if DoubleIndicator == 0
            % In the case of a stencil that crosses only one interface,
            % we'll have to compute a "translation" of Taylor terms across
            % that interface using continuity conditions for u and f across
            % that interface.
            
            % Immediately below, we're designating "left" (L) and "right"
            % (R) sides of the interface currently being treated, and
            % defining parameters accordingly (determined by whether we're
            % at the left or right interface bounding the heterogeneous
            % region).
            cL = c1; cR = c2; rhoL = rho1; rhoR = rho2;
            
            leftEvalFlag = 1;
            interfaceX = 0;
            
            if (evalX > 0) && (evalX < interfaceWidth)
                if (evalX > interfaceWidth/2)
                    leftEvalFlag = 1;
                    cL = c2; cR = c1;
                    rhoL = rho2; rhoR = rho1;
                    interfaceX = interfaceWidth;
                else
                    leftEvalFlag = 0;
                end
            else
                if evalX > 0
                    leftEvalFlag = 0;
                    interfaceX = interfaceWidth;
                    cL = c2; cR = c1; rhoL = rho2; rhoR = rho1;
                end
            end
            
            % As in the traditional FD stencils, we recenter our
            % expansions...this time to the interface itself as our new
            % origin.
            
            stencilX = stencilX-interfaceX;
            stencilScale = 1/max(abs(stencilX));
            scaledX = stencilX*stencilScale;
            
            % Below, we're creating a discrete version of the PDE we're
            % solving (the 1st-order form of the 2-way wave equation). This
            % discrete form is defined by how the PDE acts on vectors of
            % polynomial coefficients up to the degree specified by FDorder
            % (for example, with 4th-order stencils, we're requiring
            % accurate treatment of 1,x,x^2,x^3,x^4 terms).
            
            dxblock = diag(1:FDorder,1);
            zeroblock = zeros(size(dxblock));
            
            bigOpL = [zeroblock dxblock/rhoL; cL^2*rhoL*dxblock zeroblock];
            bigOpR = [zeroblock dxblock/rhoR; cR^2*rhoR*dxblock zeroblock];
            
            % the continuity matrices below will be used to translate
            % Taylor expansion data between one side of the interface and
            % the other (leveraging continuity of u and f at the
            % interface).
            
            continuityMatrixU = [];
            continuityMatrixF = [];
            
            numPolys = 1+FDorder;
            
            for counter2 = 0:FDorder
                preContMat = [bigOpL^counter2 -bigOpR^counter2];
                if mod(counter2,2) == 0
                    continuityMatrixU = [continuityMatrixU; ...
                        preContMat(1,1:numPolys) ...
                        preContMat(1,1+2*numPolys:3*numPolys)];
                    continuityMatrixF = [continuityMatrixF; ...
                        preContMat(1+numPolys,1+numPolys:2*numPolys) ...
                        preContMat(1+numPolys,1+3*numPolys:4*numPolys)];
                else
                    continuityMatrixU = [continuityMatrixU; ...
                        preContMat(1+numPolys,1:numPolys) ...
                        preContMat(1+numPolys,1+2*numPolys:3*numPolys)];
                    continuityMatrixF = [continuityMatrixF; ...
                        preContMat(1,1+numPolys:2*numPolys) ...
                        preContMat(1,1+3*numPolys:4*numPolys)];
                end
            end
            
            % One way to "fill out" the space of piecewise polynomials for
            % support of stencils that cross interfaces is to find a basis
            % for the null spaces of the continuity matrices (as done
            % below).  There are other ways of accomplishing this as well
            % (as you may see in other example code).
            
            supportPolysU = null(continuityMatrixU);
            supportPolysF = null(continuityMatrixF);
            
            % With the initialization and loop below, we'll fill out
            % collocation matrices and RHS's to set up the linear problem
            % that we can solve to find stencil weights across an
            % interface.
            
            uWeightsMatrix = zeros(numPolys);
            uRHS = zeros(numPolys,1);
            
            fWeightsMatrix = zeros(numPolys);
            fRHS = zeros(numPolys,1);
            
            for counterI = 1:numPolys
                
                % RHS setup
                
                if leftEvalFlag == 1
                    for counterBF = 2:numPolys
                        uRHS(counterI,1) = ...
                            uRHS(counterI,1) + ...
                            supportPolysU(counterBF,counterI)*...
                            (counterBF-1)*...
                            scaledX((numPolys+1)/2,1)^(counterBF-2);
                        fRHS(counterI,1) = ...
                            fRHS(counterI,1) + ...
                            supportPolysF(counterBF,counterI)*...
                            (counterBF-1)*...
                            scaledX((numPolys+1)/2,1)^(counterBF-2);
                    end
                else
                    for counterBF = 2:numPolys
                        uRHS(counterI,1) = ...
                            uRHS(counterI,1) + ...
                            supportPolysU(counterBF+numPolys,counterI)*...
                            (counterBF-1)*...
                            scaledX((numPolys+1)/2,1)^(counterBF-2);
                        fRHS(counterI,1) = ...
                            fRHS(counterI,1) + ...
                            supportPolysF(counterBF+numPolys,counterI)*...
                            (counterBF-1)*...
                            scaledX((numPolys+1)/2,1)^(counterBF-2);
                    end
                end
                
                % Weights matrix setup
                
                for counterJ = 1:numPolys
                    if scaledX(counterJ,1) < 0
                        for counterBF = 1:numPolys
                            uWeightsMatrix(counterI,counterJ) = ...
                                uWeightsMatrix(counterI,counterJ) + ...
                                supportPolysU(counterBF,counterI)*...
                                scaledX(counterJ,1)^(counterBF-1);
                            fWeightsMatrix(counterI,counterJ) = ...
                                fWeightsMatrix(counterI,counterJ) + ...
                                supportPolysF(counterBF,counterI)*...
                                scaledX(counterJ,1)^(counterBF-1);
                        end
                    else
                        for counterBF = 1:numPolys
                            uWeightsMatrix(counterI,counterJ) = ...
                                uWeightsMatrix(counterI,counterJ) + ...
                                supportPolysU(counterBF+numPolys,counterI)*...
                                scaledX(counterJ,1)^(counterBF-1);
                            fWeightsMatrix(counterI,counterJ) = ...
                                fWeightsMatrix(counterI,counterJ) + ...
                                supportPolysF(counterBF+numPolys,counterI)*...
                                scaledX(counterJ,1)^(counterBF-1);
                        end
                    end
                end
            end
            
            uWeights = (uWeightsMatrix\uRHS)*stencilScale;
            fWeights = (fWeightsMatrix\fRHS)*stencilScale;
            
            uxmatrix(counter1,trueStencilIndices) = uWeights;
            fxmatrix(counter1,trueStencilIndices) = fWeights;
        else
            % Double-overlap AC method: here we have to account for 2
            % separate interfaces crossed by the same stencil (in the case
            % of a very thin heterogeneous region).
            
            cL = c1; cR = c2; rhoL = rho1; rhoR = rho2;
            
            EvalFlag = 0; % 0 = left; 1 = mid; 2 = right
            
            if (evalX >= 0) && (evalX < interfaceWidth)
                EvalFlag = 1;
            else
                if (evalX >= interfaceWidth/2)
                    EvalFlag = 2;
                end
            end
            
            % As with the traditional FD stencils and "single-crossings,"
            % we recenter our stencil.  Here we've made the "executive
            % decision" to recenter at the left interface of a
            % double-crossing, since our left interface is always at x = 0
            % in this simple test problem.
            
            interfaceX = 0;
            stencilX = stencilX-interfaceX;
            stencilScale = 1/max(abs(stencilX));
            scaledX = stencilX*stencilScale;
            scaledWidth = interfaceWidth*stencilScale;
            
            % The continuity matrices and discrete operators seen below
            % serve the same function here in the "double-cross" case that
            % they did for stencils that cross only one interface.
            
            dxblock = diag(1:FDorder,1);
            zeroblock = zeros(size(dxblock));
            
            bigOpL = [zeroblock dxblock/rhoL; cL^2*rhoL*dxblock zeroblock];
            bigOpR = [zeroblock dxblock/rhoR; cR^2*rhoR*dxblock zeroblock];
            
            continuityMatrixU = [];
            continuityMatrixF = [];
            
            numPolys = 1+FDorder;
            
            for counter2 = 0:FDorder
                preContMat = [bigOpL^counter2 -bigOpR^counter2];
                if mod(counter2,2) == 0
                    continuityMatrixU = [continuityMatrixU; ...
                        preContMat(1,1:numPolys) ...
                        preContMat(1,1+2*numPolys:3*numPolys)];
                    continuityMatrixF = [continuityMatrixF; ...
                        preContMat(1+numPolys,1+numPolys:2*numPolys) ...
                        preContMat(1+numPolys,1+3*numPolys:4*numPolys)];
                else
                    continuityMatrixU = [continuityMatrixU; ...
                        preContMat(1+numPolys,1:numPolys) ...
                        preContMat(1+numPolys,1+2*numPolys:3*numPolys)];
                    continuityMatrixF = [continuityMatrixF; ...
                        preContMat(1,1+numPolys:2*numPolys) ...
                        preContMat(1,1+3*numPolys:4*numPolys)];
                end
            end
            
            cM1U = continuityMatrixU(1:numPolys,1:numPolys);
            cM2U = -continuityMatrixU(1:numPolys,1+numPolys:2*numPolys);
            
            cM1F = continuityMatrixF(1:numPolys,1:numPolys);
            cM2F = -continuityMatrixF(1:numPolys,1+numPolys:2*numPolys);
            
            % One of the big challenges of the "double-cross" is that we've
            % got to deal with two separate expansions: one each across the
            % left and right interfaces.  The matrix inversions below will
            % help accomplish this.

            P1givenP2U = cM1U\cM2U;
            P1givenP2F = cM1F\cM2F;
            
            supportPolysU = null(continuityMatrixU);
            supportPolysF = null(continuityMatrixF);
            
            supportPolysUThird = ...
                zeros(size(supportPolysU,1)/2,size(supportPolysU,2));
            supportPolysFThird = ...
                zeros(size(supportPolysF,1)/2,size(supportPolysF,2));
            
            supportPolyUMidWidth = zeros(numPolys,1);
            supportPolyFMidWidth = zeros(numPolys,1);
            
            % Below, after we've accomplished the translation/population of
            % piecewise Taylor information from the left to the middle of a
            % thin heterogeneous region, we have to recenter the middle
            % information into a new expansion about the right interface,
            % and then compute ANOTHER translation to the right side of the
            % region to complete the "double-cross."
            
            for columnCounter = 1:numPolys
                
                supportPolyUMidTemp = fliplr(...
                    supportPolysU(numPolys+1:2*numPolys,columnCounter)');
                supportPolyFMidTemp = fliplr(...
                    supportPolysF(numPolys+1:2*numPolys,columnCounter)');
                for orderCounter = 1:numPolys
                    supportPolyUMidWidth(orderCounter,1) = polyval(...
                        supportPolyUMidTemp,scaledWidth)/factorial(orderCounter-1);
                    supportPolyUMidTemp = polyder(supportPolyUMidTemp);
                    supportPolyFMidWidth(orderCounter,1) = polyval(...
                        supportPolyFMidTemp,scaledWidth)/factorial(orderCounter-1);
                    supportPolyFMidTemp = polyder(supportPolyFMidTemp);
                end
                supportPolysUThird(1:numPolys,columnCounter) = P1givenP2U*supportPolyUMidWidth;
                supportPolysFThird(1:numPolys,columnCounter) = P1givenP2F*supportPolyFMidWidth;
                
            end
            
            % Below: prep for and execution of collocation matrix and RHS
            % formation for "double-cross" stencil determination
            
            supportPolysU = [supportPolysU; supportPolysUThird];
            supportPolysF = [supportPolysF; supportPolysFThird];
            
            uWeightsMatrix = zeros(numPolys);
            uRHS = zeros(numPolys,1);
            
            fWeightsMatrix = zeros(numPolys);
            fRHS = zeros(numPolys,1);
            
            for counterI = 1:numPolys
                
                % RHS setup
                
                if EvalFlag == 0
                    for counterBF = 2:numPolys
                        uRHS(counterI,1) = ...
                            uRHS(counterI,1) + ...
                            supportPolysU(counterBF,counterI)*...
                            (counterBF-1)*...
                            scaledX((numPolys+1)/2,1)^(counterBF-2);
                        fRHS(counterI,1) = ...
                            fRHS(counterI,1) + ...
                            supportPolysF(counterBF,counterI)*...
                            (counterBF-1)*...
                            scaledX((numPolys+1)/2,1)^(counterBF-2);
                    end
                end
                if EvalFlag == 1
                    for counterBF = 2:numPolys
                        uRHS(counterI,1) = ...
                            uRHS(counterI,1) + ...
                            supportPolysU(counterBF+numPolys,counterI)*...
                            (counterBF-1)*...
                            scaledX((numPolys+1)/2,1)^(counterBF-2);
                        fRHS(counterI,1) = ...
                            fRHS(counterI,1) + ...
                            supportPolysF(counterBF+numPolys,counterI)*...
                            (counterBF-1)*...
                            scaledX((numPolys+1)/2,1)^(counterBF-2);
                    end
                end
                if EvalFlag == 2
                    for counterBF = 2:numPolys
                        uRHS(counterI,1) = ...
                            uRHS(counterI,1) + ...
                            supportPolysU(counterBF+2*numPolys,counterI)*...
                            (counterBF-1)*...
                            (scaledX((numPolys+1)/2,1)-scaledWidth)^(counterBF-2);
                        fRHS(counterI,1) = ...
                            fRHS(counterI,1) + ...
                            supportPolysF(counterBF+2*numPolys,counterI)*...
                            (counterBF-1)*...
                            (scaledX((numPolys+1)/2,1)-scaledWidth)^(counterBF-2);
                    end
                end
                
                % Weights matrix setup
                
                for counterJ = 1:numPolys
                    if scaledX(counterJ,1) < 0
                        for counterBF = 1:numPolys
                            uWeightsMatrix(counterI,counterJ) = ...
                                uWeightsMatrix(counterI,counterJ) + ...
                                supportPolysU(counterBF,counterI)*...
                                scaledX(counterJ,1)^(counterBF-1);
                            fWeightsMatrix(counterI,counterJ) = ...
                                fWeightsMatrix(counterI,counterJ) + ...
                                supportPolysF(counterBF,counterI)*...
                                scaledX(counterJ,1)^(counterBF-1);
                        end
                    end
                    if ((scaledX(counterJ,1) >= 0) && (scaledX(counterJ,1) < scaledWidth))
                        for counterBF = 1:numPolys
                            uWeightsMatrix(counterI,counterJ) = ...
                                uWeightsMatrix(counterI,counterJ) + ...
                                supportPolysU(counterBF+numPolys,counterI)*...
                                scaledX(counterJ,1)^(counterBF-1);
                            fWeightsMatrix(counterI,counterJ) = ...
                                fWeightsMatrix(counterI,counterJ) + ...
                                supportPolysF(counterBF+numPolys,counterI)*...
                                scaledX(counterJ,1)^(counterBF-1);
                        end
                    end
                    if scaledX(counterJ,1) >= scaledWidth
                        for counterBF = 1:numPolys
                            uWeightsMatrix(counterI,counterJ) = ...
                                uWeightsMatrix(counterI,counterJ) + ...
                                supportPolysU(counterBF+2*numPolys,counterI)*...
                                (scaledX(counterJ,1)-scaledWidth)^(counterBF-1);
                            fWeightsMatrix(counterI,counterJ) = ...
                                fWeightsMatrix(counterI,counterJ) + ...
                                supportPolysF(counterBF+2*numPolys,counterI)*...
                                (scaledX(counterJ,1)-scaledWidth)^(counterBF-1);
                        end
                    end
                end
            end
            
            uWeights = (uWeightsMatrix\uRHS)*stencilScale;
            fWeights = (fWeightsMatrix\fRHS)*stencilScale;
            
            uxmatrix(counter1,trueStencilIndices) = uWeights;
            fxmatrix(counter1,trueStencilIndices) = fWeights;
        end
    end
    
end

uxmatrix = sparse(uxmatrix);
fxmatrix = sparse(fxmatrix);

%RK4 timestep

for n = 2:timesteps
    
    d1u = fxmatrix*fxt(:,n-1)./rhomult;
    step2u = uxt(:,n-1)+tstep/2*d1u;
    
    d1f = uxmatrix*uxt(:,n-1).*c2mult.*rhomult;
    step2f = fxt(:,n-1)+tstep/2*d1f;
    
    d2u = fxmatrix*step2f./rhomult;
    step3u = uxt(:,n-1)+tstep/2*d2u;
    
    d2f = uxmatrix*step2u.*c2mult.*rhomult;
    step3f = fxt(:,n-1)+tstep/2*d2f;
    
    d3u = fxmatrix*step3f./rhomult;
    step4u = uxt(:,n-1)+tstep*d3u;
    
    d3f = uxmatrix*step3u.*c2mult.*rhomult;
    step4f = fxt(:,n-1)+tstep*d3f;
    
    d4u = fxmatrix*step4f./rhomult;
    d4f = uxmatrix*step4u.*c2mult.*rhomult;
    
    uxt(:,n) = uxt(:,n-1)+tstep/6*(d1u+2*d2u+2*d3u+d4u);
    fxt(:,n) = fxt(:,n-1)+tstep/6*(d1f+2*d2f+2*d3f+d4f);
    timeVec(1,n) = tstep*(n-1);

end

