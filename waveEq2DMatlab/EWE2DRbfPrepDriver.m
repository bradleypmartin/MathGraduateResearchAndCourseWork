%%% nodeset creation and domain parameters

yLI = 0.25; % average y-value of lower interface
yUI = 0.5;  % "                                "

const1 = 1; % these constants control variability of lambda/mu on either side
const2 = 2; % of the interfaces
const3 = 2;

N = 10000;              % number of nodes desired in periodic domain
nearestneighbors = 10;   % number of nearest neighbors used in node set creation (repulsion)
totalmovements = 100;    % number of movements in node creation ("repulsions")
plotflag = 0;            % set to 0 if you don't want to plot node set; 1 if you do
nodemakingflag = 1;      % set to 1 if you're changing node sets; 0 if you're just
                         % making new operators on an existing node set
numIntNodes = round(N^0.5);
intWrapFlag = 1;
naiveFlag = 0;
allRBFflag = 0;

if nodemakingflag == 1
    xnodes = 1;
    ynodes = 1;
end

xlb = 0;     % periodic domain boundaries (at left, setup is for the unit square)
xub = 1;
ylb = 0;
yub = 1;

lambdaMu1 = 1;
rho1 = 1;

lambdaMu2 = 4;
rho2 = 2;


%%% RBF properties for dx/dy/hyperviscosity operators

IMQshp = 0.2;       % (relative) shape parameter used in IMQ RBFs (dx/dy)
GAshp = 0.4;        % (relative) shape parameter used in GA RBFs  (hyper - accuracy -> 0.3)
stencilsize = 30;   % stencil size for operators away from interface
 
polydegree = 4;     % highest degree of added polynomials for dx/dy/hyper operators
                    % (away from interface)

polydegreeInt = 3;  % degree of polynomials near/overlapping interfaces
stencilsizeInt = 19;% stencil size near/overlapping interfaces
queryFactor = 4;

% power of the Laplacian used in hyperviscosity operator
dampLk = 3;

[bigLmatrix bigDampMatrix initdatavec xnodes ynodes lambdaMuVec rhoVec tripleVec FDflagMat] = ...
    EWE2DRbfPrep(N,GAshp,stencilsize,polydegree,...
    xlb,xub,ylb,yub,lambdaMu1,lambdaMu2,rho1,rho2,xnodes,ynodes,nearestneighbors,...
    totalmovements,plotflag,numIntNodes,intWrapFlag,...
    polydegreeInt,stencilsizeInt,naiveFlag,dampLk,...
    const1,const2,const3,yLI,yUI,queryFactor,allRBFflag);
