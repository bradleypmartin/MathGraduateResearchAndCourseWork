% Script EXEPREPRBF6 carries out all the steps necessary to prepare for
% simulating wave transmission through a periodic square (periodic UNIT
% square by default).

%%% nodeset creation and domain parameters

yLI = 0.6;
yUI = 0.8;

closedFlag = 0;
thinFlag = 0;
curvedFlag = 0;

const1 = 0;
const2 = 4;
const3 = 4;

N = 1250;               % number of nodes desired in periodic domain
nearestneighbors = 10;   % number of nearest neighbors used in node set creation (repulsion)
totalmovements = 100;    % number of movements in node creation ("repulsions")
plotflag = 0;            % set to 0 if you don't want to plot node set forming; 1 if you do
nodemakingflag = 1;      % set to 1 if you're changing node sets; 0 if you're just
RBFwarpFlag = 1;

numIntNodes = round(0.95*N^0.5);
intWrapFlag = 1;
naiveFlag = 0;

if nodemakingflag == 1
    xnodes = 1;
    ynodes = 1;
end

xlb = 0;     % periodic domain boundaries (at left, setup is for the unit square)
xub = 1;
ylb = 0;
yub = 1;

c1 = 1;       % value of parameter mu below the interface
rho1 = 1;

c2 = 1;
rho2 = 1/5;


%%% RBF properties for dx/dy/hyperviscosity operators

IMQshp = 0.2;       % (relative) shape parameter used in IMQ RBFs (dx/dy)
GAshp = 0.4;        % (relative) shape parameter used in GA RBFs  (hyper - accuracy -> 0.3)

stencilsize = 19;
polydegree = 3;

stencilsizeBound = 19;
polydegreeBound = 3;

polydegreeInt = 2;     % 3/6/30/30 worked okay with Lk = 4
overlapTrigger = 6;
stencilsizeInt = 19;
queryFactor = 5;

dampLk = 3;

if closedFlag == 0
    %
    [bigLmatrix bigDampMatrix initdatavec xnodes ynodes rhoVec tripleVec] = ...
        ExeprepRBFHeatLaplace4(N,GAshp,stencilsize,polydegree,...
        xlb,xub,ylb,yub,c1,c2,rho1,rho2,xnodes,ynodes,nearestneighbors,...
        totalmovements,plotflag,numIntNodes,intWrapFlag,...
        polydegreeInt,overlapTrigger,stencilsizeInt,naiveFlag,dampLk,...
        const1,const2,const3,yLI,yUI,queryFactor,RBFwarpFlag,thinFlag,curvedFlag,...
        stencilsizeBound,polydegreeBound);
    %}
    %{
    [bigLmatrix bigDampMatrix initdatavec xnodes ynodes rhoVec tripleVec] = ...
        ExeprepRBFHeatLaplace3(N,GAshp,stencilsize,polydegree,...
        xlb,xub,ylb,yub,c1,c2,rho1,rho2,xnodes,ynodes,nearestneighbors,...
        totalmovements,plotflag,numIntNodes,intWrapFlag,...
        polydegreeInt,overlapTrigger,stencilsizeInt,naiveFlag,dampLk,...
        const1,const2,const3,yLI,yUI,queryFactor,RBFwarpFlag,thinFlag,curvedFlag);
    %}
else
    [bigLmatrix bigDampMatrix initdatavec xnodes ynodes rhoVec tripleVec] = ...
        ExeprepRBFHeatLaplace2(N,GAshp,stencilsize,polydegree,...
        xlb,xub,ylb,yub,c1,c2,rho1,rho2,xnodes,ynodes,nearestneighbors,...
        totalmovements,plotflag,numIntNodes,intWrapFlag,...
        polydegreeInt,overlapTrigger,stencilsizeInt,naiveFlag,dampLk,...
        const1,const2,const3,yLI,yUI,queryFactor,RBFwarpFlag);
end
