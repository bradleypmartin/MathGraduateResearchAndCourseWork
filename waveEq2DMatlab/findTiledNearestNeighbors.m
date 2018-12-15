function [origidxstack distances] = ...
    neighborstilinginterp(oldxnodes, oldynodes, newxnodes, newynodes, ...
    xlb, xub, ylb, yub, stencilsize)

% NEIGHBORSTILINGINTERP is a "tiling" method for finding periodic nearest
% neighbors.  It is a variant of NEIGHBORSTILING to be used for finding
% nearest neighbors not only of each RBF node within its OWN set, but more
% generally to find nearest neighbors of each RBF node (NEWNODES) within 
% SOME set of RBF nodes (OLDNODES) - which makes it useful for local RBF 
% interpolation during dynamic node refinement.

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
        
        
        
        
        
        
        