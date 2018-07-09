function [sparsedzeromatrix origidxstack] = ...
    dzeroIMQpluspolyaltpoint(xnodes,ynodes,xeval,yeval,xlb,xub,ylb,yub,...
    shp,stencilsize,polydegree)

% This function handles the RBF-FD-based interpolation from one node set to
% the next. No special interface handling is in place here.

rows = size(xeval,1);
cols = size(xnodes,1);

% We first call NEIGHBORSTILING to find periodic nearest neighbors for all 
% RBF-FD nodes in the domain.  ORIGIDXSTACK holds the indices for the
% STENCILSIZE nearest neighbors (columns) of each of the N RBF-FD nodes
% (rows).

origidxstack = neighborstilinginterp(xnodes,ynodes,xeval,yeval,xlb,xub,ylb,yub,stencilsize);

% the KNNSEARCH below outputs the squared euclidean distance between each
% RBF-FD node and its nearest neighbor in the stack MINDISTANCES2.  This
% will be used to normalize distances in each RBF-FD stencil.

[nearidx mindistances] = knnsearch([xnodes ynodes],[xeval yeval],'k',4);
mindistances2 = mindistances(:,4).^2;
[faridx maxdistance] = ...
    knnsearch([xnodes ynodes],[xeval(1,1) yeval(1,1)],'k',stencilsize+1);
normfactorint = maxdistance(1,stencilsize+1);

sparsei = zeros(stencilsize,size(yeval,1));
sparsej = zeros(stencilsize,size(yeval,1));
sparsedzero = zeros(stencilsize,size(yeval,1));

xspan = xub-xlb;                  
yspan = yub-ylb;                  

onesVec = ones(stencilsize,1);
                                    
parfor m = 1:rows
    
    xpositions = zeros(1,stencilsize);  % within a periodic domain, we can't
    ypositions = zeros(1,stencilsize);
    
    A = ones(stencilsize,stencilsize); % For each RBF-FD
    RHSdzero = zeros(stencilsize,1);                    % stencil, we 
            
    for n = 1:stencilsize
        
        xpositions(1,n)=xnodes(origidxstack(m,n),1);
        ypositions(1,n)=ynodes(origidxstack(m,n),1);
        
        % Below, each node in every stencil is checked for the need to
        % "tile" or "ghost" that node to another part of the domain.  In
        % the future, this step could probably be incorporated into
        % NEIGHBORSTILING for speed.
        
        if abs(xeval(m,1)-(xpositions(1,n)+xspan))<...
                abs(xeval(m,1)-xpositions(1,n))
            xpositions(1,n)=xpositions(1,n)+xspan;
        end
        
        if abs(xeval(m,1)-(xpositions(1,n)-xspan))<...
                abs(xeval(m,1)-xpositions(1,n))
            xpositions(1,n)=xpositions(1,n)-xspan;
        end
        
        if abs(yeval(m,1)-(ypositions(1,n)+yspan))<...
                abs(yeval(m,1)-ypositions(1,n))
            ypositions(1,n)=ypositions(1,n)+yspan;
        end
        
        if abs(yeval(m,1)-(ypositions(1,n)-yspan))<...
                abs(yeval(m,1)-ypositions(1,n))
            ypositions(1,n)=ypositions(1,n)-yspan;
        end
        
    end
    
    distmin2 = mindistances2(m,1);   % DISTMIN2 normalizes distances for
                                     % each RBF-FD stencil so we don't have
    % to choose optimal shape parameters for each stencil separately - with
    % this done, a uniform SHP of 0.4 works quite well.
    
    for i = 1:stencilsize               % creating the RBF (A) matrix...
        for j = 1:stencilsize
            A(i,j) = ...
                1/...
                (1+shp^2*((xpositions(1,j)-xpositions(1,i))^2+...
                (ypositions(1,j)-ypositions(1,i))^2)/distmin2)^0.5;
        end
        
        % ... and the right-hand sides for each stencil of the
        % dx and dy matrices ...
        
        RHSdzero(i,1) = 1/...
                (1+shp^2*((xeval(m,1)-xpositions(1,i))^2+...
                (yeval(m,1)-ypositions(1,i))^2)/distmin2)^0.5;
        
    end
    
    RBFFDweightsdzero = pluspolyweightsdzeroaltpoint(A,RHSdzero,xpositions,ypositions,...
        xeval(m,1),yeval(m,1),polydegree,1/normfactorint);
    
    % Below, weights are put into the condensed dx and dy differentiation
    % matrices in a manner that will match up exactly with indices for all
    % RBF-FD nodes as presented in ORIGIDXSTACK.  This will allow quick
    % matrix multiplication with SPARSEMULTIPLY during time stepping.
     
    sparsei(:,m) = m*onesVec;
    sparsej(:,m) = origidxstack(m,:)';
    sparsedzero(:,m) = RBFFDweightsdzero(1,1:stencilsize)';
    
end

sparsei = reshape(sparsei,size(xeval,1)*stencilsize,1);
sparsej = reshape(sparsej,size(xeval,1)*stencilsize,1);
sparsedzero = reshape(sparsedzero,size(xeval,1)*stencilsize,1);

sparsedzeromatrix = sparse(sparsei,sparsej,sparsedzero,size(xeval,1),size(xnodes,1));

