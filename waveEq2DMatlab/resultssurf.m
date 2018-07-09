function void = ...
    resultssurf(xnodes,ynodes,datavec,xlb,xub,ylb,yub,safedist,res)

% function RESULTSSURF performs a simple interpolation of RBF-FD data
% values into a meshgrid so the results can be visualized with Matlab's
% SURF function.

% INPUTS:

% xnodes: (N,1) vector of x-coordinates for RBF-FD nodes
% ynodes: "             " y-coordinates "              "
% xlb, etc.: lower and upper boundaries on the vortex problem domain
% safedist: for periodic domains, this is a distance for "tiling" parts of
%           the domain near boundaries so that data interpolation near the
%           (periodic) edges proceeds smoothly.
% res: number of points in x and y for creating the interpolated meshgrid
%      (resolution)

% OUTPUTS:

% graph: dummy variable - not used (I don't know how to make a function 
%        behave without at least one output variable)

%%% Begin function resultssurf %%%

xoverlap = [];
yoverlap = [];
origidx = [];

xspan = xub-xlb;
yspan = yub-ylb;

counter = 1;

% Below, if nodes are located within SAFEDIST of the domain boundary, they 
% are "ghosted" or "tiled" on the other side of the periodic domain so that
% the subsequent KNNSEARCH produces the correct identity of nearest 
% neighbors for each meshgrid point near the domain boundary.

for m = 1:size(xnodes,1)
    
    if xub-xnodes(m,1)<=safedist
        xoverlap(counter,1)=xnodes(m,1)-xspan;
        yoverlap(counter,1)=ynodes(m,1);
        origidx(counter,1) = m;
        counter = counter+1;
    end
    
    if xnodes(m,1)-xlb<=safedist
        xoverlap(counter,1)=xnodes(m,1)+xspan;
        yoverlap(counter,1)=ynodes(m,1);
        origidx(counter,1) = m;
        counter = counter+1;
    end
    
    if yub-ynodes(m,1)<=safedist
        yoverlap(counter,1)=ynodes(m,1)-yspan;
        xoverlap(counter,1)=xnodes(m,1);
        origidx(counter,1) = m;
        counter = counter+1;
    end
    
    if ynodes(m,1)-ylb<=safedist
        yoverlap(counter,1)=ynodes(m,1)+yspan;
        xoverlap(counter,1)=xnodes(m,1);
        origidx(counter,1) = m;
        counter = counter+1;
    end
    
    if xub-xnodes(m,1)<=safedist
        if yub-ynodes(m,1)<=safedist
        xoverlap(counter,1)=xnodes(m,1)-xspan;
        yoverlap(counter,1)=ynodes(m,1)-yspan;
        origidx(counter,1) = m;
        counter = counter+1;
        end
    end
    
    if xub-xnodes(m,1)<=safedist
        if ynodes(m,1)-ylb<=safedist
        xoverlap(counter,1)=xnodes(m,1)-xspan;
        yoverlap(counter,1)=ynodes(m,1)+yspan;
        origidx(counter,1) = m;
        counter = counter+1;
        end
    end
    
    if xnodes(m,1)-xlb<=safedist
        if yub-ynodes(m,1)<=safedist
        xoverlap(counter,1)=xnodes(m,1)+xspan;
        yoverlap(counter,1)=ynodes(m,1)+yspan;
        origidx(counter,1) = m;
        counter = counter+1;
        end
    end
    
    if xnodes(m,1)-xlb<=safedist
        if ynodes(m,1)-ylb<=safedist
        xoverlap(counter,1)=xnodes(m,1)+xspan;
        yoverlap(counter,1)=ynodes(m,1)-yspan;
        origidx(counter,1) = m;
        counter = counter+1;
        end
    end
end


xaug = [xnodes; xoverlap];                  % XAUG and YAUG hold the x-
yaug = [ynodes; yoverlap];                  % and y-coordinates of both the
origidx = [(1:1:size(xnodes,1))'; origidx]; % original and "tiled" RBF-FD
                                            % nodes...
datavecaug = zeros(size(xaug));
datavecaug(1:size(xnodes,1)) = datavec;     % ...and DATAVECAUG holds
                                            % their data values.
for m = size(xnodes,1):size(datavecaug,1)
    datavecaug(m,1)=datavec(origidx(m,1),1);
end

F = TriScatteredInterp(xaug,yaug,datavecaug); % now we can use Matlab's
                                              % TRISCATTEREDINTERP
% function to get a decently good and fast interpolation to a meshgrid
% (below) for visualization of the data via SURF.

graphpointsx = linspace(xlb,xub,res);
graphpointsy = linspace(ylb,yub,res);
[Xgraph Ygraph] = meshgrid(graphpointsx,graphpointsy);
Zgraph = F(Xgraph,Ygraph);

surf(Xgraph,Ygraph,Zgraph)
shading('interp')

%%% End function resultssurf %%%

end

